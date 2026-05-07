using WaterLily,BiotSavartBCs,Plots,CUDA,StaticArrays,WriteVTK

# dummy overwrite
import WaterLily: @loop,scale_u!,conv_diff!,udf!,accelerate!,BDIM!,CFL

# good idea if no backgroud flow
WaterLily.CFL(a::Flow) = WaterLily.CFL(a;Δt_max=0.5f0)

# Biot-Savart momentum step with U and acceleration prescribed
import BiotSavartBCs: biot_mom_step!,biot_project!
function biot_mom_step_fall!(a::Flow{N},b,ω...;λ=quick,udf=nothing,fmm=true,U,kwargs...) where N
    a.u⁰ .= a.u; scale_u!(a,0); t₁ = sum(a.Δt); t₀ = t₁-a.Δt[end]
    # predictor u → u'
    @log "p"
    conv_diff!(a.f,a.u⁰,a.σ,λ,ν=a.ν)
    udf!(a,udf,t₀; kwargs...)
    BDIM!(a);
    biot_project!(a,b,ω...,U;fmm) # new
    # corrector u → u¹
    @log "c"
    conv_diff!(a.f,a.u,a.σ,λ,ν=a.ν)
    udf!(a,udf,t₁; kwargs...)
    BDIM!(a); scale_u!(a,0.5)
    biot_project!(a,b,ω...,U;fmm,w=0.5) # new
    push!(a.Δt,CFL(a))
end

import WaterLily: @loop
# falling body acceleration term
fall!(flow,t;acceleration) = for i ∈ 1:ndims(flow.p)
    @loop flow.f[I,i] += acceleration[i] over I ∈ CartesianIndices(flow.p)
end

# define an ellipse
@inline @fastmath ellipse(x,t;radius=1,Λ=1) = √sum(abs2, SA[x[1]/Λ,x[2]])-radius/Λ

#helper to rotate a vector
@inline @fastmath rotate(v,θ::T) where T = SA{T}[cos(θ) -sin(θ); sin(θ) cos(θ)]*v

function run(L=128;T=Float32,U=1,Λ=4.f0,radius=T(L/2),mem=Array,Re=1e4)

    # we know have to define these functions, for example:
    vtk_vorticity(a::AbstractSimulation) = (@inside a.flow.σ[I] = WaterLily.curl(3,I,a.flow.u)*a.L/a.U; a.flow.σ |> Array)
    vtk_body(a::AbstractSimulation) = (measure_sdf!(a.flow.σ, a.body, WaterLily.time(a.flow)); a.flow.σ |> Array)
    vtk_position(a::AbstractSimulation) = (f=zeros(size(a.flow.σ)...,3); f[:,:,1] .= X₀[1]; f[:,:,2] .= X₀[2]; f |> Array)

    # we might want to write some custom stuff, for this we have to create a `Dictionary`  of custom attributes
    # with the keys being the names of the field in the vtk file and the values being the functions that return the data to an `Array`.
    custom_write_attributes = Dict("ω" => vtk_vorticity,"d" => vtk_body, "pos"=>vtk_position)

    # now we prepare the vtk writer
    wr = vtkWriter("falling_ellipse"; attrib=custom_write_attributes)

    # fsi parameters
    ρ = 10.f0                               # buoyancy corrected density
    mₐ = SA{T}[π*radius^2,π*radius^2/Λ^2]   # added-mass coefficient ellipse
    m,g = ρ*π*radius^2/Λ,SA{T}[0,-U^2/L]    # mass and gravity scale
    X₀,vel,a0 = SA{T}[3L,1.5L],SA{T}[0,0],SA{T}[0,0]
    # rotation variables
    Iₘ = 0.25f0*m*(radius^2+radius^2/Λ^2)    # mass ellipse
    Iₐ = 0.125f0*π*(radius^2-radius^2/Λ^2)^2 # added mass ellipse m₆₆
    θ,ω,dω = 0.10f0,0.0f0,0.0f0

    # make the sim
    body = AutoBody((x,t)->ellipse(x,t;radius=L/2.f0,Λ=Λ), RigidMap(X₀,θ))
    sim = BiotSimulation((6L,6L),(0,0),L/2.f0;U,ν=U*L/2Re,body,T,mem,Δt=0.05f0)
    Xₘ = copy(X₀) # the moment point is constant in lab frame

    @gif for tᵢ in range(0,20.0;step=0.1)
        # update
        while sim_time(sim) < tᵢ
            # the step we are doing and the initial angle
            Δt,θ = sim.flow.Δt[end],sim.body.map.θ
            # compute pressure force and moment in lab frame
            force = -WaterLily.pressure_force(sim)
            moment = -WaterLily.pressure_moment(Xₘ,sim)[1]
            # transform to body frame
            force,a0 = rotate(force+m.*g ,θ),rotate(a0,θ)
            # update linear motion in body frame, and then back to lab frame
            a0 = rotate((force + mₐ.*a0)./(m .+ mₐ), -θ)
            vel += Δt*a0; X₀ += Δt*vel
            # update rotation ODE
            dω = (moment + dω*Iₐ)/(Iₘ + Iₐ)
            ω += Δt*dω; θ += Δt*ω
            # update the body
            sim.sim.body = setmap(sim.sim.body;θ=T(θ),ω=T(ω))
            # remeasure the sim and update the flow
            measure!(sim);
            biot_mom_step_fall!(sim.flow,sim.pois,sim.ω,sim.x₀,sim.tar,sim.ftar;
                                fmm=sim.fmm,U=-vel,udf=fall!,acceleration=-a0)
        end
        # plot vorticity
        @inside sim.flow.σ[I] = WaterLily.curl(3,I,sim.flow.u)*sim.L/sim.U
        @inside sim.flow.σ[I] = ifelse(sdf(sim.body,loc(0,I))<0,0.0,sim.flow.σ[I])
        flood(sim.flow.σ|>Array,shift=(-1.5,-1.5),clims=(-5,5), axis=([], false),
              cfill=:seismic,legend=false,border=:none,size=(1080,1080)); body_plot!(sim)
        println("tU/L=",round(tᵢ,digits=4),", Δt=",round(sim.flow.Δt[end],digits=3),
                " X₂=", round(X₀[2]/sim.L,digits=3), " θ=", round(rad2deg(θ),digits=3))
        save!(wr, sim)
    end
    close(wr)
    return sim
end

using CUDA
sim = run(;mem=CuArray) # run the sim