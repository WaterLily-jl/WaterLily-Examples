using WaterLily,Plots,CUDA,StaticArrays

# good idea if no backgroud flow
WaterLily.CFL(a::Flow) = WaterLily.CFL(a;Δt_max=0.5f0)

# define an ellipse
@inline @fastmath ellipse(x,t;radius=1,Λ=1) = √sum(abs2, SA[x[1]/Λ,x[2]])-radius/Λ

#helper to rotate a vector
@inline @fastmath rotate(v,θ::T) where T = SA{T}[cos(θ) -sin(θ); sin(θ) cos(θ)]*v

function run(L=128;T=Float32,U=1,Λ=4.f0,radius=T(L/2),mem=Array,Re=1e4)

    # fsi parameters
    ρ = 10.f0                               # buoyancy corrected density
    mₐ = SA{T}[π*radius^2,π*radius^2/Λ^2]   # added-mass coefficient ellipse
    m,g = ρ*π*radius^2/Λ,SA{T}[0,-U^2/L]    # mass and gravity scale
    X₀,vel,a0 = SA{T}[5L,9.5L],SA{T}[0,0],SA{T}[0,0]
    # rotation variables
    Iₘ = 0.25f0*m*(radius^2+radius^2/Λ^2)    # mass ellipse
    Iₐ = 0.125f0*π*(radius^2-radius^2/Λ^2)^2 # added mass ellipse m₆₆
    θ,ω,dω = 0.10f0,0.0f0,0.0f0

    # make the sim
    body = AutoBody((x,t)->ellipse(x,t;radius=L/2.f0,Λ=Λ), RigidMap(X₀,θ))
    sim = Simulation((10L,10L),(0,0),L/2.f0;U,ν=U*L/2Re,body,T,mem,Δt=0.05f0)

    @gif for tᵢ in range(0,16.0;step=0.1)
        # update
        while sim_time(sim) < tᵢ
            # the step we are doing and the initial angle
            Δt,θ = sim.flow.Δt[end],sim.body.map.θ
            # compute pressure force and moment in lab frame
            force = -WaterLily.pressure_force(sim)
            moment = -WaterLily.pressure_moment(X₀,sim)[1]
            # transform to body frame
            force,a0 = rotate(force+m.*g,θ),rotate(a0,θ)
            # update linear motion in body frame, and then back to lab frame
            a0 = rotate((force + mₐ.*a0)./(m .+ mₐ), -θ)
            vel += Δt*a0; X₀ += Δt*vel
            # update rotation ODE
            dω = (moment + dω*Iₐ)/(Iₘ + Iₐ)
            ω += Δt*dω; θ += Δt*ω
            # update the body
            sim.body = setmap(sim.body;x₀=SVector{2,T}(X₀),V=SVector{2,T}(vel),θ=T(θ),ω=T(ω))
            # measure and update flow
            sim_step!(sim;remeasure=true)
        end
        # plot vorticity
        @inside sim.flow.σ[I] = WaterLily.curl(3,I,sim.flow.u)*sim.L/sim.U
        @inside sim.flow.σ[I] = ifelse(sdf(sim.body,loc(0,I))<0,0.0,sim.flow.σ[I])
        flood(sim.flow.σ|>Array,shift=(-1.5,-1.5),clims=(-5,5), axis=([], false),
              cfill=:seismic,legend=false,border=:none,size=(1080,1080)); body_plot!(sim)
        println("tU/L=",round(tᵢ,digits=4),", Δt=",round(sim.flow.Δt[end],digits=3),
                " X₂=", round(X₀[2]/sim.L,digits=3), " θ=", round(rad2deg(θ),digits=3))
    end
    return sim
end

using CUDA
sim = run(;mem=CuArray) # run the sim