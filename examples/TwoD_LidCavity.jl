using WaterLily,GLMakie

# velocity magnitude
mag(I,u) = √sum(ntuple(i->0.25*(u[I,i]+u[I+δ(i,I),i])^2,length(I)))

import WaterLily: size_u,@loop,slice
# import explicitly BC function and overwrite
function BC_lid!(a)
    N,n = WaterLily.size_u(a)
    for j ∈ 1:n, i ∈ 1:n
        if i==1 && j==2 # lid, Dirichlet cannot be imposed, must interpolate u[i,j+1,1]+u[i,j]/2 = uBC
            @loop a[I,i] = 2.00-a[I-δ(j,I),i] over I ∈ slice(N,N[j],j)
            @loop a[I,i] = -a[I+δ(j,I),i] over I ∈ slice(N,1,j)
        elseif i==j # Normal direction, homoheneous Dirichlet
            @loop a[I,i] = 0.0 over I ∈ slice(N,1,j)
            @loop a[I,i] = 0.0 over I ∈ slice(N,N[j],j)
        else  # Tangential directions, interpolate ghost cell to homogeneous Dirichlet
            @loop a[I,i] = -a[I+δ(j,I),i] over I ∈ slice(N,1,j)
            @loop a[I,i] = -a[I-δ(j,I),i] over I ∈ slice(N,N[j],j)
        end
    end
end
import WaterLily: mom_step!,scale_u!,conv_diff!,quick,BDIM!,project!,CFL
# overwrite the mom_step! function
@fastmath function mom_step!(a::Flow{N},b::AbstractPoisson;udf=nothing,kwargs...) where N
    a.u⁰ .= a.u; scale_u!(a,0)
    # predictor u → u'
    conv_diff!(a.f,a.u⁰,a.σ,quick,ν=a.ν,perdir=a.perdir)
    BDIM!(a); BC_lid!(a.u)
    project!(a,b); BC_lid!(a.u)
    # corrector u → u¹
    conv_diff!(a.f,a.u,a.σ,quick,ν=a.ν,perdir=a.perdir)
    BDIM!(a); scale_u!(a,0.5); BC_lid!(a.u)
    project!(a,b,0.5); BC_lid!(a.u)
    push!(a.Δt,CFL(a))
end

function Lid_cavity(;L=2^6,U=1.0,Re=100,T=Float32,mem=Array)
    # make the simulations
    Simulation((L,L),(0,0),L;U=U,ν=U*L/Re,T,mem)
end

# using CUDA
sim = Lid_cavity(L=2^6)#,mem=CuArray)

# get start time
function umag(arr, sim)
    a = sim.flow.σ
    @inside a[I] = mag(I,sim.flow.u)
    copyto!(arr ,a[inside(a)])
end

duration, step = 20, 0.1
viz!(sim; f=umag, duration, step, clims=(0,1), levels=20) # add: video="lid.mp4" to store the video

## Alternative visualization with Plots
# using Plots
# t₀ = round(sim_time(sim))
# @gif for tᵢ in range(t₀,t₀+duration;step)

#     # update until time tᵢ in the background
#     sim_step!(sim,tᵢ)

#     # flood plot
#     @inside sim.flow.σ[I] = mag(I,sim.flow.u)
#     flood(sim.flow.σ|>Array; shift=(-0.5,-0.5),clims=(0,1))

#     # print time step
#     println("tU/L=",round(tᵢ,digits=4),", Δt=",round(sim.flow.Δt[end],digits=3))
# end