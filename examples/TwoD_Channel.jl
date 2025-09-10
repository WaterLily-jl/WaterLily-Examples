using WaterLily,Plots

# velocity magnitude
mag(I,u) = √sum(ntuple(i->0.25*(u[I,i]+u[I+δ(i,I),i])^2,length(I)))

import WaterLily: size_u,@loop,slice,CIj
# import explicitly BC function and overwrite
function BC_CH!(a,perdir=())
    N,n = size_u(a)
    for j ∈ 1:n, i ∈ 1:n
        if j in perdir
            @loop a[I,i] = a[CIj(j,I,N[j]-1),i] over I ∈ slice(N,1,j)
            @loop a[I,i] = a[CIj(j,I,2),i] over I ∈ slice(N,N[j],j)
        else
            if i==1 && j==2 # u-velocity on top and bottom boundaries
                @loop a[I,i] =  a[I-δ(j,I),i] over I ∈ slice(N,N[j],j) # d(ux)/dy = 0
                @loop a[I,i] = -a[I+δ(j,I),i] over I ∈ slice(N,1,j)    # ux = 0 on bottom (no slip), must interpolate
            else # v-velocity on top and bottom boundaries
                @loop a[I,i] = 0.0 over I ∈ slice(N,1,j)    # v = 0
                @loop a[I,i] = 0.0 over I ∈ slice(N,N[j],j) # v = 0, no slip
            end
        end
    end
end
import WaterLily: mom_step!,scale_u!,conv_diff!,quick,accelerate!,BDIM!,project!,CFL
# overwrite the mom_step! function
@fastmath function mom_step!(a::Flow{N},b::AbstractPoisson;udf=nothing,kwargs...) where N
    a.u⁰ .= a.u; scale_u!(a,0); t₁ = sum(a.Δt); t₀ = t₁ - a.Δt[end]
    # predictor u → u'
    conv_diff!(a.f,a.u⁰,a.σ,quick,ν=a.ν,perdir=a.perdir)
    accelerate!(a.f,t₀,a.g,a.uBC)
    BDIM!(a); BC_CH!(a.u,a.perdir)
    project!(a,b); BC_CH!(a.u,a.perdir)
    # corrector u → u¹
    conv_diff!(a.f,a.u,a.σ,quick,ν=a.ν,perdir=a.perdir)
    accelerate!(a.f,t₁,a.g,a.uBC)
    BDIM!(a); scale_u!(a,0.5); BC_CH!(a.u,a.perdir)
    project!(a,b,0.5); BC_CH!(a.u,a.perdir)
    push!(a.Δt,CFL(a))
end

function channel(;L=2^6,U=1.0,Re=100,T=Float32,mem=Array)

    # pressure gradient required to drive the flow to u~1
    g(i,x,t) = convert(typeof(t), i == 1 ? U^2/(L/2)^2 : 0)

    # make the simulations
    Simulation((4L,L),(0,0),L;U=U,ν=U*L/Re,g,perdir=(1,),T,mem)
end

# using CUDA
# make the sim
sim = channel(L=2^7)#;mem=CuArray)

# get start time
duration = 50; step = 0.1

@gif for tᵢ in range(0,duration;step)

    # update until time tᵢ in the background
    sim_step!(sim,tᵢ)

    # flood plot
    @inside sim.flow.σ[I] = mag(I,sim.flow.u)
    flood(sim.flow.σ[inside(sim.flow.σ)]|>Array; shift=(-0.5,-0.5),clims=(0,1))

    # print time step
    println("tU/L=",round(tᵢ,digits=4),", Δt=",round(sim.flow.Δt[end],digits=3))
end

# plot the profile
u_num = sim.flow.u[:,2:end,1]'./maximum(sim.flow.u[:,:,1])|>Array
Plots.plot(u_num,collect(1:sim.L+1)./sim.L,
     label=:none,xlabel="u/U",ylabel="y/L")
# analytical solution laminar chanel flow u/U ~ C*(y/L)-(y/L)^2 ∀ y ∈ [0,L/2]
Plots.plot!(4*(collect(0:0.01:0.5).-collect(0:0.01:0.5).^2),2collect(0:0.01:0.5),label="analytical")
savefig("chanel_flow.png")