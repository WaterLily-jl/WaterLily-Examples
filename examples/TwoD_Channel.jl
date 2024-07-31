using WaterLily
include("../src/TwoD_plots.jl")

# velocity magnitude
mag(I,u) = √sum(ntuple(i->0.25*(u[I,i]+u[I+δ(i,I),i])^2,length(I)))

# import explicitly BC function and overwrite
function BC_CH!(a,perdir=())
    N,n = WaterLily.size_u(a)
    for j ∈ 1:n, i ∈ 1:n
        if j in perdir
            @WaterLily.loop a[I,i] = a[WaterLily.CIj(j,I,N[j]-1),i] over I ∈ WaterLily.slice(N,1,j)
            @WaterLily.loop a[I,i] = a[WaterLily.CIj(j,I,2),i] over I ∈ WaterLily.slice(N,N[j],j)
        else
            if i==1 && j==2 # u-velocity on top and bottom boundaries
                @WaterLily.loop a[I,i] =  a[I-δ(j,I),i] over I ∈ WaterLily.slice(N,N[j],j) # d(ux)/dy = 0
                @WaterLily.loop a[I,i] = -a[I+δ(j,I),i] over I ∈ WaterLily.slice(N,1,j)    # ux = 0 on bottom (no slip), must interpolate
            else # v-velocity on top and bottom boundaries
                @WaterLily.loop a[I,i] = 0.0 over I ∈ WaterLily.slice(N,1,j)    # v = 0
                @WaterLily.loop a[I,i] = 0.0 over I ∈ WaterLily.slice(N,N[j],j) # v = 0, no slip
            end
        end
    end
end
@fastmath function WaterLily.mom_step!(a::Flow{N},b::AbstractPoisson) where N
    a.u⁰ .= a.u; WaterLily.scale_u!(a,0)
    # predictor u → u'
    WaterLily.conv_diff!(a.f,a.u⁰,a.σ,ν=a.ν,perdir=a.perdir)
    WaterLily.accelerate!(a.f,@view(a.Δt[1:end-1]),a.g,a.U)
    WaterLily.BDIM!(a); BC_CH!(a.u,a.perdir)
    WaterLily.project!(a,b); BC_CH!(a.u,a.perdir)
    # corrector u → u¹
    WaterLily.conv_diff!(a.f,a.u,a.σ,ν=a.ν,perdir=a.perdir)
    WaterLily.accelerate!(a.f,@view(a.Δt[1:end-1]),a.g,a.U)
    WaterLily.BDIM!(a); WaterLily.scale_u!(a,0.5); BC_CH!(a.u,a.perdir)
    WaterLily.project!(a,b,0.5); BC_CH!(a.u,a.perdir)
    push!(a.Δt,WaterLily.CFL(a))
end

# Initialize simulation
L = 2^7
U = 1.0
Re = 100
# pressure gradient required to drive the flow to u~1
g(i, t) = i == 1 ? U^2/(L/2)^2 : 0

# using CUDA
sim = Simulation((4L,L),(0.0,0.0),L;U=U,ν=U*L/Re,g,perdir=(1,),mem=Array)

# get start time
t₀ = round(sim_time(sim))
duration = 50; step = 0.1

anim = @animate for tᵢ in range(t₀,t₀+duration;step)

    # update until time tᵢ in the background
    sim_step!(sim,tᵢ)

    # flood plot
    @inside sim.flow.σ[I] = mag(I,sim.flow.u)
    flood(sim.flow.σ[inside(sim.flow.σ)]|>Array; shift=(-0.5,-0.5),clims=(0,1))

    # print time step
    println("tU/L=",round(tᵢ,digits=4),", Δt=",round(sim.flow.Δt[end],digits=3))
end
gif(anim,"chanel_flow.gif",fps=15)

# # plot the profile
# plot(sim.flow.u[:,2:end,1]'./maximum(sim.flow.u[:,:,1]),collect(1:129)./129,
#     label=:none,xlabel="u/U",ylabel="y/L")
# # analytical solution laminar chanel flow u/U ~ C*(y/L)-(y/L)^2 ∀ y ∈ [0,L/2]
# plot!(4*(collect(0:0.01:0.5).-collect(0:0.01:0.5).^2),2collect(0:0.01:0.5),label="analytical")
# savefig("chanel_flow.png")