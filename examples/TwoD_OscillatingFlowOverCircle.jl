using WaterLily
using StaticArrays
using Plots
using CUDA

function circle(n,m;κ=1.5,Re=250,U=1,T=Float32,mem=Array)
    # define a circle at the domain center
    radius = T(m/8)
    U, κ = T(U), T(κ)
    body = AutoBody((x,t)->√sum(abs2, x .- (n/2,m/2)) - radius)

    # define time-varying body force `g` and periodic direction `perdir`
    accelScale, timeScale = U^2/2radius, κ*radius/U
    g(i,x,t::T) where {T} = i==1 ? -2accelScale*sin(t/timeScale) : zero(T)
    Simulation((n,m), (U,zero(T)), radius; ν=U*radius/Re, body, g, perdir=(1,), T, mem)
end

function run_oscillating_flow(n=392, stop=20)
    sim = circle(n,n)
    sim_step!(sim,0.1)

    @time @gif for tᵢ in range(0.,stop;step=0.2)
        println("tU/L=",round(tᵢ,digits=4))
        sim_step!(sim,tᵢ)
        @inside sim.flow.σ[I] = WaterLily.curl(3,I,sim.flow.u)*sim.L/sim.U
        @inside sim.flow.σ[I] = ifelse(abs(sim.flow.σ[I])<0.001,0.0,sim.flow.σ[I])
        # It's important to have `|>Array` during GPU simulation as `flood` only accept CPU Array input
        flood(sim.flow.σ|>Array,shift=(-2,-1.5),clims=(-8,8), axis=([], false),
              cfill=:seismic,legend=false,border=:none,size=(n,n))
        body_plot!(sim)
    end
end

run_oscillating_flow()