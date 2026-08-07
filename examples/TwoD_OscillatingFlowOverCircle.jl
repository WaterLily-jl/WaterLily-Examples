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
    sim_gif!(sim; duration=stop, step=0.1, clims=(-8,8), plotbody=true)
end

run_oscillating_flow()