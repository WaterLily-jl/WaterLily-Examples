using WaterLily,StaticArrays

function circle(p=4;Re=250,mem=Array,U=1,T=Float32)
    # Define simulation size, geometry dimensions, viscosity
    L=2^p
    center, r, zeroT = SA{T}[3L,0.5L], T(L), zero(T)
    ν = U*L/Re

    # functions for the body
    norm2(x) = √sum(abs2,x)
    function sdf(x,t)
        norm2(SA[x[1]-center[1],mod(x[2]-center[2]+3L,6L)-3L])-r
    end
    function map(x,t)
        x.-SA[zeroT,U*t/2]
    end
    # make a body
    body = AutoBody(sdf,map)

    # return sim
    Simulation((8L,6L),(U,0),L;ν,body,mem,perdir=(2,),exitBC=true)
end

# using CUDA
sim = circle(5)#;mem=CuArray) to run on GPU
t₀ = sim_time(sim)
duration = 40.0
step = 0.1

using GLMakie
viz!(sim;duration,step,video="2DCirclePeriodicBC.mp4")

# Alternative visualization using Plots
# using Plots
# sim_gif!(sim;duration,step,video="2DCirclePeriodicBC.mp4",plotbody=true,remeasure=true)

