using WaterLily,StaticArrays,GLMakie

function circle(p=4;Re=250,mem=Array,U=1,T=Float32)
    # Define simulation size, geometry dimensions, viscosity
    L=2^p
    center, r, zeroT = SA{T}[3L,3L], T(L), zero(T)
    ν = U*L/Re

    # functions for the body
    norm2(x) = √sum(abs2,x)
    function sdf(x,t)
        norm2(SA[x[1]-center[1],mod(x[2]-6L,6L)-center[2]])-r
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

viz!(sim;duration,step,video="2DCirclePeriodicBC.mp4")

## Alternative visualization using Plots
# using Plots
# @time @gif for tᵢ in range(t₀,t₀+duration;step)
#     # update until time tᵢ in the background
#     sim_step!(sim,tᵢ,remeasure=true)

#     # print time step
#     @inside sim.flow.σ[I] = WaterLily.curl(3,I,sim.flow.u)*sim.L/sim.U
#     flood(sim.flow.σ|>Array,clims=(-10,10),shift=(-0.5,-0.5)); body_plot!(sim)
#     println("tU/L=",round(tᵢ,digits=4),", Δt=",round(sim.flow.Δt[end],digits=3))
# end

