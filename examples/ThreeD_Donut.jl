using WaterLily,StaticArrays,GLMakie

function donut(L;Re=1e3,mem=Array,U=1)
    # Define simulation size, geometry dimensions, viscosity
    center,R,r = SA[L/2,L/2,L/2], L/4, L/16
    ν = U*R/Re

    # Apply signed distance function for a torus
    norm2(x) = √sum(abs2,x)
    body = AutoBody() do xyz,t
        x,y,z = xyz - center
        norm2(SA[x,norm2(SA[y,z])-R])-r
    end

    # Initialize simulation
    Simulation((2L,L,L),(U,0,0),R;ν,body,mem)
end

function ω_θ!(arr, sim)
    dt,a = sim.L/sim.U, sim.flow.σ
    center = SA{eltype(sim.flow.σ)}[2sim.L,2sim.L,2sim.L]
    @inside a[I] = WaterLily.ω_θ(I,(1,0,0),center,sim.flow.u)*dt
    copyto!(arr, a[inside(a)])
end

# make sim and run
# using CUDA
sim = donut(2^5)#;mem=CUDA.CuArray);
t₀ = sim_time(sim)
duration = 10.0
step = 0.25

viz!(sim;f=ω_θ!,duration,step,video="donut.mp4",algorithm=:iso,isovalue=0.5) # remove video="donut.mp4" for co-visualization during runtime