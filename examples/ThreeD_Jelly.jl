using WaterLily,StaticArrays,GLMakie

function jelly(L=2^5;Re=5e2,mem=Array,U=1,T=Float32)
    # Define simulation size, geometry dimensions, & viscosity
    R = T(2L/3); h = T(4L-2R); ν = U*R/Re
    # Motion functions
    ω = T(2U/R)
    @fastmath @inline A(t) = 1 .- SA[1,1,0]*cos(ω*t)/10
    @fastmath @inline B(t) = SA[0,0,1]*((cos(ω*t)-1)*R/4-h)
    @fastmath @inline C(t) = SA[0,0,1]*sin(ω*t)*R/4

    # Build jelly from a mapped sphere and plane
    sphere = AutoBody((x,t)->abs(√sum(abs2,x)-R)-1, # sdf
                      (x,t)->A(t).*x+B(t)+C(t))     # map
    plane = AutoBody((x,t)->x[3]-h,(x,t)->x .+ C(t))
    body =  sphere-plane

    # Return initialized simulation
    Simulation((L,L,4L),(0,0,-U),R;ν,body,mem,T)
end

function ω!(arr, sim)
    a = sim.flow.σ
    WaterLily.@inside a[I] = WaterLily.ω_mag(I,sim.flow.u)
    copyto!(arr, a[inside(a)]) # copy to CPU
end

# using CUDA
# make sim and run
sim = jelly(2^5)#; mem=CuArray)
t₀ = sim_time(sim)
duration = 5.0
step = 0.1

viz!(sim;f=ω!,duration,step,video="jelly.mp4",algorithm=:mip,colormap=:algae)