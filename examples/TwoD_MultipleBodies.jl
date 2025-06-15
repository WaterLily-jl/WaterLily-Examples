using WaterLily,StaticArrays,Plots

function circle(n,m;Re=550,U=1,mem=Array,T=Float32)
    R, x0 = m/18, m/2+1
    # random position x,y ∈ [-2.5,2.5] and circle diamater r ∈ [0.75,1.5]
    body = sum(zip(eachrow(5rand(6,2).-2.5),0.75rand(6).+0.75)) do (center,radius)
        AutoBody((x,t)->√sum(abs2, x .- x0 .- 2center.*R) - radius*R)
    end
    # make a simulation
    Simulation((n,m), (U,0), R; ν=U*R/Re, body, mem, T),bodies
end

# make a simulation and run it
sim,bodies = circle(3*2^7,2^8,mem=Array);
sim_gif!(sim,duration=1,clims=(-5,5),remeasure=false,plotbody=true,axis=([], false),
         cfill=:seismic,legend=false,border=:none)

# get net force on all bodies
f_total = WaterLily.pressure_force(sim)

# force on _each_ body
each_force(sim,body) = WaterLily.pressure_force(sim)
each_force(sim,body::SetBody) = WaterLily.pressure_force(sim)

f_each = mapreduce(body->WaterLily.pressure_force(sim.flow,body) for body in bodies]

# check that this is only non-zero near a single body
WaterLily.pressure_force(sim.flow,bodies[1]);
flood(sim.flow.f[:,:,1],clims=(-1,1))