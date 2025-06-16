using WaterLily,StaticArrays,Plots

function circle(n,m;Re=550,U=1,mem=Array,T=Float32)
    R, x0 = m÷18, m÷2
    # 6 circles with random centers x,y ∈ [5,5] and radius r ∈ [0.75,1.5]
    body = sum(1:6) do _
        center,radius = 10rand(SVector{2,T}) .- 5, (3rand(T) + 3)/4
        AutoBody((x,t)->√sum(abs2, x .- x0 - center*R) - radius*R)
    end
    # make a simulation
    Simulation((n,m), (U,0), R; ν=U*R/Re, body, mem, T)
end

# make a simulation and run it
# using CUDA
sim = circle(3*2^7,2^8);#mem=CUDA.CuArray);
sim_gif!(sim,duration=30,clims=(-5,5),remeasure=false,plotbody=true,axis=([], false),
         cfill=:seismic,legend=false,border=:none)

# get net force on all bodies
f_total = WaterLily.pressure_force(sim)

# force on _each_ body
each_force(flow,body) = WaterLily.pressure_force(flow,body)
each_force(flow,body::WaterLily.SetBody{typeof(min)}) = mapreduce(bod->each_force(flow,bod),hcat,(body.a,body.b))
each_force(sim.flow,sim.body)