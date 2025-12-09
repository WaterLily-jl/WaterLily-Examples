using WaterLily, JLD2, GLMakie, CUDA, Printf

function circle(n, m; Re=250, U=1, mem=Array, T=Float32)
    radius, center = m/8, m/2
    body = AutoBody((x,t)->√sum(abs2, x .- center) - radius)
    Simulation((n,m), (U,0), radius; ν=U*radius/Re, body, mem, T)
end

sim = circle(3*2^6, 2^7; mem=CuArray)
sim_step!(sim, 40; verbose=true) # warm up
perturb!(sim)
sim_step!(sim, sim_time(sim) + 40; verbose=true) # warm up
viz!(sim; levels=60) # visualize

function run(time_max)
    next_save = sim_time(sim) + save_interval
    while sim_time(sim) < time_max
        sim_step!(sim, sim_time(sim)+stats_interval; remeasure=false, verbose=false)
        sim_info(sim)
        WaterLily.update!(meanflow, sim.flow)

        if WaterLily.sim_time(sim) > next_save || sim_time(sim) > time_max
            save!("circle_t$(@sprintf "%.2f" sim_time(sim)).jld2", sim)
            save!("circle_mean_t$(@sprintf "%.2f" sim_time(sim)).jld2", meanflow)
            next_save = sim_time(sim) + save_interval
            println("Saved simulation and mean flow statistics.")
        end
    end
end

# Init mean flow stats and run for 20 CTU
stats_interval = 0.1 # in CTU
save_interval = 30 # in CTU
meanflow = MeanFlow(sim.flow; uu_stats=true)
run(sim_time(sim) + 20)

# Visualise mean flow (not yet converged in time!)
ω_mean = zeros(size(meanflow.P))|>CuArray
@inside ω_mean[I] = WaterLily.curl(3,I,meanflow.U)
viz!(sim, ω_mean; clims=(-0.2,0.2), levels=60)

# Create a new simulation and restore from last saved snapshot
sim2 = circle(3*2^6, 2^7; mem=CuArray)
load!(sim2; fname="circle_t$(@sprintf "%.2f" sim_time(sim)).jld2")
meanflow2 = MeanFlow(sim2.flow; uu_stats=true)
load!(meanflow2; fname="circle_mean_t$(@sprintf "%.2f" sim_time(sim)).jld2")

# Check that the restart works
@assert all(sim2.flow.u .≈ sim.flow.u)
@assert all(meanflow2.U .≈ meanflow.U)

# Average flow for 100 CTU more
run(sim_time(sim) + 80)
@inside ω_mean[I] = WaterLily.curl(3,I,meanflow.U)
viz!(sim, ω_mean; clims=(-0.2,0.2), levels=60) # the mean flow should look better now
