using Printf, StaticArrays, CUDA, JLD2, WaterLily, GLMakie, BiotSavartBCs
using WaterLily: dot, sgs!

smagorinsky(I::CartesianIndex{m} where m; S, Cs, Δ) = @views (Cs*Δ)^2*sqrt(2dot(S[I,:,:],S[I,:,:])) # define the Smagorinsky-Lilly model
function save_sim(sim, meanflow)
    t_str = @sprintf("%i", sim_time(sim))
    save!(fname_save * "_t$(t_str).jld2", sim)
    save!(fname_save * "_t$(t_str)_meanflow.jld2", meanflow)
end
function sphere(m; n=5m÷2, R=m÷3, U=1, Re=3700, T=Float32, mem=CuArray)
    body = AutoBody((x,t)->√sum(abs2,x .- m÷2)-R)
    BiotSimulation((n,m,m), (U,0,0), 2R; ν=U*2R/Re, body, T, mem)
end

function run_sim(p; m=3*2^p, Re=3700, mem=CuArray, T=Float32, restart=nothing, restart_stats=nothing)
    sim = sphere(m; Re, T, mem)
    S = zeros(T, size(sim.flow.p)..., ndims(sim.flow.p), ndims(sim.flow.p)) |> mem # working array holding a tensor for each cell

    if !isnothing(restart)
        load!(sim; fname=restart)
        println("Loaded: $restart")
    end
    sim_step!(sim, stats_init; remeasure=false, verbose=true, λ, udf, νₜ=smagorinsky, S, Cs, Δ)

    meanflow = MeanFlow(sim.flow; uu_stats=true)
    if !isnothing(restart_stats)
        load!(meanflow; fname=restart_stats)
        println("Loaded: $restart_stats")
    end

    println("\nComputing mean flow statistics from: T=$(sim_time(sim))\n\
    Total accumulated: T=$(WaterLily.time(meanflow)*sim.U/sim.L)\n\
    Remaining: T=$(time_max-stats_init-WaterLily.time(meanflow)*sim.U/sim.L)\n")

    next_save = sim_time(sim) + save_interval
    while sim_time(sim) < time_max
        sim_step!(sim, sim_time(sim)+stats_interval; remeasure=false, verbose=false, λ, udf, νₜ=smagorinsky, S, Cs, Δ)
        sim_info(sim)
        WaterLily.update!(meanflow, sim.flow)

        if WaterLily.sim_time(sim) > next_save || sim_time(sim) > time_max
            save_sim(sim, meanflow)
            next_save = sim_time(sim) + save_interval
            println("Saved simulation and mean flow statistics.")
        end

    end
    return sim, meanflow
end

mem = CuArray
T = Float32
Re = 3700
time_max = 150 # in CTU
stats_init = 50 # in CTU
stats_interval = 0.1 # in CTU
save_interval = 25 # in CTU

p = 4
udf = sgs!
λ = cds
Cs = T(0.17) # Smagorinsky constant
Δ = sqrt(1^2 + 1^2 + 1^2) |> T # Filter width

datadir = "data"
fname_save = joinpath(datadir,"p$(p)")
restart = nothing #  "data/p4_t100.jld2"
restart_stats = nothing # "data/p4_t100_meanflow.jld2"

function main()
    mkpath(datadir)
    println("Running: p=$p, udf=$udf, λ=$λ")
    return run_sim(p; mem, Re, T, restart, restart_stats)
end

sim, meanflow = main()

## Visualization 3D
# viz!(sim) # 3D
# viz!(sim; d=2) # 2D
## Run and visualize
# S = zeros(T, size(sim.flow.p)..., ndims(sim.flow.p), ndims(sim.flow.p)) |> mem
# viz!(sim;duration=10,remeasure=false,λ,udf,udf_kwargs=(:S=>S,:νₜ=>smagorinsky,:Cs=>0.17,:Δ=>Δ))