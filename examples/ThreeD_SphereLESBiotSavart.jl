using Printf, StaticArrays, CUDA, JLD2, WaterLily, Pkg
Pkg.add(url="https://github.com/weymouth/BiotSavartBCs.jl")
using BiotSavartBCs
using WaterLily: dot, sgs!

# utils for command line arguments
iarg(arg) = occursin.(arg, ARGS) |> findfirst
iarg(arg, args) = occursin.(arg, args) |> findfirst
arg_value(arg) = split(ARGS[iarg(arg)], "=")[end]
arg_value(arg, args) = split(args[iarg(arg, args)], "=")[end]
metaparse(x) = eval(Meta.parse(x))
# SGS model
smagorinsky(I::CartesianIndex{m} where m; S, Cs, Δ) = @views (Cs*Δ)^2*sqrt(2dot(S[I,:,:],S[I,:,:])) # Smagorinsky-Lilly SGS model
# I/O
function save_sim(sim, meanflow, force, t)
    t_str = @sprintf("%i", sim_time(sim))
    save!(fname_save * "_t$(t_str).jld2", sim)
    save!(fname_save * "_t$(t_str)_meanflow.jld2", meanflow)
    jldsave(fname_save * "_t$(t_str)_force.jld2"; force, t)
end
function read_forces(fname)
    obj = jldopen(fname)
    return obj["force"], obj["t"]
end
# Sphere
function sphere(m; n=5m÷2, R=m÷3, U=1, Re=3700, T=Float32, mem=CuArray)
    body = AutoBody((x,t)->√sum(abs2,x .- m÷2)-R)
    BiotSimulation((n,m,m), (U,0,0), 2R; ν=U*2R/Re, body, T, mem)
end
# Metrics
recirculation_length(U) = findfirst(axes(U,1)) do i
    (U[i,m÷2,m÷2,1] < 0) && (U[i+1,m÷2,m÷2,1] > 10eps(eltype(U)))
end |> i -> (i- m÷2)/2R

function run_sim(p; m=3*2^p, R=m÷3, U=1, Re=3700, mem=CuArray, T=Float32, restart=nothing, restart_stats=nothing, restart_force=nothing)
    sim = sphere(m; R, U, Re, T, mem)
    S = zeros(T, size(sim.flow.p)..., ndims(sim.flow.p), ndims(sim.flow.p)) |> mem # working array holding a tensor for each cell
    force, t = Vector{T}[] ,T[] # force coefficients, time

    if !isnothing(restart)
        load!(sim; fname=restart)
        println("Loaded: $restart")
    end
    sim_step!(sim, stats_init; remeasure=false, verbose=true, λ, udf, νₜ=smagorinsky, S, Cs, Δ)

    meanflow = MeanFlow(sim.flow; uu_stats=true)
    if !isnothing(restart_stats)
        load!(meanflow; fname=restart_stats)
        println("Loaded: $restart_stats")
        force, t = read_forces(restart_force)
        println("Loaded: $restart_force")
    end

    println("\nComputing mean flow statistics from: T=$(sim_time(sim))\n\
    Total accumulated: T=$(WaterLily.time(meanflow)*sim.U/sim.L)\n\
    Remaining: T=$(time_max-stats_init-WaterLily.time(meanflow)*sim.U/sim.L)\n")

    next_save = sim_time(sim) + save_interval
    while sim_time(sim) < time_max
        sim_step!(sim, sim_time(sim)+stats_interval; remeasure=false, verbose=false, λ, udf, νₜ=smagorinsky, S, Cs, Δ)
        sim_info(sim)

        WaterLily.update!(meanflow, sim.flow)
        push!(force, WaterLily.total_force(sim)/(0.5*sim.U^2*sim.L^2))
        push!(t, sim_time(sim))

        if WaterLily.sim_time(sim) > next_save || sim_time(sim) > time_max
            save_sim(sim, meanflow, force, t)
            next_save = sim_time(sim) + save_interval
            println("Saved simulation and mean flow statistics.")
        end

    end
    return sim, meanflow, mapreduce(permutedims, vcat, force), t
end

mem = CuArray
T = Float32
Re = 3700
time_max = 300 # in CTU
stats_init = 100 # in CTU
stats_interval = 0.1 # in CTU
save_interval = 50 # in CTU

p = !isnothing(iarg("p", ARGS)) ? arg_value("p", ARGS) |> metaparse : 4
explicit_sgs = !isnothing(iarg("smag", ARGS)) ? arg_value("smag", ARGS) |> metaparse : true
m = 3*2^p
R = m÷3
U = 1
udf = explicit_sgs ? sgs! : nothing
λ = explicit_sgs ? cds : quick
Cs = T(0.17) # Smagorinsky constant
Δ = sqrt(1^2 + 1^2 + 1^2) |> T # Filter width

datadir = !isnothing(iarg("datadir", ARGS)) ? arg_value("datadir", ARGS) |> String : explicit_sgs ? "data_p$(p)_sgs" : "data_p$(p)"
fname_save = joinpath(datadir,"p$(p)")
restart = nothing #  "data/p4_t100.jld2"
restart_stats = nothing # "data/p4_t100_meanflow.jld2"
restart_force = nothing # "data/p4_t100_force.jld2"

function main(;load_time=nothing)
    if isnothing(load_time)
        mkpath(datadir)
        println("Running: p=$p, udf=$udf, λ=$λ, datadir=$datadir")
        return run_sim(p; m, R, U, mem, Re, T, restart, restart_stats, restart_force)
    else
        t_str = @sprintf("%i", load_time)
        sim = sphere(m; R, U, Re, T, mem)
        fname = "$(fname_save)_t$(t_str).jld2"
        load!(sim; fname)
        println("Loaded: $fname")
        meanflow = MeanFlow(sim.flow; uu_stats=true)
        fname = "$(fname_save)_t$(t_str)_meanflow.jld2"
        load!(meanflow; fname)
        println("Loaded: $fname")
        fname = "$(fname_save)_t$(t_str)_force.jld2"
        force, t = read_forces(fname)
        println("Loaded: $fname")
        return sim, meanflow,  mapreduce(permutedims, vcat, force), t
    end
end

sim, meanflow, force, t = main() # load_time=200
println("▷ ΔT [CTU] = "*@sprintf("%.2f", WaterLily.time(meanflow)/2R*U))
println("▷ CD_mean = "*@sprintf("%.2f", sum(force[2:end,1].*diff(t))/sum(diff(t))))
println("▷ L/D = "*@sprintf("%.2f", recirculation_length(meanflow.U|>Array)))

## Visualization -------
using GLMakie
viz!(sim) # 3D
viz!(sim; d=2) # 2D
# Run and visualize
S = zeros(T, size(sim.flow.p)..., ndims(sim.flow.p), ndims(sim.flow.p)) |> mem
viz!(sim;duration=10,remeasure=false,λ,udf,udf_kwargs=(:S=>S,:νₜ=>smagorinsky,:Cs=>Cs,:Δ=>Δ))

# Visualize the meanflow
sim_meanflow = deepcopy(sim)
WaterLily.copy!(sim_meanflow.flow,meanflow)
sim_meanflow.flow.μ₁ .= uu(meanflow)
function U_viz(a, sim)
    b = sim.flow.σ
    @inside b[I] = sim.flow.u[I,1]
    copyto!(a,b[inside(b)])
end
function uu_viz(a, sim)
    b = sim.flow.σ
    @inside b[I] = sim.flow.μ₁[I,1,1]
    copyto!(a,b[inside(b)])
end
viz!(sim_meanflow;f=U_viz,d=2,levels=20)
viz!(sim_meanflow;f=uu_viz,d=2,clims=(0.00552,0.0554),threshhold=1e-5,levels=12) # same contour levels as Rodriguez et al. 2011 https://doi.org/10.1017/jfm.2011.136