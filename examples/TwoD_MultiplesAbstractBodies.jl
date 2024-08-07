using WaterLily
using StaticArrays
using ParametricBodies
include("../src/TwoD_plots.jl")

function circle(L;Re=550,U=1,mem=Array,T=Float32)
    bodies = AbstractBody[]
    # a first circle using autoBody 
    push!(bodies,AutoBody((x,t)->√sum(abs2, x .- SA_F32[3.75L,2.5L]) - L/2))
    # another one from parametricbodies
    cps = SA_F32[1 1 0 -1 -1 -1  0  1 1
                0 1 1  1  0 -1 -1 -1 0]*L/2 .+ [2L,1.5L]
    weights = SA_F32[1.,√2/2,1.,√2/2,1.,√2/2,1.,√2/2,1.]
    knots =   SA_F32[0,0,0,1/4,1/4,1/2,1/2,3/4,3/4,1,1,1]
    # make a nurbs curve and a body for the simulation
    push!(bodies,ParametricBody(NurbsCurve(cps,knots,weights)))
    # combine into one body
    body = Bodies(bodies, [+])
    # make a simulation
    Simulation((10L,4L), (U,0), L/2; ν=U*L/Re, body, mem, T)
end

# make a simulation and run it
sim = circle(2^6,mem=Array)
sim_gif!(sim,duration=30,clims=(-5,5),plotbody=true,axis=([], false),
         cfill=:seismic,legend=false,border=:none)

# # get force on first body
# f1 = WaterLily.pressure_force(sim.flow,sim.body.bodies[1])
# flood(sim.flow.f[:,:,1],clims=(-1,1)) # check that this is only non-zero near the first body

# # force on the second
# f2 = WaterLily.pressure_force(sim.flow,sim.body.bodies[2])
