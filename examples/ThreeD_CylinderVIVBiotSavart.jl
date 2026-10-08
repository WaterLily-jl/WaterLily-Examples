# Vortex-induced vibration of a cylinder on a transverse spring, with and without helical strakes.
using WaterLily,StaticArrays,BiotSavartBCs
using WaterLily: pressure_force
using CUDA,GLMakie,Meshing

# Create the signed-distance function
function strakes_sdf((x,y,z),D,pitch)
    angle_from_strake = mod(atan(y,x)-2π*z/pitch+π/3,2π/3)-π/3
    across,along = hypot(x,y) .* sincos(angle_from_strake)
    hypot(along-clamp(along,D/2-2,D/2+D/4),across)-3/2
end

# Biot-Savart BCs on the x,y faces allow a snug domain, and perdir=(3,) makes the span infinite.
function viv_cylinder(D;Re,strakes=false,U=1,T=Float32,mem=Array)
    cylinder_sdf(x) = hypot(x[1],x[2])-D/2
    sdf(x,t) = (x = x-SA[2D,2D,0]; strakes ? min(cylinder_sdf(x),strakes_sdf(x,D,6D)) : cylinder_sdf(x))
    body = AutoBody(sdf,RigidMap(SA{T}[0,0,0],SA{T}[0,0,0]))
    BiotSimulation((6D,4D,6D),(U,0,0),D;body,ν=U*D/Re,perdir=(3,),T,mem)
end

# Move the cylinder on its transverse spring over one time step:
# The mass ratio m*=m/(ρπD²/4 span), reduced velocity U*=U/(fₙD) and damping ratio ζ set the spring,
# where fₙ is the natural frequency in still fluid, so the stiffness includes the added mass.
span(sim) = size(sim.flow.p,3)-2
function spring!(sim;mstar=4,Ustar=5,ζ=0)
    D,U,T = sim.L,sim.U,eltype(sim.flow.p)
    dt,x₀,V = sim.flow.Δt[end],sim.body.map.x₀,sim.body.map.V
    mass,ωₙ = mstar*π*D^2/4*span(sim), 2π*U/(Ustar*D)
    stiffness,damping = (1+1/mstar)*mass*ωₙ^2,2mass*ζ*ωₙ
    lift = -pressure_force(sim)[2]
    acceleration = (lift-stiffness*x₀[2]-damping*V[2])/mass
    V += SA{T}[0,dt*acceleration,0]; x₀ += dt*V
    sim.body = setmap(sim.body;x₀,V)
end

# Add the dynamics and displacement logging to Simulation.measure!
const history = Tuple{Float64,Float64}[]
function WaterLily.measure!(sim::Simulation,t=sum(sim.flow.Δt))
    spring!(sim)
    push!(history,(sim_time(sim),sim.body.map.x₀[2]/sim.L))
    measure!(sim.flow,sim.body;t,ϵ=sim.ϵ)
    WaterLily.update!(sim.pois)
end

# Run the simulations with and without strakes on CUDA if available
D,Re,iso_level,duration = 32,300,8,60
mem = CUDA.functional() ? CuArray : Array
histories = map((false,true)) do strakes
    sim = viv_cylinder(D;Re,strakes,mem); perturb!(sim;noise=0.05); empty!(history)
    viz!(sim;duration,video=strakes ? "viv_strakes.mp4" : "viv.mp4", isomesh=iso_level/D,
         color=RGBAf(0.35,0.55,0.8,0.6),transparency=true,body2mesh=true,body_color=:grey40,hidedecorations=true)
    copy(history)
end

# Plot motion history
fig = Figure(size=(1200,400))
ax = Axis(fig[1,1];xlabel="Convective time tU/D",ylabel="y/D")
for (history,label,color) in zip(histories,("plain","strakes"),(:black,:dodgerblue))
    lines!(ax,first.(history),last.(history);label,color)
end
axislegend(ax;position=:lt)
save(file,"viv_response.png")