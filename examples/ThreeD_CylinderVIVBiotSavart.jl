using WaterLily,StaticArrays,BiotSavartBCs
using WaterLily: pressure_force

# Vortex-induced vibration of a cylinder on a transverse spring. Biot-Savart BCs on the x,y faces
# allow a snug domain, and perdir=(3,) makes the span infinite.
function viv_cylinder(D;Re=300,U=1,T=Float32,mem=Array)
    body = AutoBody((x,t)->√(x[1]^2+x[2]^2)-D/2,RigidMap(SA{T}[D,3D÷2,0],SA{T}[0,0,0]))
    BiotSimulation((4D,3D,3D),(U,0,0),D;body,ν=U*D/Re,perdir=(3,),T,mem)
end
span(sim) = size(sim.flow.p,3)-2

# Integrate m ÿ + c ẏ + k y = F_y over one time step and update the body map, where the
# mass ratio m*=m/(ρπD²/4 span), reduced velocity U*=U/(fₙD) and damping ratio ζ set m, k and c.
# fₙ is the natural frequency in still fluid, including the added mass m/m*.
function spring!(sim,y₀;mstar=4,Ustar=5,ζ=0)
    D,U,T = sim.L,sim.U,eltype(sim.flow.p)
    m = mstar*π*D^2/4*span(sim); ωₙ = 2π*U/(Ustar*D)
    dt,x₀,V = sim.flow.Δt[end],sim.body.map.x₀,sim.body.map.V
    Fy = -pressure_force(sim)[2]
    a = T((Fy-(1+1/mstar)*m*ωₙ^2*(x₀[2]-y₀)-2m*ζ*ωₙ*V[2])/m)
    V += SA{T}[0,dt*a,0]; x₀ += dt*V
    sim.body = setmap(sim.body;x₀,V)
    return Fy/(sim.L*span(sim)*U^2/2) # lift coefficient
end

sim = viv_cylinder(16); # D=16 runs on a laptop CPU in ~10 minutes. 
# using CUDA; sim = viv_cylinder(64;mem=CuArray); # For a resolved wake use D=64 on a GPU
y₀ = sim.body.map.x₀[2]
perturb!(sim;noise=0.05) # seed the shedding and spanwise instabilities
history = map(0.1:0.1:80) do t
    CL = 0.
    while sim_time(sim) < t
        CL = spring!(sim,y₀); sim_step!(sim) # sim_step! remeasures the moved body
    end
    (t,(sim.body.map.x₀[2]-y₀)/sim.L,CL)
end

using Plots
t,y,CL = first.(history),getindex.(history,2),last.(history)
response = plot(t,[y CL],label=["y/D" "C_L"],xlabel="Convective time")
@inside sim.flow.σ[I] = WaterLily.curl(3,I,sim.flow.u)*sim.L/sim.U
plane = CartesianIndices((2:size(sim.flow.p,1)-1,2:size(sim.flow.p,2)-1,span(sim)÷2:span(sim)÷2)) # mid-span cut
wake = flood(WaterLily.squeeze(Array(sim.flow.σ[plane])),legend=false,axis=false)
body_plot!(sim;CIs=plane)
plot(response,wake,layout=(2,1),size=(600,700))
