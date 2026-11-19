# Interactive flow past a spinning cylinder: hold ←/→ to change the spin and ↑/↓ to move it.
using WaterLily,StaticArrays,GLMakie

# A RigidMap with constant ω spins the body: the surface velocity ω×r comes from ω, and since the
# geometry is axisymmetric θ is never advanced. ω>0 is counter-clockwise (right-handed about z).
function rotor(D=32;λ=2,Re=200,U=1,T=Float32,mem=Array)
    map = RigidMap(SA{T}[2D,2D],zero(T);ω=T(2λ*U/D)) # λ=ωD/2U: surface speed ratio
    Simulation((8D,4D),(U,0),D;ν=U*D/Re,body=AutoBody((x,t)->√(x'x)-D/2,map),T,mem)
end
spin_ratio(sim) = sim.body.map.ω*sim.L/2sim.U
lift_coeff(sim) = -2WaterLily.total_force(sim)[2]/(sim.U^2*sim.L) # total_force is the force on the fluid

# A steady spin needs no remeasuring. Counter-clockwise spin speeds the flow under the cylinder,
# so the Magnus lift is downward.
sim = rotor()
sim_step!(sim,10;remeasure=false)
println("λ=$(spin_ratio(sim)), C_L=$(round(lift_coeff(sim),digits=2))")

# Interactive: changing ω or moving x₀ changes the body velocity, so remeasure (viz! default).
fig,ax = viz!(sim;clims=(-10,10))
while events(fig).window_open[]
    (;x₀,ω) = sim.body.map; T = eltype(x₀)
    spin = ispressed(fig,Keyboard.left)-ispressed(fig,Keyboard.right)
    V = SA{T}[0,(ispressed(fig,Keyboard.up)-ispressed(fig,Keyboard.down))*sim.U/2]
    dt = T(sim.L/sim.U/10) # viz_step! default step: 0.1 L/U
    sim.body = setmap(sim.body;x₀=x₀+dt*V,V,ω=ω+T(spin*sim.U/sim.L/5))
    viz_step!(fig,sim)
    ax.title = "λ=$(round(spin_ratio(sim),digits=1)), C_L=$(round(lift_coeff(sim),digits=2))"
end
