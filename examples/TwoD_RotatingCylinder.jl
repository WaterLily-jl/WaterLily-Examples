using StaticArrays
using WaterLily
using Plots

function cylinder(n, m; Re=250, U=1, ω=0, mem=Array, T=Float32)
    # Define a circular cylinder at the center of the domain.
    radius = T(m / 8)
    center = SVector{2,T}(T(m / 2), T(m / 2))

    # Define the signed-distance function in the cylinder's body frame
    sdf(x, t) = sqrt(sum(abs2, x)) - radius

    # Define the rigid-body map and angular velocity
    map = RigidMap(center, zero(T); ω=T(ω))
    body = AutoBody(sdf, map)

    # Make the simulation
    Simulation(
        (n, m),
        (U, 0),
        2radius;
        T,
        ν=U * 2radius / Re,
        body,
        mem
    )
end

# Make stationary and rotating simulations
stationary = cylinder(
    3 * 2^5,
    2^6;
    Re=100,
    U=1,
    ω=0
)

# Positive ω sets the cylinder's rotation direction
rotating = cylinder(
    3 * 2^5,
    2^6;
    Re=100,
    U=1,
    ω=0.5
)

# Run the stationary-cylinder simulation
sim_gif!(
    stationary;
    duration=10,
    step=0.2,
    clims=(-5, 5),
    plotbody=true
)

# Run the rotating-cylinder simulation
sim_gif!(
    rotating;
    duration=10,
    step=0.2,
    clims=(-5, 5),
    plotbody=true
)