using WaterLily, GLMakie, StaticArrays, ParametricBodies

# function of asymmetrical wing geometry
function a_wing(; L=32, aoa=15, Ay=0, taper=0.75, sweep=15, span=1.0, center=(0.20, 0.50), Re=1e4, domain=(8,4,2.0), U=1, T=Float32, mem=Array)
    α = T(deg2rad(aoa))
    Λ = T(deg2rad(sweep))
    xc, yc = T.(center)
    n, m, l = T.(domain)
    ampy, ω = T(Ay), T(U/L)
    b = T(span * L)               

    # NACA 6412 points from (http://airfoiltools.com/airfoil/details?airfoil=naca6412-il)
    pnts = L * SA[1.0 0.958 0.906 0.84 0.8 0.7 0.6 0.5 0.4 0.3 0.2 0.12 0.085 0.035 0.0 0.031 0.077 0.1 0.21 0.3 0.4 0.5 0.6 0.7 0.8 0.868 0.9 0.955 1.0;
                    0.0   0.0154 0.0313 0.0506 0.0608 0.0808 0.0983 0.111 0.118 0.116 0.103 0.0823 0.07  0.0445 0.0  -0.0178 -0.0205 -0.0197 -0.0112 -0.0038 0.00182 0.0054 0.00787 0.00853 0.00747 0.0052 0.00386 0.00129 0.0]
    foil = interpNurbs(pnts; p=3) # interpolate a closed airfoil curve

    function extrude(x::SVector{3}, t)
        # If outside the span then no body
        if x[3] > b
            return SA[T(10), T(10)]
        end

        # Scaling due to taper ratio
        chord_scale = 1 / ((taper - 1) * (x[3]/b) + 1)
        return chord_scale * SA[x[1] - (xc * L * n + tan(Λ) * (x[3])), x[2] - (yc * L * m) + ampy * sin(ω*t)] # perform transformation of taper and sweep as wing is extruded
    end

    body = ParametricBody(foil; map=extrude, ndims=3)
    Simulation(round.(Int, (n * L, m * L, l * L)), (U*cos(α), U*sin(α), 0), L; 
         U, exitBC=true, ν=T(U*L/Re), body, T, mem)
end

# Parameters
T = Float32
mem = Array
Re = 5e4 # Reynolds number
L = 16 # length of chord in points
aoa = 10 # angle of attack in degrees
taper = 0.5 # taper ratio (tip_chord/root_chord)
sweep = 15 # leading edge wing sweep in degrees
span = 1.5 # wing span in units of L
U = 1 # velocity
Ay = 4 # amplitude of "heaving" motion in points
domain = (8, 4, 2) # domain size in units of L
center = (0.20, 0.50) # wing origin location in units of domain
azimuth = -3*pi/4 # for GLMakie viz 0.35
elevation = pi/4 # for GLMakie viz 0.5

# run and visualize the simulation
sim = a_wing(; U, aoa, Re, L, taper, sweep, Ay, domain, center, span, T, mem)
fig, ax = viz!(sim; elevation, azimuth) # visualize the geometry

viz!(sim; duration=sim_time(sim)+10, video="ThreeD_AsymmetricalWing.mp4",
    algorithm=:iso, isovalue=0.25, colormap=[:lightblue], alpha=0.9, transparency=true,
    body_color=(:grey, 0.9), body2mesh=false, hidedecorations=true, azimuth, elevation) # remove video="ThreeD_AsymmetricalWing.mp4" for co-visualization during runtime
