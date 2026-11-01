#
# This example depends on a WaterLily#pathline-viz! and the unregistered packages LilyPad and Pathlines:
#
# pkg> activate --temp  # or a persistent environment you want to keep
# pkg> add WaterLily#pathline-viz!
# pkg> add https://github.com/WaterLily-jl/LilyPad.jl, https://github.com/WaterLily-jl/Pathlines.jl
# julia> include("examples/TwoD_InteractiveCircle.jl")
#
# Once it's running, click on the Makie window and use the arrow keys to move the circle around.
#
using WaterLily,GLMakie,StaticArrays,Pathlines
function circle(n,m;Re=100,U=1,T=Float32)
    # signed distance function to circle and Rigid-body mapping
    radius, center = m/8, SVector{2,T}(m/2, m/2)
    sdf(x,t) = √(x'x) - radius
    map = RigidMap(center,zero(T))

    Simulation((n,m),   # domain size
               (U,0),   # domain velocity (& velocity scale)
               2radius; # length scale
               T,       # float type
               ν=U*2radius/Re,         # fluid viscosity
               body=AutoBody(sdf,map)) # geometry
end
circ = circle(3*2^5,2^6);
fig,ax = viz!(circ,remeasure=true); # remeasure=true is needed for interative motion of the body

begin
    # Move circle's center with arrow keys
    center = Ref(circ.body.map.x₀)
    on(events(fig).keyboardbutton) do event
        event.action == Keyboard.press || return
        if     event.key == Keyboard.up;    center[] += SVector{2,Float32}(0,0.1)
        elseif event.key == Keyboard.down;  center[] -= SVector{2,Float32}(0,0.1)
        elseif event.key == Keyboard.right; center[] += SVector{2,Float32}(0.1,0)
        elseif event.key == Keyboard.left;  center[] -= SVector{2,Float32}(0.1,0)
        else; return; end
        circ.body = setmap(circ.body; x₀=center[])
    end

    # Advance viz
    while events(fig).window_open[]
        viz_step!(fig, circ) # default step size is 0.1, can be overridden with viz_step!(fig, circ, t)
    end
end