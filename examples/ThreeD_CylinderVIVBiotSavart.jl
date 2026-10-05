using WaterLily,StaticArrays,BiotSavartBCs
using WaterLily: pressure_force

# Vortex-induced vibration of a cylinder on a transverse spring. Biot-Savart BCs on the x,y faces
# allow a snug domain, and perdir=(3,) makes the span infinite.

# Three helical strakes: D/8 segments, 3 cells thick, sticking out of the cylinder and turning once over the span λ
strake((x,y,z),D,λ) = (r=hypot(x,y); θ=mod(atan(y,x)-2π*z/λ+π/3,2π/3)-π/3; hypot(r*cos(θ)-clamp(r*cos(θ),D/2-2,D/2+D/8),r*sin(θ))-3/2)
function viv_cylinder(D;Re=300,U=1,T=Float32,mem=Array,strakes=false)
    sdf(x) = √(x[1]^2+x[2]^2)-D/2
    body = AutoBody((x,t)->strakes ? min(sdf(x),strake(x,D,3D)) : sdf(x),RigidMap(SA{T}[D,3D÷2,0],SA{T}[0,0,0]))
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

using CUDA,GLMakie,Meshing # headless: xvfb-run -s '-screen 0 1920x1080x24' julia --project -t auto <this file>
using GeometryBasics: Cylinder,Tessellation,normal_mesh
D,Re = 64,1000                  # resolved wake on a GPU
t_warm,t_video,t_end = 10,50,80 # video from t_warm to t_video, response history from 0 to t_end

# Vorticity isosurface |ω|D/U=10 as a depth-sorted mesh, in a free 3D scene looking down the span
ωD!(arr,sim) = (a=sim.flow.σ; @inside a[I] = WaterLily.ω_mag(I,sim.flow.u)*sim.L/sim.U; copyto!(arr,@view a[inside(a)]))
function scene(sim)
    fig = Figure(size=(1600,1200)); ax = LScene(fig[1,1]; show_axis=false)
    viz!(sim;fig,ax,f=ωD!,isomesh=10,color=RGBAf(0.35,0.55,0.8,0.6),transparency=true,body2mesh=true,body_color=:grey40)
    # body2mesh leaves the ends open where the body crosses the periodic z faces: cap them with a closed cylinder
    n = size(inside(sim.flow.p)); c = Vec3f(n ./ 2)
    caps(sim) = (x=sim.body.map.x₀ .+ 0.5f0; normal_mesh(Tessellation(Cylinder(Point3f(x[1],x[2],0.5f0),Point3f(x[1],x[2],n[3]+0.5f0),Float32(sim.L/2-0.5)),96)))
    cap = Observable(caps(sim)); mesh!(ax,cap;color=:grey40)
    setcam!() = update_cam!(ax.scene,c+1.6f0*maximum(n)*Vec3f(0.28,0.28,0.92),c,Vec3f(0,1,0))
    refresh!(sim,t) = (viz_step!(fig,sim,t); cap[] = caps(sim); setcam!()) # Makie re-centres on new meshes, so reset the camera
    fig,refresh!
end

# Run the coupled sim to t_end, recording the video between t_warm and t_video. viz! steps the sim itself,
# which would skip spring!, so the sim is stepped here and viz_step! only refreshes the figure.
function run_viv(strakes)
    sim = viv_cylinder(D;Re,strakes,mem=CuArray)
    y₀ = sim.body.map.x₀[2]
    perturb!(sim;noise=0.05) # seed the shedding and spanwise instabilities
    history = Tuple{Float64,Float64,Float64}[]
    function advance!(t)
        CL = 0.
        while sim_time(sim) < t
            CL = spring!(sim,y₀); sim_step!(sim) # sim_step! remeasures the moved body
        end
        push!(history,(t,(sim.body.map.x₀[2]-y₀)/sim.L,CL))
    end
    foreach(advance!,0.1:0.1:t_warm)
    fig,refresh! = scene(sim)
    Makie.record(fig,strakes ? "viv_strakes.mp4" : "viv.mp4",t_warm+0.1:0.1:t_video;framerate=30) do t
        advance!(t); refresh!(sim,t)
    end
    foreach(advance!,t_video+0.1:0.1:t_end)
    return history
end
histories = [run_viv(false),run_viv(true)]

# Compare the responses
fig = Figure(size=(1200,700))
ax1 = Axis(fig[1,1];ylabel="y/D"); ax2 = Axis(fig[2,1];xlabel="Convective time tU/D",ylabel="C_L")
for (h,label,color) in zip(histories,("plain","strakes"),(:black,:dodgerblue))
    lines!(ax1,first.(h),getindex.(h,2);label,color); lines!(ax2,first.(h),last.(h);color)
end
axislegend(ax1;position=:lt)
save("viv_response.png",fig)
