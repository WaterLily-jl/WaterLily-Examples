using WaterLily,StaticArrays,BiotSavartBCs
using WaterLily: pressure_force

# Vortex-induced vibration of a cylinder on a transverse spring. Biot-Savart BCs on the x,y faces
# allow a snug domain, and perdir=(3,) makes the span infinite.

# Three helical strakes: D/4 segments, 3 cells thick, sticking out of the cylinder and turning once over the span λ
strake((x,y,z),D,λ) = (r=hypot(x,y); θ=mod(atan(y,x)-2π*z/λ+π/3,2π/3)-π/3; hypot(r*cos(θ)-clamp(r*cos(θ),D/2-2,D/2+D/4),r*sin(θ))-3/2)
function viv_cylinder(D;Re=300,U=1,T=Float32,mem=Array,strakes=false)
    sdf(x) = √(x[1]^2+x[2]^2)-D/2
    body = AutoBody((x,t)->strakes ? min(sdf(x),strake(x,D,6D)) : sdf(x),RigidMap(SA{T}[2D,2D,0],SA{T}[0,0,0]))
    BiotSimulation((6D,4D,6D),(U,0,0),D;body,ν=U*D/Re,perdir=(3,),T,mem)
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

using CUDA,GLMakie,Meshing # headless: xvfb-run -s '-screen 0 2560x1600x24' julia --project -t auto <this file>
using GeometryBasics: Cylinder,Tessellation,normal_mesh
using JLD2 # save!/load! restart files
D,Re,level = 32,300,8           # resolved wake on a GPU; level is the |ω|D/U of the isosurface
t_warm,t_video,t_end = 10,50,80 # video from t_warm to t_video, response history from 0 to t_end

# Vorticity isosurface as a depth-sorted mesh, in a free 3D scene looking down the span
ωD!(arr,sim) = (a=sim.flow.σ; @inside a[I] = WaterLily.ω_mag(I,sim.flow.u)*sim.L/sim.U; copyto!(arr,@view a[inside(a)]))
function scene(sim)
    fig = Figure(size=(1600,1200)); ax = LScene(fig[1,1]; show_axis=false)
    viz!(sim;fig,ax,f=ωD!,isomesh=level,color=RGBAf(0.35,0.55,0.8,0.6),transparency=true,body2mesh=true,body_color=:grey40,verbose=false)
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
# The history (<name>.csv) and a restart file (<name>.jld2) are written as soon as the run ends.
function run_viv(strakes)
    name = strakes ? "viv_strakes" : "viv"
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
        isinteger(round(t;digits=1)) && (println(name,": t=",round(t;digits=1),", y/D=",round(history[end][2];digits=3)); flush(stdout)) # stdout is buffered in a log file
    end
    foreach(advance!,0.1:0.1:t_warm)
    fig,refresh! = scene(sim)
    Makie.record(fig,name*".mp4",t_warm+0.1:0.1:t_video;framerate=30) do t
        advance!(t); refresh!(sim,t)
    end
    GLMakie.closeall() # release the GL screen before the next run
    foreach(advance!,t_video+0.1:0.1:t_end)
    open(name*".csv","w") do io
        println(io,"t,y/D,C_L"); foreach(h->println(io,join(h,",")),history)
    end
    save!(name*".jld2",sim) # restart file: p, u and Δt (load! with fname). The body's position and velocity are saved alongside
    jldsave(name*"_body.jld2";x₀=sim.body.map.x₀,V=sim.body.map.V)
    return history
end
histories = [run_viv(false),run_viv(true)]

# Compare the responses with Plots (GR, no display needed) rather than another GL screen
ENV["GKSwstype"] = "100"
import Plots
p = Plots.plot(layout=(2,1),size=(1200,700),xlabel=["" "Convective time tU/D"],ylabel=["y/D" "C_L"])
for (h,label,color) in zip(histories,("plain","strakes"),(:black,:dodgerblue))
    Plots.plot!(p[1],first.(h),getindex.(h,2);label,color)
    Plots.plot!(p[2],first.(h),last.(h);label="",color)
end
Plots.savefig(p,"viv_response.png")
