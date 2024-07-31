using GLMakie
GLMakie.activate!()
function makie_video!(makie_plot,sim,dat,obs_update!;remeasure=false,name="file.mp4",duration=1,step=0.1,framerate=30,compression=20)
    # Set up viz data and figure
    obs = obs_update!(dat,sim) |> Observable;
    f = makie_plot(obs)
    
    # Run simulation and update figure data
    t₀ = round(sim_time(sim))
    t = range(t₀,t₀+duration;step)
    Makie.record(f, name, t; framerate, compression) do tᵢ
        sim_step!(sim,tᵢ;remeasure)
        obs[] = obs_update!(dat,sim)
        println("simulation ",round(Int,(tᵢ-t₀)/duration*100),"% complete")
    end
    return f
end

function flow_λ₂!(dat,sim)
    a = sim.flow.σ
    @inside a[I] = max(0,log10(-min(-1e-6,WaterLily.λ₂(I,sim.flow.u)*(sim.L/sim.U)^2))+.25)
    copyto!(dat,a[inside(a)])                  # copy to CPU
end
function flow_λ₂(sim)
    dat = sim.flow.σ[inside(sim.flow.σ)] |> Array
    flow_λ₂!(dat,sim)
    dat
end

using Meshing
function geom!(md,d,sim,t=WaterLily.time(sim))
    a = sim.flow.σ
    WaterLily.measure_sdf!(a,sim.body,t)
    copyto!(d,a[inside(a)]) # copy to CPU
    mirrorto!(md,d)         # mirror quadrant
    alg = Meshing.MarchingCubes()
    ranges = range.((0, 0, 0), size(md))
    points, faces = Meshing.isosurface(md, alg, ranges...)
    p3f = Point3f.(points)
    gltriangles = GLMakie.GLTriangleFace.(faces)
    return GLMakie.normal_mesh(p3f, gltriangles)
end

function ω!(md,d,sim)
    a,dt = sim.flow.σ,sim.L/sim.U
    @inside a[I] = WaterLily.ω_mag(I,sim.flow.u)*dt
    copyto!(d,a[inside(a)]) # copy to CPU
    mirrorto!(md,d)         # mirror quadrant
end