# WaterLily.jl — Agent Quick Reference

**Canonical examples repo:** https://github.com/WaterLily-jl/WaterLily-Examples  
**Main solver repo:** https://github.com/WaterLily-jl/WaterLily.jl  
**Paper:** Weymouth & Font, *Computer Physics Communications* 2025, doi:10.1016/j.cpc.2025.109748

Before writing any WaterLily script, read the relevant example(s) in `examples/` or `notebooks/`. They are always up to date. This file is a map and a list of gotchas — not a substitute for the examples.

---

## Simulation constructor

```julia
Simulation(dims, u_BC, L; ν, body, U=1, T=Float64, mem=Array, uλ=nothing,
           g=nothing, perdir=(), exitBC=false, ϵ=1)
```

- `dims` — Int tuple defining the (interior) array size, e.g. `(3*2^6, 2^7)` for 2D, `(L,L,L)` for 3D. **Use powers of 2 for the multigrid solver.** Actual arrays have a ghost layer, so scalar array size (i.e. pressure) is `ndims .+ 2`, vector array (i.e. velocity) size is `(ndims .+ 2 ..., length(ndims))`.
- `u_BC` — tuple of boundary velocities, one per dimension, e.g. `(U,0)` or `(U,0,0)`. Can also be a function `Ut(i,x,t)` for time-varying BCs.
- `L` — reference length in cells (sets Re via `ν=U*L/Re`)
- `ν` — kinematic viscosity (Greek nu, U+03BD). **Do not use Latin "v".**
- `body` — an `AutoBody` or `AbstractBody`; omit for body-free flow
- `T` — float type; use `T=Float32` for GPU efficiency
- `mem` — array backend: `Array` (CPU), `CuArray` (NVIDIA), `ROCArray` (AMD)
- `uλ(i, xyz)` — optional initial velocity field function; `i` is component index, `xyz` is position vector
- `g(i, x, t)` — optional spatially-uniform body force function
- `perdir` — tuple of periodic directions, e.g. `perdir=(1,)` or `perdir=(1,2)`
- `exitBC=true` — convective outlet BC on x-direction exit face

---

## Defining bodies with AutoBody

```julia
body = AutoBody(sdf)               # static body
body = AutoBody(sdf, map)          # moving/deforming body
```

- `sdf(x, t)` — signed distance function; negative inside body, positive outside
- `map(x, t)` — coordinate transform applied before sdf; returns transformed `x`
- **Always use `StaticArrays.SA[...]` inside sdf/map** to avoid per-cell allocations
- AutoDiff computes normals and velocities automatically from sdf and map

```julia
# Circle (static)
body = AutoBody((x,t) -> √sum(abs2, x .- center) - radius)

# Moving body: map shifts/rotates x before passing to sdf
using StaticArrays
function map(x, t)
    α = amp*cos(t*U/L)
    R = SA[cos(α) sin(α); -sin(α) cos(α)]
    R * (x - SA[3L - L*sin(t*U/L), 4L])
end
```

### Set operations on bodies

```julia
body = bodyA ∩ bodyB     # intersection (\cap<tab>)
body = bodyA ∪ bodyB     # union (\cup<tab>), also: bodyA + bodyB
body = sum(iter) do ...; AutoBody(...); end   # union of many bodies
```

---

## Running a simulation

```julia
sim_step!(sim)                 # take a single time step
sim_step!(sim, t_end)          # advance to time t_end (convective time units)
sim_time(sim)                  # current simulation time
```

**Standard time-loop pattern:**

```julia
t₀ = round(sim_time(sim))
for tᵢ in range(t₀, t₀+duration; step=0.1)
    sim_step!(sim, tᵢ)
    # ... postprocess and plot ...
end
```

---

## Accessing flow fields

```julia
sim.flow.u          # velocity field, shape (n,m,2) for 2D or (n,m,p,3) for 3D
                    # component-last indexing: u[I,i] is component i at CartesianIndex I
sim.flow.p          # pressure field
sim.flow.σ          # scratch scalar field (reuse for derived quantities)
sim.body            # body object
sim.L               # reference length
sim.U               # reference velocity
sim.flow.ν          # viscosity
sim.flow.Δt[end]    # current timestep
```

**Ghost cells:** arrays include one ghost cell layer on each face. Use `inside(a)` to get the interior index range, or use `@inside` for assignments.

```julia
@inside sim.flow.σ[I] = WaterLily.curl(3, I, sim.flow.u)   # 2D vorticity
@inside sim.flow.σ[I] = WaterLily.λ₂(I, sim.flow.u)        # 3D λ₂ criterion
```

For GPU simulations, pipe to CPU before plotting or saving:

```julia
sim.flow.u |> Array
sim.flow.p |> Array
sim.flow.σ |> Array
```

**Staggered grid:** velocity components are face-centred, pressure is cell-centred. To get the velocity magnitude at a cell centre you must average adjacent face values. See `TwoD_LidCavity.jl` for an example `mag(I, u)` helper.

---

## Force and metrics functions

```julia
using WaterLily
total_force(sim)      # total force on all bodies (pressure + viscous)
pressure_force(sim)   # pressure contribution only
viscous_force(sim)    # viscous contribution only
```

`total_force` returns the combined force on **all bodies**. To isolate forces on one body in a multi-body simulation, see `TwoD_MultipleBodies.jl`.

---

## Visualization

**With Plots.jl (2D, gif output):**

```julia
using WaterLily, Plots
sim_gif!(sim; duration=10, clims=(-5,5), plotbody=true)
```

**With GLMakie.jl (2D or 3D, live rendering or save a movie):**

```julia
using WaterLily, GLMakie
function λ₂!(arr, sim)
    a = sim.flow.σ
    @inside a[I] = log10(max(1e-6, -WaterLily.λ₂(I, sim.flow.u)*sim.L/sim.U))
    copyto!(arr, a[inside(a)])
end
viz!(sim; f=λ₂!, duration=15, step=0.05, sym=(1,1,1),
     algorithm=:absorption, colormap=:Reds)
```

`viz!` takes a function `f(arr, sim)` that fills a pre-allocated CPU array `arr` with the scalar field to render. See docstring for details

---

## GPU usage

```julia
# NVIDIA
using CUDA, WaterLily
sim = circle(3*2^6, 2^7; T=Float32, mem=CuArray)

# AMD (requires Julia 1.9+)
using AMDGPU, WaterLily
sim = circle(3*2^6, 2^7; T=Float32, mem=ROCArray)

# Multi-threaded CPU
# Start Julia with: julia --threads auto
# KernelAbstractions detects threads automatically; no code change needed.
```

---

## VTK output (3D post-processing, ParaView)

```julia
using WaterLily, WriteVTK

writer = vtkWriter("my_sim")         # opens my_sim.pvd
save!(writer, sim)                   # writes current state
close(writer)

# Custom fields:
vtk_velocity(a) = a.flow.u |> Array
vtk_pressure(a) = a.flow.p |> Array
vtk_body(a) = (measure_sdf!(a.flow.σ, a.body, WaterLily.time(a)); a.flow.σ |> Array)
writer = vtkWriter("my_sim"; attrib=Dict("Velocity"=>vtk_velocity,
                                          "Pressure"=>vtk_pressure,
                                          "Body"=>vtk_body))
```

**Paraview tip:** arrays include ghost cells. Use ExtractSubset → VOI to trim first and last cell in each direction.

**Tensor fields** need one manual `permutedims` before passing to the writer: `permutedims(T, (4,1,2,3))` for a 3D tensor.

---

## Restart: JLD2

```julia
using WaterLily, JLD2
sim_step!(sim, 20)
save!("checkpoint.jld2", sim)

sim2 = make_sim(...)   # identical constructor call
load!(sim2; fname="checkpoint.jld2")
```

## Restart: VTK

```julia
using WaterLily, ReadVTK
writer = restart_sim!(sim; fname="file_restart.pvd")
save!(writer, sim)    # appends to existing file
close(writer)
```

---

## Mean flow statistics

```julia
meanflow = MeanFlow(sim.flow; uu_stats=true)

while sim_time(sim) < t_max
    sim_step!(sim, sim_time(sim) + 0.1)
    WaterLily.update!(meanflow, sim.flow)
end

meanflow.U   # time-averaged velocity
meanflow.P   # time-averaged pressure
τ = uu(meanflow)   # Reynolds stresses (requires uu_stats=true)
```

---

## Initial velocity field

```julia
function TGV(L; Re=1600, U=1, T=Float32, mem=Array)
    κ, U = T(π/L), T(U)
    function uλ(i, xyz)
        x, y, z = @. xyz * κ
        i==1 && return -U*sin(x)*cos(y)*cos(z)
        i==2 && return  U*cos(x)*sin(y)*cos(z)
        return zero(U)
    end
    Simulation((L,L,L), (0,0,0), L; U, uλ, ν=U*L/Re, T, mem)
end
```

---

## Pressure solver logging

```julia
WaterLily.logger("my_sim")          # start logging residuals to my_sim.log
using Logging; disable_logging(Logging.Debug)   # suppress GPU debug noise
# ... run sim ...
using Plots; plot_logger("my_sim.log")
```

---

## Common mistakes

- **Wrong `ν` symbol:** `ν` is U+03BD (Greek nu), not Latin `v`. Copy from examples.
- **Array sizes not powers of 2:** The geometric multigrid requires this. Use `2^n` or multiples.
- **Plotting GPU arrays directly:** Always do `|> Array` before any CPU-side plot or file write.
- **Forgetting `StaticArrays` inside sdf/map:** Causes massive allocations. Always `using StaticArrays` and use `SA[...]`.
- **`total_force` returns combined force:** It sums over all bodies. See `TwoD_MultipleBodies.jl` for per-body force extraction.
- **Ghost cell indexing:** Use `@inside` or `inside(a)` to restrict to physical cells. Whole-array operations include ghost layers.
- **`sim_step!(sim, t_end)` is absolute, not relative:** It advances to time `t_end`, not by `t_end`.
- **CFL is diffusion-limited at low Re:** `CFL` includes a `5ν` diffusive term, so `Δt ∝ 1/(U/Δx + 5ν/Δx²)`. At low Re (high ν) the timestep is set by diffusion, not convection. 

## Adding no-slip domain BCs

The default domain BCs are reflection with Neumann conditions for the tangential velocity. There are two examples applying no-slip tangential BCs, TwoD_LidCavity.jl and TwoD_Channel.jl. Both use a simple helper function to update the tangential ghost velocities and overwrite `WaterLily.mom_step!` to call this function after each flow update.

---

## Example index

| File | What it demonstrates |
|---|---|
| `TwoD_Circle.jl` | Minimal 2D external flow, pressure logger |
| `TwoD_CirclePeriodicBC.jl` | Periodic BCs, moving body with `mod` in sdf |
| `TwoD_CircleVIV.jl` | FSI with OrdinaryDiffEq.jl, 1-DOF VIV |
| `TwoD_FreeRotatingEllipse.jl` | RigidMap FSI, free rotation |
| `TwoD_MeanCircleJLD2.jl` | JLD2 restart, MeanFlow statistics |
| `TwoD_Hover.jl` | Moving+rotating body via map, manual gif loop |
| `TwoD_Julia.jl` | Complex SDF (Julia logo) |
| `TwoD_LidCavity.jl` | Custom BC!, overriding mom_step! |
| `TwoD_Square.jl` | AutoBody set operations (∩) |
| `TwoD_MultipleBodies.jl` | Union of many AutoBodies, per-body forces |
| `TwoD_MultipleAbstractBodies.jl` | ParametricBodies + AutoBody combined |
| `TwoD_OscillatingFlowOverCircle.jl` | Body force `g`, periodic BCs |
| `TwoD_SlowStartCircle.jl` | Time-varying `u_BC` function |
| `TwoD_TandemFoilOptim.jl` | AutoDiff optimisation through sim |
| `TwoD_Triangle.jl` | Custom analytic SDF |
| `TwoD_UnicodePlots.jl` | Terminal-only visualisation |
| `TwoD_Channel.jl` | Channel flow, periodic and body-force in x, no-slip in y |
| `ThreeD_TaylorGreenVortex.jl` | 3D, `uλ` IC, `T=Float32`, GPU, GLMakie viz |
| `ThreeD_Donut.jl` | 3D torus SDF, GPU, live Makie |
| `ThreeD_Jelly.jl` | 3D moving body, GPU |
| `ThreeD_HoverWriteVTK.jl` | 3D extruded sdf, VTK output |
| `ThreeD_CylinderVTKRestart.jl` | VTK write + restart |
| `ThreeD_SphereLESBiotSavart.jl` | LES, BiotSavart BCs, turbulence statistics |
| `notebooks/Shark.jl` | Pluto notebook: fish SDF + traveling wave map |
