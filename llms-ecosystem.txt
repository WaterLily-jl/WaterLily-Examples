# WaterLily Ecosystem — Agent Repository Guide

This document helps agents select the correct repository in the WaterLily ecosystem.

For API usage, always consult the examples and README of the selected repository.

---

## WaterLily.jl

Repository: https://github.com/WaterLily-jl/WaterLily.jl

Purpose: Core immersed-boundary incompressible Navier–Stokes solver.

Extras: metric & utility functions, RigidMap, visualization extensions

---

## WaterLily-Examples

Repository: https://github.com/WaterLily-jl/WaterLily-Examples

Purpose: Canonical usage examples for the WaterLily ecosystem.

Important: Find closest matching examples before generating code.

---

## WaterLilyMeshBodies.jl

Repository: https://github.com/WaterLily-jl/WaterLilyMeshBodies.jl

Purpose: Import STL and triangle-mesh geometry directly into WaterLily simulations.

Main type: MeshBody

Key features:

* Signed-distance queries from arbitrary triangle meshes.
* Bounding Volume Hierarchy (BVH) acceleration.
* Static and deforming meshes.
* Motion interpolation from mesh snapshots.

Important:

* WaterLily uses unit-cell coordinates; meshes usually require scaling and mapping.
* Closed watertight meshes should use `boundary=true` instead of `boundary=false` default.

---

## ParametricBodies.jl

Repository: https://github.com/WaterLily-jl/ParametricBodies.jl

Purpose: Construct bodies from smooth parametric curves, splines, NURBS, and mappings.

Main types: ParametricBody, PlanarBody, NurbsCurve, DynamicNurbsBody

Key features:

* NURBS support.
* Dynamic control-point motion.
* Extrusion and revolution mappings.

Important:

* Parametric curves should return `SVector`.
* 3D swept and revolved bodies require `ndims=3`.

Additional WaterLily-example: `TwoD_MultipleAbstractBodies.jl`

---

## BiotSavartBCs.jl

Repository: https://github.com/WaterLily-jl/BiotSavartBCs.jl

Purpose: External-flow boundary conditions based on the Biot–Savart equation. 

Main types: BiotSimulation, BiotSavartPoisson

Key features:

* Specialized projection & domain update for the BiotSavartPoisson type
* Reduced dependence on domain size means tiny O(L) domains can be used
* Fast multipole / multilevel acceleration.

Important:

* Mixed Biot–Savart and periodic boundary conditions are not supported.
* Individual domain faces can disable Biot–Savart updates using `nonbiotfaces`.
* Enforcing symmetry conditions on `nonbiotface` requires overwriting `symmetry` function. See README.md, examples/square_sym.jl, and gallery/jelly.jl.

Additional Examples: 

* MeshBodies.jl/examples/*.jl
* WaterLily-examples/ThreeD_SphereLESBiotSavart.jl

References: 
* Weymouth & Lauber (2024), "Using Biot-Savart boundary conditions for unbounded external flow on Eulerian meshes".
* Weymouth & Lauber (2026), "Stability of Kirigami parachutes in effectively infinite numerical domains".

---

## LilyPad.jl

Repository: https://github.com/WaterLily-jl/LilyPad.jl

Purpose: Semi-Lagrangian variant of WaterLily.

Main types: LilyPadSim, LilyFlow

Key features:

* Specialize Semi-Lagrangian predictor-corrector methods for the LilyFlow type
* Stable for any time step size (dt=1.5 default) allowing huge speed-up

Important:

* Currently under active development
* No explicit viscous damping model yet

---

## CoupledSimulations.jl

Repository: https://github.com/WaterLily-jl/CoupledSimulations.jl

Purpose: Framework for coupling multiple simulations and physics models.

See repository README and examples for supported coupling patterns.

---

## Repository Selection Guide

| Task                                | Repository             |
| ----------------------------------- | ---------------------- |
| Standard CFD simulation             | WaterLily.jl           |
| Learn API usage                     | WaterLily-Examples     |
| STL or CAD geometry                 | WaterLilyMeshBodies.jl |
| NURBS, splines, procedural geometry | ParametricBodies.jl    |
| Compact external-flow domains       | BiotSavartBCs.jl       |
| Semi-Lagrangian solver              | LilyPad.jl             |
| Multiphysics or solver coupling     | CoupledSimulations.jl  |
