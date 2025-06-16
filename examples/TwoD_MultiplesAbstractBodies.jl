using WaterLily, StaticArrays, ParametricBodies, Plots

function circle_and_foil(L;Re=550,U=1,mem=Array,T=Float32)
    # A moving circle using AutoBody
    circle = AutoBody((x,t)->√sum(abs2,x)-L÷2, (x,t)->x-SA[L+U*t,3L÷2-10])

    # A foil defined by the difference between an upper and lower surface
    upper = ParametricBody(BSplineCurve(SA{T}[5L 4L 3L 3L; 2L 5L÷2 5L÷2 2L],degree=3))
    lower = ParametricBody(BSplineCurve(SA{T}[3L 5L; 2L 2L],degree=1))
    foil = upper ∩ lower

    # Simulate the flow past the two general bodies
    Simulation((10L,4L), (U,0), L; ν=U*L/Re, body=foil+circle, mem, T)
end

# make a simulation and run it
#using CUDA
sim = circle_and_foil(2^6);#mem=CUDA.CuArray);
sim_gif!(sim,duration=8,clims=(-5,5),remeasure=true,plotbody=true,axis=([], false),
         cfill=:seismic,legend=false,border=:none)
