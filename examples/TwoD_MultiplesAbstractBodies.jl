using WaterLily, StaticArrays, ParametricBodies, GLMakie, CUDA

function circle_and_foil(L=2^6;Re=550,U=1,mem=CUDA.CuArray,T=Float32)
    # A moving circle using AutoBody
    circle = AutoBody((x,t)->√sum(abs2,x)-L÷2, (x,t)->x-SA[L+U*t,3L÷2-10])

    # A foil defined by the intersection of the upper and lower surface
    upper = ParametricBody(BSplineCurve(SA{T}[5L 4L 3L 3L; 2L 5L÷2 5L÷2 2L],degree=3))
    lower = ParametricBody(BSplineCurve(SA{T}[3L 5L; 2L 2L],degree=1))
    foil = upper ∩ lower

    # Simulate the flow past the two general bodies
    Simulation((10L,4L), (U,0), L; ν=U*L/Re, body=foil+circle, mem, T)
end

# make a simulation and run it
viz!(circle_and_foil(),duration=9,clims=(-0.1,0.1))
