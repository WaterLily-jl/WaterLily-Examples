using WaterLily,ParametricBodies,StaticArrays,Plots

function main(L=64;Re=1000,U=1,T=Float32,pivot=T(-0.16),thk=T(0.5),mem=Array,tend=50)
    # Define the ellipse geometry, initial angle and rotation-rate
    a,b,x₀,xₚ,θ,ω = T(L/2),T(thk*L/2),SA{T}[2L,2L],SA{T}[pivot*L,0],T(0.1),T(0)
    cps = SA{T}[a a 0 -a -a -a  0  a a
                0 b b  b  0 -b -b -b 0]
    weights = SA[1.,√2/2,1.,√2/2,1.,√2/2,1.,√2/2,1.]
    knots = SA[0,0,0,1/4,1/4,1/2,1/2,3/4,3/4,1,1,1]
    body = ParametricBody(NurbsCurve(cps,knots,weights);map=RigidMap(x₀,θ;xₚ,ω))

    # make sim
    sim = Simulation((6L,4L),(U,0),L;ν=U*L/Re,body,T,mem)

    # Set dynamical properties
    m,ma = π*a*b,π*a^2                   # mass (density = 1) and added mass
    Is = π*a^3*b/4+m*T(pivot*L)^2        # solid moment of inertia
    Ia = π/8*(a^2-b^2)^2+ma*T(pivot*L)^2 # added moment of inertia
    α = zero(T)                          # initial angular acceleration

    # Integrate and plot
    @time @gif for tᵢ in range(0,tend;step=0.1)
        # update until time tᵢ in the background
        while sim_time(sim) < tᵢ
            # the step we are doing
            Δt = sim.flow.Δt[end]
            # get moment
            moment = -WaterLily.pressure_moment(x₀+xₚ, sim.flow, sim.body)[1]
            # update rotation ODE
            α = (moment + Ia * α) / (Is + Ia)
            ω += Δt * α; θ += Δt * ω # Verlet
            # update the body
            sim.body = setmap(sim.body; θ=T(θ), ω=T(ω))
            # measure and update flow
            sim_step!(sim;remeasure=true)
        end
        @inside sim.flow.σ[I] = WaterLily.curl(3, I, sim.flow.u) * sim.L / sim.U
        @inside sim.flow.σ[I] = ifelse(abs(sim.flow.σ[I])<0.001,0.0,sim.flow.σ[I])
        flood(sim.flow.σ|>Array,shift=(-1.5,-1.5),clims=(-5,5),axis=([],false),
              cfill=:seismic,legend=false,border=:none,size=(6sim.L,4sim.L)); body_plot!(sim)
        plot!(title="tU/L $tᵢ, θ=$(round(rad2deg(θ), digits=1))°")
        println("tU/L=", round(tᵢ, digits=4), ", θ=", round(rad2deg(θ), digits=1), "°, ω=", round(ω*sim.L/sim.U, digits=3))
    end
    return sim
end

# run the main function
# using CUDA
chaotic = main();#mem=CuArray);
stable_invert = main(Re=250,thk=0.12,pivot=0.25,tend=15);#,mem=CuArray);