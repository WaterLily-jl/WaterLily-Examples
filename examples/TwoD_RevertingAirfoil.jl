using WaterLily,ParametricBodies,StaticArrays,Plots

function main(L=64;Re=250,U=1,T=Float32,pivot=T(L/4),a=T(L/2),mem=Array)
    # Define foil as a thickened line segment with an initial angle/rotation-rate
    x₀,xₚ,θ,ω = SA{T}[2L, 2L],SA{T}[pivot,0],T(0.1),T(0)
    thk(s) = 1.2f0L*(0.2969f0*√s - 0.126f0*(s) - 0.3516f0*(s)^2 + 0.2843f0*(s)^3 - 0.1036f0*(s)^4)
    body = ParametricBody(BSplineCurve(SA[-a a;0 0]);map=RigidMap(x₀,θ;xₚ,ω),thk,boundary=false)

    # make sim
    sim = Simulation((6L,4L),(U,0),L;ν=U*L/Re,body,T,mem)

    # Set dynamical properties
    m,ma = T(0.12)*L^2,π*a^2   # mass (density = 1) and added mass
    Is = m*a^2/3 + m*pivot^2   # solid moment of inertia
    Ia = ma*a^2/8 + ma*pivot^2 # added moment of inertia
    α = zero(T)                # initial angular acceleration

    # Integrate and plot
    @time @gif for tᵢ in range(0,20;step=0.1)
        # update until time tᵢ in the background
        while sim_time(sim) < tᵢ
            # the step we are doing
            Δt = sim.flow.Δt[end]
            # get moment
            moment = -WaterLily.pressure_moment(x₀+xₚ, sim.flow, sim.body)[1]
            # update rotation ODE
            α = (moment + Ia * α) / (Is + Ia)
            θ += Δt * (ω + Δt * α / 2)
            ω += Δt * α
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
        println("tU/L=", round(tᵢ, digits=4), ", θ=", round(rad2deg(θ), digits=1), "°, ω=", round(ω, digits=3))
    end
    return sim
end

# run the main function
# using CUDA
sim = main(); #mem=CuArray);