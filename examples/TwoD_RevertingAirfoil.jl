using WaterLily,StaticArrays,Plots

function main(L=64;Re=250,U=1,thick=L/16.f0,T=Float32,mem=Array)

    # Pivot at quarter-chord point
    center,pivot = SA[2.f0L, 2.f0L],SA[L/4.f0,0.f0]

    # Moment of inertia for a thin plate, rotating around pivot
    # using parallel axis theorem: I_pivot = I_center + m*d^2
    mass = L * thick # assuming density = 1
    Ia = (mass*(L^2+thick^2)/12.f0 + mass*(L/4.f0)^2)*10.f0
    moment,α,θ,ω = 0.f0,0.f0,0.1f0,0.f0

    # SDF for a foil (thick plate) centered at `center`, aligned with x-axis
    thickness(s) = 0.6f0*(0.2969f0*√s - 0.126f0*(s) - 0.3516f0*(s)^2 + 0.2843f0*(s)^3 - 0.1036f0*(s)^4)*L
    function sdf(x, t)
        ξ = clamp(x[1], -L/2.f0, L/2.f0);
        s_norm = (ξ + L/2.f0) / L;
        p = x .- SA[ξ, 0.f0];
        return √sum(abs2, p) - thickness(s_norm)
    end
    body = AutoBody(sdf, RigidMap(center, θ; xₚ=pivot))

    # make sim
    sim = Simulation((6L,4L),(U,0),L;ν=U*L/Re,body,T,mem)

    # get start time
    duration=20.0; step=0.1

    # run
    @time @gif for tᵢ in range(0,duration;step)
        # update until time tᵢ in the background
        while sim_time(sim) < tᵢ
            # the step we are doing
            Δt = sim.flow.Δt[end]
            # get moment
            moment = -WaterLily.pressure_moment(center+pivot, sim.flow, sim.body)[1]
            # update rotation ODE
            α = (moment + Ia * α) / 2Ia # assuming some added moment of inertia I_a
            θ += Δt * (ω + Δt * α / 2.f0)
            ω += Δt * α
            # update the body
            sim.body = update!(sim.body; θ=θ, ω=ω)
            # measure and update flow
            sim_step!(sim;remeasure=true)
        end
        @inside sim.flow.σ[I] = WaterLily.curl(3, I, sim.flow.u) * sim.L / sim.U
        @inside sim.flow.σ[I] = ifelse(abs(sim.flow.σ[I])<0.001,0.0,sim.flow.σ[I])
        flood(sim.flow.σ|>Array,shift=(-1.5,-1.5),clims=(-5,5), axis=([], false),
              cfill=:seismic,legend=false,border=:none,size=(6sim.L,4sim.L)); body_plot!(sim)
        plot!(title="tU/L $tᵢ, θ=$(round(rad2deg(θ), digits=1))°")
        println("tU/L=", round(tᵢ, digits=4), ", θ=", round(rad2deg(θ), digits=1), "°, M=",
                round(moment, digits=3), ", ω=", round(ω, digits=3))
    end
    return sim
end

# run the main function
sim = main(64; mem=Array);