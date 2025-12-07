using WaterLily,StaticArrays,Plots,OrdinaryDiffEq

function main(D;Vᵣ=4.f0,Mᵣ=4.6,ξ=0.01,U=1,Re=500,T=Float32,mem=Array)
    # fsi parameters
    mₐ = π*(D/2.f0)^2           # added-mass coefficient circle
    m = Mᵣ*mₐ                   # mass as a fraction of the added-mass
    k = (U/(2π*Vᵣ*D))^2*(m+mₐ)  # spring stiffness from reduced velocity
    c = 2ξ*m*√(k/(m+mₐ))        # structural damping
    Tₙ = inv(2π*c/(2ξ*m))
    println("Period of oscillation TₙU/D: ", Tₙ/D)

    # the FSI ODE
    function fsi!(du, u, p, t)
        # unpack parameters and functions, here we can
        # pass a0 and Fy to the ODE solver
        y₁,v₁,a₁,Fy = u
        # unpack parameters, never change
        m,mₐ,k,c = p
        # rates, second-order ODE as first-order system, we need to fill the acceleration
        # in the solution vector to use it in the solver, the rate is d[3:4] = 0.
        du[1] = v₁
        u[3] = du[2] = (Fy - k*y₁ - c*v₁ + mₐ*a₁)/(m + mₐ)
        du[3:4] .= 0
    end

    # initial condition FSI, pull on the spring a bit
    y0,v0,a0,Fy = u₀ = [0.15f0D,0.f0,0.f0,0.f0]
    # parameters and time span
    params = (m,mₐ,k,c)
    tspan = (0,1000)

    # fsi problem and solver
    fsi = init(ODEProblem(fsi!, u₀, tspan, params), Tsit5(),
               reltol=1e-6, abstol=1e-6, save_everystep=false)

    # make a body
    body = AutoBody((x,t)->√sum(abs2,x)-D/2.f0, RigidMap(SA[2.f0D,2.f0D+y0],0.f0))

    # generate sim
    sim = Simulation((6D,4D), (U,0), D; ν=U*D/Re, body, T, mem)

    # get start time
    duration=200.0; step=0.1
    data = []

    @time @gif for tᵢ in range(0,duration;step)
        # update until time tᵢ in the background
        while sim_time(sim) < tᵢ
            # pressure force
            Fy = -WaterLily.pressure_force(sim)[2]
            # compute motion and acceleration 1DOF
            SciMLBase.set_u!(fsi, [y0,v0,a0,Fy])
            OrdinaryDiffEq.step!(fsi, sim.flow.Δt[end], true)
            y0,v0,a0 = fsi.u[1:3]
            # update the body, pass new position and velocity
            sim.body = update!(sim.body;x₀=SA[2.f0sim.L,2.f0sim.L+y0],V=SA[0.f0,v0])
            # measure body, update flow
            sim_step!(sim,remeasure=true)
            # store state
            push!(data,[sim_time(sim), y0, v0, a0, Fy])
        end
        @inside sim.flow.σ[I] = WaterLily.curl(3,I,sim.flow.u)*sim.L/sim.U
        @inside sim.flow.σ[I] = ifelse(abs(sim.flow.σ[I])<0.001,0.0,sim.flow.σ[I])
        flood(sim.flow.σ|>Array,shift=(-1.5,-1.5),clims=(-5,5), axis=([], false),
              cfill=:seismic,legend=false,border=:none,size=(10sim.L,8sim.L)); body_plot!(sim)
        plot!(title="tU/D $tᵢ, y/D=$(round(y0/sim.L, digits=2))")
        println("tU/L=",round(tᵢ,digits=4),", y/D=",round(y0/sim.L, digits=2),", Fy=",
                round(2Fy/sim.L, digits=3), ", a=", round(a0, digits=3))
    end
    return sim,data
end

# run
# using CUDA
sim,data = main(64) #;mem=CuArray);

# postprocess and plot
using FFTW
p1=plot(getindex.(data,1), getindex.(data,2)./sim.L,
        xlabel="tU/D", ylabel="y/D", label=:none, legend=:topright)
F = rfft(getindex.(data,2))
F = F ./ sum(abs.(F))
freqs = 2rfftfreq(length(getindex.(data,1)), sim.L)
p2=plot(freqs, abs.(F), xlims=(0,1), xlabel="fD/U", ylabel="PSD(y/D)",label=:none)
plot(p1,p2, layout=(2,1))
savefig("VIV_response.png")