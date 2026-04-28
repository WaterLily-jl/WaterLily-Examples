using WaterLily,ParametricBodies,Plots,StaticArrays

function TUDflame(;L=2^5,Re=1_000,U=1,mem=Array,T=Float32)
    # points for the flame
    pnts = hcat([[0.8312, 0.8299],[0.8266, 0.8317],[0.7852, 0.7831],[0.6958, 0.7351],
                 [0.5562, 0.6954],[0.4084, 0.6432],[0.2853, 0.5667],[0.2143, 0.4897],
                 [0.1693, 0.3725],[0.1716, 0.2742],[0.1852, 0.2185],[0.1970, 0.2185],
                 [0.2119, 0.3275],[0.2681, 0.4039],[0.3391, 0.4602],[0.4036, 0.5153],
                 [0.4414, 0.5656],[0.4491, 0.5639],[0.4207, 0.4881],[0.3810, 0.4276],
                 [0.3662, 0.3773],[0.3952, 0.3311],[0.4626, 0.3401],[0.5070, 0.3792],
                 [0.5247, 0.4076],[0.5324, 0.4058],[0.5223, 0.3342],[0.4951, 0.2643],
                 [0.4430, 0.2057],[0.3933, 0.1778],[0.3939, 0.1671],[0.4418, 0.1707],
                 [0.5104, 0.1944],[0.5962, 0.2543],[0.6613, 0.3248],[0.7051, 0.3899],
                 [0.7388, 0.4575],[0.7572, 0.5155],[0.7460, 0.5190],[0.6785, 0.4722],
                 [0.6105, 0.4645],[0.5733, 0.5018],[0.5905, 0.5817],[0.6538, 0.6410],
                 [0.7348, 0.6990],[0.8070, 0.7589],[0.8360, 0.8282],[0.8312, 0.8299]]...)
    # scale and move to center of domain
    points = SMatrix{size(pnts)...,T}(pnts.*2L.+SA[2L,2L])
    # make the body
    body = ParametricBody(BSplineCurve(points;degree=2))
    # make the sim and return
    Simulation((10L,6L),(U,0),L;U,ν=U*L/Re,body,T=T,mem=mem)
end

# make a sim
# using CUDA
sim = TUDflame(;L=2^7) #,mem=CuArray)

# set -up simulations time and time-step for ploting
t₀,duration,tstep = round(sim_time(sim)), 50.0, 0.05

# run
@gif for tᵢ in range(t₀,t₀+duration;step=tstep)

    # update until time tᵢ in the background
    sim_step!(sim,tᵢ,remeasure=false)

    # flood plot of the masked vorticity
    @inside sim.flow.σ[I] = WaterLily.curl(3,I,sim.flow.u)*sim.L/sim.U
    @inside sim.flow.σ[I] = ifelse(abs(sim.flow.σ[I])<0.001,0.0,sim.flow.σ[I])
    flood(sim.flow.σ,shift=(0,0),clims=(-10,10),
          cfill=:starrynight,legend=false,border=:none,size=(1200,800))
    # binary array for the body
    @inside sim.flow.σ[I] = ifelse(sdf(sim.body,loc(0,I),0.0)<=0,1.0,NaN)
    a = Array(sim.flow.σ[inside(sim.flow.p)])
    contourf!(axes(a,1).+2,axes(a,2).+1.5,a',color="#00A6D6")
    # plot curve as line
    c = [sim.body.curve(s,0.0) for s ∈ 0:1/100:1]
    plot!(getindex.(c,1).+2,getindex.(c,2).+1.5,color="#00A6D6",lw=1)
    # print time step
    println("tU/L=",round(tᵢ,digits=4),", Δt=",round(sim.flow.Δt[end],digits=3))
end