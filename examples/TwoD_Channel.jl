using WaterLily,GLMakie

import WaterLily: size_u,@loop,slice
# ghost-cell correction: enforce no-slip on bottom wall via interpolation
function channel_BC!(u)
    N,_ = size_u(u)
    @loop u[I,1] = -u[I+δ(2,I),1] over I ∈ slice(N,1,2)
end

import WaterLily: mom_step!,mom_predict!,mom_correct!,mom_project!,scale_u!,CFL,AbstractFlow
# overwrite mom_step! to inject channel BC after each BC! call
@fastmath function mom_step!(a::AbstractFlow,b::AbstractPoisson;udf=nothing,kwargs...)
    a.u⁰ .= a.u; scale_u!(a,0); t₁ = sum(a.Δt); t₀ = t₁-a.Δt[end]
    # predictor u → u'
    mom_predict!(a,t₀,t₁;udf,kwargs...); channel_BC!(a.u)
    mom_project!(a,b,1,t₁); channel_BC!(a.u)
    # corrector u → u¹
    mom_correct!(a,t₁;udf,kwargs...); channel_BC!(a.u)
    mom_project!(a,b,0.5,t₁); channel_BC!(a.u)
    push!(a.Δt,CFL(a))
end

function channel(;L=2^6,U=1.0,Re=200,T=Float32,mem=Array)
    # pressure gradient required to drive the flow to u~1
    g(i,x,t) = convert(typeof(t), i == 1 ? U^2/(L/2)^2 : 0)
    # make the simulations
    Simulation((L,L),(U,0),L;U=U,ν=U*L/Re,g,perdir=(1,),T,mem)
end

# using CUDA
sim = channel(L=2^7)#;mem=CuArray)

function umag(arr, sim)
    a = sim.flow.σ
    @inside a[I] = √WaterLily.ke(I,sim.flow.u)
    copyto!(arr ,a[inside(a)])
end

viz!(sim; f=umag, duration=12, step=0.1, clims=(0,1), levels=20) # add: video="channel.mp4" to store the video