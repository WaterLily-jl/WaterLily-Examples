using WaterLily,GLMakie

import WaterLily: size_u,@loop,slice
# ghost-cell correction: enforce lid tanjential velocity U=1 at top face and 0 on other faces via interpolation
function lid_BC!(u)
    N,_ = size_u(u)
    @loop u[I,1] = 2-u[I-δ(2,I),1] over I ∈ slice(N,N[2],2)
    @loop u[I,1] = -u[I+δ(2,I),1] over I ∈ slice(N,1,2)
    @loop u[I,2] = -u[I-δ(1,I),2] over I ∈ slice(N,N[1],1)
    @loop u[I,2] = -u[I+δ(1,I),2] over I ∈ slice(N,1,1)
end

import WaterLily: mom_step!,mom_predict!,mom_correct!,mom_project!,scale_u!,CFL,AbstractFlow,quick
# overwrite mom_step! to inject lid BC after each BC! call
@fastmath function mom_step!(a::AbstractFlow,b::AbstractPoisson;λ=quick,udf=nothing,kwargs...)
    a.u⁰ .= a.u; scale_u!(a,0); t₁ = sum(a.Δt); t₀ = t₁-a.Δt[end]
    # predictor u → u'
    mom_predict!(a,t₀,t₁;λ,udf,kwargs...); lid_BC!(a.u)
    mom_project!(a,b,1,t₁); lid_BC!(a.u)
    # corrector u → u¹
    mom_correct!(a,t₁;λ,udf,kwargs...); lid_BC!(a.u)
    mom_project!(a,b,0.5,t₁); lid_BC!(a.u)
    push!(a.Δt,CFL(a))
end

function Lid_cavity(;L=2^6,U=1,Re=100,T=Float32,mem=Array)
    # make the simulations
    Simulation((L,L),(0,0),L;U=U,ν=U*L/Re,T,mem)
end

# using CUDA
sim = Lid_cavity(L=2^6)#,mem=CuArray)

function umag(arr, sim)
    a = sim.flow.σ
    @inside a[I] = √WaterLily.ke(I,sim.flow.u)
    copyto!(arr ,a[inside(a)])
end

viz!(sim; f=umag, duration=12, step=0.1, clims=(0,1), levels=20) # add: video="lid.mp4" to store the video
