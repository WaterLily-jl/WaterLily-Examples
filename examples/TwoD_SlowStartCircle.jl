using WaterLily,Plots

function circle(m,n;a0=0.5,Re=250,U=1,mem=Array)
    # define a circle at the domain center
    R = n/8
    body = AutoBody((x,t)->√sum(abs2, x .- (m/4,n/2)) - R)
    # define time-varying velocity boundary conditions
    Ut(i,x,t::T) where T = i==1 ? convert(T,a0*t/R+(1.0+tanh(31.4*(t/R-1.0/a0)))/2.0*(1-a0*t/R)) : zero(T)
    Simulation((m,n), Ut, R; U, ν=U*R/Re, body, mem)
end
# using CUDA
sim = circle(2*196,196)#;mem=CuArray)
sim_gif!(sim,duration=20,clims=(-8,8),plotbody=true,axis=([], false),
         cfill=:seismic,legend=false,border=:none)