using WaterLily,GLMakie

function channel(;L=2^6,U=1.0,Re=200,T=Float32,mem=Array)
    # pressure gradient required to drive the flow to u~1
    g(i,x,t) = convert(typeof(t), i == 1 ? U^2/(L/2)^2 : 0)
    # periodic in x with a no-slip floor at rest (uBC=0 on the wall), starting from u=U
    Simulation((L,L),(0,0),L;U,u0=(U,0),ν=U*L/Re,g,perdir=(1,),walls=(-2,),T,mem)
end

# using CUDA
sim = channel(L=2^7)#;mem=CuArray)

function umag(arr, sim)
    a = sim.flow.σ
    @inside a[I] = √WaterLily.ke(I,sim.flow.u)
    copyto!(arr ,a[inside(a)])
end

viz!(sim; f=umag, duration=12, step=0.1, clims=(0,1), levels=20) # add: video="channel.mp4" to store the video