using WaterLily,GLMakie

function Lid_cavity(;L=2^6,U=1,Re=100,T=Float32,mem=Array)
    # no-slip walls on all faces, the lid (top wall) moves with uBC=U
    lid(i,x,t) = convert(eltype(x), i==1 && x[2]≥L ? U : 0)
    Simulation((L,L),lid,L;U,ν=U*L/Re,walls=(-1,1,-2,2),T,mem)
end

# using CUDA
sim = Lid_cavity(L=2^6)#,mem=CuArray)

function umag(arr, sim)
    a = sim.flow.σ
    @inside a[I] = √WaterLily.ke(I,sim.flow.u)
    copyto!(arr ,a[inside(a)])
end

viz!(sim; f=umag, duration=12, step=0.1, clims=(0,1), levels=20) # add: video="lid.mp4" to store the video
