using WaterLily,StaticArrays,Plots

function make_sim(;L=2^5,U=1,Re=250,mem=Array)

    # plane sdf
    function plane(x,t,center,normal)
        normal = normal/√sum(abs2,normal)
        sum((x .- center).*normal)
    end
    # map to center of domain
    map(x,t) = x-SA[2L,2L]

    # square is intersection of four planes
    body = AutoBody((x,t)->plane(x,t,SA[-L/2,0],SA[-1, 0]),map) ∩
           AutoBody((x,t)->plane(x,t,SA[ 0,L/2],SA[ 0, 1]),map) ∩
           AutoBody((x,t)->plane(x,t,SA[L/2, 0],SA[ 1, 0]),map) ∩
           AutoBody((x,t)->plane(x,t,SA[0,-L/2],SA[ 0,-1]),map)

    # make a simulation
    Simulation((8L,4L),(U,0),L;U,ν=U*L/Re,body,T=Float64,mem)
end
# using CUDA
# intialize and run
sim = make_sim()#;mem=CuArray)
sim_gif!(sim,duration=10,clims=(-10,10),plotbody=true)