using WaterLily,StaticArrays,Plots

# helper
norm(x) = √sum(abs2,x)

function make_sim(;L=2^5,U=1,Re=250,mem=Array)

    # triangle sdf
    function triangle(p,t)
        r,k = L/2,sqrt(3.0)
        x = abs(p[1]) - r
        y = p[2] + r/k
        p = SA[clamp(x,-2r,2),y]
        if x+k*y>0.0
            p = SA[clamp(x-k*y,-2r,2),-k*x-y]./2.0
        end
        return -norm(p)*sign(p[2])
    end
    # map to center of domain
    map(x,t) = x-SA[2L,2L]

    # construct the body
    body = AutoBody(triangle,map)

    # make a simulation
    Simulation((8L,4L),(U,0),L;U,ν=U*L/Re,body,T=Float64,mem)
end
# using CUDA
# intialize
sim = make_sim()#mem=CuArray)
sim_gif!(sim,duration=10,clims=(-5,5),plotbody=true,dpi=300)