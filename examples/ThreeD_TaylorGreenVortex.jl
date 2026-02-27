using WaterLily,GLMakie

function TGV(L; Re=1600, U=1, T=Float32, mem=Array)
    # wavenumber, velocity
    κ, U = T(π/L), T(U)
    # Taylor-Green-Vortex initial velocity field
    function uλ(i,xyz)
        x,y,z = @. xyz*κ                       # scaled coordinates
        i==1 && return -U*sin(x)*cos(y)*cos(z) # u_x
        i==2 && return  U*cos(x)*sin(y)*cos(z) # u_y
        return zero(U)                         # u_z
    end
    # Initialize simulation
    return Simulation((L, L, L), (0, 0, 0), L; U, uλ, ν = U*L/Re, T, mem)
end

function λ₂!(arr, sim)                          # compute log10(-λ₂)
    a = sim.flow.σ
    @inside a[I] = log10(max(1e-6,-WaterLily.λ₂(I,sim.flow.u)*sim.L/sim.U))
    copyto!(arr ,a[inside(a)])                  # copy to CPU
end

# Initialize CUDA simulation
# using CUDA
sim = TGV(2^6; T=Float32)#, mem=CuArray);
t₀ = sim_time(sim)
duration = 15.0
step = 0.05

viz!(sim;f=λ₂!,duration,step,algorithm=:absorption,colormap=:Reds)