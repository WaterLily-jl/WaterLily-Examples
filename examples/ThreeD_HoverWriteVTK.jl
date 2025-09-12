using WaterLily,StaticArrays,WriteVTK

function hover(L=2^6;Re=250,U=1,amp=π/4,ϵ=0.5,thk=2ϵ+√3,mem=Array)
    # z-periodic plate SDF
    function sdf(x,t)
        y = x .- SA[0,clamp(x[2],-L/2,L/2),x[3]]
        √sum(abs2,y)-thk/2
    end
    # Oscillating motion and rotation
    function map(x,t)
        α = amp*cos(t*U/L); R = SA[cos(α) sin(α) 0; -sin(α) cos(α) 0; 0 0 1]
        R * (x - SA[3L-L*sin(t*U/L),4L,0])
    end
    Simulation((6L,6L,L),(0,0,0),L;U,ν=U*L/Re,body=AutoBody(sdf,map),ϵ,perdir=(3,),mem)
end

# using CUDA
# make the sim
sim = hover(64)#;mem=CuArray)

# we know have to define these functions, for example:
vtk_velocity(a::AbstractSimulation) = a.flow.u |> Array
vtk_pressure(a::AbstractSimulation) = a.flow.p |> Array
vtk_vorticity(a::AbstractSimulation) = (@inside a.flow.σ[I] = WaterLily.curl(3,I,a.flow.u)*a.L/a.U; a.flow.σ |> Array)
vtk_body(a::AbstractSimulation) = (measure_sdf!(a.flow.σ, a.body, WaterLily.time(a.flow)); a.flow.σ |> Array)
function vtk_laplacian(a::AbstractSimulation)
    L = copy(a.flow.μ₁); N = length(size(a.flow.μ₁))
    WaterLily.@loop L[I,:,:] .= WaterLily.S(I,a.flow.u) over I in WaterLily.inside_u(a.flow.u)
    return permutedims(L,[N,1:N-1...]) |> Array # permute dims once, the writer will do it another time
end

# we might want to write some custom stuff, for this we have to create a `Dictionary`  of custom attributes
# with the keys being the names of the field in the vtk file and the values being the functions that return the data to an `Array`.
custom_write_attributes = Dict("u" => vtk_velocity,
                               "p" => vtk_pressure,
                               "ω" => vtk_vorticity,
                               "∇²u" => vtk_laplacian,
                               "d" => vtk_body)

# now we prepare the vtk writer
wr = vtkWriter("ThreeD_hover"; attrib=custom_write_attributes)

# run the simulation for 10 convective cycles and write every 0.1 convective cycle
@time for tᵢ in range(0.,10;step=0.1)
    println("tU/L=",round(tᵢ,digits=4))
    sim_step!(sim,tᵢ;remeasure=true)
    save!(wr, sim)
end
close(wr)