
include("SEM.jl")
using DataFrames, XLSX
using Gridap

σ = 0.1 #eddy dimensions, the same in all the directions
b = 5.0
a = 0.0
y = a:0.1:b
z = y
Vboxinfo = Virtual_Box(y,z,σ)

N = Vboxinfo.N #you can override it 
t = 0
dt = 0.001

U₀ = 1.0 #Convective Velocity


TI = 0.1 #turbulence intensity

#Isotropic turbulence
u_p = (U₀ * TI)^2

Re_stress = [u_p 0.0 0.0; 
            0.0 u_p 0.0;
            0.0 0.0 u_p]



A = cholesky_decomposition(Re_stress)
t_old = 1.5
cache = (Vboxinfo, dt, U₀, σ,  A, t_old, nothing)


mutable struct SEM_cache
    Vboxinfo::Virtual_Box
    dt::Float64
    U₀::Float64
    σ::Float64
    A::Matrix{Float64}
    t_old::Float64
    Eddies::Union{Vector{SEM_EDDY},Nothing}
end

S_cache = SEM_cache(Vboxinfo, dt, U₀, σ, A, t_old, nothing)



function SEM_u_!(x,t, SEM_cache)


N =  SEM_cache.Vboxinfo.N

if SEM_cache.Eddies === nothing
    SEM_cache.Eddies = initialize_eddies(N, SEM_cache.σ, SEM_cache.Vboxinfo)
end


if t>SEM_cache.t_old
    #convecting eddies
    SEM_cache.Eddies = map(xmap -> convect_eddy(SEM_cache.dt, xmap, SEM_cache.U₀, SEM_cache.σ, SEM_cache.Vboxinfo), SEM_cache.Eddies)
    println("Convecting eddies")

    #Update cache: time and Eddies
    SEM_cache.t_old = t
elseif t<SEM_cache.t_old
    error("cannot compute eddies at a previous time step")
else
    println("Eddy at the same time")
end




u_fluct = compute_uᵢₚ([x[1], x[2], x[3]], SEM_cache.Eddies, SEM_cache.U₀, SEM_cache.Vboxinfo)

U =  compute_U(u_fluct, SEM_cache.A, SEM_cache.U₀)

return U

end

xp = VectorValue(0.0,2.5,3.0)
tp = 1.8
U =  SEM_u_!(xp,tp , S_cache)

S_cache.t_old
