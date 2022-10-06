using Plots
using Trapz
using PlotlyJS, CSV, HTTP, DataFrames
using LinearAlgebra
using PlotlyBase
using BenchmarkTools
using Statistics
using FFTW #For fast fourier transformation

"Inital and inlet turbulent conditions, from paper 10.1016/j.ijheatfluidflow.2006.02.006"


"Cholesky Decomposition of Reynolds Stress Tensor"
function cholesky_decomposition(R::Matrix{Float64})
    if length(R) == 4
        a11 = (R[1, 1])^0.5
        a21 = (R[2, 1]) / a11
        a22 = (R[2, 2] - a21)^0.5

        A = [a11 0.0
            a21 a22]

    elseif length(R) == 9
        a11 = (R[1, 1])^0.5
        a21 = (R[2, 1]) / a11
        a22 = (R[2, 2] - a21)^0.5
        a31 = (R[3, 1]) / a11
        a32 = (R[3, 2] - a21 * a31) / a22
        a33 = (R[3, 3] - a31^2 - a32^2)^0.5
        A = [a11 0.0 0.0
            a21 a22 0.0
            a31 a32 a33]
    end

    return A
end




"""
function gauss_fun(x, σ, μ)
    1 / ((2π) * σ) * exp(-0.5 * ((x - μ)^2 / (σ^2)))
end


function fσ_old(x, σ)
    x = norm(x)
    if x > σ || x < -σ
        return 0
    else
        return exp(-x^2 / (σ^2 - (x^2)))

    end

end


function fσ(x, σ::Float64, C::Float64)
    x = norm(x)
    if x > σ || x < -σ
        return 0
    else
        return C * exp(-9 .* x^2 / (2 * σ^2))

    end
end



"Compute che normalization constant of the function"
function I(σ::Float64)
    xt = -σ:0.001:σ
    Cm = trapz(xt, fσ.(xt, σ, 1.0) .^ 2)
    return (1 / Cm)^0.5
end


"""

"Shape function"
function fσ(x)
   
    if abs(x[1]) <= 1 && abs(x[2]) <= 1 && abs(x[3]) <= 1

        return sqrt(1.5) * (1 - abs(x[1])) *
        sqrt(1.5) * (1 - abs(x[2])) *
        sqrt(1.5) * (1 - abs(x[3]))
    else
        return 0
        

    end
end









"Structure with the property of a single eddy"
mutable struct SEM_EDDY
    eddy_num::Int64     # Eddy specification number
    σ::Float64  # Eddy length scale
    xᵢ::Vector{Float64} # Eddy's position in the computational box [x,y,z]
    ϵᵢ::Vector{Float64}  # Eddy's intensity (+1 or -1) in [x,y,z]

end


"Volume box where the eddies are created"
struct Virtual_Box
    Y::Vector{Float64}
    Z::Vector{Float64}
    σ::Float64

    N::Int64
    V_b::Float64
    Y_start::Float64
    Y_end::Float64
    Z_start::Float64
    Z_end::Float64

    function Virtual_Box(Y, Z, σ)
        Y_start = Y[1] - σ
        Y_end = Y[end] + σ
        Z_start = Z[1] - σ
        Z_end = Z[end] + σ
  

        Sₚ = (Y_end - Y_start) * (Z_end - Z_start)
        Sₛ = σ * σ
        N = Int(round(Sₚ / Sₛ))
        V_b = 2*σ * (Y_end - Y_start) * (Z_end - Z_start)

        new(Y, Z, σ, N, V_b, Y_start, Y_end, Z_start, Z_end)
    end

end




"Initialize Eddy position and intensity"
function initialize_eddies(N::Int64, σ::Float64, Vbinfo::Virtual_Box)
    SEM_Eddy = SEM_EDDY[]
    for i =1:1:N
        ϵᵢ = rand((-1,1), 3)
        xᵢ = new_rand_position(Vbinfo)
        push!(SEM_Eddy, SEM_EDDY(i, σ,  xᵢ,  ϵᵢ))
    end
    return SEM_Eddy
end


"Random position of an eddy inside the Virtual Box"
function new_rand_position(Vbinfo::Virtual_Box)

    xx = (rand() .- 0.5) .* 2 .*Vbinfo.σ
    yy = (rand() .- 0.5) .* (Vbinfo.Y_end - Vbinfo.Y_start) .+ (Vbinfo.Y_end + Vbinfo.Y_start) ./ 2
    zz =(rand() .- 0.5) .* (Vbinfo.Z_end - Vbinfo.Z_start) .+ (Vbinfo.Z_end + Vbinfo.Z_start) ./ 2

    return[xx,yy,zz]
end




function uᵢ(vec_points::Vector{Vector{Float64}}, ϵᵢ::Float64, xᵢ::Vector{Float64}, σ::Float64)
    map(x -> ϵᵢ .* fσ((x .- xᵢ)./σ ), vec_points)
end


function uᵢ(vec_points::Vector{Float64}, ϵᵢ::Float64, xᵢ::Vector{Float64}, σ::Float64)
    ϵᵢ .* fσ((vec_points .- xᵢ)./σ )
end




"Compute the new position of all the Eddies. We consider only the convective velocity along x axis. If outside the Virtual Box, a new eddy is randomly generated inside the Virtual Box"
function convect_eddy(dt, Eddy, U₀, σ, Vbinfo)
    x_tmp = Eddy.xᵢ[1] + dt * U₀
    if x_tmp < σ
        Eddy.xᵢ = [x_tmp, Eddy.xᵢ[2], Eddy.xᵢ[3]]
    else
        Eddy.xᵢ = new_rand_position(Vbinfo)
        Eddy.ϵᵢ = rand((-1,1), 3)
    end
    return Eddy
end

"The velocity in the 3 directions is computed in each point provided in x"
function compute_uᵢₚ(x::Vector{Vector{Float64}}, dt::Float64, Eddies::Vector{SEM_EDDY}, U₀::Float64, Vbinfo::Virtual_Box)
    contribution = zeros(length(x),3)
    for j = 1:1:length(Eddies)
        Eddies[j] = convect_eddy(dt, Eddies[j], U₀, σ,Vbinfo)
        contribution[:,1] .+= uᵢ(x, Eddies[j].ϵᵢ[1], Eddies[j].xᵢ, σ)
        contribution[:,2] .+= uᵢ(x, Eddies[j].ϵᵢ[2], Eddies[j].xᵢ, σ)
        contribution[:,3] .+= uᵢ(x, Eddies[j].ϵᵢ[3], Eddies[j].xᵢ, σ)

    end

    return sqrt(Vbinfo.V_b/(Vbinfo.σ^3)) ./ (Vbinfo.N)^0.5 .* contribution, Eddies
end

"Non-convected version of compute_ui"
function compute_uᵢₚ(x::Vector{Float64}, Eddies::Vector{SEM_EDDY}, U₀::Float64, Vbinfo::Virtual_Box)
    contribution = zeros(length(x))
    for j = 1:1:length(Eddies)
        contribution[1] += uᵢ(x, Eddies[j].ϵᵢ[1], Eddies[j].xᵢ, σ)
        contribution[2] += uᵢ(x, Eddies[j].ϵᵢ[2], Eddies[j].xᵢ, σ)
        contribution[3] += uᵢ(x, Eddies[j].ϵᵢ[3], Eddies[j].xᵢ, σ)

    end

    return sqrt(Vbinfo.V_b/(Vbinfo.σ^3)) ./ (Vbinfo.N)^0.5 .* contribution
end


function create_vector_points(x, y, z)
    vector_points = Vector{Float64}[]
    for i = 1:1:length(x)
        for j = 1:1:length(y)
            for k = 1:1:length(z)
                push!(vector_points, [x[i], y[j], z[k]])
            end
        end
        
    end
    return vector_points
end

"Compute the acutual velocity and the turbulent kinetic energy. The convective velocity is just in the x direction"
function compute_U_k(q::Matrix{Float64}, A::Matrix{Float64}, U₀::Float64)
    U = A * q'
    U = U'
    U[1] = U[1] .+ U₀
    k = zeros(size(U)[1])
    for i = 1:1:size(U)[1]
        k[i] = 0.5.* (U[i,1].^2 .+ U[i,2].^2 .+ U[i,3].^2)

    end
    return U, k    
end

"Non*vectorized version"
function compute_U(q::Vector{Float64}, A::Matrix{Float64}, U₀::Float64)
    U = A * q
    U = U
    U[1] = U[1] .+ U₀
    
    return U  
end



#Spectral Tools


function fft_from_signal(q,dt)
    nt=length(q)
    fhat=fft(q)
    
    PSD = fhat.*conj(fhat)/(nt)
    PSD = real(fftshift(PSD))
    freqs = fftshift(fftfreq(nt,1/dt))
    idx = findall(x -> x>0, freqs)

return PSD[idx], freqs[idx]
end