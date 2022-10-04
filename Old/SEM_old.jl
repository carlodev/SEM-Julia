using Plots
using Trapz
using PlotlyJS, CSV, HTTP, DataFrames
using LinearAlgebra
using PlotlyBase
using BenchmarkTools
using Statistics

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






"Compute a random position, between a and b, in one dimension, for an eddy"
function random_position(b::Float64, a::Float64)
    (rand() .- 0.5) .* (b - a) .+ (b + a) ./ 2
end



mutable struct SEM_EDDY
    eddy_num::Int64     # Eddy specification number
    σ::Float64  # Eddy length scale
    X_pos::Float64      # Eddy's X position
    Y_pos::Float64      # Eddy's Y position
    Z_pos::Float64      # Eddy's Z position
    X_int::Float64      # Eddy's X intensity
    Y_int::Float64      # Eddy's Y intensity
    Z_int::Float64      # Eddy's Z intensity
end









function uᵢ(vec_points::Vector{Vector{Float64}}, ϵᵢ::Float64, xᵢ::Vector{Float64}, σ::Float64)
    map(x -> ϵᵢ .* fσ((x .- xᵢ)./σ ), vec_points)
end







function xp(dt::Float64, xᵢ::Vector{Float64}, U₀::Float64, σ::Float64, b::Float64, a::Float64)
    x_tmp = xᵢ[1] + dt * U₀
    if x_tmp < σ
        return [x_tmp, xᵢ[2], xᵢ[3]]
    else
        xx = random_position(σ, -σ)
        yy = random_position(b, a)
        zz = random_position(b, a)
        return [xx, yy, zz]
    end

end


function eddy_number(b::Float64, a::Float64, σ::Float64)
    Sₚ = (b - a) * (b - a)
    Sₛ = σ * σ
    N = Int(round(Sₚ / Sₛ))
    Vb = 2*σ * (b - a) * (b - a)
    return N, Vb
end


function initialize_eddy(b::Float64, a::Float64, σ::Float64, N::Int)
    ϵᵢ = Float64[]
    xᵢ₀ = Vector{Float64}[]

    for i = 1:1:N
        push!(ϵᵢ, rand((-1, 1)))

        xx = random_position(σ, -σ)
        yy = random_position(b, a)
        zz = random_position(b, a)

        push!(xᵢ₀, [xx, yy, zz])

    end

    return ϵᵢ, xᵢ₀
end



function compute_uᵢₚ(x::Vector{Vector{Float64}}, dt::Float64, xᵢ::Vector{Vector{Float64}}, ϵᵢ::Vector{Float64}, U₀::Float64, σ, N::Int, b::Float64, a::Float64)
    contribution = zeros(length(x))
    for j = 1:1:N
        xᵢ[j] = xp(dt, xᵢ[j], U₀, σ, b, a)
        contribution .+= uᵢ(x, ϵᵢ[j], xᵢ[j], σ)
    end

    return 1 ./ (N)^0.5 .* contribution, xᵢ
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


function compute_U_k(u1::Vector{Float64}, u2::Vector{Float64}, u3::Vector{Float64}, A::Matrix{Float64}, U₀::Float64)
    U = Re_stress * [u1, u2, u3]
    U[1] = U[1] .+ U₀
    k = 0.5.* (u1.^2 .+ u2.^2 .+ u3.^2)
    return U, k    
end




