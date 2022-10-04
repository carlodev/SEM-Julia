"""Translation and adaptation from Fortran 90
https://github.com/blackcata/SEM
PURPOSE : To make inflow generation by using SEM (Synthetic Eddy Method)
"""

using Plots

include("SEM_module.jl")

include("SEM_setup.jl")
include("SEM_read.jl")


include("SEM_eddy_settings.jl")
include("SEM_fluctuation_generation.jl")
include("SEM_combine_slices.jl")
include("SEM_convection.jl")


TKE = Float64[]

for it = 1:1:Nt
    time = it * dt
    println(time)
    SEM_Eddy = EDDY_SETTING(N, Ny, Nz, SIGMA, Y, Z, SEM_Eddy)

    U_inlet, V_inlet, W_inlet = FLUCT_GEN(N, Ny, Nz, V_b, Y, Z, SEM_Eddy)

    U_comb, V_comb, W_comb = comb_SLICE(Ny, Nz, U_inlet, V_inlet, W_inlet, A)
    U_comb = add_mean_flow(U_comb, U0)
    TKE_tmp = compute_TKE(Ny, Nz, U_comb, V_comb, W_comb)
    push!(TKE, TKE_tmp[40,100])
       
    SEM_Eddy = CONVECT_EDDY(N, Ny, Nz, SIGMA, dt, Y, Z, SEM_Eddy, U_comb, V_comb, W_comb)

end



## Power Spectral Density Analysis
using Statistics, FFTW

function fft_from_signal(q,dt)
    nt=length(q)
    fhat=fft(q)
    
    PSD = fhat.*conj(fhat)/(nt)
    PSD = real(fftshift(PSD))
    freqs = fftshift(fftfreq(nt,1/dt))
    idx = findall(x -> x>0, freqs)

return PSD[idx], freqs[idx]
end

PSD, freqs = fft_from_signal(TKE, dt)
Plots.plot(xaxis=:log, yaxis=:log, xlim = [0.5, 10000], ylims =[1e-7, 100], xlabel="k", ylabel="E(k)")
Plots.plot!(freqs, PSD, label = "SEM")
