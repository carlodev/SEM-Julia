
include("SEM.jl")

σ = 0.1 #eddy dimensions, the same in all the directions
b = 5.0
a = 0.0

x = -σ:0.1:+σ
y = a:0.1:b
z = a:0.1:b


Vboxinfo = Virtual_Box(y,z,σ)

N = Vboxinfo.N #you can override it 
Eddies = initialize_eddies(N, σ, Vboxinfo)



t = 0
dt = 0.001

U₀ = 1.0


TI = 0.6 #turbulence intensity


#Isotropic turbulence
u_p = (U₀ * TI)^2

Re_stress = [u_p 0.0 0.0; 
            0.0 u_p 0.0;
            0.0 0.0 u_p]



A = cholesky_decomposition(Re_stress)









vector_points = [[0.0, b/2, b/2]]



Nt = 10000
q = zeros(Nt, 3)

for i = 1:1:Nt
    q[i,:] = compute_uᵢₚ(vector_points, dt, Eddies, U₀, Vboxinfo)[1]

end

U, Ek =  compute_U_k(q, A, U₀)

Plots.plot(1:Nt, Ek)

Statistics.std(U[:, 2])



X, Y, Z = mgrid(x, y, z)

vector_points = create_vector_points(x, y, z)
value = compute_uᵢₚ(vector_points, dt, Eddies, U₀, Vboxinfo)[1]

plotlyjs()
iso_surfaces = isosurface(
    x=X[:],
    y=Y[:],
    z=Z[:],
    value=value[:,1],
    isomin=0.1,
    isomax=1,
    surface_count=3,
    opacity=0.5,
    caps=attr(x_show=false, y_show=false)
)

layout=Layout(yaxis=attr(scaleanchor="x", scaleratio=1), zaxis=attr(scaleanchor="x", scaleratio=1))
io = PlotlyJS.plot(iso_surfaces, Layout(yaxis=attr(scaleanchor="x", scaleratio=1)))

 


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



PSD_mean = 0.0
freqs_mean = 0.0


k = 0.1:100
E = (k).^(-5/3)



N_restart = 20
freqs = 0.0   
Nt = 2000
PSD = 0.0
vector_points = [[0.0, b/2, b/2]]
for i=1:1:N_restart
    q = zeros(Nt, 3)
for i = 1:1:Nt
    q[i,:] = compute_uᵢₚ(vector_points, dt, Eddies, U₀, Vboxinfo)[1]
end


    U, Ek =  compute_U_k(q, A, U₀)
    PSD_tmp, freqs = fft_from_signal(Ek, dt)
    PSD = PSD .+ PSD_tmp./N_restart

end

N_rand = 1000
PSD_rand_tot = 0.0
freqs_rand = 0.0
for i = 1:1:N_rand
    rand_signal = randn(3000)
    PSD_rand, freqs_rand = fft_from_signal(rand_signal.^2 ,dt)
    PSD_rand_tot = PSD_rand_tot .+ 1/N_rand .*PSD_rand
end




Plots.plot(xaxis=:log, yaxis=:log, xlim = [0.5, 100000], ylims =[1e-7, 1000], xlabel="k", ylabel="E(k)")
Plots.plot!(freqs, PSD, label = "SEM")

#Plots.plot!(freqs_mean, PSD_mean, label = "SEM mean")
Plots.plot!(freqs_rand, PSD_rand_tot, label = "RAND")
Plots.plot!(k, E, linestyle=:dash, label = "E(k)∝k^-5/3")

Plots.savefig("SEM_vs_RAND.pdf")