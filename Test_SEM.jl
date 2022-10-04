
include("SEM.jl")

σ = 0.1 #eddy dimensions, the same in all the directions
b = 5.0
a = 0.0

#defining the domamain
x = -σ:0.1:+σ 
y = a:0.1:b
z = a:0.1:b


Vboxinfo = Virtual_Box(y,z,σ)

N = Vboxinfo.N #you can override it 
t = 0
dt = 0.001

U₀ = 1.0 #Convective Velocity


TI = 0.5 #turbulence intensity

#Isotropic turbulence
u_p = (U₀ * TI)^2

Re_stress = [u_p 0.0 0.0; 
            0.0 u_p 0.0;
            0.0 0.0 u_p]



A = cholesky_decomposition(Re_stress)

Eddies = initialize_eddies(N, σ, Vboxinfo)


vector_points = [[0.0, b/2, b/2]]

Nt = 20000
q = zeros(Nt, 3)

for i = 1:1:Nt
    q[i,:] = compute_uᵢₚ(vector_points, dt, Eddies, U₀, Vboxinfo)[1]

end

U, Ek =  compute_U_k(q, A, U₀)


Statistics.std(U)



Statistics.std(Ek)





## Plotting 3D iso curves (good for visualizing the distribution and evolution of the eddies)

X, Y, Z = mgrid(x, y, z)
vector_points = create_vector_points(x, y, z)

value = compute_uᵢₚ(vector_points, dt, Eddies, U₀, Vboxinfo)[1]
Eddies

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



k = 0.1:1000
E = (k).^(-5/3) .*100 #multiplied by 100 for shifting the curve in the top part



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
    PSD_rand, freqs_rand = fft_from_signal(3/2 .* rand_signal.^2 ,dt)
    PSD_rand_tot = PSD_rand_tot .+ 1/N_rand .*PSD_rand
end



plotlyjs()

gr()
Plots.plot(xaxis=:log, yaxis=:log, xlim = [0.5, 1e3], ylims =[1e-7, 1e2], xlabel="k", ylabel="E(k)", legend=:bottomright)
Plots.plot!(freqs, PSD, label = "SEM")
Plots.plot!(freqs_rand, PSD_rand_tot, label = "RAND")
Plots.plot!(k, E, linestyle=:dash, label = "E(k)∝k^-5/3")

#Plots.plot!(freqs_mean, PSD_mean, label = "SEM mean")
#Plots.savefig("SEM_vs_RAND.pdf")