

Y = LinRange(0,Ly, Ny)
Z = LinRange(0, Lz, Nz)

V_b, N = Eddy_properties(Y, Z, SIGMA) #V_b: computation box, N: number of Eddies

A = cholesky_decomposition(Re_stress)

SEM_Eddy = SEM_EDDY[]
for i = 1:1:N
    push!(SEM_Eddy, SEM_EDDY(0, 0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0 ))

end