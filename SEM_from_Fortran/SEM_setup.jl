#------------------------------------------------------------------
#                         Constants for SEM                        
#------------------------------------------------------------------

U0 = 5.0

N  = 1000

Ny = 96
Nz = 256

Ly = 5.0
Lz = 5.0

Nt    = 1000
dt    = 5e-3
SIGMA = 0.20



TI = 0.1 #turbulence intensity

#Isotropic turbulence
u_p = U0 * TI
Re_stress = [u_p 0.0 0.0; 
            0.0 u_p 0.0;
            0.0 0.0 u_p]

#------------------------------------------------------------------#
#                         Initial Conditions                       #
#------------------------------------------------------------------#
Y = zeros(Ny)
Z = zeros(Nz)

U = zeros(Ny, Nz)
V = zeros(Ny, Nz)
W = zeros(Ny, Nz)

RS = zeros(6, Ny, Nz)
THS = zeros(4, Ny, Nz)

U_inlet = zeros(Ny, Nz)
V_inlet = zeros(Ny, Nz)
W_inlet = zeros(Ny, Nz)

U_comb = zeros(Ny, Nz)
V_comb = zeros(Ny, Nz)
W_comb = zeros(Ny, Nz)

U_c = zeros(Ny, Nz)

U_pr = zeros(0, Ny)
rms_pr= zeros(8,Ny)



