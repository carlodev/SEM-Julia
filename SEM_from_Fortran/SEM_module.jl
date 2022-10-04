"""
   VARIABLES : dt    : Time step                                              #
               N     : The number of eddies                                   #
               SIGMA : Eddy length scale                                      #
               V_b   : Volume of box including eddies                         #
               Nt    : The number of iterations                               #
                                                                              #
               Y,Z      : Y,Z coordinates                                     #
               U,V,W    : Mean velocity arrays                                #
               T        : Mean temperature arrays                             #
               RS       : Reynolds stress                                     #
               SEM_EDDY : Each eddies properties including positions,         #
                          intensities, length scales.                         #
               U,V,W,T_INLET : Stochastic components of inflow surface        #
               U,V,W,T_COMB  : Reconstructed components of inflow              #
                                                                              #
               U_c    : Local convection velocities                           #
               U_pr   : Mean profiles (U,V,W,T)                               #
               rms_pr : Reynolds stress profiles (uu,vv,ww,tt,uv,ut,vt,wt)    #
                       
"""

"Eddy Structure"
mutable struct SEM_EDDY
    eddy_num::Int64     # Eddy specification number
    eddy_len::Float64  # Eddy length scale
    X_pos::Float64      # Eddy's X position
    Y_pos::Float64      # Eddy's Y position
    Z_pos::Float64      # Eddy's Z position
    X_int::Float64      # Eddy's X intensity
    Y_int::Float64      # Eddy's Y intensity
    Z_int::Float64      # Eddy's Z intensity
end

"ϵ: intensity, it can be +1 or -1"
function intensity_det(x_int)
intensity = 0

if ( x_int > 0 ) 
     intensity = 1
else
    intensity = -1
end
    return intensity
end

"From the Reynolds Stress 3x3 matrix return the Cholesky decomposition"
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

"Compute the volume of the 'generation box' and the number of Eddy"
function Eddy_properties(Y, Z, SIGMA)
    S = ( Y[end] - Y[1] + 2*SIGMA ) * ( Z[end] - Z[1] + 2*SIGMA )
    V_b = S*2*SIGMA

    N = Int(round( ( Y[end] - Y[1] ) * ( Z[end] - Z[1]  ) / SIGMA^2))

    return V_b, N
end
