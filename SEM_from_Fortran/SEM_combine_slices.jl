function comb_SLICE(Ny, Nz, U_inlet, V_inlet, W_inlet, A)
  U_comb = zeros(Ny,Nz)
  V_comb = zeros(Ny,Nz)
  W_comb = zeros(Ny,Nz)

for k = 1:1:Nz
  for j = 1:1:Ny


    u_tmp  = [U_inlet[j,k],V_inlet[j,k], W_inlet[j,k]]

    """
    R_loc(1,1:4)  = (/RS(1,j,k),RS(4,j,k),RS(5,j,k),THS(2,j,k)/)
    R_loc(2,1:4)  = (/RS(4,j,k),RS(2,j,k),RS(6,j,k),THS(3,j,k)/)
    R_loc(3,1:4)  = (/RS(5,j,k),RS(6,j,k),RS(3,j,k),THS(4,j,k)/)
    R_loc(4,1:4)  = (/THS(2,j,k),THS(3,j,k),THS(4,j,k),THS(1,j,k)/)
    """

    u_fluc = A * u_tmp

    U_comb[j,k] = u_fluc[1]
    V_comb[j,k] = u_fluc[2]
    W_comb[j,k] = u_fluc[3]

  end
end
   
return U_comb, V_comb, W_comb
end


function add_mean_flow(U_comb, U0)
  U_comb .+ U0

end


function compute_TKE(Ny, Nz, U_comb, V_comb, W_comb)
  TKE = zeros(Ny,Nz)

  TKE = (U_comb .^2 .+ V_comb .^2 .+ W_comb .^2) .*0.5

  return TKE
end
    