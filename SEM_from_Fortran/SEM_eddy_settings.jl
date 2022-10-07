function EDDY_SETTING(N, Ny, Nz, SIGMA, Y, Z, SEM_Eddy)


Y_start = Y[1] - SIGMA
Y_end   = Y[Ny] + SIGMA

Z_start = Z[1] - SIGMA
Z_end   = Z[Nz] + SIGMA

Int_X = rand(N)
Int_Y = rand(N)
Int_Z = rand(N)


for it = 1:1:N
  SEM_Eddy[it].eddy_num = it
  SEM_Eddy[it].eddy_len = SIGMA

  tmp = rand(3)
  SEM_Eddy[it].X_pos = -SIGMA + 2*SIGMA.*tmp[1]
  SEM_Eddy[it].Y_pos = Y_start + (Y_end-Y_start)*tmp[2]
  SEM_Eddy[it].Z_pos = Z_start + (Z_end-Z_start)*tmp[3]

  SEM_Eddy[it].X_int = intensity_det(Int_X[it]-0.5)
  SEM_Eddy[it].Y_int = intensity_det(Int_Y[it]-0.5)
  SEM_Eddy[it].Z_int = intensity_det(Int_Z[it]-0.5)
end
    return SEM_Eddy
end
