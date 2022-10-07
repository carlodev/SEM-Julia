
function FLUCT_GEN(N, Ny, Nz, V_b, Y, Z, SEM_Eddy)
    
    U_inlet = zeros(Ny, Nz)
    V_inlet = zeros(Ny, Nz)
    W_inlet = zeros(Ny, Nz)


    #$OMP PARALLEL DO private(k,j,it,x0,y0,z0,f)
    for k = 1:1:Nz
        for j = 1:1:Ny

            for it = 1:1:N
                x0 = (0 - SEM_Eddy[it].X_pos) / SEM_Eddy[it].eddy_len
                y0 = (Y[j] - SEM_Eddy[it].Y_pos) / SEM_Eddy[it].eddy_len
                z0 = (Z[k] - SEM_Eddy[it].Z_pos) / SEM_Eddy[it].eddy_len

                #------------------------------------------------------------#
                #                        Shape function                      #
                #------------------------------------------------------------#
                if (abs(x0) <= 1 && abs(y0) <= 1 && abs(z0) <= 1)
                    f = sqrt(1.5) * (1 - abs(x0)) *
                        sqrt(1.5) * (1 - abs(y0)) *
                        sqrt(1.5) * (1 - abs(z0))

                    U_inlet[j, k] = U_inlet[j, k] .+
                                    sqrt(V_b / SEM_Eddy[it].eddy_len .^ 3) .*
                                    SEM_Eddy[it].X_int .* f

                    V_inlet[j, k] = V_inlet[j, k] .+
                                    sqrt(V_b / SEM_Eddy[it].eddy_len .^ 3) .*
                                    SEM_Eddy[it].Y_int .* f

                    W_inlet[j, k] = W_inlet[j, k] .+
                                    sqrt(V_b / SEM_Eddy[it].eddy_len .^ 3) .*
                                    SEM_Eddy[it].Z_int .* f
                end

            end
        end
    end
    #OMP END PARALLEL

    U_inlet = U_inlet ./ sqrt(N)
    V_inlet = V_inlet ./ sqrt(N)
    W_inlet = W_inlet ./ sqrt(N)

    return U_inlet, V_inlet, W_inlet
end #end function




    
    