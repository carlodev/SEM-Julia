 #------------------------------------------------------------------------------ #
 #                                                                               #
 #   PROGRAM : SEM_convection.f90                                                #
 #                                                                               #
 #   PURPOSE : Convect each eddies by convective velocity                        #
 #                                                                               #
 #                                                             2017.03.03 K.Noh  #
 #                                                                               #
 #------------------------------------------------------------------------------ #

        function CONVECT_EDDY(N, Ny, Nz, SIGMA, dt, Y, Z,  SEM_Eddy, U_comb, V_comb, W_comb)

      
            Y_start = Y[1] - SIGMA
            Y_end   = Y[Ny] + SIGMA

            Z_start = Z[1] - SIGMA
            Z_end   = Z[Nz] + SIGMA


            U_conv = sum(U_comb)/(Ny*Nz)
            V_conv = sum(V_comb)/(Ny*Nz)
            W_conv = sum(W_comb)/(Ny*Nz)

            for it = 1:1:N
               SEM_Eddy[it].X_pos =  SEM_Eddy[it].X_pos + U_conv*dt
               SEM_Eddy[it].Y_pos =  SEM_Eddy[it].Y_pos + V_conv*dt
               SEM_Eddy[it].Z_pos =  SEM_Eddy[it].Z_pos + W_conv*dt

              if ( ( SEM_Eddy[it].X_pos-(-SIGMA))*                               
                   ( SEM_Eddy[it].X_pos-(SIGMA)) > 0 ||                       
                   ( SEM_Eddy[it].Y_pos - Y_start)*                              
                   ( SEM_Eddy[it].Y_pos - Y_end) > 0 ) 

                    SEM_Eddy[it].eddy_len = SIGMA
                    SEM_Eddy[it].X_pos = - SIGMA
                    
                    tmp = rand(5)
                    SEM_Eddy[it].X_pos = 0 - SIGMA
                    SEM_Eddy[it].Y_pos = Y_start + (Y_end-Y_start)*tmp[1]
                    SEM_Eddy[it].Z_pos = Z_start + (Z_end-Z_start)*tmp[2]

                    SEM_Eddy[it].X_int = intensity_det(tmp[3]-0.5)
                    SEM_Eddy[it].Y_int = intensity_det(tmp[4]-0.5)
                    SEM_Eddy[it].Z_int = intensity_det(tmp[5]-0.5)

              end #if

               #---------------------------------------------------------------- #
               #                 Periodic boundary conditions                    #
               #---------------------------------------------------------------- #
              if (  SEM_Eddy[it].Z_pos < Z_start ) 
                tmp_z = Z_start -  SEM_Eddy[it].Z_pos
                 SEM_Eddy[it].Z_pos = Z_end - tmp_z
              end #end if 

              if (  SEM_Eddy[it].Z_pos > Z_end ) 
                tmp_z =  SEM_Eddy[it].Z_pos - Z_end
                 SEM_Eddy[it].Z_pos = Z_start + tmp_z
              end #end if 

            end #end for 

            return SEM_Eddy
        end #function
