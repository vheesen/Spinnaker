#include "dec.h"
#include "cosmic_rays.h"
#include "runge_kutta.h"
#include "read_parameters.h"
#include "synchrotron.h"
#include "jet.h"

/****************************************************************************/
/* Command line argument that reads the name of the parameter file*/
int main(int argc,char* argv[])
{
    int choice_rk, rk_energy, counter;

    if(argc == 1)
        strcpy (parameter_file_name,"parameters");
    else
        strcpy (parameter_file_name,argv[1]);

    
/* Read the parameters from a file */
    read_parameters ();

/* Setup of the 2D grid in spatial and energy coordinates */    
    setup_initial_grid ();

/* Use the most accurate Runge-Kutta alogrithm of 4th order */    
    choice_rk = 4;

/* Experimental of Runge-Kutta also in energy direction */    
    rk_energy = -1;
    
    for (i=0; i <= grid_size; i++)
    {
            
        if (choice_rk == 1)
        {
                
            for (j=0; j <= nu_channel+1; j++)
            {
                
                if (j == 0)
                {
                    if ( ( (cr[i][j+1].N > 0.) && (cr[i][j].N > 0) ) && (cr[i][j+1].N  < cr[i][j].N > 0) )
                    {
                        gamma_cr();
                        dN_dE();
                        dN_dz (cr[i][j].z, cr[i][j].N, cr[i][j].E, cr[i][j].gamma, cr[i][j].dN_dE);
                        cr[i][j].alpha = (cr[i][j].gamma - 1.0) / 2.0;
                        cr[i+1][j].N = cr[0][j].N * pow((1.0 - b * cr[i+1][0].z / cr[0][j].E / v_z[i]), (gamma_in - 2.0));
                    }

                    else
                    {
                        cr[i+1][j].N = cr[i][j].N;
                        cr[i][j].integrate_true = -1;
                    }
                }

                else
                {
                    if ( ( (cr[i][j].N > 0.) && (cr[i][j-1].N > 0) ) && (cr[i][j].N  < cr[i][j-1].N ) )
                    {
                        gamma_cr();
                        dN_dE();
                        dN_dz (cr[i][j].z, cr[i][j].N, cr[i][j].E, cr[i][j].gamma, cr[i][j].dN_dE);
                        cr[i][j].alpha = (cr[i][j].gamma - 1.0) / 2.0;
                        cr[i+1][j].N = cr[0][j].N * pow((1.0 - b * cr[i+1][0].z / cr[0][j].E / v_z[i]), (gamma_in - 2.0));
                    }
                    
                    else
                    {
                        cr[i+1][j].N = cr[i][j].N;
                        cr[i][j].integrate_true = -1;
                    }
                }

                if ( cr[i+1][j].N < 0. )
                    cr[i+1][j].integrate_true = -1;
            }
        }
        

        else if (choice_rk == 2)
        {
               
            for (j=0; j <= nu_channel+1; j++)
            {
                    
                gamma_cr();
                    
                dN_dE();
                cr[i][j].alpha = (cr[i][j].gamma - 1.0) / 2.0;

//                if (mode == 1)
//                    rk_2_conv (delta_z);
//                else
                rk_2 (delta_z);
//                    if (j==nu_channel+1)
//                        cr[i+1][j].N = cr[1][j].N * pow((1.0 - b * cr[i+1][j].z / cr[i+1][j].E / v_convect), (gamma_in - 2.0));

            }
        }
            
                
            

        else if ((choice_rk == 4) && (rk_energy != 1))
        {
               
            for (j=0; j <= nu_channel+1; j++)
            {
                    
                if (j == 0)
                {
                    
                    if ( ( (cr[i][j+1].N > 0.) && (cr[i][j].N > 0) ) && (cr[i][j+1].N  < cr[i][j].N > 0) )
                    {
                        gamma_cr();
                        dN_dE();
                        cr[i][j].alpha = (cr[i][j].gamma - 1.0) / 2.0;
                    
                        if (mode == 1)
                            rk_4_conv (delta_z);
                        else if (mode == 2)
                            rk_4 (delta_z);
                        else if (mode == 3)
                            rk_4_cylindrical (delta_z);
                        else
                            rk_4_radial (delta_z);
                        
                    }
                
                    else
                        
                    {
                        cr[i+1][j].N = cr[i][j].N;
                        cr[i][j].integrate_true = -1;
                    }

                }

                else
                {
                    
                    if ( ( (cr[i][j].N > 0.) && (cr[i][j-1].N > 0) ) && (cr[i][j].N  < cr[i][j-1].N ) )
                    {

                        gamma_cr();
                        dN_dE();
                        cr[i][j].alpha = (cr[i][j].gamma - 1.0) / 2.0;
                    
                        if (mode == 1)
                            rk_4_conv (delta_z);
                        else if (mode == 2)
                            rk_4 (delta_z);
                        else if (mode == 3)
                            rk_4_cylindrical (delta_z);
                        else
                            rk_4_radial (delta_z);
                                            
                    }

                    else
                    {
                        cr[i+1][j].N = cr[i][j].N;
                        cr[i][j].integrate_true = -1;
                    }

                }
                
                if ( cr[i+1][j].N < 0. )
                    cr[i+1][j].integrate_true = -1;
                
            }

        }
        
/*This is the version with a Runge-Kutta also in the energy direction. But does not have improved solutions.                */
    
        else if ((choice_rk == 4) && (rk_energy == 1) && (mode == 1))
        {
                    
            for (j=0; j <= nu_channel+1; j++)
            {
                if (j == 0)
                {
                    if ( ( (cr[i][j+1].N > 0.) && (cr[i][j].N > 0) ) && (cr[i][j+1].N  < cr[i][j].N > 0) )
                    {
                        dN_dE();
                        gamma_cr();
                        cr[i][j].alpha = (cr[i][j].gamma - 1.0) / 2.0;
                        rk_4_conv_1 (delta_z);

                    }

                }

                else
                {
                    if ( ( (cr[i][j].N > 0.) && (cr[i][j-1].N > 0) ) && (cr[i][j].N  < cr[i][j-1].N ) )
                    {
                        dN_dE();
                        gamma_cr();
                        cr[i][j].alpha = (cr[i][j].gamma - 1.0) / 2.0;
                        rk_4_conv_1 (delta_z);
                    }
                    
                }
                
            }

            for (j=1; j <= nu_channel; j++)
            {
                 if ( ( (cr[i][j+1].N > 0.) && (cr[i][j].N > 0) ) && (cr[i][j+1].N  < cr[i][j].N > 0) )
                     rk_4_conv_2 (delta_z);
            }
            
            for (j=0; j <= nu_channel+1; j++)
            {
                if (j == 0)
                {
                    if ( ( (cr[i][j+1].N > 0.) && (cr[i][j].N > 0) ) && (cr[i][j+1].N  < cr[i][j].N > 0) )
                        rk_4_conv_2 (delta_z);
                }

                else
                {
                   if ( ( (cr[i][j].N > 0.) && (cr[i][j-1].N > 0) ) && (cr[i][j].N  < cr[i][j-1].N ) )
                       rk_4_conv_2 (delta_z);
                }
                
            }

            for (j=1; j <= nu_channel; j++)
            {
                 if ( ( (cr[i][j+1].N > 0.) && (cr[i][j].N > 0) ) && (cr[i][j+1].N  < cr[i][j].N > 0) )
                     rk_4_conv_3 (delta_z);
            }
            
            for (j=0; j <= nu_channel+1; j++)
            {
                if (j == 0)
                {
                    if ( ( (cr[i][j+1].N > 0.) && (cr[i][j].N > 0) ) && (cr[i][j+1].N  < cr[i][j].N > 0) )
                        rk_4_conv_3 (delta_z);
                }

                else
                {
                   if ( ( (cr[i][j].N > 0.) && (cr[i][j-1].N > 0) ) && (cr[i][j].N  < cr[i][j-1].N ) )
                       rk_4_conv_3 (delta_z);
                }
            
            }
            
            for (j=1; j <= nu_channel; j++)
            {
                 if ( ( (cr[i][j+1].N > 0.) && (cr[i][j].N > 0) ) && (cr[i][j+1].N  < cr[i][j].N > 0) )
                    rk_4_conv_4 (delta_z); 
            }

            for (j=0; j <= nu_channel+1; j++)
            {
                if (j == 0)
                {
                    if ( ( (cr[i][j+1].N > 0.) && (cr[i][j].N > 0) ) && (cr[i][j+1].N  < cr[i][j].N > 0) )
                        rk_4_conv_4 (delta_z);

                    else
                    {
                        cr[i+1][j].N = cr[i][j].N;
                        cr[i][j].integrate_true = -1;
                    }

                }

                else
                {
                   if ( ( (cr[i][j].N > 0.) && (cr[i][j-1].N > 0) ) && (cr[i][j].N  < cr[i][j-1].N ) )
                       rk_4_conv_4 (delta_z);

                   else
                   {
                       cr[i+1][j].N = cr[i][j].N;
                       cr[i][j].integrate_true = -1;
                   }

                }
            
            }

            for (j=0; j <= nu_channel+1; j++)
            {
                if ( cr[i+1][j].N < 0. )
                    cr[i+1][j].integrate_true = -1;
            }
            
        }

    }
    
      
    
            
/*     if (mode == 2) */
/*     { */
/*         for (i=0; i <= grid_size; i++) */
/*         { */
/*             if (choice_rk == 2) */
/*             { */
/*                 for (j=1; j <= nu_channel+1; j++) */
/*                 { */
/*                     gamma_cr(); */
/*                     dN_dE(); */
/*                     cr[i][j].alpha = (cr[i][j].gamma - 1.0) / 2.0; */
/*                     rk_2 (delta_z); */
/*                 } */

/*             } */

/*             if (choice_rk == 4) */
/*             { */
/*                 for (j=1; j <= nu_channel+1; j++) */
/*                 { */
/*                     gamma_cr(); */
/*                     dN_dE(); */
/*                     cr[i][j].alpha = (cr[i][j].gamma - 1.0) / 2.0; */
/*                     rk_4 (delta_z); */
/*                 } */
/* /\* 	 */
/*                 if (cr[i+1][j].N < 0.0) */
/*                 { */
/*                     printf("Diffusion to large. N < 0!. Stop\n"); */
/*                     printf("i = %d j = %d\n", i, j); */
/*                     printf ("gamma = %g alpha = %g\n", cr[i][j].gamma, */
/*                             cr[i][j].alpha); */
/*                     printf("N = %g y = %g\n", cr[i+1][j].N, cr[i+1][j].y); */
/*                     cr[i+1][j].N = 10.0; */
/*                 } */
/* *\/	     */

/*             } */
/*         } */
        
    
    
    
    output_file (grid_size);    
   

    return 0;
}
