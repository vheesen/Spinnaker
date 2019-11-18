#include "dec.h"
#include "cosmic_rays.h"
#include "runge_kutta.h"
#include "read_parameters.h"
#include "synchrotron.h"
#include "jet.h"
#include "read_parameters2.h"
#include <time.h>
#include "output.h"

/****************************************************************************/

/* Command line argument that reads the name of the parameter file*/

int main(int argc,char* argv[])
{
    int choice_rk, rk_energy, counter;
	time_t seconds_0;
	time_t seconds_1;
	struct timespec ts_start;
	struct timespec ts_end;
    if(argc == 1)
	{
        strcpy (parameter_file_name,"parameters");
//	clock_gettime(CLOCK_MONOTONIC, &ts_start);
	read_parameters ();
//	clock_gettime(CLOCK_MONOTONIC, &ts_end);
//	printf("Old Reading took %ld microseconds\n",(ts_end.tv_nsec-ts_start.tv_nsec)/1000);
	}
    if(argc == 2)
	{
        strcpy (parameter_file_name,argv[1]);
//	clock_gettime(CLOCK_MONOTONIC, &ts_start);
	read_parameters ();
//	clock_gettime(CLOCK_MONOTONIC, &ts_end);
//	printf("Old Reading took %ld seconds\n",(ts_end.tv_nsec-ts_start.tv_nsec)/1000);
	}

// Read the parameters from argv[2]. No need for a file 
    if(argc == 3)
	{
	if(strcmp(argv[1],"-x")==0)
		{
//			clock_gettime(CLOCK_MONOTONIC, &ts_start);
			read_parameters2(argv[2]);
		//	clock_gettime(CLOCK_MONOTONIC, &ts_end);
		//	printf("New Reading took %ld microseconds\n",(ts_end.tv_nsec-ts_start.tv_nsec)/1000);
		}
	}
    
// print variables from parameters.
int show=0;	//for debugging purposes.
if(show==1)
{
//	printf("%s\n",argv[2]);
	printf("grid_size = %d\n",grid_size);
	printf("nu_channel = %d\n",nu_channel);
	printf("grid_delta = %d\n",grid_delta);
	printf("z_halo = %lg\n",z_halo/kpc);
	printf("first_data_point_at_0kpc = %d\n",first_data_point_at_0kpc);
	printf("normalize_intensities = %d\n",normalize_intensities);
	printf("nu_1 = %lg\n",nu_1);
	printf("nu_2 = %lg\n",nu_2);
	printf("nu_3 = %lg\n",nu_3);
	printf("nu_4 = %lg\n",nu_4);
	printf("mode = %d\n",mode);
	printf("epsilon = %d\n",epsilon);
	printf("FWHM_effective_beam = %lg\n",FWHM_effective_beam);
	printf("gamma_in = %lg\n",gamma_in);
	printf("rad_field = %lg\n",rad_field);
	printf("V0 = %lg\n",V0);
	printf("Velocity_field = %d\n",velocity_field);
	printf("h_V = %lg\n",h_V);
	printf("adiabatic_losses = %d\n",adiabatic_losses);
	printf("D0 = %lg\n",D0); 
	printf("mu_diff = %lg\n",mu_diff);
	printf("galaxy_mode = %d\n",galaxy_mode);
	printf("z1 = %lg\n",z1);
	printf("B0 = %lg\n",B0);
	printf("B1 = %lg\n",B1);
	printf("h_B1 = %lg\n",h_B1);
	printf("h_B2 = %lg\n",h_B2);
	printf("model = %d\n",model);
	printf("initialize_model = %d\n",initialize_model);
	printf("model_north = %d\n",model_north);
	printf("update_model = %d\n",update_model); 
	printf("xi = %lg\n",xi); 
	printf("beta = %lg\n",beta); 
	printf("R0 = %lg\n",R0);
}
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
                        rk_4_conv_1 (delta_z);

                    }

                }

                else
                {
                    if ( ( (cr[i][j].N > 0.) && (cr[i][j-1].N > 0) ) && (cr[i][j].N  < cr[i][j-1].N ) )
                    {
                        dN_dE();
                        gamma_cr();
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
/*                     rk_2 (delta_z); */
/*                 } */

/*             } */

/*             if (choice_rk == 4) */
/*             { */
/*                 for (j=1; j <= nu_channel+1; j++) */
/*                 { */
/*                     gamma_cr(); */
/*                     dN_dE(); */
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
        
    
	if(argc == 3)
	{
		if(strcmp(argv[1],"-x")==0)
			{
				output_stdout( (int) (grid_size / 2.0) );
			}
	}
	else
	{
  	  	output_file ( (int) (grid_size / 2.0) );        
	}

    return 0;
}
