#include "dec.h"
#include "cosmic_rays.h"
#include "runge_kutta.h"
#include "read_parameters.h"
#include "synchrotron.h"


/***********************************************************************/
int main()
{
    int ii;
    
    read_parameters ();

    setup_initial_grid ();
    /* printf ("Choose propagation mode:\n"); */
    /* printf ("1 = convection\n"); */
    /* printf ("2 = diffusion\n"); */
    /* scanf("%d",&mode); getchar(); */

    choice_rk = 4;
    

    /* for (ii=1; ii<=3; ii++) */
	/* { */
	/*     printf("Choose order of Runge Kutta:   1, 2, 4\n"); */
	/*     scanf("%d",&choice_rk); getchar(); */
	    
	    
	/*     if (choice_rk == 1 || choice_rk == 2 || choice_rk == 4) */
	/*     { */
    /*         printf("Compute Runge-Kuttta %d-th order.\n", choice_rk); */
    /*         break; */
	/*     } */
	/*     printf("Not allowed!\n"); */
	
	/*     if (ii==3) */
	/*     { */
    /*         printf("Not possible.\n"); */
    /*         printf ("Stop.\n"); */
    /*         exit(0); */
	/*     } */
	/* } */


    
   
    for (i=1; i <= grid_size; i++)
    {

        if (cr[i][1].z / parsec / 1.e3 < z0)
            mode = mode_0;
        if ((cr[i][1].z / parsec / 1.e3 >= z0) && (z / parsec / 1.e3 < z1 ))
            mode = mode_1;
        if (cr[i][1].z / parsec / 1.e3 >= z1)
            mode = mode_2;

//        printf ("z = %g, mode = %i\n", cr[i][1].z / parsec, mode);
        

        //        mode = 1;
        
            
        if (choice_rk == 1)
        {
                
            for (j=1; j <= nu_channel; j++)
            {
                gamma_cr();
                dN_dE();
                dN_dz (cr[i][j].z, cr[i][j].N, cr[i][j].E, cr[i][j].gamma, cr[i][j].dN_dE);
                cr[i][j].alpha = (cr[i][j].gamma - 1.0) / 2.0;
                cr[i+1][j].N = cr[1][j].N * pow((1.0 - b * cr[i+1][j].z / cr[i+1][j].E / v1), (gamma_in - 2.0));
                  
            }
        }

        if (choice_rk == 2)
        {
               
            for (j=1; j <= nu_channel+1; j++)
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
            
                
            

        if (choice_rk == 4)
        {
               
            for (j=1; j <= nu_channel+1; j++)
            {
                    
                gamma_cr();
                dN_dE();
                cr[i][j].alpha = (cr[i][j].gamma - 1.0) / 2.0;

                if (mode == 1)
                    rk_4_conv (delta_z);
                else
                    rk_4 (delta_z);

                    
//                    cr[i+1][j].N = Newton_conv();

                        
                    
                    
            }

//This is the version with a Runge-Kutta also in the energy direction. But does not have improved solutions.                
/*
                
                for (j=0; j <= nu_channel+1; j++)
                {
                    dN_dE();
                    gamma_ext();
                    cr[i][j].alpha = (cr[i][j].gamma - 1.0) / 2.0;
                    rk_4_conv_1 (delta_z);
                }
                for (j=1; j <= nu_channel; j++)
                    rk_4_conv_2 (delta_z);
                for (j=0; j <= nu_channel+1; j++)
                    rk_4_conv_2 (delta_z);
                for (j=1; j <= nu_channel; j++)
                    rk_4_conv_3 (delta_z);
                for (j=0; j <= nu_channel+1; j++)
                    rk_4_conv_3 (delta_z);
                for (j=1; j <= nu_channel; j++)
                    rk_4_conv_4 (delta_z);
                for (j=0; j <= nu_channel+1; j++)
                    rk_4_conv_4 (delta_z);
*/
                    
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
