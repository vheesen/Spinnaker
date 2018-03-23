#include "dec.h"
#include "cosmic_rays.h"
#include "runge_kutta.h"
#include "read_parameters.h"
#include "synchrotron.h"
#include "jet.h"


/***********************************************************************/
/* Runge Kutta 2nd order, see Numerical Recipies §16*/
/* Auxiliary variable y2=dy1_dz*/
/* Convert second order equation into to ordinary diff eq.*/

void rk_2 (double h)

{

    double k1_y, k2_y, k1_N, k2_N;


    k1_y = h * d2N_dz2 (cr[i][j].z, cr[i][j].N, cr[i][j].y, 
		cr[i][j].E, cr[i][j].gamma, cr[i][j].dN_dE);

    
 /*    printf ("k1_y = %g", k1_y); */
/*     exit(0); */

    k1_N = h * cr[i][j].y;

    k2_y = h * d2N_dz2 (cr[i][j].z + h / 2.0, cr[i][j].N + k1_N / 2.0, 
			cr[i][j].y + k1_y / 2.0, 
			cr[i][j].E, cr[i][j].gamma, cr[i][j].dN_dE);

    

    k2_N = h * (cr[i][j].y + k1_y / 2.0);

    cr[i+1][j].y = cr[i][j].y + k2_y;
    cr[i+1][j].N = cr[i][j].N + k2_N;

  
//    printf("rk_2: %d %g %g\n",i, x[i], rk[i+1].y1);
   
}
/***********************************************************************/
/* Runge Kutta 2nd order, see Numerical Recipies §16*/
/* Auxiliary variable y2=dy1_dz*/
/* Convert second order equation into to ordinary diff eq.*/

void rk_2_conv (double h)

{

    double k1, k2;


    k1 = h * dN_dz (cr[i][j].z, cr[i][j].N, cr[i][j].E, cr[i][j].gamma, cr[i][j].dN_dE);
    

    k2 = h * dN_dz (cr[i][j].z + h / 2.0, cr[i][j].N + k1 / 2.0, 
			cr[i][j].E, cr[i][j].gamma, cr[i][j].dN_dE);

    

    cr[i+1][j].N = cr[i][j].N + k2;

  
}

/***********************************************************************/
/* Runge Kutta 4th order, see Numerical Recipies §16*/
/*void rk_4(double h)
{
    double k1_y, k2_y, k3_y, k4_y;
    double k1_N, k2_N, k3_N, k4_N;

    k1_y = h * d2N_dz2 (cr[i][j].z, cr[i][j].N, cr[i][j].y, 
		cr[i][j].E, cr[i][j].gamma, cr[i][j].alpha,
		cr[i][j].Q);
    
 
    k1_N = h * (cr[i][j].y + k1_y);

    k2_y = h * d2N_dz2 (cr[i][j].z + h / 2.0, cr[i][j].N + k1_N / 2.0, 
			cr[i][j].y + k1_y / 2.0, 
			cr[i][j].E, cr[i][j].gamma, cr[i][j].alpha,
			cr[i][j].Q);

    k2_N = h * (cr[i][j].y + k2_y / 2.0);

   

    k3_y = h * d2N_dz2 (cr[i][j].z + h / 2.0, cr[i][j].N + k2_N / 2.0, 
			cr[i][j].y + k2_y / 2.0, 
			cr[i][j].E, cr[i][j].gamma, cr[i][j].alpha,
			cr[i][j].Q);

    k3_N = h * (cr[i][j].y + k2_y / 2.0);

    k4_y = h * d2N_dz2 (cr[i][j].z + h, cr[i][j].N + k3_N, 
			cr[i][j].y + k3_y, 
			cr[i][j].E, cr[i][j].gamma, cr[i][j].alpha,
			cr[i][j].Q);

    k4_N = h * (cr[i][j].y + k3_y);

    cr[i+1][j].y = cr[i][j].y + k1_y / 6.0 + k2_y / 3.0 + 
	k3_y / 3.0 + k4_y / 6.0;
    cr[i+1][j].N = cr[i][j].N + k1_N / 6.0 + k2_N / 3.0 +
	k3_N / 3.0 + k4_N / 6.0;

  
    }*/
/**********************************************************************/
/***********************************************************************/
/* Runge Kutta 4th order, see Numerical Recipies §16*/
void rk_4(double h)
{
    double k1_y, k2_y, k3_y, k4_y;
    double k1_N, k2_N, k3_N, k4_N;

    k1_y = h * d2N_dz2 (cr[i][j].z, cr[i][j].N, cr[i][j].y, 
		cr[i][j].E, cr[i][j].gamma, cr[i][j].dN_dE);
 
    k1_N = h * (cr[i][j].y + k1_y);

    k2_y = h * d2N_dz2 (cr[i][j].z + h / 2.0, cr[i][j].N + k1_N / 2.0, 
			cr[i][j].y + k1_y / 2.0, 
			cr[i][j].E, cr[i][j].gamma, cr[i][j].dN_dE);

    k2_N = h * (cr[i][j].y + k2_y / 2.0);

    k3_y = h * d2N_dz2 (cr[i][j].z + h / 2.0, cr[i][j].N + k2_N / 2.0, 
			cr[i][j].y + k2_y / 2.0, 
			cr[i][j].E, cr[i][j].gamma, cr[i][j].dN_dE);

    k3_N = h * (cr[i][j].y + k2_y / 2.0);

    k4_y = h * d2N_dz2 (cr[i][j].z + h, cr[i][j].N + k3_N, 
			cr[i][j].y + k3_y, 
			cr[i][j].E, cr[i][j].gamma, cr[i][j].dN_dE);

    k4_N = h * (cr[i][j].y + k3_y);

    cr[i+1][j].y = cr[i][j].y + k1_y / 6.0 + k2_y / 3.0 + 
        k3_y / 3.0 + k4_y / 6.0;
    cr[i+1][j].N = cr[i][j].N + k1_N / 6.0 + k2_N / 3.0 +
        k3_N / 3.0 + k4_N / 6.0;

  
}
/**********************************************************************/
/***********************************************************************/
void rk_4_conv(double h)
{
    double k1_N, k2_N, k3_N, k4_N;

    k1_N = h * dN_dz (cr[i][j].z, cr[i][j].N, cr[i][j].E,
                            cr[i][j].gamma, cr[i][j].dN_dE);
   
    k2_N = h * dN_dz (cr[i][j].z + h / 2.0, cr[i][j].N + k1_N / 2.0, 
                            cr[i][j].E, cr[i][j].gamma, cr[i][j].dN_dE);
   

    k3_N = h * dN_dz (cr[i][j].z + h / 2.0, cr[i][j].N + k2_N / 2.0, 
                            cr[i][j].E, cr[i][j].gamma, cr[i][j].dN_dE);


    k4_N = h * dN_dz (cr[i][j].z + h, cr[i][j].N + k3_N, 
                                cr[i][j].E, cr[i][j].gamma, cr[i][j].dN_dE);

    cr[i+1][j].N = cr[i][j].N + k1_N / 6.0 + k2_N / 3.0 + k3_N / 3.0 + k4_N / 6.0;

}
/***********************************************************************/
void rk_4_conv_1(double h)
{
    cr[i][j].k1_N = h * dN_dz (cr[i][j].z, cr[i][j].N, cr[i][j].E,
                            cr[i][j].gamma, cr[i][j].dN_dE);

 
}
/***********************************************************************/
void rk_4_conv_2(double h)
{

    double dN_dE_rk;
    

    cr[i][j].gamma_rk = - log ( (cr[i][j+1].N + cr[i][j+1].k1_N / 2.0) / (cr[i][j-1].N + cr[i][j-1].k1_N / 2.0) ) / 
        log (cr[i][j+1].E / cr[i][j-1].E);
    cr[i][j].dN_dE_rk = ( (cr[i][j+1].N + cr[i][j+1].k1_N / 2.0) - (cr[i][j-1].N + cr[i][j-1].k1_N / 2.0) ) / 
        (cr[i][j+1].E - cr[i][j-1].E);

    if (j==nu_channel+1)
    {
        cr[i][j].gamma_rk = 2.0 * cr[i][j-1].gamma_rk - cr[i][j-2].gamma_rk;
        cr[i][j].dN_dE_rk = ( (cr[i][j].N + cr[i][j].k1_N / 2.0) - (cr[i][j-1].N + cr[i][j-1].k1_N / 2.0) ) / 
        (cr[i][j].E - cr[i][j-1].E);
    }

    if (j==0)
    {
        cr[i][j].gamma_rk = 2.0 * cr[i][1].gamma_rk - cr[i][2].gamma_rk;
        cr[i][j].dN_dE_rk = ( (cr[i][j+1].N + cr[i][j+1].k1_N / 2.0) - (cr[i][j].N + cr[i][j].k1_N / 2.0) ) / 
        (cr[i][j+1].E - cr[i][j].E);
        
    }

//    if (j==0)
//        printf ("Conv2: i = %d gamma=%g\n", i, cr[i][1].gamma_rk);

    
//    if (j==nu_channel)
//        printf ("conv2: gamma_rk = %g, gamma = %g \n", cr[i][j].gamma_rk, cr[i][j].gamma);
    
        
    cr[i][j].k2_N = h * dN_dz (cr[i][j].z + h / 2.0, cr[i][j].N + cr[i][j].k1_N / 2.0, 
                                     cr[i][j].E, cr[i][j].gamma_rk,  cr[i][j].dN_dE_rk);

//    if (j==50)
//        printf ("i = %d conv2: N = %g\n", i, cr[i][j].k2_N);

    
}
/***********************************************************************/
void rk_4_conv_3(double h)
{
    double dN_dE_rk;
    
    
    cr[i][j].gamma_rk = - log ( (cr[i][j+1].N + cr[i][j+1].k2_N / 2.0) / (cr[i][j-1].N + cr[i][j-1].k2_N / 2.0) ) / 
        log (cr[i][j+1].E / cr[i][j-1].E);
    cr[i][j].dN_dE_rk = ( (cr[i][j+1].N + cr[i][j+1].k2_N / 2.0) - (cr[i][j-1].N + cr[i][j-1].k2_N / 2.0) ) / 
        (cr[i][j+1].E - cr[i][j-1].E);
    if (j==nu_channel+1)
    {
        cr[i][j].gamma_rk = 2.0 * cr[i][j-1].gamma_rk - cr[i][j-2].gamma_rk;
        cr[i][j].dN_dE_rk = ( (cr[i][j].N + cr[i][j].k2_N / 2.0) - (cr[i][j-1].N + cr[i][j-1].k2_N / 2.0) ) / 
        (cr[i][j].E - cr[i][j-1].E);
    }
    
    if (j==0)
    {
        cr[i][j].gamma_rk = 2.0 * cr[i][1].gamma_rk - cr[i][2].gamma_rk;
        cr[i][j].dN_dE_rk = ( (cr[i][j+1].N + cr[i][j+1].k2_N / 2.0) - (cr[i][j].N + cr[i][j].k2_N / 2.0) ) / 
        (cr[i][j+1].E - cr[i][j].E);
    }


    
    
//    if (j==nu_channel)
//        printf ("conv3: gamma_rk = %g, gamma = %g \n", cr[i][j].gamma_rk, cr[i][j].gamma);


    cr[i][j].k3_N = h * dN_dz (cr[i][j].z + h / 2.0, cr[i][j].N + cr[i][j].k2_N / 2.0, 
                                     cr[i][j].E, cr[i][j].gamma_rk, cr[i][j].dN_dE_rk);

//    if (j==50)
//        printf ("i = %d conv3: N = %g\n", i, cr[i][j].k3_N);

    
}
/***********************************************************************/
void rk_4_conv_4(double h)
{
    double dN_dE_rk;
        
    cr[i][j].gamma_rk = - log ( (cr[i][j+1].N + cr[i][j+1].k3_N) / (cr[i][j-1].N + cr[i][j-1].k3_N) ) / 
        log (cr[i][j+1].E / cr[i][j-1].E);
    cr[i][j].dN_dE_rk = ( (cr[i][j+1].N + cr[i][j+1].k3_N) - (cr[i][j-1].N + cr[i][j-1].k3_N) ) / 
        (cr[i][j+1].E - cr[i][j-1].E);

    if (j==nu_channel+1)
    {
        cr[i][j].gamma_rk = 2.0 * cr[i][j-1].gamma_rk - cr[i][j-2].gamma_rk;
        cr[i][j].dN_dE_rk = ( (cr[i][j].N + cr[i][j].k3_N) - (cr[i][j-1].N + cr[i][j-1].k3_N) ) / 
        (cr[i][j].E - cr[i][j-1].E);
    }
    
    if (j==0)
    {
        cr[i][j].gamma_rk = 2.0 * cr[i][1].gamma_rk - cr[i][2].gamma_rk;
        cr[i][j].dN_dE_rk = ( (cr[i][j+1].N + cr[i][j+1].k3_N) - (cr[i][j].N + cr[i][j].k3_N) ) / 
        (cr[i][j+1].E - cr[i][j].E);
    }
    
//    if (j==nu_channel)
//        printf ("conv4: gamma_rk = %g, gamma = %g \n", cr[i][j].gamma_rk, cr[i][j].gamma);

    cr[i][j].k4_N = h * dN_dz (cr[i][j].z + h / 2.0, cr[i][j].N + cr[i][j].k3_N, 
                                     cr[i][j].E, cr[i][j].gamma_rk, cr[i][j].dN_dE_rk);

//    if (j==50)
//        printf ("i = %d conv4: N = %g\n", i, cr[i][j].k4_N);
    

    cr[i+1][j].N = cr[i][j].N + cr[i][j].k1_N / 6.0 + cr[i][j].k2_N / 3.0 + cr[i][j].k3_N / 3.0 + cr[i][j].k4_N / 6.0;
    

//    if (j==50)
//        printf ("i = %d conv4: N+1 = %g\n", i, cr[i+1][j].N);
}
