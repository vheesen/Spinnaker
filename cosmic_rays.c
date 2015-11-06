#include"dec.h"
#include "cosmic_rays.h"
#include "runge_kutta.h"
#include "read_parameters.h"
#include "synchrotron.h"
/**********************************************************************/


double dN_dz (double z, double N, double E, double gamma, double dN_dE)
{
 
    double db_dE;
//    double t_ad = z_halo/v0 * 3.0;
//    double t_ad = 100.0e6 * 3.14e7;
    double v_z;

            
    b = 4.0 / 3.0 * sigma_t * c_light * 
        pow ((E / (m_electron * pow(c_light, 2))), 2) * (u_B[i] * (1.0 + rad_field) + u_CMB);
    bz[i][j]= 4.0 / 3.0 * sigma_t * c_light * 
        pow ((E / (m_electron * pow(c_light, 2))), 2) * u_B[i];
    db_dE = 2.0 * b / E;

// This is the derivative of the power-law distribution. Can be used as an alternative.
//    dN_dE = -gamma * N / E;
    
//    v_z = v0 * (1.0 + z/z_halo);
    if (z / parsec / 1.e3 < z0)
        v_z = v0 * (1. - z / parsec / 1.e3 / 180.);
    if ((z / parsec / 1.e3 >= z0) && (z / parsec / 1.e3 < z1 ))
        v_z = v1;
    if (z / parsec / 1.e3 >= z1)
        v_z = v2;

    /* if (j==122) */
    /*     printf("i=%i z=%g E=%g t=%g b*t/E=%g\n", i, z, E, z/v1/3.14e13, b*z/v1/E); */
    

 
    return ((db_dE * N + b * dN_dE) / v_z);
    

}

/***************************************************************************/
/*second derivative of N with respect to z*/
/*Diffusion equation for cosmic rays*/
double d2N_dz2 (double z, double N, double y, double E, double gamma, double dN_dE)
{
 
    
    double diff;
    double db_dE;
    b = 4.0 / 3.0 * sigma_t * c_light * 
        pow ((E / (m_electron * pow(c_light, 2))), 2) * (u_B[i] * (1.0 + rad_field) + u_CMB);
    bz[i][j]= 4.0 / 3.0 * sigma_t * c_light * 
        pow ((E / (m_electron * pow(c_light, 2))), 2) * u_B[i];
    db_dE = 2.0 * b / E;
//    dN_dE = -gamma * N / E;

    if (z / parsec / 1.e3 < z0)
        diff_0 = D0;
    if ((z / parsec / 1.e3 >= z0) && (z / parsec / 1.e3 < z1 ))
        diff_0 = D1;
    if (z / parsec / 1.e3 >= z1)
        diff_0 = D2;
    
    diff = diff_0 * pow (E / 1.6e-3, mu_diff);

    return ((db_dE * N + b * dN_dE) / diff);

}
/***************************************************************************/

/* Compute the spectral index: compare nonthermal emission at two separate frequencies */
void gamma_cr (void)
{

    if (j>0 && j<= nu_channel)
    {
        if (cr[i][j-1].N == 0 || cr[i][j-1].E == 0)
        {
            printf("Problem in Subroutine gamma_cr\n");
            printf("cr[i][j-1].N = %g cr[i][j-1].E = %g\n", 
                   cr[i][j-1].N, cr[i][j-1].E); 
            printf("gamma= %g\n", - log (cr[i][j+1].N / cr[i][j-1].N) / 
                   log (cr[i][j+1].E / cr[i][j-1].E));
            printf("Stop at i = %ld j = %ld\n", i, j);
            exit(0);
        }
        cr[i][j].gamma = - log (cr[i][j+1].N / cr[i][j-1].N) / 
            log (cr[i][j+1].E / cr[i][j-1].E);
    }
    
    if (j==nu_channel+1)
        cr[i][j].gamma = - log (cr[i][j].N / cr[i][j-1].N) / 
            log (cr[i][j].E / cr[i][j-1].E);
//    if (j==nu_channel+1)
//        cr[i][j].gamma = 2.0 * cr[i][j-1].gamma - cr[i][j-2].gamma;

    if (j==0)
        cr[i][j].gamma = gamma_in;
    
/*    if (j==0)
        cr[i][j].gamma = - log (cr[i][j+1].N / cr[i][j].N) / 
        log (cr[i][j+1].E / cr[i][j].E);*/
//        cr[i][j].gamma = 2.0 * cr[i][1].gamma - cr[i][2].gamma;
}

/***************************************************************************/
void dN_dE (void)
{

    if (j>0 && j<= nu_channel)
        cr[i][j].dN_dE = (cr[i][j+1].N - cr[i][j-1].N) / (cr[i][j+1].E - cr[i][j-1].E);

    if (j==0)
        cr[i][j].dN_dE = -gamma_in * cr[i][j].N / cr[i][j].E;
    
/*    if (j==0)
        cr[i][j].dN_dE = (cr[i][j+1].N - cr[i][j].N) / (cr[i][j+1].E - cr[i][j].E);
//        cr[i][j].dN_dE = 2.0 * cr[i][1].dN_dE - cr[i][2].dN_dE;*/
    
    if (j==nu_channel+1)
        cr[i][j].dN_dE = (cr[i][j].N - cr[i][j-1].N) / (cr[i][j].E - cr[i][j-1].E);
//        cr[i][j].dN_dE = 2.0 * cr[i][j-1].dN_dE - cr[i][j-2].dN_dE;

//    if (j==nu_channel+1)
//        printf ("j = %d dN_dE = %g dN_dE=%g\n", j, cr[i][j].dN_dE, cr[i][j-1].dN_dE);

   


}

/****************************************************************************/
struct grid_1d setup_initial_grid (void)
{
    double B_CMB, R;
   
    
    B_CMB = 3.2e-6;

    
    
    u_CMB = 1.0 / 8.0 / pi * pow (B_CMB, 2.0);
    delta_z = z_halo / ((double) grid_size);

     
        
/* equal distance grid*/
    for (j=0; j <= nu_channel + 1; j++)
    {
        for (i=0; i <= grid_size + 1; i++)
            cr[i][j].delta_z = delta_z;
    }
     
    
    for (j=0; j <= nu_channel + 1; j++)
        cr[0][j].z = 0.0;
    

    for (j=0; j <= nu_channel + 1; j++)
    {
        for (i=0; i <= grid_size; i++)
        {
            cr[i+1][j].z = cr[i][j].z + cr[i][j].delta_z;
        }
    }

    
/*Handle a set of frequqency channels*/
//    xxxdelta_nu = (nu_high - nu_low) / ((double) nu_channel - 1.0);
    delta_nu = (1000.e9 - 0.001e9) / ((double) nu_channel - 1.0);
    delta_nu_factor = exp(log(1000.e9 / 0.001e9) / ((double) nu_channel));
        
    for (i=0; i <= grid_size + 1; i++)
        cr[i][0].nu = 0.001e9;
//        cr[i][0].nu = nu_low - delta_nu;

    for (j=1; j <= nu_channel + 1; j++)
    {
        for (i=0; i <= grid_size + 1; i++)
        {
//            cr[i][j].nu = cr[i][j-1].nu + delta_nu;
            cr[i][j].nu = cr[i][j-1].nu * delta_nu_factor;
        }
    }


    for (i=0; i <= grid_size + 1; i++)
    {

        
        if (cr[i][1].z / parsec / 1.e3 < z0)
            B_field[i] = B0 * exp(-cr[i][1].z / parsec / 1.e3 / 190.);

        if ((cr[i][1].z / parsec / 1.e3 >= z0) && (cr[i][1].z / parsec / 1.e3 <= z1))
        {
            
            B_field[i] = B0 * exp(-z0 / 190.) * exp(-(cr[i][1].z / parsec / 1.e3 - z0) / h_B1);

            if (power_law == 1)

            {
                R = R0 + (R1 - R0) * (cr[i][1].z / parsec / 1.e3 - z0) / (z1 - z0);
                                    
                B_field[i] = B0 * exp(-z0 / 190.) * pow(R0 / R, 2.);
            }
            
            
        }
        
        if (cr[i][1].z / parsec / 1.e3 >= z1)
        {
            
            B_field[i] = B0 * exp(-z0 / 190.) * exp(-(z1-z0) / h_B1) * exp(-(cr[i][1].z / parsec / 1.e3 - z1) / h_B2);

            if (power_law == 1)

            {

                R = R1 + (R2 - R1) * (cr[i][1].z / parsec / 1.e3 - z1) / (z_halo_parsec / 1.e3- z1);
                
                B_field[i] = B0 * exp(-z0 / 190.) * pow(R0 / R1, 2.) * pow(R1 / R, 2.);
                
            }
            

        }
        

        u_B[i] = 1.0 / 8.0 / pi * pow (B_field[i], 2.0);
    }

    
/* Compute energy from observed frequency*/
/* make the assumption that all energy of a electron at energy E is emitted at a single frequency*/
    for (j=0; j <= nu_channel + 1; j++)
    {
        for (i=0; i <= grid_size + 1; i++)
        {
            cr[i][j].E = 1.6e-3 * pow (cr[i][j].nu / 16.1e6 / (B0 / 1.e-6), 0.5);
                        
        }
    }

/*Setup the energy distribution of the cosmic ray electrons 
  in the galactic disk*/
/*use the known injection spectral index gamma_in*/

    cr[1][1].N = 1.0; //arbitrary units
    

    for (j=0; j <= nu_channel + 1; j++)
    {
        for (i=0; i <= grid_size + 1; i++)
            cr[i][j].N = cr[1][1].N * 
                pow (cr[1][j].E / cr[1][1].E, -gamma_in);
    }

    printf ("E(1)=%g\n", cr[1][1].E);
    
    return cr[i][j];
    

}

/*************************************************************************/
void output_file (long i_max)
{
    FILE *f1;
    FILE *f2;
    FILE *f3;
    FILE *f4;
    int ii, jj;
    int nu_1, nu_2, nu_3, nu_spec_ref;
    double intensity_nu1[402], intensity_nu2[402], i_syn1[402], i_syn2[402];
    double x1[402][402], x2[402][402];
    double integrate_true;
    
    
    nu_1 = (int) ((nu_low - 0.001e9) / delta_nu);
    nu_2 = (int) ((nu_high - 0.001e9) / delta_nu);

    nu_1 = (int) (log(nu_low / 0.001e9) / log(delta_nu_factor));
    nu_2 = (int) (log(nu_high / 0.001e9) / log(delta_nu_factor));
    nu_3 = (int) (log(1.0e9 / 0.001e9) / log(delta_nu_factor));
    nu_spec_ref = (int) (log(0.03e9 / 0.001e9) / log(delta_nu_factor));
    
    printf ("nu_1=%i %g\n", nu_1, cr[1][nu_1].nu);
    printf ("nu_2=%i %g\n", nu_2, cr[1][nu_2].nu);
    printf ("log(2.7) = %g %g\n", delta_nu_factor, nu_low);
    

    double A1;
    double A2;
    double syn_z, nu_crit_corr, gamma_mean;

/*
//    Northern radio tail
    long spec_1=8;
    long spec_2=16;
    long spec_3=30;
    long spec_4=43;
    long spec_5=57;
    long spec_6=72;
*/
//    JP model test
    long spec_1=0;
    long spec_2=20;
    long spec_3=200;
    long spec_4=200;
    long spec_5=200;
    long spec_6=200;

    
    double alpha1, alpha2, alpha3, alpha4, alpha5, alpha6;
    long integrate_true_1, integrate_true_2, integrate_true_3, integrate_true_4, integrate_true_5, integrate_true_6;
//    double x_spec_1[402][402], x_spec_2[402][402], x_spec_3[402][402], x_spec_4[402][402], x_spec_5[402][402], x_spec_6[402][402];
    double nu_crit_corr_1, nu_crit_corr_2, nu_crit_corr_3, nu_crit_corr_4, nu_crit_corr_5, nu_crit_corr_6;
//    double i_syn_spec_1[402], i_syn_spec_2[402], i_syn_spec_3[402], i_syn_spec_4[402], i_syn_spec_5[402], i_syn_spec_6[402];
    double A_spec_1, A_spec_2, A_spec_3, A_spec_4, A_spec_5, A_spec_6;
    
//    printf("z(spec_4) =%g\n", cr[spec_4][1].z / parsec / 1.e3);
    
    
   f1=fopen("./n1.dat", "w");
   f2=fopen("./b.dat", "w");
   f3=fopen("./int.dat", "w");
   f4=fopen("./spec.dat", "w");

    if (f1 == NULL)
	printf("Could not open 'n1.dat'.\n");

    A1 = (cr[1][nu_1].nu / 1.0e9)/ (0.0168 * B_field[1] / 1.0e-10);
    A2 = (cr[1][nu_2].nu / 1.0e9)/ (0.0168 * B_field[1] / 1.0e-10);

    /* printf ("A1= %g\n", A1); */
    /* printf ("A2= %g\n", A2); */


    for (ii=1; ii <= i_max; ii++)
    {
        for (jj=1; jj <= nu_channel+1; jj++)
        {
            nu_crit_corr = pow(B0 / B_field[ii], 1.0);
            x1[ii][jj] = cr[1][nu_1].nu * nu_crit_corr / cr[ii][jj].nu;
            x2[ii][jj] = cr[1][nu_2].nu * nu_crit_corr / cr[ii][jj].nu;
            if (i==1)
                printf("x1=%g\n", x1[ii][jj]);
            
        }

    }

    nu_crit_corr_1 = pow(B0 / B_field[spec_1], 1.0);
    nu_crit_corr_2 = pow(B0 / B_field[spec_2], 1.0);
    nu_crit_corr_3 = pow(B0 / B_field[spec_3], 1.0);
    nu_crit_corr_4 = pow(B0 / B_field[spec_4], 1.0);
    nu_crit_corr_5 = pow(B0 / B_field[spec_5], 1.0);
    nu_crit_corr_6 = pow(B0 / B_field[spec_6], 1.0);
            
    /* for (ii=1; ii <= nu_channel+1; ii++) */
    /* { */
    /*     for (jj=1; jj <= nu_channel+1; jj++) */
    /*     { */
         
    /*         x_spec_1[ii][jj] = cr[spec_1][jj].nu * nu_crit_corr_1 / cr[spec_1][ii].nu; */
    /*         x_spec_2[ii][jj] = cr[spec_2][jj].nu * nu_crit_corr_2 / cr[spec_2][ii].nu; */
    /*         x_spec_3[ii][jj] = cr[spec_3][jj].nu * nu_crit_corr_3 / cr[spec_3][ii].nu; */
    /*         x_spec_4[ii][jj] = cr[spec_4][jj].nu * nu_crit_corr_4 / cr[spec_4][ii].nu; */
    /*         x_spec_5[ii][jj] = cr[spec_5][jj].nu * nu_crit_corr_5 / cr[spec_5][ii].nu; */
    /*         x_spec_6[ii][jj] = cr[spec_6][jj].nu * nu_crit_corr_6 / cr[spec_6][ii].nu; */
                     
    /*     } */
    /* } */
    


    
    for (ii=1; ii <= i_max; ii++)
    {
         
        intensity_nu1[ii] = 0.0;
        intensity_nu2[ii] = 0.0;

    
        /* if (ii==100) */
        /* { */
            
        /*     printf ("A1= %g\n", A1); */
        /*     printf ("A2= %g\n", A2); */
            
        /* } */

        integrate_true = 1;

        
        for (jj=1; jj <= nu_channel; jj++)
        {
            nu_crit_corr = pow(B0 / B_field[ii], 1.0);
            A1 = x1[ii][jj] * pow(cr[ii][jj].E, 2.0) / nu_crit_corr;
            A2 = x2[ii][jj] * pow(cr[ii][jj].E, 2.0) / nu_crit_corr;



            if (ii==100 && jj==100)
                printf ("A1= %g\n", A1);

            if (cr[ii][jj].N < 0.0)
                    integrate_true = -1;
            
            if ( (x1[ii][jj] < 100.0) && (x1[ii][jj] > 0.001) )
            {
                
                if (integrate_true == 1)
                    intensity_nu1[ii] = intensity_nu1[ii] + sqrt(A1) * cr[ii][jj].N *synchrotron(x1[ii][jj]) * pow(nu_crit_corr * cr[1][nu_1].nu, -0.5) * pow(nu_crit_corr * cr[1][jj].nu, -0.5) * (cr[ii][jj+1].nu - cr[ii][jj].nu);
                else
                    intensity_nu1[ii] = intensity_nu1[ii];
                


                /* if (ii==50) */
                /*     printf("jj=%i, N=%g, x1=%g, syn=%g, delta_int=%g, int_nu1=%g\n", jj, cr[ii][jj].N, x1[ii][jj], synchrotron(x1[ii][jj]), cr[ii][jj].N * synchrotron(x1[ii][jj]), intensity_nu1[ii]); */
                
            }

            if ( (x2[ii][jj] < 100.0) && (x2[ii][jj] > 0.001))
            {
                if (integrate_true == 1)
                    intensity_nu2[ii] = intensity_nu2[ii] + sqrt(A2) * cr[ii][jj].N *synchrotron(x2[ii][jj]) * pow(nu_crit_corr * cr[1][nu_2].nu, -0.5) * pow(nu_crit_corr * cr[1][jj].nu, -0.5) * (cr[ii][jj+1].nu - cr[ii][jj].nu);
                else
                    intensity_nu2[ii] = intensity_nu2[ii];



                                    /* intensity_nu2[ii] = intensity_nu2[ii] - -0.5 * sqrt(A2) * cr[ii][jj].N *synchrotron(x2[ii][jj]) * pow(x2[ii][jj], -1.5) * (x2[ii][jj+1] - x2[ii][jj]); */
 

                /* if (ii==50) */
                /*     printf("jj=%i, N=%g, x2=%g, syn=%g, delta_int=%g, int_nu1=%g, switch=%i\n", jj, cr[ii][jj].N, x2[ii][jj], synchrotron(x2[ii][jj]), sqrt(A2) * cr[ii][jj].N * synchrotron(x2[ii][jj]) * pow(nu_crit_corr * cr[1][nu_2].nu, -0.5) * pow(nu_crit_corr * cr[1][jj].nu, -0.5) * (cr[ii][jj+1].nu - cr[ii][jj].nu), intensity_nu2[ii], integrate_true); */

            }
        }
    }


/*     for (jj=1; jj <= nu_channel; jj++) */
/*     { */

/*         i_syn_spec_1[jj] = 0.0; */
/*         i_syn_spec_2[jj] = 0.0; */
/*         i_syn_spec_3[jj] = 0.0; */
/*         i_syn_spec_4[jj] = 0.0; */
/*         i_syn_spec_5[jj] = 0.0; */
/*         i_syn_spec_6[jj] = 0.0; */

/*         integrate_true_1 = 1; */
/*         integrate_true_2 = 1; */
/*         integrate_true_3 = 1; */
/*         integrate_true_4 = 1; */
/*         integrate_true_5 = 1; */
/*         integrate_true_6 = 1; */


/*         for (ii=1; ii <= nu_channel; ii++) */
/*         { */

            
/* //            if (jj==100) */
/* //                printf("x=%g f(x)=%g\n", x_spec_1[ii][jj], synchrotron(x_spec_1[ii][jj])); */
 
/*             A_spec_1 = x_spec_1[ii][jj] * pow(cr[spec_1][ii].E, 2.0) / nu_crit_corr_1; */
/*             A_spec_2 = x_spec_2[ii][jj] * pow(cr[spec_2][ii].E, 2.0) / nu_crit_corr_2; */
/*             A_spec_3 = x_spec_3[ii][jj] * pow(cr[spec_3][ii].E, 2.0) / nu_crit_corr_3; */
/*             A_spec_4 = x_spec_4[ii][jj] * pow(cr[spec_4][ii].E, 2.0) / nu_crit_corr_4; */
/*             A_spec_5 = x_spec_5[ii][jj] * pow(cr[spec_5][ii].E, 2.0) / nu_crit_corr_5; */
/*             A_spec_6 = x_spec_6[ii][jj] * pow(cr[spec_6][ii].E, 2.0) / nu_crit_corr_6; */



/*             if (cr[spec_1][ii].N < 0.0) */
/*                 integrate_true_1 = -1; */
/*             if (cr[spec_2][ii].N < 0.0) */
/*                 integrate_true_2 = -1; */
/*             if (cr[spec_3][ii].N < 0.0) */
/*                 integrate_true_3 = -1; */
/*             if (cr[spec_4][ii].N < 0.0) */
/*                 integrate_true_4 = -1; */
/*             if (cr[spec_5][ii].N < 0.0) */
/*                 integrate_true_5 = -1; */
/*             if (cr[spec_6][ii].N < 0.0) */
/*                 integrate_true_6 = -1; */


            
/*             if ( (x_spec_1[ii][jj] < 100.0) && (x_spec_1[ii][jj] > 0.001) ) */
/*             { */
/*                 if (integrate_true_1 == 1) */
/*                     i_syn_spec_1[jj] = i_syn_spec_1[jj] + sqrt(A_spec_1) * cr[spec_1][ii].N *synchrotron(x_spec_1[ii][jj]) * pow(nu_crit_corr_1 * cr[spec_1][jj].nu, -0.5) * pow(nu_crit_corr_1 * cr[spec_1][ii].nu, -0.5) * (cr[spec_1][ii+1].nu - cr[spec_1][ii].nu); */
/*                 else */
/*                     i_syn_spec_1[jj] = i_syn_spec_1[jj]; */
/*             } */

/*             /\* if (jj==50) *\/ */
/*             /\*     printf("ii=%i, nu=%g, N=%g, x1=%g, syn=%g, delta_int=%g, int_nu1=%g, switch=%i\n", ii, cr[spec_1][ii].nu, cr[spec_1][ii].N, x_spec_1[ii][jj], synchrotron(x_spec_1[ii][jj]), sqrt(A_spec_1) * cr[spec_1][ii].N *synchrotron(x_spec_1[ii][jj]) * pow(nu_crit_corr_1 * cr[spec_1][ii].nu, -0.5) * pow(nu_crit_corr_1 * cr[spec_1][ii].nu, -0.5) * (cr[spec_1][ii+1].nu - cr[spec_1][ii].nu), i_syn_spec_1[ii], integrate_true_1); *\/ */


/*             if ( (x_spec_2[ii][jj] < 100.0) && (x_spec_2[ii][jj] > 0.001) ) */
/*             { */
/*                 if (integrate_true_2 == 1) */
/*                     i_syn_spec_2[jj] = i_syn_spec_2[jj] + sqrt(A_spec_2) * cr[spec_2][ii].N *synchrotron(x_spec_2[ii][jj]) * pow(nu_crit_corr_2 * cr[spec_2][jj].nu, -0.5) * pow(nu_crit_corr_2 * cr[spec_2][ii].nu, -0.5) * (cr[spec_2][ii+1].nu - cr[spec_2][ii].nu); */
/*                 else */
/*                     i_syn_spec_2[jj] = i_syn_spec_2[jj]; */

/*                 /\* if (jj==114) *\/ */
/*                 /\*     printf("ii=%i, nu=%g, N=%g, x2=%g, syn=%g, delta_int=%g, int_nu1=%g, switch=%i\n", ii, cr[spec_2][ii].nu, cr[spec_2][ii].N, x_spec_2[ii][jj], synchrotron(x_spec_2[ii][jj]), sqrt(A_spec_2) * cr[spec_2][ii].N *synchrotron(x_spec_2[ii][jj]) * pow(nu_crit_corr_2 * cr[spec_2][ii].nu, -0.5) * pow(nu_crit_corr_2 * cr[spec_2][ii].nu, -0.5) * (cr[spec_2][ii+1].nu - cr[spec_2][ii].nu), i_syn_spec_2[ii], integrate_true_2); *\/ */

/*             } */


/*             /\* if (jj==114) *\/ */
/*             /\*     printf("ii=%i, nu=%g, N=%g, x2=%g, syn=%g, switch=%i\n", ii, cr[spec_2][ii].nu, cr[spec_2][ii].N, x_spec_2[ii][jj], synchrotron(0.007), integrate_true_2); *\/ */


/*             if ( (x_spec_3[ii][jj] < 100.0) && (x_spec_3[ii][jj] > 0.001) ) */
/*             { */
/*                 if (integrate_true_3 == 1) */
/*                     i_syn_spec_3[jj] = i_syn_spec_3[jj] + sqrt(A_spec_3) * cr[spec_3][ii].N *synchrotron(x_spec_3[ii][jj]) * pow(nu_crit_corr_3 * cr[spec_3][jj].nu, -0.5) * pow(nu_crit_corr_3 * cr[spec_3][ii].nu, -0.5) * (cr[spec_3][ii+1].nu - cr[spec_3][ii].nu); */
/*                 else */
/*                     i_syn_spec_3[jj] = i_syn_spec_3[jj]; */

/*                 /\* if (jj==133) *\/ */
/*                 /\*     printf("ii=%i, nu=%g, N=%g, x3=%g, syn=%g, delta_int=%g, int_nu1=%g, switch=%i\n", ii, cr[spec_3][ii].nu, cr[spec_3][ii].N, x_spec_3[ii][jj], synchrotron(x_spec_3[ii][jj]), sqrt(A_spec_3) * cr[spec_3][ii].N *synchrotron(x_spec_3[ii][jj]) * pow(nu_crit_corr_3 * cr[spec_3][ii].nu, -0.5) * pow(nu_crit_corr_3 * cr[spec_3][ii].nu, -0.5) * (cr[spec_3][ii+1].nu - cr[spec_3][ii].nu), i_syn_spec_3[jj], integrate_true_3); *\/ */

                
/*             } */

/*             if ( (x_spec_4[ii][jj] < 100.0) && (x_spec_4[ii][jj] > 0.001) ) */
/*             { */
/*                 if (integrate_true_4 == 1) */
/*                     i_syn_spec_4[jj] = i_syn_spec_4[jj] + sqrt(A_spec_4) * cr[spec_4][ii].N *synchrotron(x_spec_4[ii][jj]) * pow(nu_crit_corr_4 * cr[spec_4][jj].nu, -0.5) * pow(nu_crit_corr_4 * cr[spec_4][ii].nu, -0.5) * (cr[spec_4][ii+1].nu - cr[spec_4][ii].nu); */
/*                 else */
/*                     i_syn_spec_4[jj] = i_syn_spec_4[jj]; */
/*             } */

/*             if ( (x_spec_5[ii][jj] < 100.0) && (x_spec_5[ii][jj] > 0.001) ) */
/*             { */
/*                 if (integrate_true_5 == 1) */
/*                     i_syn_spec_5[jj] = i_syn_spec_5[jj] + sqrt(A_spec_5) * cr[spec_5][ii].N *synchrotron(x_spec_5[ii][jj]) * pow(nu_crit_corr_5 * cr[spec_5][jj].nu, -0.5) * pow(nu_crit_corr_5 * cr[spec_5][ii].nu, -0.5) * (cr[spec_5][ii+1].nu - cr[spec_5][ii].nu); */
/*                 else */
/*                     i_syn_spec_5[jj] = i_syn_spec_5[jj]; */
/*             } */

/*             if ( (x_spec_6[ii][jj] < 100.0) && (x_spec_6[ii][jj] > 0.001) ) */
/*             { */
/*                 if (integrate_true_6 == 1) */
/*                     i_syn_spec_6[jj] = i_syn_spec_6[jj] + sqrt(A_spec_6) * cr[spec_6][ii].N *synchrotron(x_spec_6[ii][jj]) * pow(nu_crit_corr_6 * cr[spec_6][jj].nu, -0.5) * pow(nu_crit_corr_6 * cr[spec_6][ii].nu, -0.5) * (cr[spec_6][ii+1].nu - cr[spec_6][ii].nu); */
/*                 else */
/*                     i_syn_spec_6[jj] = i_syn_spec_6[jj]; */
/*             } */

/*         } */
        
/*     } */
   

     fprintf(f1, "# z[kpc], N1, N2, N1_norm, N2_norm, alpha\n");
     fprintf(f3, "# z[kpc], Isyn1, Isyn2, alpha_syn, I1, I2, alpha\n");
    
    for (ii=1; ii <= i_max; ii++)
    {
        
        syn_z = pow(B_field[ii], (gamma_in + 1.0) / 2.0);
        

        i_syn1[ii] = cr[ii][nu_1].N / cr[1][nu_1].N * syn_z;
        i_syn2[ii] = cr[ii][nu_2].N / cr[1][nu_2].N * syn_z;
        
        
        gamma_mean = (- log (cr[ii][nu_2].N / cr[ii][nu_1].N) /
                      log (cr[ii][nu_2].E / cr[ii][nu_1].E) );

       
        
                    
        if((ii-(int)(grid_delta/2.0))%grid_delta == 0)
        {


            fprintf(f1, "%10g %10g %10g %10g %10g %10g\n",
                    cr[ii][1].z / parsec / 1000.0, cr[ii][nu_1].N, cr[ii][nu_2].N,
                    cr[ii][nu_1].N / cr[1][nu_1].N,
                    cr[ii][nu_2].N / cr[1][nu_2].N,
                    cr[ii][(int)((nu_1+nu_2)/2.0)].alpha);
            
            fprintf(f2, "%10g %10g %10g\n",
                    cr[ii][1].z / parsec / 1000.0, B_field[ii], B_field[ii] / B_field[1]);
            fprintf(f3, "%10g %10g %10g %10g %10g %10g %10g\n",
                    cr[ii][1].z / parsec / 1000.0, intensity_nu1[ii] / intensity_nu1[1], intensity_nu2[ii] / intensity_nu1[1], -log(intensity_nu1[ii]/intensity_nu2[ii])/log(cr[1][nu_1].nu/cr[1][nu_2].nu),  i_syn1[ii] / i_syn1[1],  i_syn2[ii] / i_syn2[1], -log(i_syn1[ii]/i_syn2[ii]) / log(nu_low/nu_high) + (gamma_in-1.0)/2.0);

            
            
        }
	
    }

    
    for (ii=1; ii <= nu_channel; ii++)
    {

       

        /* i_syn1_int[ii] = cr[spec1][ii].N * pow(B_field[spec1], (gamma_in + 1.0) / 2.0); */
        /* i_syn2_int[ii] = cr[spec2][ii].N * pow(B_field[spec2], (gamma_in + 1.0) / 2.0); */
        /* i_syn3_int[ii] = cr[spec3][ii].N * pow(B_field[spec3], (gamma_in + 1.0) / 2.0); */
        /* i_syn4_int[ii] = cr[spec4][ii].N * pow(B_field[spec4], (gamma_in + 1.0) / 2.0); */
        /* i_syn5_int[ii] = cr[spec5][ii].N * pow(B_field[spec5], (gamma_in + 1.0) / 2.0); */
        /* i_syn6_int[ii] = cr[spec6][ii].N * pow(B_field[spec6], (gamma_in + 1.0) / 2.0); */
    }

    for (ii=1; ii <= nu_channel; ii++)
    {


        
        /* alpha1 = -log(i_syn1_int[ii]/i_syn1_int[1]) / log(cr[1][ii].nu / cr[1][1].nu) - (gamma_in-1.0)/2.0; */
        /* alpha2 = -log(i_syn2_int[ii]/i_syn2_int[1]) / log(cr[1][ii].nu / cr[1][1].nu) - (gamma_in-1.0)/2.0; */
        /* alpha3 = -log(i_syn3_int[ii]/i_syn3_int[1]) / log(cr[1][ii].nu / cr[1][1].nu) - (gamma_in-1.0)/2.0; */
        /* alpha4 = -log(i_syn4_int[ii]/i_syn4_int[1]) / log(cr[1][ii].nu / cr[1][1].nu) - (gamma_in-1.0)/2.0; */
        /* alpha5 = -log(i_syn5_int[ii]/i_syn5_int[1]) / log(cr[1][ii].nu / cr[1][1].nu) - (gamma_in-1.0)/2.0; */
        /* alpha6 = -log(i_syn6_int[ii]/i_syn6_int[1]) / log(cr[1][ii].nu / cr[1][1].nu) - (gamma_in-1.0)/2.0; */


        
        


        /* i_syn1[ii] = pow(cr[1][ii].nu / cr[1][1].nu, -alpha1); */
        /* i_syn2[ii] = pow(cr[1][ii].nu / cr[1][1].nu, -alpha2); */
        /* i_syn3[ii] = pow(cr[1][ii].nu / cr[1][1].nu, -alpha3); */
        /* i_syn4[ii] = pow(cr[1][ii].nu / cr[1][1].nu, -alpha4); */
        /* i_syn5[ii] = pow(cr[1][ii].nu / cr[1][1].nu, -alpha5); */
        /* i_syn6[ii] = pow(cr[1][ii].nu / cr[1][1].nu, -alpha6); */
        
/*
        fprintf(f4, "%10g %10g %10g %10g %10g %10g %10g\n",
                cr[1][ii].nu, i_syn_spec_1[ii]/i_syn_spec_1[nu_spec_ref], i_syn_spec_2[ii]/i_syn_spec_2[nu_spec_ref], i_syn_spec_3[ii]/i_syn_spec_3[nu_spec_ref], i_syn_spec_4[ii]/i_syn_spec_4[nu_spec_ref], i_syn_spec_5[ii]/i_syn_spec_5[nu_spec_ref], i_syn_spec_6[ii]/i_syn_spec_6[nu_spec_ref]);
*/        
        
    }
    
   

    fclose(f1);
    fclose(f2);
    fclose(f3);
    fclose(f4);
}


/**********************************************************************/
