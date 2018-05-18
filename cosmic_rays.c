#include "dec.h"
#include "cosmic_rays.h"
#include "runge_kutta.h"
#include "read_parameters.h"
#include "synchrotron.h"
#include "jet.h"

/**********************************************************************/
/*First derivative of N with respect to z*/
/*Advection equation for cosmic rays*/


double dN_dz (double z, double N, double E, double gamma, double dN_dE)
{
 
    double db_dE;

    
   if (adiabatic_losses == 1)
       adiabatic();
   
 
    if (adiabatic_losses == 1)
        b = 4.0 / 3.0 * sigma_t * c_light * 
            pow ((E / (m_electron * pow(c_light, 2))), 2) *
            (u_B[i] * (1.0 + rad_field) + u_CMB) + E / t_ad[i];
    else
        b = 4.0 / 3.0 * sigma_t * c_light * 
            pow ((E / (m_electron * pow(c_light, 2))), 2) *
            (u_B[i] * (1.0 + rad_field) + u_CMB);
    
    
    if (adiabatic_losses == 1)
        db_dE = 8.0 / 3.0 * sigma_t * c_light * 
            pow ((E / (m_electron * pow(c_light, 2))), 2) *
            (u_B[i] * (1.0 + rad_field) + u_CMB)  / E + 1./t_ad[i];
    else
        db_dE = 2.0 * b / E;

  
    return ((db_dE * N + b * dN_dE) / v_z[i]);
    

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
    db_dE = 2.0 * b / E;
    
    diff = D0 * pow (E / 1.6e-3, mu_diff);

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
            printf("Stop at i = %i j = %i\n", i, j);
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
    {
        if ( cr[i][j+1].N > 0. && cr[i][j-1].N > 0. )
            cr[i][j].dN_dE = (cr[i][j+1].N - cr[i][j-1].N) / (cr[i][j+1].E - cr[i][j-1].E);
        else if ( cr[i][j].N > 0 && cr[i][j-1].N > 0. )
            cr[i][j].dN_dE = (cr[i][j].N - cr[i][j-1].N) / (cr[i][j].E - cr[i][j-1].E);
        else
        {
            printf("Should not happen, i=%i, j=%i, N=%g\n", i, j, cr[i][j].N);
            
            cr[i][j].dN_dE = (cr[i][j+1].N - cr[i][j-1].N) / (cr[i][j+1].E - cr[i][j-1].E);

        }
        
    }
    

    if (j==0)
        cr[i][j].dN_dE = -gamma_in * cr[i][j].N / cr[i][j].E;
    
/*    if (j==0)
        cr[i][j].dN_dE = (cr[i][j+1].N - cr[i][j].N) / (cr[i][j+1].E - cr[i][j].E);
//        cr[i][j].dN_dE = 2.0 * cr[i][1].dN_dE - cr[i][2].dN_dE;*/
    
    if (j==nu_channel+1)
    {
        if ( cr[i][j].N > 0. && cr[i][j-1].N > 0. )
            cr[i][j].dN_dE = (cr[i][j].N - cr[i][j-1].N) / (cr[i][j].E - cr[i][j-1].E);
        else
        {
             printf("Should not happen, boundary\n");
             cr[i][j].dN_dE = (cr[i][j].N - cr[i][j-1].N) / (cr[i][j].E - cr[i][j-1].E);
        }
        

    }
    



    
//        cr[i][j].dN_dE = 2.0 * cr[i][j-1].dN_dE - cr[i][j-2].dN_dE;

//    if (j==nu_channel+1)
//        printf ("j = %d dN_dE = %g dN_dE=%g\n", j, cr[i][j].dN_dE, cr[i][j-1].dN_dE);

   


}


/***********************************************************************/


void adiabatic (void)

{
    double dV_dz, dr_dz;
    

    if ((i > 0) && (i < grid_size + 1))
        dV_dz = (v_z[i+1] - v_z[i-1]) / (cr[i+1][0].z - cr[i-1][0].z);
    if (i == 0)
        dV_dz = (v_z[i+1] - v_z[i]) / (cr[i+1][0].z - cr[i][0].z);
    if (i == grid_size + 1)
        dV_dz = (v_z[i] - v_z[i-1]) / (cr[i][0].z - cr[i-1][0].z);

    if ( velocity_field  == -1 )
    {
        
        if ((i > 0) && (i < grid_size + 1))
            dr_dz = (radius(cr[i+1][0].z / kpc) - radius(cr[i-1][0].z / kpc) ) / (cr[i+1][0].z - cr[i-1][0].z);
        if (i == 0)
            dr_dz = (radius(cr[i+1][0].z / kpc) - radius(cr[i][0].z / kpc) ) / (cr[i+1][0].z - cr[i][0].z);
        if (i == grid_size + 1)
            dr_dz = (radius(cr[i][0].z / kpc) - radius(cr[i-1][0].z / kpc) ) / (cr[i][0].z - cr[i-1][0].z);
    }

    else
        dr_dz = 0.0;
    
    

    if (adiabatic_losses == 1)
    {
            
/* Adiabatic losses/gains due to lateral expansion/contraction */
        if (fabs(dr_dz) > 0.)
            t_ad_r[i] = pow(2./3. * dr_dz * v_z[i] / radius(cr[i][1].z / kpc), -1.0);
        else
            t_ad_r[i] = 3.e17;

/* Adiabatic losses/gains due to longitudinal expansion/contraction */        
        if (fabs(dV_dz) > 0.)
            t_ad_v[i] = pow(dV_dz / 3., -1.0);
        else
            t_ad_v[i] = 3.e17;

        t_ad[i] = pow(1. / t_ad_r[i] + 1. / t_ad_v[i], -1.);
    }
    
}


/****************************************************************************/
struct grid_1d setup_initial_grid (void)
{
    double B_CMB;

    B_CMB = 3.2e-6;
    
    u_CMB = 1.0 / 8.0 / pi * pow (B_CMB, 2.0);
    delta_z = z_halo / ((double) grid_size);
     
        
/* Equal distance grid in space */
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

    
/* Logarithmic grid in frequency */

    delta_nu = (1000.e9 - 0.001e9) / ((double) nu_channel - 1.0);
    delta_nu_factor = exp(log(1000.e9 / 0.001e9) / ((double) nu_channel));
        
    for (i=0; i <= grid_size + 1; i++)
        cr[i][0].nu = 0.001e9;

    for (j=1; j <= nu_channel + 1; j++)
    {
        for (i=0; i <= grid_size + 1; i++)
        {
            cr[i][j].nu = cr[i][j-1].nu * delta_nu_factor;
        }
    }

/* Magnetic field strength setup */

    if (model == 1 &&  initialize_model != 1)
    {
        read_magnetic_field_model();
        for (i=0; i <= grid_size + 1; i++)
            B_field[i] = magnetic_field (cr[i][0].z / kpc);
    }

    else
    {
        

/*For galaxy_mode=1 this is a superposition of thin and thick disc.
  For galaxy_mode = -1 it is a piecewise exponential function
  For galaxy_mode = 2 the magnetic field is constant (useful for a diffusion kernel) */
        
        for (i=0; i <= grid_size + 1; i++)
        {

            if (galaxy_mode == 1)
                B_field[i] = B1 * exp(-(cr[i][0].z / kpc) / h_B1) + (B0-B1) * exp(-(cr[i][0].z / kpc) / h_B2);
            else if (galaxy_mode == -1)
            {
                if (cr[i][1].z / kpc < z1)
                    B_field[i] = B0 * exp(-cr[i][0].z / kpc / h_B1);
                else
                    B_field[i] = B0 * exp(-z1 / h_B1) * exp(-(cr[i][0].z / kpc - z1) / h_B2);
            }
            else
                B_field[i] = B0;
            
        }
            
    }


/* Magnetic field energy density */
    
    for (i=0; i <= grid_size + 1; i++)
    {
        
        u_B[i] = 1.0 / 8.0 / pi * pow (B_field[i], 2.0);
    

/* Set up the velocity distribution */
        
        if (velocity_field == 0)
            v_z[i] = V0;

        else if (velocity_field == -1)
            v_z[i] = V0 * pow(radius(cr[i][0].z / kpc) / R0, beta);
        
        
//Old version with radius not set
//v_z[i] = V0 * pow( (R0 / kpc + cr[i][0].z / kpc / h_V) / R0 * kpc , beta);            

        else if (velocity_field == 1)
            v_z[i] = V0 * exp(cr[i][0].z / kpc / h_V);
        else if (velocity_field == 2)
            v_z[i] = V0 * pow(1.0 + cr[i][0].z / kpc / h_V, beta);
        else
        {
            printf("Wrong velocity field.\n");
            printf("Stop.\n");
            exit(0);
        }
    }


    
/* Compute CRE energy in units of ergs from observed frequency
   If calculating edge-on galaxies or jets, the magnetic field
   is assumed to be ordered and perpendicular to the line of sight
   If the magnetic field is constant (galaxy_mode = 2),
   the magnetic field is assumed to be turbulent and isotropic */
    
    for (j=0; j <= nu_channel + 1; j++)
    {
        for (i=0; i <= grid_size + 1; i++)
        {
            if ( galaxy_mode != 2)
                cr[i][j].E = 1.6e-3 * pow (cr[i][j].nu / 16.1e6 / (B0 / 1.e-6), 0.5);
            else
                cr[i][j].E = 1.6e-3 * pow (cr[i][j].nu / 16.1e6 / (sqrt(2./3.) * B0 / 1.e-6), 0.5);
                        
        }
    }

/*Setup the energy distribution of the cosmic ray electrons in the galactic disk */
/* Use the injection spectral index gamma_in */

    cr[0][0].N = 1.0; //arbitrary units
    

    for (j=0; j <= nu_channel + 1; j++)
    {
        for (i=0; i <= grid_size + 1; i++)
            cr[i][j].N = cr[0][0].N * 
                pow (cr[0][j].E / cr[0][0].E, -gamma_in);
    }

    
/* For diffusion, dN/dz = 0 at z = 0 */
/* This boundary condition ensures diffusion dominates over advection at z = 0 */
    
    for ( j = 0; j <= nu_channel + 1; j++ )
        cr[0][j].y = 0.0;
    

/* The integrate true switch is that the */
    for (j=0; j <= nu_channel + 1; j++)
    {
        for (i=0; i <= grid_size + 1; i++)
            cr[i][j].integrate_true = 1;
    }

    if (model == 1 || initialize_model == 1)
    {
        
        for (i=1; i <= grid_size; i++)
            set_interpolate_values (cr[i][1].z / kpc, i);
        
    }
    

    
    return cr[i][j];
    
}

/*************************************************************************/
void output_file (int i_max)
{
    FILE *f1;
    FILE *f2;
    FILE *f3;
    FILE *f4;
    FILE *f5;
    FILE *f6;
    FILE *f7;
    FILE *f8;

    int ii, jj, kk, ii_mod;
    int nu1, nu2, nu3, nu4, nu_spec_ref, nu1_crit[402], nu2_crit[402];
    double intensity_nu1[402], intensity_nu2[402], intensity_nu3[402], intensity_nu4[402];
    double x1[402][402], x2[402][402];
    double intensity_interp1, intensity_interp2, intensity_interp3, intensity_interp4;
    double t_adv, v_z_interp;
    double grid_print_modulo;
    
        
    
    nu1 = (int) (log(nu_1 / 0.001e9) / log(delta_nu_factor));
    nu2 = (int) (log(nu_2 / 0.001e9) / log(delta_nu_factor));
    nu3 = (int) (log(nu_3 / 0.001e9) / log(delta_nu_factor));
    nu4 = (int) (log(nu_4 / 0.001e9) / log(delta_nu_factor));

/* This is for the spectral output */
    nu_spec_ref = (int) (log(0.03e9 / 0.001e9) / log(delta_nu_factor));

/* The frequency position the code finds is always a bit low */
    /* printf ("nu1=%i nu_1=%g nu_1(+1)=%g nu_1(-1)=%g\n", nu1, cr[0][nu1].nu, cr[0][nu1-1].nu, cr[0][nu1+1].nu); */
    /* printf ("nu2=%i nu_2=%g nu_2(+1)=%g nu_2(-1)=%g\n", nu2, cr[0][nu2].nu, cr[0][nu2-1].nu, cr[0][nu2+1].nu); */
    /* printf ("nu3=%i nu_3=%g nu_3(+1)=%g nu_3(-1)=%g\n", nu3, cr[0][nu3].nu, cr[0][nu3-1].nu, cr[0][nu3+1].nu); */
    /* printf ("nu4=%i nu_4=%g nu_4(+1)=%g nu_4(-1)=%g\n", nu4, cr[0][nu4].nu, cr[0][nu4-1].nu, cr[0][nu4+1].nu); */

    printf ("nu1=%i nu_1=%g Hz E1=%g GeV\n", nu1, cr[0][nu1].nu, cr[0][nu1].E / 1.6e-3);
    printf ("nu2=%i nu_2=%g Hz E2=%g GeV\n", nu2, cr[0][nu2].nu, cr[0][nu2].E / 1.6e-3);
    printf ("nu3=%i nu_3=%g Hz E3=%g GeV\n", nu3, cr[0][nu3].nu, cr[0][nu3].E / 1.6e-3);
    printf ("nu4=%i nu_4=%g Hz E4=%g GeV\n", nu4, cr[0][nu4].nu, cr[0][nu4].E / 1.6e-3);

    double gamma_mean, nu_crit_corr;

/*    JP model test */

    int spec_1 = 12;
    int spec_2 = 51;
    int spec_3 = 104;
    int spec_4 = 153;
    int spec_5 = 215;
    int spec_6 = 300;
    
/* 
    if (epsilon == 1)
        f3=fopen("./epsilon.dat", "w");
    else
*/
    
    f1=fopen("./ne.dat", "w");
    f2=fopen("./b.dat", "w");
    f3=fopen("./int.dat", "w");
    f4=fopen("./spec.dat", "w");
    f5=fopen("./ne_spec.dat", "w");
    if (model == 1 || initialize_model == 1)
    {
        
        f6=fopen("./epsilon_interp.dat", "w");
        f7=fopen("./epsilon_interp2.dat", "w");
        f8=fopen("./b2.dat", "w");
    }
   
    if (f1 == NULL)
        printf("Could not open 'n1.dat'.\n");

    for (ii=0; ii <= i_max; ii++)
    {
        for (jj=0; jj <= nu_channel+1; jj++)
        {
            nu_crit_corr = pow(B0 / B_field[ii], 1.0);
            x1[ii][jj] = nu_1 * nu_crit_corr / cr[ii][jj].nu;
            x2[ii][jj] = nu_2 * nu_crit_corr / cr[ii][jj].nu;
            nu1_crit[ii] = (int) (log(nu_crit_corr * nu_1 / 0.001e9) / log(delta_nu_factor));
            nu2_crit[ii] = (int) (log(nu_crit_corr * nu_2 / 0.001e9) / log(delta_nu_factor));
        }
    }

    for (ii=0; ii <= i_max; ii++)
    {
        intensity_nu1[ii] = synchrotron_intensity (nu_1, ii);
        intensity_nu2[ii] = synchrotron_intensity (nu_2, ii);
        intensity_nu3[ii] = synchrotron_intensity (nu_3, ii);
        intensity_nu4[ii] = synchrotron_intensity (nu_4, ii);
    }

         

    for (ii=0; ii <= i_max; ii++)
    {
        
        if (ii == spec_1)
        {
            kk = 1;
            synchrotron_spectrum (kk, ii);
        }

         if (ii == spec_2)
        {
            kk = 2;
            synchrotron_spectrum (kk, ii);
        }

         if (ii == spec_3)
        {
            kk = 3;
            synchrotron_spectrum (kk, ii);
        }

         if (ii == spec_4)
        {
            kk = 4;
            synchrotron_spectrum (kk, ii);
        }

         if (ii == spec_5)
        {
            kk = 5;
            synchrotron_spectrum (kk, ii);
        }

         if (ii == spec_6)
        {
            kk = 6;
            synchrotron_spectrum (kk, ii);
        }

  
    }

 
    fprintf(f1, "# z[kpc], n(nu_1), n(nu_2), n(nu_1_crit), n(nu_2_crit), R/R0, V/V0\n");

    
    fprintf(f2, "# z[kpc], B [G], V [cm s^-1], R [cm], t_ad [s], t_ad_R [s], t_ad_V [s], t_adv [s]\n");
    if (epsilon == 1)
        fprintf(f3, "# z[kpc], epsilon(nu_1), epsilon(nu_2), epsilon(nu_3), epsilon(nu_4), alpha(nu_1-nu_2), alpha(nu_2-nu_3), alpha(nu_3-nu_4)\n");
    else
        fprintf(f3, "# z[kpc], I(nu_1), I(nu_2), I(nu_3), I(nu_4), alpha(nu_1-nu_2), alpha(nu_2-nu_3), alpha(nu_3-nu_4)\n");
    fprintf(f4, "# nu[Hz], I(z_1), I(z_2), I(z_3), I(z_4), I(z_5), I(z_6)\n");
    fprintf(f5, "# nu[Hz], N(z_1), N(z_2), N(z_3), N(z_4), N(z_5), N(z_6)\n");

    t_adv = 0.;
    
    
    for (ii=0; ii <= i_max; ii++)
    {
       


        gamma_mean = (- log (cr[ii][nu2].N / cr[ii][nu1].N) /
                      log (cr[ii][nu2].E / cr[ii][nu1].E) );

       
        
//        if((ii-(int)(grid_delta/2.0))%grid_delta == 0)

        if (first_data_point_at_0kpc == 1)
            grid_print_modulo = 1.0;
        else
            grid_print_modulo = 2.0;
                    
        if((ii-(int)(grid_delta/grid_print_modulo))%grid_delta == 0)
        {
            

/* Output file: ne.dat */
            fprintf(f1, "% 10e % 10e % 10e % 10e % 10e \n",
                    cr[ii][0].z / kpc,
                    pow (R0 / radius(cr[ii][0].z / kpc), 2.0) * V0 / v_z[ii] *
                    interpolate_frequency(cr[ii][nu1].N, cr[ii][nu1-1].N, cr[ii][nu1+1].N, nu_1, nu1) /
                    interpolate_frequency(cr[0][nu1].N, cr[0][nu1-1].N, cr[0][nu1+1].N, nu_1, nu1),
                    pow (R0 / radius(cr[ii][0].z / kpc), 2.0) * V0 / v_z[ii] *
                    interpolate_frequency(cr[ii][nu2].N, cr[ii][nu2-1].N, cr[ii][nu2+1].N, nu_2, nu2) /
                    interpolate_frequency(cr[0][nu1].N, cr[0][nu1-1].N, cr[0][nu1+1].N, nu_1, nu1),
                    pow (R0 / radius(cr[ii][0].z / kpc), 2.0) * V0 / v_z[ii] *
                    interpolate_frequency(cr[ii][nu1_crit[ii]].N, cr[ii][nu1_crit[ii]-1].N, cr[ii][nu1_crit[ii]+1].N,
                                          pow(B0 / B_field[ii], 1.0) * nu_1, nu1_crit[ii]) /
                    interpolate_frequency(cr[0][nu1_crit[ii]].N, cr[0][nu1_crit[ii]-1].N,
                                          cr[0][nu1_crit[ii]+1].N, pow(B0 / B_field[ii], 1.0) * nu_1, nu1_crit[ii]),
                    pow (R0 / radius(cr[ii][0].z / kpc), 2.0) * V0 / v_z[ii] *
                    interpolate_frequency(cr[ii][nu2_crit[ii]].N, cr[ii][nu2_crit[ii]-1].N, cr[ii][nu2_crit[ii]+1].N,
                                          pow(B0 / B_field[ii], 1.0) * nu_2, nu2_crit[ii]) /
                    interpolate_frequency(cr[0][nu2_crit[ii]].N, cr[0][nu2_crit[ii]-1].N, cr[0][nu2_crit[ii]+1].N,
                                          pow(B0 / B_field[ii], 1.0) * nu_2, nu2_crit[ii]) );




            /*           if (model == 1 || initialize_model == 1)
            {
                
                fprintf(f1, "% 10e % 10e % 10e % 10e % 10e \n",
                        cr[ii][0].z / kpc,
                        pow (R0 / radius(cr[ii][0].z / kpc), 2.0) * V0 / v_z[ii] *
                        interpolate_frequency(cr[ii][nu1].N, cr[ii][nu1-1].N, cr[ii][nu1+1].N, nu_1, nu1) / cr[0][nu1].N,
                        pow (R0 / radius(cr[ii][0].z / kpc), 2.0) * V0 / v_z[ii] *
                        interpolate_frequency(cr[ii][nu2].N, cr[ii][nu2-1].N, cr[ii][nu2+1].N, nu_2, nu2) / cr[0][nu1].N,
                        pow (R0 / radius(cr[ii][0].z / kpc), 2.0) * V0 / v_z[ii] *
                        interpolate_frequency(cr[ii][nu1_crit[ii]].N, cr[ii][nu1_crit[ii]-1].N, cr[ii][nu1_crit[ii]+1].N,
                                              pow(B0 / B_field[ii], 1.0) * nu_1, nu1_crit[ii]) /
                        interpolate_frequency(cr[0][nu1_crit[ii]].N, cr[0][nu1_crit[ii]-1].N,
                                              cr[0][nu1_crit[ii]+1].N, pow(B0 / B_field[ii], 1.0) * nu_1, nu1_crit[ii]),
                        pow (R0 / radius(cr[ii][0].z / kpc), 2.0) * V0 / v_z[ii] *
                        interpolate_frequency(cr[ii][nu2_crit[ii]].N, cr[ii][nu2_crit[ii]-1].N, cr[ii][nu2_crit[ii]+1].N,
                                              pow(B0 / B_field[ii], 1.0) * nu_2, nu2_crit[ii]) /
                        interpolate_frequency(cr[0][nu2_crit[ii]].N, cr[0][nu2_crit[ii]-1].N, cr[0][nu2_crit[ii]+1].N,
                                              pow(B0 / B_field[ii], 1.0) * nu_2, nu2_crit[ii]) );
            }

            else

            {
                
                if (velocity_field == -1)
                    fprintf(f1, "% 10e % 10e % 10e % 10e % 10e % 10e % 10e \n",
                            cr[ii][0].z / kpc,
                            pow((R0 / kpc + cr[ii][0].z / kpc / h_V) / R0 * kpc, -2.0) * V0 / v_z[ii] *
                            interpolate_frequency(cr[ii][nu1].N, cr[ii][nu1-1].N, cr[ii][nu1+1].N, nu_1, nu1)
                            / interpolate_frequency(cr[0][nu1].N, cr[0][nu1-1].N, cr[0][nu1+1].N, nu_1, nu1),
                            V0 / v_z[ii] * interpolate_frequency(cr[ii][nu2].N, cr[ii][nu2-1].N,
                                                                 cr[ii][nu2+1].N, nu_2, nu2) / cr[0][nu1].N,
                            V0 / v_z[ii] * interpolate_frequency(cr[ii][nu1_crit[ii]].N, cr[ii][nu1_crit[ii]-1].N,
                                                                 cr[ii][nu1_crit[ii]+1].N,
                                                                 pow(B0 / B_field[ii], 1.0) * nu_1, nu1_crit[ii]) /
                            interpolate_frequency(cr[0][nu1_crit[ii]].N, cr[0][nu1_crit[ii]-1].N,
                                                  cr[0][nu1_crit[ii]+1].N, pow(B0 / B_field[ii], 1.0) * nu_1, nu1_crit[ii]),
                            V0 / v_z[ii] * interpolate_frequency(cr[ii][nu2_crit[ii]].N, cr[ii][nu2_crit[ii]-1].N,
                                                                 cr[ii][nu2_crit[ii]+1].N,
                                                                 pow(B0 / B_field[ii], 1.0) * nu_2, nu2_crit[ii]) /
                            interpolate_frequency(cr[0][nu2_crit[ii]].N, cr[0][nu2_crit[ii]-1].N, cr[0][nu2_crit[ii]+1].N,
                                                  pow(B0 / B_field[ii], 1.0) * nu_2, nu2_crit[ii]),
                            (R0 / kpc + cr[ii][0].z / kpc / h_V) / R0 * kpc, v_z[ii] / V0);
                else
                    fprintf(f1, "% 10e % 10e % 10e % 10e % 10e \n",
                            cr[ii][0].z / kpc,
                            V0 / v_z[ii] * interpolate_frequency(cr[ii][nu1].N, cr[ii][nu1-1].N, cr[ii][nu1+1].N, nu_1, nu1) / cr[0][nu1].N,
                            V0 / v_z[ii] * interpolate_frequency(cr[ii][nu2].N, cr[ii][nu2-1].N,
                                                                 cr[ii][nu2+1].N, nu_2, nu2) / cr[0][nu1].N,
                            V0 / v_z[ii] * interpolate_frequency(cr[ii][nu1_crit[ii]].N, cr[ii][nu1_crit[ii]-1].N,
                                                                 cr[ii][nu1_crit[ii]+1].N,
                                                                 pow(B0 / B_field[ii], 1.0) * nu_1, nu1_crit[ii]) /
                            interpolate_frequency(cr[0][nu1_crit[ii]].N, cr[0][nu1_crit[ii]-1].N,
                                                  cr[0][nu1_crit[ii]+1].N, pow(B0 / B_field[ii], 1.0) * nu_1, nu1_crit[ii]),
                            V0 / v_z[ii] * interpolate_frequency(cr[ii][nu2_crit[ii]].N, cr[ii][nu2_crit[ii]-1].N,
                                                                 cr[ii][nu2_crit[ii]+1].N,
                                                                 pow(B0 / B_field[ii], 1.0) * nu_2, nu2_crit[ii]) /
                            interpolate_frequency(cr[0][nu2_crit[ii]].N, cr[0][nu2_crit[ii]-1].N, cr[0][nu2_crit[ii]+1].N,
                                                  pow(B0 / B_field[ii], 1.0) * nu_2, nu2_crit[ii]));
                
                            
            
            }

            */
            
/*            fprintf(f1, "% 10e % 10e % 10e % 10e % 10e % 10e\n", */
/*                    cr[ii][0].z / kpc, */
/*                    cr[ii][nu1].N / cr[0][nu1].N, */
/*                    cr[ii][nu2].N / cr[0][nu1].N, */
/*                    cr[ii][nu1_crit[ii]].N / cr[0][nu1_crit[ii]].N, */
/*                    cr[ii][nu2_crit[ii]].N / cr[0][nu2_crit[ii]].N, */
/*                    cr[ii][(int)((nu1+nu2)/2.0)].alpha); */
           
/* /\* Output file: b.dat *\/ */
            fprintf(f2, "% 10e % 10e % 10e % 10e % 10e % 10e % 10e % 10e \n",
                    cr[ii][0].z / kpc, B_field[ii], v_z[ii], radius(cr[ii][0].z / kpc), t_ad[ii] , t_ad_r[ii], t_ad_v[ii], t_adv );
            
            /* if (model == 1 || initialize_model == 1) */
            /* { */
            /*     fprintf(f2, "% 10e % 10e % 10e % 10e % 10e % 10e % 10e % 10e \n", */
            /*             cr[ii][0].z / kpc, B_field[ii], v_z[ii], pi * pow(radius(cr[ii][0].z / kpc), 2.), t_ad[ii] , t_ad_r[ii], t_ad_v[ii], t_adv ); */
            /*     fprintf(f8, "%10e \n", B_field[ii] ); */
            /* } */
            
            /* else */
            /*     fprintf(f2, "% 10e % 10e % 10e % 10e\n", cr[ii][0].z / kpc, B_field[ii], v_z[ii], (R0 / kpc + cr[ii][0].z / kpc / h_V) / R0 * kpc); */
            

            
/* Output file: int.dat */
/* Export emissivities */            
            if ( epsilon == 1 )
            {

                if ( ii == 0 )
                    printf("radius = %g, R0 = %g\n", radius(cr[ii][0].z / kpc), R0);
                
                
                if (normalize_intensities == 1)
                    fprintf(f3, "% 10e % 10e % 10e % 10e % 10e % 10e % 10e % 10e \n",
                            cr[ii][0].z / kpc, R0 / radius(cr[ii][0].z / kpc) * R0 / radius(cr[ii][0].z / kpc) * V0 / v_z[ii] * intensity_nu1[ii] / intensity_nu1[0], R0 / radius(cr[ii][0].z / kpc) * R0 / radius(cr[ii][0].z / kpc) * V0 / v_z[ii] * intensity_nu2[ii] / intensity_nu2[0],  R0 / radius(cr[ii][0].z / kpc) * R0 / radius(cr[ii][0].z / kpc) * V0 / v_z[ii] * intensity_nu3[ii] / intensity_nu3[0], R0 / radius(cr[ii][0].z / kpc) *  R0 / radius(cr[ii][0].z / kpc) * V0 / v_z[ii] * intensity_nu4[ii] / intensity_nu4[0], log(intensity_nu1[ii]/intensity_nu2[ii])/log(nu_1/nu_2), log(intensity_nu2[ii]/intensity_nu3[ii])/log(nu_2/nu_3), log(intensity_nu3[ii]/intensity_nu4[ii])/log(nu_3/nu_4) );
                else
                    fprintf(f3, "% 10e % 10e % 10e % 10e % 10e % 10e % 10e % 10e \n",
                            cr[ii][0].z / kpc, kpc / radius(cr[ii][0].z / kpc) * R0 / radius(cr[ii][0].z / kpc) * V0 / v_z[ii] * intensity_nu1[ii] / intensity_nu1[0], kpc / radius(cr[ii][0].z / kpc) * R0 / radius(cr[ii][0].z / kpc) * V0 / v_z[ii] * intensity_nu2[ii] / intensity_nu1[0],  kpc / radius(cr[ii][0].z / kpc) * R0 / radius(cr[ii][0].z / kpc) * V0 / v_z[ii] * intensity_nu3[ii] / intensity_nu1[0], kpc / radius(cr[ii][0].z / kpc) *  R0 / radius(cr[ii][0].z / kpc) * V0 / v_z[ii] * intensity_nu4[ii] / intensity_nu1[0], log(intensity_nu1[ii]/intensity_nu2[ii])/log(nu_1/nu_2), log(intensity_nu2[ii]/intensity_nu3[ii])/log(nu_2/nu_3), log(intensity_nu3[ii]/intensity_nu4[ii])/log(nu_3/nu_4) );
                    
            }
            
/* Output file: int.dat */
/* Export intensities */              
            else
            {
                if (normalize_intensities == 1)
                    fprintf(f3, "% 10e % 10e % 10e % 10e % 10e % 10e % 10e % 10e \n",
                            cr[ii][0].z / kpc, V0 / v_z[ii] * intensity_nu1[ii] / intensity_nu1[0], V0 / v_z[ii] * intensity_nu2[ii] / intensity_nu2[0],  V0 / v_z[ii] * intensity_nu3[ii] / intensity_nu3[0],  V0 / v_z[ii] * intensity_nu4[ii] / intensity_nu4[0], log(intensity_nu1[ii]/intensity_nu2[ii])/log(nu_1/nu_2), log(intensity_nu2[ii]/intensity_nu3[ii])/log(nu_2/nu_3), log(intensity_nu3[ii]/intensity_nu4[ii])/log(nu_3/nu_4) );
                else          
                    fprintf(f3, "% 10e % 10e % 10e % 10e % 10e % 10e % 10e % 10e \n",
                            cr[ii][0].z / kpc, V0 / v_z[ii] * intensity_nu1[ii] / intensity_nu1[0], V0 / v_z[ii] * intensity_nu2[ii] / intensity_nu1[0],  V0 / v_z[ii] * intensity_nu3[ii] / intensity_nu1[0],  V0 / v_z[ii] * intensity_nu4[ii] / intensity_nu1[0], log(intensity_nu1[ii]/intensity_nu2[ii])/log(nu_1/nu_2), log(intensity_nu2[ii]/intensity_nu3[ii])/log(nu_2/nu_3), log(intensity_nu3[ii]/intensity_nu4[ii])/log(nu_3/nu_4) );
            }
            
            t_adv = t_adv + (cr[ii+1][0].z - cr[ii][0].z) / v_z[ii];


        }
	
    }

/* Output file: spec.dat */

    for (ii=1; ii <= nu_channel; ii++)
    {

        fprintf(f4, "% 10e % 10e % 10e % 10e % 10e % 10e % 10e\n",
                cr[0][ii].nu, i_syn_spec[1][ii]/i_syn_spec[1][nu_spec_ref], i_syn_spec[2][ii]/i_syn_spec[2][nu_spec_ref], i_syn_spec[3][ii]/i_syn_spec[3][nu_spec_ref], i_syn_spec[4][ii]/i_syn_spec[4][nu_spec_ref], i_syn_spec[5][ii]/i_syn_spec[5][nu_spec_ref], i_syn_spec[6][ii]/i_syn_spec[6][nu_spec_ref]);        
        
    }
    
/* Output file: ne_spec.dat */

    for (ii=1; ii <= nu_channel; ii++)
    {

        fprintf(f5, "% 10e % 10e % 10e % 10e % 10e % 10e % 10e\n",
                cr[0][ii].nu, cr[spec_1][ii].N, cr[spec_2][ii].N, cr[spec_3][ii].N, cr[spec_4][ii].N, cr[spec_5][ii].N, cr[spec_6][ii].N);        
        
    }


    if (model == 1 || initialize_model == 1)
    {
        
        for (ii=0; ii <= i_max; ii++)
            set_interpolate_values (cr[ii][0].z / kpc, ii);

//    printf("Number = %i\n", number_of_data_points);
    

        for (ii=0; ii <= number_of_data_points; ii++)
        {
            ii_mod = mod[ii].ii;
            v_z_interp = interpolated_value (v_z[ii_mod], v_z[ii_mod-1], v_z[ii_mod+1], ii, ii_mod);
            if (epsilon == 1)
            {
                
            intensity_interp1 = 1. / mod[ii].radius * R0 / mod[ii].radius / kpc * V0 / v_z_interp * interpolated_value (intensity_nu1[ii_mod], intensity_nu1[ii_mod-1], intensity_nu1[ii_mod+1], ii, ii_mod);
            intensity_interp2 = 1. / mod[ii].radius * R0 / mod[ii].radius / kpc * V0 / v_z_interp * interpolated_value (intensity_nu2[ii_mod], intensity_nu2[ii_mod-1], intensity_nu2[ii_mod+1], ii, ii_mod);
            intensity_interp3 = 1. / mod[ii].radius * R0 / mod[ii].radius / kpc * V0 / v_z_interp * interpolated_value (intensity_nu3[ii_mod], intensity_nu3[ii_mod-1], intensity_nu3[ii_mod+1], ii, ii_mod);
            intensity_interp4 = 1. / mod[ii].radius * R0 / mod[ii].radius / kpc * V0 / v_z_interp * interpolated_value (intensity_nu4[ii_mod], intensity_nu4[ii_mod-1], intensity_nu4[ii_mod+1], ii, ii_mod);

            }

            else
            {
                intensity_interp1 = V0 / v_z_interp * interpolated_value (intensity_nu1[ii_mod], intensity_nu1[ii_mod-1], intensity_nu1[ii_mod+1], ii, ii_mod);
            intensity_interp2 = V0 / v_z_interp * interpolated_value (intensity_nu2[ii_mod], intensity_nu2[ii_mod-1], intensity_nu2[ii_mod+1], ii, ii_mod);
            intensity_interp3 = V0 / v_z_interp * interpolated_value (intensity_nu3[ii_mod], intensity_nu3[ii_mod-1], intensity_nu3[ii_mod+1], ii, ii_mod);
            intensity_interp4 = V0 / v_z_interp * interpolated_value (intensity_nu4[ii_mod], intensity_nu4[ii_mod-1], intensity_nu4[ii_mod+1], ii, ii_mod);
            }
            

            fprintf(f6, "%10e %10e %10e %10e %10e %10e %10e\n",
                    mod[ii].z, intensity_interp1 / intensity_nu1[0], intensity_interp2 / intensity_nu1[0], intensity_interp3 / intensity_nu1[0], intensity_interp4 / intensity_nu1[0], -log(intensity_interp1/intensity_interp2)/log(nu_1/nu_2), v_z_interp);
            fprintf(f7, "%10e\n", intensity_interp2 / intensity_nu1[0]);
        }

    }
    

    fclose(f1);
    fclose(f2);
    fclose(f3);
    fclose(f4);
    fclose(f5);

    if (model == 1 || initialize_model == 1)
    {
        fclose(f6);
        fclose(f7);
        fclose(f8);
    }
    
}


/**********************************************************************/
double interpolate_frequency (double value, double value_low, double value_high, double nu, int jj)
{

    double interp;

    if ( nu >= cr[0][jj].nu )
        interp = value + (value_high  - value) / ( cr[0][jj+1].nu  - cr[0][jj].nu ) * ( nu - cr[0][jj].nu);
    else
        interp = value + (value  - value_low) / ( cr[0][jj].nu - cr[0][jj-1].nu ) * ( nu - cr[0][jj].nu);
    

//    printf("value=%g, value_low = %g, value_high = %g, interp=%g\n", value, value_low, value_high, interp);
    
       

    return (interp);
    

    
}

/**********************************************************************/
