#include "dec.h"
#include "cosmic_rays.h"
#include "runge_kutta.h"
#include "read_parameters.h"
#include "synchrotron.h"
#include "jet.h"
#include "output.h"

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

    if ( velocity_field  == -1 || velocity_field == 3 )
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
    int nu_break, ii, warning_flag;
    double fnewton, fnewton_prime, v_wind, rhs, z_crit, c_const;
    

    B_CMB = 3.2e-6;
    
    u_CMB = 1.0 / 8.0 / pi * pow (B_CMB, 2.0);
    delta_z = 2.0 * z_halo / ((double) grid_size);
     
        
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

/* Set up the velocity distribution */
    if (velocity_field == 3)
    {

        z_crit = 1.0;
        
        for (ii = 0; ii <= 100; ii++)
                
        {
            fnewton =  pow(z_crit, beta) - 2.0 * beta * pow(V0, 2.0) / pow(V_rot, 2.0) * (R0/kpc) * pow(z_crit, beta -1.0) *
                exp(z_crit / h_grav) + pow(h_V, beta);
            fnewton_prime = beta * pow(z_crit, beta - 1.0) - 2.0 * beta * pow(V0 / V_rot, 2.0) * (R0/kpc) * exp(z_crit / h_grav) *
                (pow(z_crit, beta -2.0) + pow(z_crit, beta -1.0) / h_grav);

            z_crit = z_crit - fnewton / fnewton_prime;

            
                

            if (fabs(fnewton / fnewton_prime) < 0.001)
            {
                
                printf("z_crit = %g kpc\n", z_crit);
                break;
            }
            
            
                
//            printf("fnewton = %g fnewton_prime = %g z_crit =%g \n", fnewton, fnewton_prime, z_crit);
            if (ii == 100)
                printf("Warning: solution of the critical point did not converge\n");
                
        }


        i_crit = 0;
        for (i=0; i <= grid_size + 1; i++)
            if ( fabs(cr[i][0].z - z_crit * kpc) < fabs(cr[i_crit][0].z - z_crit * kpc) )
                i_crit = i;
                
//        printf("z_crit = %g z_crit_approx = %g\n", z_crit, cr[i_crit][0].z/kpc );
        

        
        
        
        
    }
    
    
    for (i=0; i <= grid_size + 1; i++)
    {
        if (velocity_field == 0)
            v_z[i] = V0;

/* Velocity field for 3C31 */         
        else if (velocity_field == -1)
            v_z[i] = V0 * pow(radius(cr[i][0].z / kpc) / R0, beta);

/* Exponential velocoty field */
        else if (velocity_field == 1)
            v_z[i] = V0 * exp(cr[i][0].z / kpc / h_V);

/* Power-law velocity field */
        else if (velocity_field == 2)
            v_z[i] = V0 * (1.0 + pow(cr[i][0].z / kpc / h_V, beta) );

/* Wind solution */        
        else if (velocity_field == 3)
        {
            
            v_z[i] = sqrt (pow(V0 * (pow(1.0 + 4. * log(radius(cr[i][0].z / kpc) / R0 ), 1./2. ) ), 2.0)
                               - pow(25.e5 * cr[i][0].z / kpc, 2.0) );

            if (i == 0)
                v_wind = 1.0;
            else
                v_wind = v_z[i] / V0;

            if (cr[i][0].z / kpc < z_crit)
                v_wind = 0.2;
            else
                v_wind = 2.0;
            
            warning_flag = 0;
            
            for (ii = 0; ii <= 10; ii++)
                
            {

                
                c_const = 1.0 - 2.0 * log(1.0 + pow(z_crit / h_V, beta)) - pow(V_rot / V0, 2.0) * h_grav * kpc / R0 * exp(-z_crit / h_grav);
                                
                rhs = 2.0 * log(1.0 + pow(cr[i][0].z / kpc / h_V, beta)) + pow(V_rot / V0, 2.0) * h_grav * kpc / R0 * exp(-cr[i][0].z / kpc / h_grav) + c_const;
                fnewton = pow(v_wind, 2.0) - 2.0 * log(v_wind) - rhs;
                fnewton_prime = 2.0 * v_wind - 2.0 / v_wind;
                
                v_wind = v_wind - fnewton / fnewton_prime;
/*                if (ii ==0)
                    printf("v_wind = %g z = %g j =%i\n", v_wind, fnewton, j);
                if (v_wind < fnewton / fnewton_prime)
                printf("v_wind = %g z = %g j =%i\n", v_wind, fnewton, j);*/

                if ((ii == 10) && (fnewton >= 1.0e-6))
                {
                    
                    printf("Problem in the solution of the wind equation.\n");
                    printf("v_wind = %g fnewton = %g rhs = %g\n", v_wind, fnewton, rhs);
                    printf("Reset velocity at i = %i z = %g kpc\n", i, cr[i][0].z / kpc);
                    warning_flag = 1;
                    
//                    exit(0);
                }
                
                
            }

            if (warning_flag == 0)
                v_z[i] = v_wind * V0;
            else
                v_z[i] = V0;
            
            
        }
        
        else
        {
            printf("Wrong velocity field.\n");
            printf("Stop.\n");
            exit(0);
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
            if (velocity_field == 3)
            {
/* Magnetic field models from Baum et al. (1997) */                
/* Adiabatic longitudinal magnetic field */                
//                B_field[i] = B0 * R0 * R0 / pow(radius(cr[i][0].z / kpc), 2.0); 
/* Adiabatic radial and toroidal magnetic field */
                B_field[i] = B0 * R0 / radius(cr[i][0].z / kpc) * v_z[0] / v_z[i];
/* Flux freezing */                
//                B_field[i] = B0 * R0 / radius(cr[i][0].z / kpc) * pow(v_z[0] / v_z[i], 0.5);
                
            }
            
            else
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

    }

/* Magnetic field energy density */
    
    for (i=0; i <= grid_size + 1; i++)
    {
        
        u_B[i] = 1.0 / 8.0 / pi * pow (B_field[i], 2.0);
    
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

/* Option for a broken power-law injection spectrum */
    nu_break = (int) (log(10.0e9 / 0.001e9) / log(delta_nu_factor));

    for (j=0; j <= nu_channel + 1; j++)
    {
        for (i=0; i <= grid_size + 1; i++)
            cr[i][j].N = cr[0][0].N *
                pow (cr[0][j].E / cr[0][0].E, -gamma_in);
    }

    /* for (j=0; j <= nu_channel + 1; j++) */
    /* { */
    /*     if (j <= nu_break) */
    /*     { */
            
    /*         for (i=0; i <= grid_size + 1; i++) */
    /*             cr[i][j].N = cr[0][0].N * */
    /*                 pow (cr[0][j].E / cr[0][0].E, -gamma_in); */

    /*     } */

    /*     else */

    /*     { */
    /*         for (i=0; i <= grid_size + 1; i++) */
    /*             cr[i][j].N = cr[0][nu_break].N * */
    /*                 pow (cr[0][j].E / cr[0][nu_break].E, -gamma_in-0.3); */

    /*     } */
        
    /* } */

/* Advection time */
    
    for (j=0; j <= nu_channel + 1; j++)
    {
        cr[0][j].t_adv = 0.0;
        
        for (i=0; i <= grid_size; i++)
            cr[i+1][j].t_adv = cr[i][j].t_adv + (cr[i+1][0].z - cr[i][0].z) / v_z[i];
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

/****************************************************************************/
