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




double synchrotron (double x)

{

    double x_diff, x_min, intensity;
    int i_syn, i_syn_min;

    syn[0].x =  0.00;
    syn[1].x =  0.01; 
    syn[2].x =  0.02; 
    syn[3].x =  0.03; 
    syn[4].x =  0.04; 
    syn[5].x =  0.05; 
    syn[6].x =  0.06; 
    syn[7].x =  0.07; 
    syn[8].x =  0.08; 
    syn[9].x =  0.09; 
    syn[10].x = 0.10; 
    syn[11].x = 0.12; 
    syn[12].x = 0.14; 
    syn[13].x = 0.16; 
    syn[14].x = 0.18; 
    syn[15].x = 0.20; 
    syn[16].x = 0.22; 
    syn[17].x = 0.24; 
    syn[18].x = 0.26; 
    syn[19].x = 0.28; 
    syn[20].x = 0.30; 
    syn[21].x = 0.40; 
    syn[22].x = 0.50; 
    syn[23].x = 0.60; 
    syn[24].x = 0.70; 
    syn[25].x = 0.80; 
    syn[26].x = 0.90; 
    syn[27].x = 1.0; 
    syn[28].x = 1.1; 
    syn[29].x = 1.2; 
    syn[30].x = 1.3; 
    syn[31].x = 1.4; 
    syn[32].x = 1.5; 
    syn[33].x = 1.6; 
    syn[34].x = 1.7; 
    syn[35].x = 1.8; 
    syn[36].x = 1.9; 
    syn[37].x = 2.0; 
    syn[38].x = 2.5; 
    syn[39].x = 3.0; 
    syn[40].x = 3.5; 
    syn[41].x = 4.0; 
    syn[42].x = 4.5; 
    syn[43].x = 5.0; 
    syn[44].x = 5.5; 
    syn[45].x = 6.0; 
    syn[46].x = 6.5; 
    syn[47].x = 7.0; 
    syn[48].x = 7.5; 
    syn[49].x = 8.0; 
    syn[50].x = 8.5; 
    syn[51].x = 9.0; 
    syn[52].x = 9.5; 
    syn[53].x = 10.0; 
  
    syn[0].f =  0.0;
    syn[1].f =  0.4450; 
    syn[2].f =  0.5472; 
    syn[3].f =  0.6136; 
    syn[4].f =  0.6628; 
    syn[5].f =  0.7016; 
    syn[6].f =  0.7332; 
    syn[7].f =  0.7597; 
    syn[8].f =  0.7822; 
    syn[9].f =  0.8015; 
    syn[10].f = 0.8182; 
    syn[11].f = 0.8454; 
    syn[12].f = 0.8662; 
    syn[13].f = 0.8822; 
    syn[14].f = 0.8943; 
    syn[15].f = 0.9034; 
    syn[16].f = 0.9099; 
    syn[17].f = 0.9143; 
    syn[18].f = 0.9169; 
    syn[19].f = 0.9179; 
    syn[20].f = 0.9177; 
    syn[21].f = 0.9019; 
    syn[22].f = 0.8708; 
    syn[23].f = 0.8315; 
    syn[24].f = 0.7879; 
    syn[25].f = 0.7424; 
    syn[26].f = 0.6966; 
    syn[27].f = 0.6514; 
    syn[28].f = 0.6075; 
    syn[29].f = 0.5653; 
    syn[30].f = 0.5250; 
    syn[31].f = 0.4867; 
    syn[32].f = 0.4506; 
    syn[33].f = 0.4167; 
    syn[34].f = 0.3849; 
    syn[35].f = 0.3551; 
    syn[36].f = 0.3274; 
    syn[37].f = 0.3016; 
    syn[38].f = 0.1981; 
    syn[39].f = 0.1286; 
    syn[40].f = 0.0827; 
    syn[41].f = 0.0528; 
    syn[42].f = 0.0336; 
    syn[43].f = 0.0213; 
    syn[44].f = 0.0134; 
    syn[45].f = 8.37e-3; 
    syn[46].f = 5.30e-3; 
    syn[47].f = 3.32e-3; 
    syn[48].f = 2.01e-3; 
    syn[49].f = 1.30e-3; 
    syn[50].f = 8.12e-4; 
    syn[51].f = 5.07e-4; 
    syn[52].f = 3.18e-4; 
    syn[53].f = 1.92e-4; 

    syn[0].g =  0.0;
    syn[1].g =  0.2310; 
    syn[2].g =  0.2900; 
    syn[3].g =  0.3305; 
    syn[4].g =  0.3621; 
    syn[5].g =  0.3881; 
    syn[6].g =  0.4102; 
    syn[7].g =  0.4295; 
    syn[8].g =  0.4465; 
    syn[9].g =  0.4617; 
    syn[10].g = 0.4753; 
    syn[11].g = 0.4988; 
    syn[12].g = 0.5184; 
    syn[13].g = 0.5348; 
    syn[14].g = 0.5486; 
    syn[15].g = 0.5604; 
    syn[16].g = 0.5703; 
    syn[17].g = 0.5786; 
    syn[18].g = 0.5855; 
    syn[19].g = 0.5913; 
    syn[20].g = 0.5960; 
    syn[21].g = 0.6069; 
    syn[22].g = 0.6069; 
    syn[23].g = 0.5897; 
    syn[24].g = 0.5703; 
    syn[25].g = 0.5471; 
    syn[26].g = 0.5214; 
    syn[27].g = 0.4945; 
    syn[28].g = 0.4669; 
    syn[29].g = 0.4394; 
    syn[30].g = 0.4123; 
    syn[31].g = 0.3859; 
    syn[32].g = 0.3604; 
    syn[33].g = 0.3359; 
    syn[34].g = 0.3125; 
    syn[35].g = 0.2904; 
    syn[36].g = 0.2694; 
    syn[37].g = 0.2502; 
    syn[38].g = 0.1682; 
    syn[39].g = 0.1112; 
    syn[40].g = 0.07256; 
    syn[41].g = 0.04692; 
    syn[42].g = 0.03012; 
    syn[43].g = 0.01922; 
    syn[44].g = 0.01221; 
    syn[45].g = 7.73e-3; 
    syn[46].g = 4.87e-3; 
    syn[47].g = 3.06e-3; 
    syn[48].g = 1.92e-3; 
    syn[49].g = 1.20e-3; 
    syn[50].g = 7.52e-4; 
    syn[51].g = 4.69e-4; 
    syn[52].g = 2.92e-4; 
    syn[53].g = 1.86e-4;


    


    x_min = 100.0;
    for (i_syn = 1; i_syn <= 53; i_syn++)
    {
        x_diff = fabs(x - syn[i_syn].x);
        
       
        if (x_diff < x_min)
        {
            x_min = x_diff;
            i_syn_min = i_syn;
            
        }
    }

    intensity = syn[i_syn_min].f;

    if (x < 0.01)
        intensity = 2.1495 * pow(x, 1./3.);
    if (x > 10.)
        intensity = 1.2533 * exp(-x) * pow(x, 0.5);

    

    return (intensity);
    
//    return (syn[i_syn_min].f);
    
    
    


}


void synchrotron_spectrum (int k, int i_spec)

{
    
    int integrate_true;
    double nu_crit_corr;
    double x_spec[402][402];
    double A_spec;

    int ii, jj;
    
    nu_crit_corr = pow(B0 / B_field[i_spec], 1.0);

    for (ii=0; ii <= nu_channel+1; ii++)
    {
        for (jj=0; jj <= nu_channel+1; jj++)
        {
            
            x_spec[ii][jj] = cr[i_spec][jj].nu * nu_crit_corr / cr[i_spec][ii].nu;
                     
        }
    }
    
    for (jj=0; jj <= nu_channel; jj++)
    {

        i_syn_spec[k][jj] = 0.0;
        integrate_true = 1;
        
        for (ii=0; ii <= nu_channel; ii++)
        {

            A_spec = x_spec[ii][jj] * pow(cr[i_spec][ii].E, 2.0) / nu_crit_corr;
            
            /* if (cr[i_spec][ii].N < 0.0) */
            /*     integrate_true = -1; */
            
            
            if ( (x_spec[ii][jj] < 100.0) && (x_spec[ii][jj] > 0.001) )
            {
                if (cr[i_spec][ii].integrate_true == 1)
                    i_syn_spec[k][jj] = i_syn_spec[k][jj] + sqrt(A_spec) * cr[i_spec][ii].N *synchrotron(x_spec[ii][jj]) * pow(nu_crit_corr * cr[i_spec][jj].nu, -0.5) * pow(nu_crit_corr * cr[i_spec][ii].nu, -0.5) * (cr[i_spec][ii+1].nu - cr[i_spec][ii].nu);
                else
                    i_syn_spec[k][jj] = i_syn_spec[k][jj];

                /* if (k==6) */
                /*     printf("z=%g, jj=%i, N=%g, x=%g, syn=%g, delta_int=%g, int_nu=%g, int_true=%i\n", cr[i_spec][jj].z/kpc, jj, cr[i_spec][ii].N, x_spec[ii][jj], synchrotron(x_spec[ii][jj]), cr[i_spec][ii].N * synchrotron(x_spec[ii][jj]), i_syn_spec[k][jj], cr[i_spec][ii].integrate_true); */

            }
        }
        
        
    }

}


double synchrotron_intensity (double nu, int ii)

{

    double intensity_nu;
    double x[402];
    double nu_crit_corr;
    double A;
    int jj;
    
    
    
/* B_field[0]?*/
    
    A = (nu / 1.0e9) / (0.0168 * B_field[1] / 1.0e-10);

    for (jj=0; jj <= nu_channel+1; jj++)
    {
        nu_crit_corr = pow(B0 / B_field[ii], 1.0);
        x[jj] = nu * nu_crit_corr / cr[ii][jj].nu;
            
    }

    intensity_nu = 0.0;

 

    for (jj=0; jj <= nu_channel; jj++)
    {
        nu_crit_corr = pow(B0 / B_field[ii], 1.0);
        A = x[jj] * pow(cr[ii][jj].E, 2.0) / nu_crit_corr;
        
        /* if ((x[jj] < 10.0) && (x[jj] > 0.001) && cr[ii][jj].N < 0) */
        /*     integrate_true = -1; */

            
        if ( (x[jj] < 100.0) && (x[jj] > 0.001) )
        {
                                    
            if (cr[ii][jj].integrate_true == 1)
                intensity_nu = intensity_nu + sqrt(A) * cr[ii][jj].N * synchrotron(x[jj]) * pow(nu_crit_corr * nu, -0.5) * pow(nu_crit_corr * cr[1][jj].nu, -0.5) * (cr[ii][jj+1].nu - cr[ii][jj].nu);
            else
                intensity_nu = intensity_nu;
                
            /* if (ii==350) */
            /*     printf("Hier,\n"); */


            /* if (ii==265) */
            /*     printf("z=%g, jj=%i, N=%g, x=%g, syn=%g, delta_int=%g, int_nu=%g, int_true=%i\n", cr[ii][jj].z/kpc, jj, cr[ii][jj].N, x[jj], synchrotron(x[jj]), cr[ii][jj].N * synchrotron(x[jj]), intensity_nu, cr[ii][jj].integrate_true); */
                    
            }

        }

    return intensity_nu;
    
    
}
/**********************************************************************/
double convolve_intensity_nu1 (int ii_ref)
{

    double int_conv, sigma_beam;
    int ii;

    sigma_beam = 0.425 * FWHM_effective_beam * kpc;
    
    int_conv = 0.0;

    if (sigma_beam != 0.0)
    
        for (ii=0; ii <= grid_size; ii++)
        {

            int_conv = int_conv + intensity_nu1[ii] * exp ( -pow ( cr[ii][0].z -  cr[ii_ref][0].z, 2.0 ) /
                                                            ( 2.0 * sigma_beam * sigma_beam ) )
                / sqrt(2.0 * pi * sigma_beam * sigma_beam)* cr[ii][0].delta_z;
            
            if (ii != 0 )
                int_conv = int_conv + intensity_nu1[ii] * exp ( -pow ( -cr[ii][0].z -  cr[ii_ref][0].z, 2.0 ) /
                                                                ( 2.0 * sigma_beam * sigma_beam ) )
                    / sqrt(2.0 * pi * sigma_beam * sigma_beam)* cr[ii][0].delta_z;
            
            
        }

    else
        int_conv = intensity_nu1[ii_ref];
    
    
    return (int_conv);
    
    
}

/**********************************************************************/
double convolve_intensity_nu2 (int ii_ref)
{

    double int_conv, sigma_beam;
    int ii;
    
    sigma_beam = 0.425 * FWHM_effective_beam * kpc;

    int_conv = 0.0;

    if (sigma_beam != 0.0)
    
        for (ii=0; ii <= grid_size; ii++)
        {

            int_conv = int_conv + intensity_nu2[ii] * exp ( -pow ( cr[ii][0].z -  cr[ii_ref][0].z, 2.0 ) /
                                                            ( 2.0 * sigma_beam * sigma_beam ) )
                / sqrt(2.0 * pi * sigma_beam * sigma_beam)* cr[ii][0].delta_z;

            if (ii != 0 )
                int_conv = int_conv + intensity_nu2[ii] * exp ( -pow ( -cr[ii][0].z -  cr[ii_ref][0].z, 2.0 ) /
                                                                ( 2.0 * sigma_beam * sigma_beam ) )
                    / sqrt(2.0 * pi * sigma_beam * sigma_beam)* cr[ii][0].delta_z;
            
        }

    else
        int_conv = intensity_nu2[ii_ref];
    
    
    return (int_conv);
    
    
}
/**********************************************************************/
double convolve_intensity_nu3 (int ii_ref)
{

    double int_conv, sigma_beam;
    int ii;
    
    sigma_beam = 0.425 * FWHM_effective_beam * kpc;

    int_conv = 0.0;

    if (sigma_beam != 0.0)
    
        for (ii=0; ii <= grid_size; ii++)
        {

            int_conv = int_conv + intensity_nu3[ii] * exp ( -pow ( cr[ii][0].z -  cr[ii_ref][0].z, 2.0 ) /
                                                            ( 2.0 * sigma_beam * sigma_beam ) )
                / sqrt(2.0 * pi * sigma_beam * sigma_beam)* cr[ii][0].delta_z;

            if (ii != 0 )
                int_conv = int_conv + intensity_nu3[ii] * exp ( -pow ( -cr[ii][0].z -  cr[ii_ref][0].z, 2.0 ) /
                                                                ( 2.0 * sigma_beam * sigma_beam ) )
                    / sqrt(2.0 * pi * sigma_beam * sigma_beam)* cr[ii][0].delta_z;
    
        }

    else
        int_conv = intensity_nu3[ii_ref];
    
    
    return (int_conv);
    
    
}
/**********************************************************************/
double convolve_intensity_nu4 (int ii_ref)
{

    double int_conv, sigma_beam;
    int ii;
    
    sigma_beam = 0.425 * FWHM_effective_beam * kpc;

    int_conv = 0.0;

    if (sigma_beam != 0.0)
    
        for (ii=0; ii <= grid_size; ii++)
        {

            int_conv = int_conv + intensity_nu4[ii] * exp ( -pow ( cr[ii][0].z -  cr[ii_ref][0].z, 2.0 ) /
                                                            ( 2.0 * sigma_beam * sigma_beam ) )
                / sqrt(2.0 * pi * sigma_beam * sigma_beam)* cr[ii][0].delta_z;

            if (ii != 0 )
                int_conv = int_conv + intensity_nu4[ii] * exp ( -pow ( -cr[ii][0].z -  cr[ii_ref][0].z, 2.0 ) /
                                                                ( 2.0 * sigma_beam * sigma_beam ) )
                    / sqrt(2.0 * pi * sigma_beam * sigma_beam)* cr[ii][0].delta_z;
            
        }

    else
        int_conv = intensity_nu4[ii_ref];
    
    
    return (int_conv);
    
    
}
/**********************************************************************/
