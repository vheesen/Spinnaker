#include "dec.h"
#include "cosmic_rays.h"
#include "runge_kutta.h"
#include "read_parameters.h"
#include "synchrotron.h"
#include "jet.h"
/***********************************************************************/
int set_radius_north ( void )

{

    int counter_max;
            
    mod[0].z =  3.7499999999999943E-002  ;
    mod[1].z =  0.11249999999999993      ; 
    mod[2].z =  0.18749999999999994      ; 
    mod[3].z =  0.26249999999999996      ; 
    mod[4].z =  0.33749999999999991      ; 
    mod[5].z =  0.41249999999999998      ; 
    mod[6].z =  0.48749999999999993      ; 
    mod[7].z =  0.56249999999999989      ; 
    mod[8].z =  0.63749999999999996      ; 
    mod[9].z =  0.71249999999999991      ; 
    mod[10].z = 0.78749999999999987      ; 
    mod[11].z = 0.86249999999999993      ; 
    mod[12].z = 0.93749999999999989      ; 
    mod[13].z = 1.0124999999999997       ; 
    mod[14].z = 1.0874999999999999       ; 

/* //Judith's radius values     */
/*     mod[0].radius =  2.9; */
/*     mod[1].radius =  2.9; */
/*     mod[2].radius =  2.9; */
/*     mod[3].radius =  2.9; */
/*     mod[4].radius =  2.9; */
/*     mod[5].radius =  2.9; */
/*     mod[6].radius =  2.9; */
/*     mod[7].radius =  2.9; */
/*     mod[8].radius =  3.4; */
/*     mod[9].radius =  3.9; */
/*     mod[10].radius = 4.4; */
/*     mod[11].radius = 4.8; */
/*     mod[12].radius = 4.9; */
/*     mod[13].radius = 4.9; */
/*     mod[14].radius = 4.9; */
/*     mod[15].radius = 4.8; */
/*     mod[16].radius = 4.8; */
/*     mod[17].radius = 4.9; */
/*     mod[18].radius = 5.0; */
/*     mod[19].radius = 5.1; */
/*     mod[20].radius = 5.1; */
/*     mod[21].radius = 5.0; */
/*     mod[22].radius = 5.0; */
/*     mod[23].radius = 4.9; */
/*     mod[24].radius = 4.9; */
/*     mod[25].radius = 5.2; */
    
//My radius values        
    mod[0].radius =  1.0;
    mod[1].radius =  1.0;
    mod[2].radius =  1.0;
    mod[3].radius =  1.0;
    mod[4].radius =  1.0;
    mod[5].radius =  1.0;
    mod[6].radius =  1.0;
    mod[7].radius =  1.0;
    mod[8].radius =  1.0;
    mod[9].radius =  1.0;
    mod[10].radius = 1.0;
    mod[11].radius = 1.0;
    mod[12].radius = 1.0;
    mod[13].radius = 1.0;
    mod[14].radius = 1.0;
 

    mod[0].intensity =  22542.218027048875;
    mod[1].intensity =  16422.703969676702;
    mod[2].intensity =  8495.0737633390581;
    mod[3].intensity =  3839.3229612047526;
    mod[4].intensity =  2066.4930081369598;
    mod[5].intensity =  1637.1505038466182;
    mod[6].intensity =  1690.4037358798128;
    mod[7].intensity =  1691.8618051409182;
    mod[8].intensity =  1558.7999215636960;
    mod[9].intensity =  1292.6620792394947;
    mod[10].intensity = 1058.6352482678687;
    mod[11].intensity = 704.85684752413738;
    mod[12].intensity = 498.05118835832411;
    mod[13].intensity = 374.08472057874548;
    mod[14].intensity = 247.62929432337728;

/* 52-145 MHz spectral index */
    mod[0].alpha =  -0.55;
    mod[1].alpha =  -0.55;
    mod[2].alpha =  -0.55;
    mod[3].alpha =  -0.55;
    mod[4].alpha =  -0.55;
    mod[5].alpha =  -0.55;
    mod[6].alpha =  -0.55;
    mod[7].alpha =  -0.55;
    mod[8].alpha =  -0.55;
    mod[9].alpha =  -0.55;
    mod[10].alpha = -0.55;
    mod[11].alpha = -0.55;
    mod[12].alpha = -0.55;
    mod[13].alpha = -0.55;
    mod[14].alpha = -0.55;

    counter_max = 14;

    return (counter_max);
    
    

}
/***********************************************************************/
int set_radius_south ( void )

{
    int counter_max;

    mod[0].z =  3.7499999999999943E-002  ;
    mod[1].z =  0.11249999999999993      ; 
    mod[2].z =  0.18749999999999994      ; 
    mod[3].z =  0.26249999999999996      ; 
    mod[4].z =  0.33749999999999991      ; 
    mod[5].z =  0.41249999999999998      ; 
    mod[6].z =  0.48749999999999993      ; 
    mod[7].z =  0.56249999999999989      ; 
    mod[8].z =  0.63749999999999996      ; 
    mod[9].z =  0.71249999999999991      ; 
    mod[10].z = 0.78749999999999987      ; 
    mod[11].z = 0.86249999999999993      ; 
    mod[12].z = 0.93749999999999989      ; 
    mod[13].z = 1.0124999999999997       ; 
    mod[14].z = 1.0874999999999999       ; 
    mod[15].z = 1.1625000000000001      ; 
    
    mod[0].radius =  1.0;
    mod[1].radius =  1.0;
    mod[2].radius =  1.0;
    mod[3].radius =  1.0;
    mod[4].radius =  1.0;
    mod[5].radius =  1.0;
    mod[6].radius =  1.0;
    mod[7].radius =  1.0;
    mod[8].radius =  1.0;
    mod[9].radius =  1.0;
    mod[10].radius = 1.0;
    mod[11].radius = 1.0;
    mod[12].radius = 1.0;
    mod[13].radius = 1.0;
    mod[14].radius = 1.0;
    mod[15].radius = 1.0;
 

    mod[0].intensity =  19429.725072258294 ; 
    mod[1].intensity =  12668.009585292770 ; 
    mod[2].intensity =  8340.0027795255373 ; 
    mod[3].intensity =  5704.5551581289319 ; 
    mod[4].intensity =  4688.9805968239034 ; 
    mod[5].intensity =  3690.1799763785261 ; 
    mod[6].intensity =  2307.8891400010884 ; 
    mod[7].intensity =  1769.6790689883874 ; 
    mod[8].intensity =  1851.0528740522714 ; 
    mod[9].intensity =  1696.0334024130152 ; 
    mod[10].intensity = 1483.6373298646195 ; 
    mod[11].intensity = 1285.1745205939187 ; 
    mod[12].intensity = 1076.4726184075657 ; 
    mod[13].intensity = 1104.6618516034471 ; 
    mod[14].intensity = 935.87244297902100 ; 
    mod[15].intensity = 409.65714396855253 ; 
                                           
                                           
/* 52-145 MHz spectral index */            
    mod[0].alpha =  -0.55;                 
    mod[1].alpha =  -0.55;
    mod[2].alpha =  -0.55;
    mod[3].alpha =  -0.55;
    mod[4].alpha =  -0.55;
    mod[5].alpha =  -0.55;
    mod[6].alpha =  -0.55;
    mod[7].alpha =  -0.55;
    mod[8].alpha =  -0.55;
    mod[9].alpha =  -0.55;
    mod[10].alpha = -0.55;
    mod[11].alpha = -0.55;
    mod[12].alpha = -0.55;
    mod[13].alpha = -0.55;
    mod[14].alpha = -0.55;
    mod[15].alpha = -0.55;
    counter_max = 15;

    return (counter_max);
    
}

/***********************************************************************/

double radius (double z)

{
    
    double z_diff, z_min, radius_jet;
    int counter, counter_min;

    z_min = 100.0;
    for (counter = 0; counter <= number_of_data_points; counter++)
    {
        z_diff = fabs(z - mod[counter].z);
       
        if (z_diff < z_min)
        {
            z_min = z_diff;
            counter_min = counter;
        }
    }

    if ( (counter_min > 0) && (counter_min < number_of_data_points) )
    {
        if ( z >= mod[counter_min].z )
            radius_jet = mod[counter_min].radius + (mod[counter_min+1].radius - mod[counter_min].radius) * ( z  - mod[counter_min].z ) / (mod[counter_min+1].z - mod[counter_min].z);
        else
            radius_jet = mod[counter_min].radius + (mod[counter_min].radius - mod[counter_min-1].radius) * ( z - mod[counter_min].z ) / ( mod[counter_min].z - mod[counter_min-1].z );
    }

    else if (counter_min == 0)
    {
        if ( z >= mod[counter_min].z )
            radius_jet = mod[counter_min].radius + (mod[counter_min+1].radius - mod[counter_min].radius) * ( z  - mod[counter_min].z ) / (mod[counter_min+1].z - mod[counter_min].z);
        else
            radius_jet = mod[counter_min].radius - (mod[counter_min+1].radius - mod[counter_min].radius) * ( z  - mod[counter_min].z ) / (mod[counter_min+1].z - mod[counter_min].z);
    }

    else if (counter_min == number_of_data_points)
    {
        if ( z >= mod[counter_min].z )
            radius_jet = mod[counter_min].radius - (mod[counter_min].radius - mod[counter_min-1].radius) * ( z - mod[counter_min].z ) / ( mod[counter_min].z - mod[counter_min-1].z );
        else
            radius_jet = mod[counter_min].radius + (mod[counter_min].radius - mod[counter_min-1].radius) * ( z - mod[counter_min].z ) / ( mod[counter_min].z - mod[counter_min-1].z );
    }
        
    else
    {
        printf("Error in function 'radius'. Stop.\n");
        exit (0);
    }
    
    /* if ( i == 0 && j == 0) */
    /*     printf("z=%g z_mod=%g, counter_min = %i  radius =%g radius-1=%g radius+1=%g interp=%g\n", z, mod[counter_min].z, counter_min, mod[counter_min].radius, mod[counter_min-1].radius, mod[counter_min+1].radius, radius_jet); */


/* Return jet radius in cm */
    return (kpc * radius_jet);

}
/***********************************************************************/

double dr_dz (double z)

{

    double z_diff, z_min, radius_jet, dr_dz;
    int counter, counter_min;

    z_min = 100.0;
    for (counter = 0; counter <= number_of_data_points; counter++)
    {
        z_diff = fabs(z - mod[counter].z);
       
        if (z_diff < z_min)
        {
            z_min = z_diff;
            counter_min = counter;
        }
    }


    if ((counter_min > 0) && (counter_min < number_of_data_points))
        dr_dz = (mod[counter_min + 1].radius - mod[counter_min - 1].radius) / (mod[counter_min + 1].z - mod[counter_min - 1].z);

    if (counter_min == 0)
        dr_dz = (mod[counter_min + 1].radius - mod[counter_min].radius) / (mod[counter_min + 1].z - mod[counter_min].z);

    if (counter_min == number_of_data_points)
        dr_dz = (mod[counter_min].radius - mod[counter_min - 1].radius) / (mod[counter_min].z - mod[counter_min - 1].z);

    return (dr_dz);

}

/***********************************************************************/


double magnetic_field (double z)
{

    double z_diff, z_min, magnetic_field_strength, intensity, intensity2, alpha;
    int counter, counter_min;

    z_min = 100.0;
    for (counter = 1; counter <= number_of_data_points; counter++)
    {
        z_diff = fabs(z - mod[counter].z);
        
        if (z_diff < z_min)
        {
            z_min = z_diff;
            counter_min = counter;
        }
    }

    if ( z >= mod[counter_min].z )
        intensity = mod[counter_min].intensity + (mod[counter_min+1].intensity - mod[counter_min].intensity) * ( z  - mod[counter_min].z ) / (mod[counter_min+1].z - mod[counter_min].z) ;
    else
        intensity = mod[counter_min].intensity + (mod[counter_min].intensity - mod[counter_min-1].intensity) * ( z - mod[counter_min].z ) / ( mod[counter_min].z - mod[counter_min-1].z );

    if ( z >= mod[counter_min].z )
        intensity2 = mod[counter_min].intensity2 + (mod[counter_min+1].intensity2 - mod[counter_min].intensity2) * ( z  - mod[counter_min].z ) / (mod[counter_min+1].z - mod[counter_min].z) ;
    else
        intensity2 = mod[counter_min].intensity2 + (mod[counter_min].intensity2 - mod[counter_min-1].intensity2) * ( z - mod[counter_min].z ) / ( mod[counter_min].z - mod[counter_min-1].z );

    if ( z >= mod[counter_min].z )
        alpha = mod[counter_min].alpha + (mod[counter_min+1].alpha - mod[counter_min].alpha) * ( z  - mod[counter_min].z ) / (mod[counter_min+1].z - mod[counter_min].z) ;
    else
        alpha = mod[counter_min].alpha + (mod[counter_min].alpha - mod[counter_min-1].alpha) * ( z - mod[counter_min].z ) / ( mod[counter_min].z - mod[counter_min-1].z );

    if (update_model == 1)
    {
//        printf("B=%g, intensity =%g, intensity2=%g, alpha=%g\n", B_field[i], intensity, intensity2, alpha);
        
        magnetic_field_strength =  B_field[i] * pow(intensity / (factor_model * intensity2), 1./(1.-alpha));
    }
    
    else
        magnetic_field_strength =  B_field[i];
    
    
    return (magnetic_field_strength);

}

/***********************************************************************/

void set_interpolate_values (double z, int ii)

{

    double z_diff, z_min, radius_jet;
    int counter, counter_min;
    
 
    z_min = 100.0;
    for (counter = 1; counter <= number_of_data_points; counter++)
    {
        z_diff = fabs(z - mod[counter].z);
       
        if (z_diff < z_min)
        {
            z_min = z_diff;
            counter_min = counter;
        }
    }


    if ( ( fabs(cr[ii][0].z / kpc -  mod[counter_min].z) <= fabs(cr[ii-1][0].z / kpc -  mod[counter_min].z) ) && ( fabs(cr[ii][0].z / kpc -  mod[counter_min].z) <= fabs(cr[ii+1][0].z / kpc -  mod[counter_min].z) ) )
        mod[counter_min].ii = ii;

    /* if (ii==35) */
    /*     printf("ii=%i, z=%g, diff1=%g, diff2=%g, diff3=%g\n", mod[counter_min].ii, mod[counter_min].z, fabs(cr[ii][0].z / kpc -  mod[counter_min].z), fabs(cr[ii-1][0].z / kpc -  mod[counter_min].z) , fabs(cr[ii+1][0].z / kpc -  mod[counter_min].z) );   */
    

}

/***********************************************************************/

double interpolated_value (double value, double value_low, double value_high, int ii, int ii_mod)
{

    double interp;
    
    if ( mod[ii].z >= cr[ii_mod][0].z / kpc )
        interp = value + (value_high  - value) * ( mod[ii].z  - cr[ii_mod][0].z / kpc ) / (cr[ii_mod+1][0].z / kpc - cr[ii_mod][0].z / kpc) ;
    else
        interp = value + (value  - value_low) * ( mod[ii].z  - cr[ii_mod][0].z / kpc ) / (cr[ii_mod][0].z / kpc - cr[ii_mod-1][0].z / kpc );


    /* if (ii_mod == 39) */
    /*     printf("z=%g value=%g, value_low = %g, value_high = %g, interp=%g\n", mod[ii].z, value, value_low, value_high, interp) */;
    
    return (interp);
    
}

/***********************************************************************/

void read_intensity_model (void)
{
  FILE *myfile;
  double myvariable;
  int ii;
  int kk;
  int height, width;
  double value[100];

  myfile=fopen("int_interp2.dat", "r");
  kk=0;

  for(ii = 0; ii <= number_of_data_points; ii++)
  {
      fscanf(myfile,"%lf",&myvariable);
      mod[kk].intensity2 = myvariable;
      kk = kk + 1;
//      printf("k=%i intensity=%g\n", kk, mod[kk].intensity2);
  }

  fclose(myfile);
}

/***********************************************************************/

void read_magnetic_field_model (void)
{
  FILE *myfile;
  double myvariable;
  int ii;
  double value[100];

  if( access( "b2.dat", F_OK ) != -1 )
  {
    myfile=fopen("b2.dat", "r");
  }
  else
  {
      printf("File b2.dat does not exist.\n");
      printf("Stop.\n");
      exit(0);
  }
  

  for(ii = 0; ii <= grid_size+1; ii++)
  {
      fscanf(myfile,"%lf",&myvariable);
      B_field[ii] = myvariable;
//      printf("ii=%li B=%g\n", ii, B_field[ii]);
  }

  fclose(myfile);
}
/*EOF************************************************************************/
