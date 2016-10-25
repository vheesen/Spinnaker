#include "dec.h"
#include "cosmic_rays.h"
#include "runge_kutta.h"
#include "read_parameters.h"
#include "synchrotron.h"
#include "magnetic_field.h"
#include "adiabatic.h"

/***********************************************************************/


double radius_north (double z)

{

    double z_diff, z_min, radius_jet;
    long counter, counter_min;

    mod[1].z =  103.00; 
    mod[2].z =  120.19; 
    mod[3].z =  134.97; 
    mod[4].z =  149.07; 
    mod[5].z =  163.17; 
    mod[6].z =  179.32; 
    mod[7].z =  195.14; 
    mod[8].z =  208.20; 
    mod[9].z =  220.58; 
    mod[10].z = 232.61; 
    mod[11].z = 245.68; 
    mod[12].z = 260.12; 
    mod[13].z = 272.15; 
    mod[14].z = 286.25; 
    mod[15].z = 300.34; 
    mod[16].z = 315.47; 
    mod[17].z = 329.56; 
    mod[18].z = 342.63; 
    mod[19].z = 353.63; 
    mod[20].z = 367.73; 
    mod[21].z = 382.51; 
    mod[22].z = 395.92; 
    mod[23].z = 408.29; 
    mod[24].z = 420.67; 
    mod[25].z = 433.74; 
    mod[26].z = 445.08; 
    mod[27].z = 455.74; 
    mod[28].z = 466.74; 
    mod[29].z = 470.18; 
    mod[30].z = 488.40; 
    mod[31].z = 499.75; 
    mod[32].z = 512.81; 
    mod[33].z = 531.03; 
    mod[34].z = 551.66; 
    mod[35].z = 570.91; 
    mod[36].z = 592.57; 
 
    mod[1].radius =  51.23;
    mod[2].radius =  49.51;
    mod[3].radius =  47.10;
    mod[4].radius =  46.76;
    mod[5].radius =  56.04;
    mod[6].radius =  74.60;
    mod[7].radius =  77.36;
    mod[8].radius =  77.70;
    mod[9].radius =  82.86;
    mod[10].radius = 85.61;
    mod[11].radius = 89.04;
    mod[12].radius = 92.14;
    mod[13].radius = 89.39;
    mod[14].radius = 94.20;
    mod[15].radius = 98.67;
    mod[16].radius = 97.64;
    mod[17].radius = 94.55; 
    mod[18].radius = 96.61; 
    mod[19].radius = 92.48; 
    mod[20].radius = 92.83; 
    mod[21].radius = 92.48; 
    mod[22].radius = 89.04; 
    mod[23].radius = 82.51; 
    mod[24].radius = 74.60; 
    mod[25].radius = 72.89; 
    mod[26].radius = 67.38; 
    mod[27].radius = 63.60; 
    mod[28].radius = 56.73; 
    mod[29].radius = 55.01; 
    mod[30].radius = 55.70; 
    mod[31].radius = 57.41; 
    mod[32].radius = 58.45; 
    mod[33].radius = 59.82;
    mod[34].radius = 58.10;
    mod[35].radius = 72.20;
    mod[36].radius = 51.57;

    z_min = 100.0;
    for (counter = 1; counter <= 36; counter++)
    {
        z_diff = fabs(z - mod[counter].z);
        
       
        if (z_diff < z_min)
        {
            z_min = z_diff;
            counter_min = counter;
            
        }
    }

    radius_jet = mod[counter_min].radius;


 /* Return jet radius in cm */
    return (3.09e21* radius_jet / 2.);


}
/***********************************************************************/


double dr_dz_north (double z)

{

    double z_diff, z_min, dr_dz;
    long counter, counter_min;

    mod[1].z =  103.00; 
    mod[2].z =  120.19; 
    mod[3].z =  134.97; 
    mod[4].z =  149.07; 
    mod[5].z =  163.17; 
    mod[6].z =  179.32; 
    mod[7].z =  195.14; 
    mod[8].z =  208.20; 
    mod[9].z =  220.58; 
    mod[10].z = 232.61; 
    mod[11].z = 245.68; 
    mod[12].z = 260.12; 
    mod[13].z = 272.15; 
    mod[14].z = 286.25; 
    mod[15].z = 300.34; 
    mod[16].z = 315.47; 
    mod[17].z = 329.56; 
    mod[18].z = 342.63; 
    mod[19].z = 353.63; 
    mod[20].z = 367.73; 
    mod[21].z = 382.51; 
    mod[22].z = 395.92; 
    mod[23].z = 408.29; 
    mod[24].z = 420.67; 
    mod[25].z = 433.74; 
    mod[26].z = 445.08; 
    mod[27].z = 455.74; 
    mod[28].z = 466.74; 
    mod[29].z = 470.18; 
    mod[30].z = 488.40; 
    mod[31].z = 499.75; 
    mod[32].z = 512.81; 
    mod[33].z = 531.03; 
    mod[34].z = 551.66; 
    mod[35].z = 570.91; 
    mod[36].z = 592.57; 
 
    mod[1].radius =  51.23;
    mod[2].radius =  49.51;
    mod[3].radius =  47.10;
    mod[4].radius =  46.76;
    mod[5].radius =  56.04;
    mod[6].radius =  74.60;
    mod[7].radius =  77.36;
    mod[8].radius =  77.70;
    mod[9].radius =  82.86;
    mod[10].radius = 85.61;
    mod[11].radius = 89.04;
    mod[12].radius = 92.14;
    mod[13].radius = 89.39;
    mod[14].radius = 94.20;
    mod[15].radius = 98.67;
    mod[16].radius = 97.64;
    mod[17].radius = 94.55; 
    mod[18].radius = 96.61; 
    mod[19].radius = 92.48; 
    mod[20].radius = 92.83; 
    mod[21].radius = 92.48; 
    mod[22].radius = 89.04; 
    mod[23].radius = 82.51; 
    mod[24].radius = 74.60; 
    mod[25].radius = 72.89; 
    mod[26].radius = 67.38; 
    mod[27].radius = 63.60; 
    mod[28].radius = 56.73; 
    mod[29].radius = 55.01; 
    mod[30].radius = 55.70; 
    mod[31].radius = 57.41; 
    mod[32].radius = 58.45; 
    mod[33].radius = 59.82;
    mod[34].radius = 58.10;
    mod[35].radius = 72.20;
    mod[36].radius = 51.57;

    z_min = 100.0;
    for (counter = 1; counter <= 36; counter++)
    {
        z_diff = fabs(z - mod[counter].z);
        
       
        if (z_diff < z_min)
        {
            z_min = z_diff;
            counter_min = counter;
            
        }
    }


    if ((counter_min > 1) && (counter_min < 36))
        dr_dz = 1. / 2. * (mod[counter_min + 1].radius - mod[counter_min - 1].radius) / (mod[counter_min + 1].z - mod[counter_min - 1].z);

    if (counter_min == 1)
        dr_dz = 1. / 2. * (mod[counter_min + 1].radius - mod[counter_min].radius) / (mod[counter_min + 1].z - mod[counter_min].z);

    if (counter_min == 36)
        dr_dz = 1. / 2. * (mod[counter_min].radius - mod[counter_min - 1].radius) / (mod[counter_min].z - mod[counter_min - 1].z);
    
        
/* Return jet radius in cm */
    return (dr_dz);


}
/***********************************************************************/


double radius_south (double z)

{

    double z_diff, z_min, radius_jet;
    long counter, counter_min;

    mod[1].z =  116.00; 
    mod[2].z =  123.56; 
    mod[3].z =  131.81; 
    mod[4].z =  141.10; 
    mod[5].z =  152.10; 
    mod[6].z =  163.10; 
    mod[7].z =  173.76; 
    mod[8].z =  186.14; 
    mod[9].z =  197.82; 
    mod[10].z = 210.20; 
    mod[11].z = 223.61; 
    mod[12].z = 237.71; 
    mod[13].z = 254.55; 
    mod[14].z = 272.09; 
    mod[15].z = 290.31; 
    mod[16].z = 309.56; 
    mod[17].z = 324.34; 
    mod[18].z = 352.53; 
    mod[19].z = 364.22; 
    mod[20].z = 382.79; 
    mod[21].z = 399.64; 
    mod[22].z = 416.14; 
    mod[23].z = 434.36; 
    mod[24].z = 448.45; 
    mod[25].z = 462.55; 
    mod[26].z = 478.02; 
    mod[27].z = 494.52; 
    mod[28].z = 509.99; 
    mod[29].z = 524.43; 
    mod[30].z = 540.94; 
    mod[31].z = 562.25; 
    mod[32].z = 583.22; 
    mod[33].z = 604.20; 
    mod[34].z = 626.20; 
    mod[35].z = 648.55; 
    mod[36].z = 669.17; 
    mod[37].z = 691.52; 
    mod[38].z = 711.81; 
    mod[39].z = 729.34; 
    mod[40].z = 743.78; 
    mod[41].z = 758.22; 
    mod[42].z = 772.66; 
    mod[43].z = 794.66; 

    
    mod[1].radius =  42.98;
    mod[2].radius =  40.22;
    mod[3].radius =  34.38;
    mod[4].radius =  33.35;
    mod[5].radius =  32.66;
    mod[6].radius =  39.88;
    mod[7].radius =  45.73;
    mod[8].radius =  62.57;
    mod[9].radius =  84.57;
    mod[10].radius = 110.36;
    mod[11].radius = 114.14;
    mod[12].radius = 105.89;
    mod[13].radius = 95.92;
    mod[14].radius = 93.86;
    mod[15].radius = 96.61;
    mod[16].radius = 107.95;
    mod[17].radius = 109.33; 
    mod[18].radius = 132.02; 
    mod[19].radius = 129.96; 
    mod[20].radius = 85.26; 
    mod[21].radius = 67.38; 
    mod[22].radius = 67.38; 
    mod[23].radius = 76.32; 
    mod[24].radius = 89.73; 
    mod[25].radius = 92.48; 
    mod[26].radius = 92.48; 
    mod[27].radius = 83.54; 
    mod[28].radius = 81.48; 
    mod[29].radius = 84.92; 
    mod[30].radius = 84.92; 
    mod[31].radius = 75.98; 
    mod[32].radius = 69.10; 
    mod[33].radius = 62.23;
    mod[34].radius = 66.70;
    mod[35].radius = 79.07;
    mod[36].radius = 69.79;
    mod[37].radius = 72.54;
    mod[38].radius = 67.04;
    mod[39].radius = 62.23;
    mod[40].radius = 62.92;
    mod[41].radius = 70.48;
    mod[42].radius = 74.95;
    mod[43].radius = 47.10;
 
    z_min = 100.0;
    for (counter = 1; counter <= 36; counter++)
    {
        z_diff = fabs(z - mod[counter].z);
        
       
        if (z_diff < z_min)
        {
            z_min = z_diff;
            counter_min = counter;
            
        }
    }

    radius_jet = mod[counter_min].radius;


 /* Return jet radius in cm */
    return (3.09e21* radius_jet / 2.);


}
/***********************************************************************/


double dr_dz_south (double z)

{

    double z_diff, z_min, dr_dz;
    long counter, counter_min;

    mod[1].z =  116.00; 
    mod[2].z =  123.56; 
    mod[3].z =  131.81; 
    mod[4].z =  141.10; 
    mod[5].z =  152.10; 
    mod[6].z =  163.10; 
    mod[7].z =  173.76; 
    mod[8].z =  186.14; 
    mod[9].z =  197.82; 
    mod[10].z = 210.20; 
    mod[11].z = 223.61; 
    mod[12].z = 237.71; 
    mod[13].z = 254.55; 
    mod[14].z = 272.09; 
    mod[15].z = 290.31; 
    mod[16].z = 309.56; 
    mod[17].z = 324.34; 
    mod[18].z = 352.53; 
    mod[19].z = 364.22; 
    mod[20].z = 382.79; 
    mod[21].z = 399.64; 
    mod[22].z = 416.14; 
    mod[23].z = 434.36; 
    mod[24].z = 448.45; 
    mod[25].z = 462.55; 
    mod[26].z = 478.02; 
    mod[27].z = 494.52; 
    mod[28].z = 509.99; 
    mod[29].z = 524.43; 
    mod[30].z = 540.94; 
    mod[31].z = 562.25; 
    mod[32].z = 583.22; 
    mod[33].z = 604.20; 
    mod[34].z = 626.20; 
    mod[35].z = 648.55; 
    mod[36].z = 669.17; 
    mod[37].z = 691.52; 
    mod[38].z = 711.81; 
    mod[39].z = 729.34; 
    mod[40].z = 743.78; 
    mod[41].z = 758.22; 
    mod[42].z = 772.66; 
    mod[43].z = 794.66; 

    
    mod[1].radius =  42.98;
    mod[2].radius =  40.22;
    mod[3].radius =  34.38;
    mod[4].radius =  33.35;
    mod[5].radius =  32.66;
    mod[6].radius =  39.88;
    mod[7].radius =  45.73;
    mod[8].radius =  62.57;
    mod[9].radius =  84.57;
    mod[10].radius = 110.36;
    mod[11].radius = 114.14;
    mod[12].radius = 105.89;
    mod[13].radius = 95.92;
    mod[14].radius = 93.86;
    mod[15].radius = 96.61;
    mod[16].radius = 107.95;
    mod[17].radius = 109.33; 
    mod[18].radius = 132.02; 
    mod[19].radius = 129.96; 
    mod[20].radius = 85.26; 
    mod[21].radius = 67.38; 
    mod[22].radius = 67.38; 
    mod[23].radius = 76.32; 
    mod[24].radius = 89.73; 
    mod[25].radius = 92.48; 
    mod[26].radius = 92.48; 
    mod[27].radius = 83.54; 
    mod[28].radius = 81.48; 
    mod[29].radius = 84.92; 
    mod[30].radius = 84.92; 
    mod[31].radius = 75.98; 
    mod[32].radius = 69.10; 
    mod[33].radius = 62.23;
    mod[34].radius = 66.70;
    mod[35].radius = 79.07;
    mod[36].radius = 69.79;
    mod[37].radius = 72.54;
    mod[38].radius = 67.04;
    mod[39].radius = 62.23;
    mod[40].radius = 62.92;
    mod[41].radius = 70.48;
    mod[42].radius = 74.95;
    mod[43].radius = 47.10;

    z_min = 100.0;
    for (counter = 1; counter <= 36; counter++)
    {
        z_diff = fabs(z - mod[counter].z);
        
       
        if (z_diff < z_min)
        {
            z_min = z_diff;
            counter_min = counter;
            
        }
    }


    if ((counter_min > 1) && (counter_min < 36))
        dr_dz = 1. / 2. * (mod[counter_min + 1].radius - mod[counter_min - 1].radius) / (mod[counter_min + 1].z - mod[counter_min - 1].z);

    if (counter_min == 1)
        dr_dz = 1. / 2. * (mod[counter_min + 1].radius - mod[counter_min].radius) / (mod[counter_min + 1].z - mod[counter_min].z);

    if (counter_min == 36)
        dr_dz = 1. / 2. * (mod[counter_min].radius - mod[counter_min - 1].radius) / (mod[counter_min].z - mod[counter_min - 1].z);
    
        
/* Return jet radius in cm */
    return (dr_dz);


}

/***********************************************************************/


void set_interpolate_values (double z, int ii)

{

    double z_diff, z_min, radius_jet;
    int counter, counter_min;
    

    mod[1].z =  103.00; 
    mod[2].z =  120.19; 
    mod[3].z =  134.97; 
    mod[4].z =  149.07; 
    mod[5].z =  163.17; 
    mod[6].z =  179.32; 
    mod[7].z =  195.14; 
    mod[8].z =  208.20; 
    mod[9].z =  220.58; 
    mod[10].z = 232.61; 
    mod[11].z = 245.68; 
    mod[12].z = 260.12; 
    mod[13].z = 272.15; 
    mod[14].z = 286.25; 
    mod[15].z = 300.34; 
    mod[16].z = 315.47; 
    mod[17].z = 329.56; 
    mod[18].z = 342.63; 
    mod[19].z = 353.63; 
    mod[20].z = 367.73; 
    mod[21].z = 382.51; 
    mod[22].z = 395.92; 
    mod[23].z = 408.29; 
    mod[24].z = 420.67; 
    mod[25].z = 433.74; 
    mod[26].z = 445.08; 
    mod[27].z = 455.74; 
    mod[28].z = 466.74; 
    mod[29].z = 470.18; 
    mod[30].z = 488.40; 
    mod[31].z = 499.75; 
    mod[32].z = 512.81; 
    mod[33].z = 531.03; 
    mod[34].z = 551.66; 
    mod[35].z = 570.91; 
    mod[36].z = 592.57; 
 
    z_min = 100.0;
    for (counter = 1; counter <= 36; counter++)
    {
        z_diff = fabs(z - mod[counter].z);
        
       
        if (z_diff < z_min)
        {
            z_min = z_diff;
            counter_min = counter;
            
        }
    }


    if ( ( fabs(cr[ii][1].z / parsec / 1.e3 -  mod[counter_min].z) <= fabs(cr[ii-1][1].z / parsec / 1.e3 -  mod[counter_min].z) ) && ( fabs(cr[ii][1].z / parsec / 1.e3 -  mod[counter_min].z) <= fabs(cr[ii+1][1].z / parsec / 1.e3 -  mod[counter_min].z) ) )
        mod[counter_min].ii = ii;

/*    if (ii==50)
        printf("ii=%i, z=%g, diff1=%g, diff2=%g, diff3=%g\n", mod[counter_min].ii, mod[counter_min].z, fabs(cr[ii][1].z / parsec / 1.e3 -  mod[counter_min].z), fabs(cr[ii-1][1].z / parsec / 1.e3 -  mod[counter_min].z) , fabs(cr[ii+1][1].z / parsec / 1.e3 -  mod[counter_min].z) );
*/  
    

    

}

double interpolated_value (double value, double value_low, double value_high, int ii, int ii_mod)
{

    double interp;
    

    if ( mod[ii].z >= cr[ii_mod][1].z / parsec / 1.e3 )
        interp = value + (value_high  - value) * ( mod[ii].z  - cr[ii_mod][1].z / parsec / 1.e3 ) / (cr[ii_mod+1][1].z / parsec / 1.e3 - cr[ii_mod][1].z / parsec / 1.e3) ;
    else
        interp = value + (value  - value_low) * ( mod[ii].z  - cr[ii_mod][1].z / parsec / 1.e3 ) / (cr[ii_mod][1].z / parsec / 1.e3 - cr[ii_mod-1][1].z / parsec / 1.e3 );

//    printf("value=%g, value_low = %g, value_high = %g, interp=%g\n", value, value_low, value_high, interp);
    
       

    return (interp);
    

    
}




