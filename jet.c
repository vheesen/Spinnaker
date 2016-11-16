#include "dec.h"
#include "cosmic_rays.h"
#include "runge_kutta.h"
#include "read_parameters.h"
#include "synchrotron.h"
#include "jet.h"

/***********************************************************************/


void set_radius_north ( void )

{

    mod[0].z =  0.00  ;
    mod[1].z =  2.06  ; 
    mod[2].z =  4.16  ; 
    mod[3].z =  6.26  ; 
    mod[4].z =  8.63  ; 
    mod[5].z =  10.62 ; 
    mod[6].z =  12.48 ; 
    mod[7].z =  14.68 ; 
    mod[8].z =  16.64 ; 
    mod[9].z =  18.63 ; 
    mod[10].z = 20.77 ; 
    mod[11].z = 23.00 ; 
    mod[12].z = 24.99 ; 
    mod[13].z = 27.02 ; 
    mod[14].z = 29.29 ; 
    mod[15].z = 30.94 ; 
    mod[16].z = 32.52 ; 
    mod[17].z = 33.86 ; 
    mod[18].z = 35.93 ; 
    mod[19].z = 37.89 ; 
    mod[20].z = 39.43 ; 
    mod[21].z = 41.50 ; 
    mod[22].z = 43.56 ; 
    mod[23].z = 45.76 ; 
    mod[24].z = 48.54 ; 
    mod[25].z = 52.12 ; 
    mod[26].z = 55.25 ; 
    mod[27].z = 58.14 ; 
    mod[28].z = 61.44 ; 
    mod[29].z = 65.56 ; 
    mod[30].z = 71.37 ; 
    mod[31].z = 76.43 ; 
    mod[32].z = 80.66 ; 
    mod[33].z = 85.19 ; 
    mod[34].z = 88.80 ; 
    mod[35].z = 92.41 ; 
    mod[36].z = 96.37 ; 
    mod[37].z = 100.94; 
    mod[38].z = 105.24; 
    mod[39].z = 113.28; 
    mod[40].z = 119.09; 
    mod[41].z = 127.27; 
    mod[42].z = 134.39; 
    mod[43].z = 141.92; 
    mod[44].z = 151.10; 
    mod[45].z = 163.58; 
    mod[46].z = 172.86; 
    mod[47].z = 180.19; 
    mod[48].z = 187.75; 
    mod[49].z = 195.66; 
    mod[50].z = 204.56; 
    mod[51].z = 216.63; 
    mod[52].z = 228.70; 
    mod[53].z = 240.14; 
    mod[54].z = 257.44; 
    mod[55].z = 272.26; 
    mod[56].z = 284.87; 
    mod[57].z = 299.62; 
    mod[58].z = 311.72; 
    mod[59].z = 324.44; 
    mod[60].z = 339.40; 
    mod[61].z = 355.70; 
    mod[62].z = 368.66; 
    mod[63].z = 385.43; 
    mod[64].z = 399.19; 
    mod[65].z = 413.08; 
    mod[66].z = 432.53; 
    mod[67].z = 451.20; 
    mod[68].z = 467.60; 
    mod[69].z = 486.65; 
    mod[70].z = 505.21; 
    mod[71].z = 523.44; 
    mod[72].z = 541.24; 
    mod[73].z = 561.29; 
    mod[74].z = 580.75; 
    mod[75].z = 614.89; 


//Judith's radius values    
    mod[0].radius =  2.9;
    mod[1].radius =  2.9;
    mod[2].radius =  2.9;
    mod[3].radius =  2.9;
    mod[4].radius =  2.9;
    mod[5].radius =  2.9;
    mod[6].radius =  2.9;
    mod[7].radius =  2.9;
    mod[8].radius =  3.4;
    mod[9].radius =  3.9;
    mod[10].radius = 4.4;
    mod[11].radius = 4.8;
    mod[12].radius = 4.9;
    mod[13].radius = 4.9;
    mod[14].radius = 4.9;
    mod[15].radius = 4.8;
    mod[16].radius = 4.8;
    mod[17].radius = 4.9;
    mod[18].radius = 5.0;
    mod[19].radius = 5.1;
    mod[20].radius = 5.1;
    mod[21].radius = 5.0;
    mod[22].radius = 5.0;
    mod[23].radius = 4.9;
    mod[24].radius = 4.9;
    mod[25].radius = 5.2;
    

//My radius values        
    mod[26].radius = 4.97 ;
    mod[27].radius = 4.47 ; 
    mod[28].radius = 7.10 ; 
    mod[29].radius = 8.82 ; 
    mod[30].radius = 11.48; 
    mod[31].radius = 12.81; 
    mod[32].radius = 12.89; 
    mod[33].radius = 13.75;
    mod[34].radius = 14.96;
    mod[35].radius = 15.59;
    mod[36].radius = 16.00;
    mod[37].radius = 15.30; 
    mod[38].radius = 15.13; 
    mod[39].radius = 23.21; 
    mod[40].radius = 22.86; 
    mod[41].radius = 22.86; 
    mod[42].radius = 22.52; 
    mod[43].radius = 22.35;
    mod[44].radius = 22.35; 
    mod[45].radius = 22.35; 
    mod[46].radius = 25.27; 
    mod[47].radius = 25.79; 
    mod[48].radius = 26.64; 
    mod[49].radius = 26.82; 
    mod[50].radius = 25.79; 
    mod[51].radius = 30.94; 
    mod[52].radius = 43.32; 
    mod[53].radius = 48.82; 
    mod[54].radius = 52.77; 
    mod[55].radius = 52.09; 
    mod[56].radius = 50.71; 
    mod[57].radius = 50.02; 
    mod[58].radius = 49.85; 
    mod[59].radius = 47.96; 
    mod[60].radius = 49.34; 
    mod[61].radius = 51.57; 
    mod[62].radius = 53.98; 
    mod[63].radius = 55.52; 
    mod[64].radius = 57.24; 
    mod[65].radius = 59.31; 
    mod[66].radius = 60.34; 
    mod[67].radius = 56.90; 
    mod[68].radius = 52.43; 
    mod[69].radius = 51.05; 
    mod[70].radius = 50.02; 
    mod[71].radius = 49.34; 
    mod[72].radius = 54.66; 
    mod[73].radius = 56.90; 
    mod[74].radius = 59.82; 
    mod[75].radius = 43.49; 


    mod[0].intensity =  2.09E+00;
    mod[1].intensity =  2.31E+00;
    mod[2].intensity =  2.42E+00;
    mod[3].intensity =  2.55E+00;
    mod[4].intensity =  2.62E+00;
    mod[5].intensity =  2.68E+00;
    mod[6].intensity =  2.56E+00;
    mod[7].intensity =  2.52E+00;
    mod[8].intensity =  2.32E+00;
    mod[9].intensity =  2.12E+00;
    mod[10].intensity = 2.10E+00;
    mod[11].intensity = 1.97E+00;
    mod[12].intensity = 1.81E+00;
    mod[13].intensity = 1.67E+00;
    mod[14].intensity = 1.48E+00;
    mod[15].intensity = 1.35E+00;
    mod[16].intensity = 1.28E+00;
    mod[17].intensity = 1.15E+00;
    mod[18].intensity = 1.02E+00;
    mod[19].intensity = 9.86E-01;
    mod[20].intensity = 8.61E-01;
    mod[21].intensity = 7.95E-01;
    mod[22].intensity = 7.50E-01;
    mod[23].intensity = 6.76E-01;
    mod[24].intensity = 7.52E-01;
    mod[25].intensity = 7.78E-01;
    mod[26].intensity = 7.17E-01;
    mod[27].intensity = 6.53E-01; 
    mod[28].intensity = 5.38E-01; 
    mod[29].intensity = 5.13E-01; 
    mod[30].intensity = 5.37E-01; 
    mod[31].intensity = 5.87E-01; 
    mod[32].intensity = 6.09E-01; 
    mod[33].intensity = 6.18E-01;
    mod[34].intensity = 6.10E-01;
    mod[35].intensity = 5.20E-01;
    mod[36].intensity = 5.10E-01;
    mod[37].intensity = 4.90E-01; 
    mod[38].intensity = 4.42E-01; 
    mod[39].intensity = 3.43E-01; 
    mod[40].intensity = 2.82E-01; 
    mod[41].intensity = 2.33E-01; 
    mod[42].intensity = 2.02E-01; 
    mod[43].intensity = 1.89E-01;
    mod[44].intensity = 1.88E-01; 
    mod[45].intensity = 1.69E-01; 
    mod[46].intensity = 1.30E-01; 
    mod[47].intensity = 1.24E-01; 
    mod[48].intensity = 1.12E-01; 
    mod[49].intensity = 1.03E-01; 
    mod[50].intensity = 9.91E-02; 
    mod[51].intensity = 7.25E-02; 
    mod[52].intensity = 6.52E-02; 
    mod[53].intensity = 5.90E-02; 
    mod[54].intensity = 5.15E-02; 
    mod[55].intensity = 4.91E-02; 
    mod[56].intensity = 4.59E-02; 
    mod[57].intensity = 4.20E-02; 
    mod[58].intensity = 3.70E-02; 
    mod[59].intensity = 3.47E-02; 
    mod[60].intensity = 3.28E-02; 
    mod[61].intensity = 3.00E-02; 
    mod[62].intensity = 2.66E-02; 
    mod[63].intensity = 2.53E-02; 
    mod[64].intensity = 2.20E-02; 
    mod[65].intensity = 1.86E-02; 
    mod[66].intensity = 1.62E-02; 
    mod[67].intensity = 1.33E-02; 
    mod[68].intensity = 1.14E-02; 
    mod[69].intensity = 1.09E-02; 
    mod[70].intensity = 7.64E-03; 
    mod[71].intensity = 5.91E-03; 
    mod[72].intensity = 5.62E-03; 
    mod[73].intensity = 5.32E-03; 
    mod[74].intensity = 5.46E-03; 
    mod[75].intensity = 5.74E-03;

    read_intensity_model ();

    /* mod[0].intensity2 =  5.676130e-01; */
    /* mod[1].intensity2 =  5.598767e-01; */
    /* mod[2].intensity2 =  5.491067e-01; */
    /* mod[3].intensity2 =  5.394363e-01; */
    /* mod[4].intensity2 =  5.277154e-01; */
    /* mod[5].intensity2 =  5.204199e-01; */
    /* mod[6].intensity2 =  5.123332e-01; */
    /* mod[7].intensity2 =  4.955614e-01; */
    /* mod[8].intensity2 =  4.549984e-01; */
    /* mod[9].intensity2 =  4.024273e-01; */
    /* mod[10].intensity2 = 3.588948e-01; */
    /* mod[11].intensity2 = 3.262482e-01; */
    /* mod[12].intensity2 = 3.098353e-01; */
    /* mod[13].intensity2 = 3.042415e-01; */
    /* mod[14].intensity2 = 2.973426e-01; */
    /* mod[15].intensity2 = 2.934466e-01; */
    /* mod[16].intensity2 = 2.889325e-01; */
    /* mod[17].intensity2 = 2.838057e-01; */
    /* mod[18].intensity2 = 2.754202e-01; */
    /* mod[19].intensity2 = 2.672270e-01; */
    /* mod[20].intensity2 = 2.621485e-01; */
    /* mod[21].intensity2 = 2.569883e-01; */
    /* mod[22].intensity2 = 2.531906e-01; */
    /* mod[23].intensity2 = 2.480645e-01; */
    /* mod[24].intensity2 = 2.411382e-01; */
    /* mod[25].intensity2 = 2.282281e-01; */
    /* mod[26].intensity2 = 2.204369e-01; */
    /* mod[27].intensity2 = 2.130653e-01;  */
    /* mod[28].intensity2 = 1.631918e-01;  */
    /* mod[29].intensity2 = 1.247250e-01;  */
    /* mod[30].intensity2 = 9.615842e-02;  */
    /* mod[31].intensity2 = 8.413710e-02;  */
    /* mod[32].intensity2 = 7.971345e-02;  */
    /* mod[33].intensity2 = 7.386794e-02; */
    /* mod[34].intensity2 = 6.728974e-02; */
    /* mod[35].intensity2 = 6.279237e-02; */
    /* mod[36].intensity2 = 5.937152e-02; */
    /* mod[37].intensity2 = 5.704494e-02;  */
    /* mod[38].intensity2 = 5.485812e-02;  */
    /* mod[39].intensity2 = 3.709192e-02;  */
    /* mod[40].intensity2 = 3.438757e-02;  */
    /* mod[41].intensity2 = 3.026353e-02;  */
    /* mod[42].intensity2 = 2.728060e-02;  */
    /* mod[43].intensity2 = 2.423384e-02; */
    /* mod[44].intensity2 = 2.108815e-02;  */
    /* mod[45].intensity2 = 1.740590e-02;  */
    /* mod[46].intensity2 = 1.380531e-02;  */
    /* mod[47].intensity2 = 1.218613e-02;  */
    /* mod[48].intensity2 = 1.126405e-02;  */
    /* mod[49].intensity2 = 1.065468e-02;  */
    /* mod[50].intensity2 = 1.003964e-02;  */
    /* mod[51].intensity2 = 7.877240e-03;  */
    /* mod[52].intensity2 = 5.160969e-03;  */
    /* mod[53].intensity2 = 4.069793e-03;  */
    /* mod[54].intensity2 = 3.283385e-03;  */
    /* mod[55].intensity2 = 2.923457e-03;  */
    /* mod[56].intensity2 = 2.632353e-03;  */
    /* mod[57].intensity2 = 2.317273e-03;  */
    /* mod[58].intensity2 = 2.085199e-03;  */
    /* mod[59].intensity2 = 1.864595e-03;  */
    /* mod[60].intensity2 = 1.630413e-03;  */
    /* mod[61].intensity2 = 1.356470e-03;  */
    /* mod[62].intensity2 = 1.157242e-03;  */
    /* mod[63].intensity2 = 9.783691e-04;  */
    /* mod[64].intensity2 = 8.396047e-04;  */
    /* mod[65].intensity2 = 6.996493e-04;  */
    /* mod[66].intensity2 = 5.793159e-04;  */
    /* mod[67].intensity2 = 4.723123e-04;  */
    /* mod[68].intensity2 = 4.037440e-04;  */
    /* mod[69].intensity2 = 3.240497e-04;  */
    /* mod[70].intensity2 = 2.648708e-04;  */
    /* mod[71].intensity2 = 2.188423e-04;  */
    /* mod[72].intensity2 = 1.434696e-04;  */
    /* mod[73].intensity2 = 1.060036e-04;  */
    /* mod[74].intensity2 = 7.477332e-05;  */
    /* mod[75].intensity2 = 4.859734e-05; */


    mod[0].alpha =  -0.55;
    mod[1].alpha =  -0.55;
    mod[2].alpha =  -0.55;
    mod[3].alpha =  -0.55;
    mod[4].alpha =  -0.55;
    mod[5].alpha =  -0.55;
    mod[6].alpha =  -0.55;
    mod[7].alpha =  -0.48;
    mod[8].alpha =  -0.55;
    mod[9].alpha =  -0.55;
    mod[10].alpha = -0.55;
    mod[11].alpha = -0.55;
    mod[12].alpha = -0.55;
    mod[13].alpha = -0.55;
    mod[14].alpha = -0.55;
    mod[15].alpha = -0.55;
    mod[16].alpha = -0.55;
    mod[17].alpha = -0.55;
    mod[18].alpha = -0.55;
    mod[19].alpha = -0.55;
    mod[20].alpha = -0.55;
    mod[21].alpha = -0.55;
    mod[22].alpha = -0.55;
    mod[23].alpha = -0.55;
    mod[24].alpha = -0.55;
    mod[25].alpha = -0.55;
    mod[26].alpha = -0.55;
    mod[27].alpha = -0.55; 
    mod[28].alpha = -0.55; 
    mod[29].alpha = -0.55; 
    mod[30].alpha = -0.55; 
    mod[31].alpha = -0.55; 
    mod[32].alpha = -0.55; 
    mod[33].alpha = -0.55;
    mod[34].alpha = -0.55;
    mod[35].alpha = -0.55;
    mod[36].alpha = -0.55;
    mod[37].alpha = -0.55; 
    mod[38].alpha = -0.55; 
    mod[39].alpha = -0.55; 
    mod[40].alpha = -0.55; 
    mod[41].alpha = -0.55; 
    mod[42].alpha = -0.55; 
    mod[43].alpha = -0.55;
    mod[44].alpha = -0.61; 
    mod[45].alpha = -0.68; 
    mod[46].alpha = -0.71; 
    mod[47].alpha = -0.73; 
    mod[48].alpha = -0.75; 
    mod[49].alpha = -0.76; 
    mod[50].alpha = -0.76; 
    mod[51].alpha = -0.82; 
    mod[52].alpha = -0.89; 
    mod[53].alpha = -0.96; 
    mod[54].alpha = -1.01; 
    mod[55].alpha = -1.07; 
    mod[56].alpha = -1.12; 
    mod[57].alpha = -1.07; 
    mod[58].alpha = -1.13; 
    mod[59].alpha = -1.16; 
    mod[60].alpha = -1.19; 
    mod[61].alpha = -1.16; 
    mod[62].alpha = -1.19; 
    mod[63].alpha = -1.27; 
    mod[64].alpha = -1.3 ; 
    mod[65].alpha = -1.4 ; 
    mod[66].alpha = -1.54; 
    mod[67].alpha = -1.95; 
    mod[68].alpha = -2.27; 
    mod[69].alpha = -2.27; 
    mod[70].alpha = -2.27; 
    mod[71].alpha = -2.27; 
    mod[72].alpha = -2.27; 
    mod[73].alpha = -2.27; 
    mod[74].alpha = -2.27; 
    mod[75].alpha = -2.27; 

}

/***********************************************************************/

double radius_north (double z)

{
    
    double z_diff, z_min, radius_jet;
    long counter, counter_min;

    z_min = 100.0;
    for (counter = 0; counter <= 75; counter++)
    {
        z_diff = fabs(z - mod[counter].z);
        
       
        if (z_diff < z_min)
        {
            z_min = z_diff;
            counter_min = counter;
            
        }
    }

//    radius_jet = mod[counter_min].radius;


    if ( z >= mod[counter_min].z )
        radius_jet = mod[counter_min].radius + (mod[counter_min+1].radius - mod[counter_min].radius) * ( z  - mod[counter_min].z ) / (mod[counter_min+1].z - mod[counter_min].z) ;
    else
        radius_jet = mod[counter_min].radius + (mod[counter_min].radius - mod[counter_min-1].radius) * ( z - mod[counter_min].z ) / ( mod[counter_min].z - mod[counter_min-1].z );

    /* if ( z < 200. ) */
    /*     radius_jet = 8.e22/3.09e21; */
    
       


   /* if ( mod[ii].z >= cr[ii_mod][1].z / parsec / 1.e3 ) */
   /*      interp = value + (value_high  - value) * ( mod[ii].z  - cr[ii_mod][1].z / parsec / 1.e3 ) / (cr[ii_mod+1][1].z / parsec / 1.e3 - cr[ii_mod][1].z / parsec / 1.e3) ; */
   /*  else */
   /*      interp = value + (value  - value_low) * ( mod[ii].z  - cr[ii_mod][1].z / parsec / 1.e3 ) / (cr[ii_mod][1].z / parsec / 1.e3 - cr[ii_mod-1][1].z / parsec / 1.e3 ); */

    

    /* if ( j == 144) */
    /*     printf("z=%g z_mod=%g, counter_min = %li  radius =%g radius-1=%g radius+1=%g interp=%g\n", z, mod[counter_min].z, counter_min, mod[counter_min].radius, mod[counter_min-1].radius, mod[counter_min+1].radius, radius_jet); */

//    radius_jet = mod[counter_min].radius;

  

//    radius_jet = pow(10.,21.9551+0.00853176*z-1.88229e-05*z*z+1.36584e-08*z*z*z);
    


 /* Return jet radius in cm */
    return (3.09e21* radius_jet);
//    return (radius_jet);


}
/***********************************************************************/


double dr_dz_north (double z)

{

    double z_diff, z_min, radius_jet, dr_dz;
    long counter, counter_min;

    z_min = 100.0;
    for (counter = 0; counter <= 75; counter++)
    {
        z_diff = fabs(z - mod[counter].z);
        
       
        if (z_diff < z_min)
        {
            z_min = z_diff;
            counter_min = counter;
            
        }
    }


    if ((counter_min > 0) && (counter_min < 75))
        dr_dz = (mod[counter_min + 1].radius - mod[counter_min - 1].radius) / (mod[counter_min + 1].z - mod[counter_min - 1].z);

    if (counter_min == 0)
        dr_dz = (mod[counter_min + 1].radius - mod[counter_min].radius) / (mod[counter_min + 1].z - mod[counter_min].z);

    if (counter_min == 75)
        dr_dz = (mod[counter_min].radius - mod[counter_min - 1].radius) / (mod[counter_min].z - mod[counter_min - 1].z);


/* Return jet radius in cm */
    return (dr_dz);


}


/***********************************************************************/


void set_radius_south ( void )

{


    mod[1].z =  15.; 
    mod[2].z =  17.; 
    mod[3].z =  19.; 
    mod[4].z =  21.; 
    mod[5].z =  23.; 
    mod[6].z =  25.; 
    mod[7].z =  27.; 
    mod[8].z =  29.; 
    mod[9].z =  31.; 
    mod[10].z = 33.; 
    mod[11].z = 35.; 
    mod[12].z = 37.; 
    mod[13].z = 39.; 
    mod[14].z = 41.; 
    mod[15].z = 43.; 
    mod[16].z = 45.; 
    mod[17].z = 47.; 
    mod[18].z = 49.; 
    mod[19].z = 51.; 
    mod[20].z = 53.; 
    mod[21].z = 55.; 
    mod[22].z = 57.; 
    mod[23].z = 59.; 
    mod[24].z = 61.; 
    mod[25].z = 63.; 
    mod[26].z = 65.; 
    mod[27].z = 67.; 
    mod[28].z = 69.; 
    mod[29].z = 71.; 
    mod[30].z = 73.; 
    mod[31].z = 75.; 
    mod[32].z = 77.; 
    mod[33].z = 79.; 
    mod[34].z = 81.; 
    mod[35].z = 83.; 
    mod[36].z = 85.; 
    mod[37].z = 87.; 
    mod[38].z = 89.; 
    mod[39].z = 91.; 
    mod[40].z = 93.; 
    mod[41].z = 95.; 
    mod[42].z = 97.; 
    mod[43].z = 99.; 
    mod[44].z = 101.; 
    mod[45].z = 103.; 
    mod[46].z = 105.; 
    mod[47].z = 107.; 
    mod[48].z = 109.; 
    mod[49].z = 111.; 
    mod[50].z = 113.; 
    mod[51].z = 115.; 

/* My values */
    mod[52].z = 116.00; 
    mod[53].z = 123.56; 
    mod[54].z = 131.81; 
    mod[55].z = 141.10; 
    mod[56].z = 152.10; 
    mod[57].z = 163.10; 
    mod[58].z = 173.76; 
    mod[59].z = 186.14; 
    mod[60].z = 197.82; 
    mod[61].z = 210.20; 
    mod[62].z = 223.61; 
    mod[63].z = 237.71; 
    mod[64].z = 254.55; 
    mod[65].z = 272.09; 
    mod[66].z = 290.31; 
    mod[67].z = 309.56; 
    mod[68].z = 324.34; 
    mod[69].z = 352.53; 
    mod[70].z = 364.22; 
    mod[71].z = 382.79; 
    mod[72].z = 399.64; 
    mod[73].z = 416.14; 
    mod[74].z = 434.36; 
    mod[75].z = 448.45; 
    mod[76].z = 462.55; 
    mod[77].z = 478.02; 
    mod[78].z = 494.52; 
    mod[79].z = 509.99; 
    mod[80].z = 524.43; 
    mod[81].z = 540.94; 
    mod[82].z = 562.25; 
    mod[83].z = 583.22; 
    mod[84].z = 604.20; 
    mod[85].z = 626.20; 
    mod[86].z = 648.55; 
    mod[87].z = 669.17; 
    mod[88].z = 691.52; 
    mod[89].z = 711.81; 
    mod[90].z = 729.34; 
    mod[91].z = 743.78; 
    mod[92].z = 758.22; 
    mod[93].z = 772.66; 
    mod[94].z = 794.66; 


    mod[1].radius =  2.46 ; 
    mod[2].radius =  2.89 ; 
    mod[3].radius =  3.52 ; 
    mod[4].radius =  4.12 ; 
    mod[5].radius =  4.63 ; 
    mod[6].radius =  4.85 ; 
    mod[7].radius =  4.88 ; 
    mod[8].radius =  4.84 ; 
    mod[9].radius =  4.84 ; 
    mod[10].radius = 4.90 ; 
    mod[11].radius = 4.99 ; 
    mod[12].radius = 5.12 ; 
    mod[13].radius = 5.24 ; 
    mod[14].radius = 5.36 ; 
    mod[15].radius = 5.44 ; 
    mod[16].radius = 5.53 ; 
    mod[17].radius = 5.64 ; 
    mod[18].radius = 5.84 ; 
    mod[19].radius = 6.11 ; 
    mod[20].radius = 6.44 ; 
    mod[21].radius = 6.79 ; 
    mod[22].radius = 7.15 ; 
    mod[23].radius = 7.47 ; 
    mod[24].radius = 7.77 ; 
    mod[25].radius = 8.02 ; 
    mod[26].radius = 8.26 ; 
    mod[27].radius = 8.52 ; 
    mod[28].radius = 8.78 ; 
    mod[29].radius = 9.10 ; 
    mod[30].radius = 9.45 ; 
    mod[31].radius = 9.88 ; 
    mod[32].radius = 10.34; 
    mod[33].radius = 10.88; 
    mod[34].radius = 11.44; 
    mod[35].radius = 12.09; 
    mod[36].radius = 12.74; 
    mod[37].radius = 13.48; 
    mod[38].radius = 14.21; 
    mod[39].radius = 15.02; 
    mod[40].radius = 15.82; 
    mod[41].radius = 16.71; 
    mod[42].radius = 17.56; 
    mod[43].radius = 18.50; 
    mod[44].radius = 19.41; 
    mod[45].radius = 20.34; 
    mod[46].radius = 21.35; 
    mod[47].radius = 22.31; 
    mod[48].radius = 23.34; 
    mod[49].radius = 24.33; 
    mod[50].radius = 25.38; 
    mod[51].radius = 26.38; 

/* My values */
    mod[52].radius = 21.49; 
    mod[53].radius = 20.11; 
    mod[54].radius = 17.19; 
    mod[55].radius = 16.67; 
    mod[56].radius = 16.33; 
    mod[57].radius = 19.94; 
    mod[58].radius = 22.86; 
    mod[59].radius = 31.29; 
    mod[60].radius = 42.29; 
    mod[61].radius = 55.18; 
    mod[62].radius = 57.07; 
    mod[63].radius = 52.95; 
    mod[64].radius = 47.96; 
    mod[65].radius = 46.93; 
    mod[66].radius = 48.30; 
    mod[67].radius = 53.98; 
    mod[68].radius = 54.66; 
    mod[69].radius = 66.01; 
    mod[70].radius = 64.98; 
    mod[71].radius = 42.63; 
    mod[72].radius = 33.69; 
    mod[73].radius = 33.69; 
    mod[74].radius = 38.16; 
    mod[75].radius = 44.87; 
    mod[76].radius = 46.24; 
    mod[77].radius = 46.24; 
    mod[78].radius = 41.77; 
    mod[79].radius = 40.74; 
    mod[80].radius = 42.46; 
    mod[81].radius = 42.46; 
    mod[82].radius = 37.99; 
    mod[83].radius = 34.55; 
    mod[84].radius = 31.11; 
    mod[85].radius = 33.35; 
    mod[86].radius = 39.54; 
    mod[87].radius = 34.90; 
    mod[88].radius = 36.27; 
    mod[89].radius = 33.52; 
    mod[90].radius = 31.11; 
    mod[91].radius = 31.46; 
    mod[92].radius = 35.24; 
    mod[93].radius = 37.47; 
    mod[94].radius = 23.55;

}


/***********************************************************************/

double radius_south (double z)

{

    double z_diff, z_min, radius_jet;
    long counter, counter_min;

    
    z_min = 100.0;
    for (counter = 1; counter <= 94; counter++)
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
    return (3.09e21* radius_jet);


}
/***********************************************************************/


double dr_dz_south (double z)

{

    double z_diff, z_min, dr_dz;
    long counter, counter_min;

    z_min = 100.0;
    for (counter = 1; counter <= 94; counter++)
    {
        z_diff = fabs(z - mod[counter].z);
        
       
        if (z_diff < z_min)
        {
            z_min = z_diff;
            counter_min = counter;
            
        }
    }


    if ((counter_min > 1) && (counter_min < 94))
        dr_dz = (mod[counter_min + 1].radius - mod[counter_min - 1].radius) / (mod[counter_min + 1].z - mod[counter_min - 1].z);

    if (counter_min == 1)
        dr_dz = (mod[counter_min + 1].radius - mod[counter_min].radius) / (mod[counter_min + 1].z - mod[counter_min].z);

    if (counter_min == 94)
        dr_dz = (mod[counter_min].radius - mod[counter_min - 1].radius) / (mod[counter_min].z - mod[counter_min - 1].z);
    
        
/* Return jet radius in cm */
    return (dr_dz);


}



/***********************************************************************/


double magnetic_field_north (double z)

{

    double z_diff, z_min, magnetic_field_tube, intensity, intensity2, alpha;
    long counter, counter_min;

    mod2[0].z =  15.00 ;
    mod2[1].z =  17.00 ; 
    mod2[2].z =  19.00 ; 
    mod2[3].z =  21.00 ; 
    mod2[4].z =  23.00 ; 
    mod2[5].z =  25.00 ; 
    mod2[6].z =  27.00 ; 
    mod2[7].z =  29.00 ; 
    mod2[8].z =  31.00 ; 
    mod2[9].z =  33.00 ; 
    mod2[10].z = 35.00 ; 
    mod2[11].z = 37.00 ; 
    mod2[12].z = 39.00 ; 
    mod2[13].z = 41.00 ; 
    mod2[14].z = 43.00 ; 
    mod2[15].z = 45.00 ; 
    mod2[16].z = 47.00 ; 
    mod2[17].z = 49.00 ; 
    mod2[18].z = 51.00 ; 
    mod2[19].z = 53.00 ; 
    mod2[20].z = 55.00 ; 
    mod2[21].z = 57.00 ; 
    mod2[22].z = 59.00 ; 
    mod2[23].z = 61.00 ; 
    mod2[24].z = 63.00 ; 
    mod2[25].z = 65.00 ; 
    mod2[26].z = 67.00 ; 
    mod2[27].z = 69.00 ; 
    mod2[28].z = 71.00 ; 
    mod2[29].z = 73.00 ; 
    mod2[30].z = 75.00 ; 
    mod2[31].z = 77.00 ; 
    mod2[32].z = 79.00 ; 
    mod2[33].z = 81.00 ; 
    mod2[34].z = 83.00 ; 
    mod2[35].z = 85.00 ; 
    mod2[36].z = 87.00 ; 
    mod2[37].z = 89.00 ; 
    mod2[38].z = 91.00 ; 
    mod2[39].z = 93.00 ; 
    mod2[40].z = 95.00 ; 
    mod2[41].z = 97.00 ; 
    mod2[42].z = 99.00 ; 
    mod2[43].z = 101.00; 
    mod2[44].z = 103.00; 
    mod2[45].z = 105.00; 
    mod2[46].z = 107.00; 
    mod2[47].z = 109.00; 
    mod2[48].z = 111.00; 
    mod2[49].z = 113.00; 
    mod2[50].z = 115.00; 
    mod2[51].z = 117.00; 
    mod2[52].z = 119.00; 
 
    mod2[0].B_field =  15.09;
    mod2[1].B_field =  14.44; 
    mod2[2].B_field =  14.02; 
    mod2[3].B_field =  13.74; 
    mod2[4].B_field =  13.58; 
    mod2[5].B_field =  13.48; 
    mod2[6].B_field =  13.35; 
    mod2[7].B_field =  13.19; 
    mod2[8].B_field =  13.03; 
    mod2[9].B_field =  12.87; 
    mod2[10].B_field = 12.71; 
    mod2[11].B_field = 12.54; 
    mod2[12].B_field = 12.38; 
    mod2[13].B_field = 12.25; 
    mod2[14].B_field = 12.14; 
    mod2[15].B_field = 12.09; 
    mod2[16].B_field = 12.08; 
    mod2[17].B_field = 12.09; 
    mod2[18].B_field = 12.10; 
    mod2[19].B_field = 12.09; 
    mod2[20].B_field = 12.06; 
    mod2[21].B_field = 11.97; 
    mod2[22].B_field = 11.84; 
    mod2[23].B_field = 11.66; 
    mod2[24].B_field = 11.45; 
    mod2[25].B_field = 11.23; 
    mod2[26].B_field = 11.02; 
    mod2[27].B_field = 10.82; 
    mod2[28].B_field = 10.64; 
    mod2[29].B_field = 10.50; 
    mod2[30].B_field = 10.37; 
    mod2[31].B_field = 10.26; 
    mod2[32].B_field = 10.17; 
    mod2[33].B_field = 10.08; 
    mod2[34].B_field = 10.00; 
    mod2[35].B_field = 9.91 ; 
    mod2[36].B_field = 9.82 ; 
    mod2[37].B_field = 9.73 ; 
    mod2[38].B_field = 9.63 ; 
    mod2[39].B_field = 9.52 ; 
    mod2[40].B_field = 9.41 ; 
    mod2[41].B_field = 9.30 ; 
    mod2[42].B_field = 9.19 ; 
    mod2[43].B_field = 9.08 ; 
    mod2[44].B_field = 8.98 ; 
    mod2[45].B_field = 8.88 ; 
    mod2[46].B_field = 8.78 ; 
    mod2[47].B_field = 8.69 ; 
    mod2[48].B_field = 8.60 ; 
    mod2[49].B_field = 8.52 ; 
    mod2[50].B_field = 8.43 ; 
    mod2[51].B_field = 8.35 ; 
    mod2[52].B_field = 8.27 ; 


    z_min = 100.0;
    for (counter = 1; counter <= 75; counter++)
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

//    intensity = mod[counter_min].intensity;

    if ( z >= mod[counter_min].z )
        intensity2 = mod[counter_min].intensity2 + (mod[counter_min+1].intensity2 - mod[counter_min].intensity2) * ( z  - mod[counter_min].z ) / (mod[counter_min+1].z - mod[counter_min].z) ;
    else
        intensity2 = mod[counter_min].intensity2 + (mod[counter_min].intensity2 - mod[counter_min-1].intensity2) * ( z - mod[counter_min].z ) / ( mod[counter_min].z - mod[counter_min-1].z );

//    intensity2 = mod[counter_min].intensity2;

    if ( z >= mod[counter_min].z )
        alpha = mod[counter_min].alpha + (mod[counter_min+1].alpha - mod[counter_min].alpha) * ( z  - mod[counter_min].z ) / (mod[counter_min+1].z - mod[counter_min].z) ;
    else
        alpha = mod[counter_min].alpha + (mod[counter_min].alpha - mod[counter_min-1].alpha) * ( z - mod[counter_min].z ) / ( mod[counter_min].z - mod[counter_min-1].z );

//    alpha = mod[counter_min].alpha;

    
/*    magnetic_field_tube =  B0 * exp(-cr[i][1].z / parsec / 1.e3 / h_B0) *
      pow(intensity / (factor_int * intensity2), 1./1.55);*/

    if (update_model == 1)
        magnetic_field_tube =  B_field[i] * pow(intensity / (factor_model * intensity2), 1./(1.-alpha));
    else
        magnetic_field_tube =  B_field[i];
    


    /* if ( z >= mod2[counter_min].z ) */
    /*     magnetic_field_tube = mod2[counter_min].B_field + (mod2[counter_min+1].B_field - mod2[counter_min].B_field) * ( z  - mod2[counter_min].z ) / (mod2[counter_min+1].z - mod2[counter_min].z) ; */
    /* else */
    /*     magnetic_field_tube = mod2[counter_min].B_field + (mod2[counter_min].B_field - mod2[counter_min-1].B_field) * ( z - mod2[counter_min].z ) / ( mod2[counter_min].z - mod2[counter_min-1].z ); */


//    magnetic_field_tube = mod2[counter_min].B_field;
    
    return (magnetic_field_tube);

    
//    return (1.e-6 * magnetic_field_tube);


}

/***********************************************************************/


double velocity_field_north (double z)

{

    double z_diff, z_min, velocity_value;
    long counter, counter_min;


    mod2[0].z =  15.00 ;
    mod2[1].z =  17.00 ; 
    mod2[2].z =  19.00 ; 
    mod2[3].z =  21.00 ; 
    mod2[4].z =  23.00 ; 
    mod2[5].z =  25.00 ; 
    mod2[6].z =  27.00 ; 
    mod2[7].z =  29.00 ; 
    mod2[8].z =  31.00 ; 
    mod2[9].z =  33.00 ; 
    mod2[10].z = 35.00 ; 
    mod2[11].z = 37.00 ; 
    mod2[12].z = 39.00 ; 
    mod2[13].z = 41.00 ; 
    mod2[14].z = 43.00 ; 
    mod2[15].z = 45.00 ; 
    mod2[16].z = 47.00 ; 
    mod2[17].z = 49.00 ; 
    mod2[18].z = 51.00 ; 
    mod2[19].z = 53.00 ; 
    mod2[20].z = 55.00 ; 
    mod2[21].z = 57.00 ; 
    mod2[22].z = 59.00 ; 
    mod2[23].z = 61.00 ; 
    mod2[24].z = 63.00 ; 
    mod2[25].z = 65.00 ; 
    mod2[26].z = 67.00 ; 
    mod2[27].z = 69.00 ; 
    mod2[28].z = 71.00 ; 
    mod2[29].z = 73.00 ; 
    mod2[30].z = 75.00 ; 
    mod2[31].z = 77.00 ; 
    mod2[32].z = 79.00 ; 
    mod2[33].z = 81.00 ; 
    mod2[34].z = 83.00 ; 
    mod2[35].z = 85.00 ; 
    mod2[36].z = 87.00 ; 
    mod2[37].z = 89.00 ; 
    mod2[38].z = 91.00 ; 
    mod2[39].z = 93.00 ; 
    mod2[40].z = 95.00 ; 
    mod2[41].z = 97.00 ; 
    mod2[42].z = 99.00 ; 
    mod2[43].z = 101.00; 
    mod2[44].z = 103.00; 
    mod2[45].z = 105.00; 
    mod2[46].z = 107.00; 
    mod2[47].z = 109.00; 
    mod2[48].z = 111.00; 
    mod2[49].z = 113.00; 
    mod2[50].z = 115.00; 
    mod2[51].z = 117.00; 
    mod2[52].z = 119.00; 
 
    mod2[0].velocity =  5.80E+09;
    mod2[1].velocity =  5.54E+09; 
    mod2[2].velocity =  5.27E+09; 
    mod2[3].velocity =  5.04E+09; 
    mod2[4].velocity =  4.90E+09; 
    mod2[5].velocity =  4.89E+09; 
    mod2[6].velocity =  4.97E+09; 
    mod2[7].velocity =  5.07E+09; 
    mod2[8].velocity =  5.16E+09; 
    mod2[9].velocity =  5.18E+09; 
    mod2[10].velocity = 5.16E+09; 
    mod2[11].velocity = 5.17E+09; 
    mod2[12].velocity = 5.18E+09; 
    mod2[13].velocity = 5.25E+09; 
    mod2[14].velocity = 5.31E+09; 
    mod2[15].velocity = 5.36E+09; 
    mod2[16].velocity = 5.36E+09; 
    mod2[17].velocity = 5.33E+09; 
    mod2[18].velocity = 5.26E+09; 
    mod2[19].velocity = 5.11E+09; 
    mod2[20].velocity = 4.88E+09; 
    mod2[21].velocity = 4.58E+09; 
    mod2[22].velocity = 4.23E+09; 
    mod2[23].velocity = 3.88E+09; 
    mod2[24].velocity = 3.57E+09; 
    mod2[25].velocity = 3.31E+09; 
    mod2[26].velocity = 3.09E+09; 
    mod2[27].velocity = 2.96E+09; 
    mod2[28].velocity = 2.89E+09; 
    mod2[29].velocity = 2.87E+09; 
    mod2[30].velocity = 2.88E+09; 
    mod2[31].velocity = 2.91E+09; 
    mod2[32].velocity = 2.94E+09; 
    mod2[33].velocity = 2.98E+09; 
    mod2[34].velocity = 2.99E+09; 
    mod2[35].velocity = 3.00E+09; 
    mod2[36].velocity = 3.00E+09; 
    mod2[37].velocity = 2.99E+09; 
    mod2[38].velocity = 2.98E+09; 
    mod2[39].velocity = 2.96E+09; 
    mod2[40].velocity = 2.95E+09; 
    mod2[41].velocity = 2.94E+09; 
    mod2[42].velocity = 2.92E+09; 
    mod2[43].velocity = 2.91E+09; 
    mod2[44].velocity = 2.90E+09; 
    mod2[45].velocity = 2.89E+09; 
    mod2[46].velocity = 2.87E+09; 
    mod2[47].velocity = 2.85E+09; 
    mod2[48].velocity = 2.84E+09; 
    mod2[49].velocity = 2.82E+09; 
    mod2[50].velocity = 2.81E+09; 
    mod2[51].velocity = 2.79E+09; 
    mod2[52].velocity = 2.76E+09; 


    z_min = 100.0;
    for (counter = 1; counter <= 52; counter++)
    {
        z_diff = fabs(z - mod2[counter].z);
        
       
        if (z_diff < z_min)
        {
            z_min = z_diff;
            counter_min = counter;
            
        }
    }


    if ( z >= mod2[counter_min].z )
        velocity_value = mod2[counter_min].velocity + (mod2[counter_min+1].velocity - mod2[counter_min].velocity) * ( z  - mod2[counter_min].z ) / (mod2[counter_min+1].z - mod2[counter_min].z) ;
    else
        velocity_value = mod2[counter_min].velocity + (mod2[counter_min].velocity - mod2[counter_min-1].velocity) * ( z - mod2[counter_min].z ) / ( mod2[counter_min].z - mod2[counter_min-1].z );


//    velocity_value = mod2[counter_min].velocity;
    

    return (velocity_value);


}




/***********************************************************************/


double magnetic_field_south (double z)

{

    double z_diff, z_min, magnetic_flux, magnetic_field_tube;
    long counter, counter_min;

    mod[1].B_field =  14.84e-6;
    mod[2].B_field =  14.23e-6;
    mod[3].B_field =  13.75e-6;
    mod[4].B_field =  13.36e-6;
    mod[5].B_field =  13.12e-6;
    mod[6].B_field =  12.97e-6;
    mod[7].B_field =  12.74e-6;
    mod[8].B_field =  12.61e-6;
    mod[9].B_field =  12.52e-6;
    mod[10].B_field = 12.26e-6;
    mod[11].B_field = 12.08e-6;
    mod[12].B_field = 12.06e-6;
    mod[13].B_field = 12.e-6;
    mod[14].B_field = 11.89e-6;
    mod[15].B_field = 11.77e-6;
    mod[16].B_field = 11.66e-6;
    mod[17].B_field = 11.57e-6; 
    mod[18].B_field = 11.52e-6; 
    mod[19].B_field = 11.47e-6; 
    mod[20].B_field = 11.39e-6; 
    mod[21].B_field = 11.27e-6; 
    mod[22].B_field = 11.14e-6; 
    mod[23].B_field = 11.05e-6; 
    mod[24].B_field = 10.98e-6; 
    mod[25].B_field = 10.9e-6; 
    mod[26].B_field = 10.83e-6; 
    mod[27].B_field = 10.74e-6; 
    mod[28].B_field = 10.62e-6; 
    mod[29].B_field = 10.47e-6; 
    mod[30].B_field = 10.34e-6; 
    mod[31].B_field = 10.21e-6; 
    mod[32].B_field = 10.09e-6; 
    mod[33].B_field = 9.97e-6;
    mod[34].B_field = 9.85e-6;
    mod[35].B_field = 9.74e-6;
    mod[36].B_field = 9.63e-6;
    mod[37].B_field = 9.51e-6;
    mod[38].B_field = 9.4e-6;
    mod[39].B_field = 9.29e-6; 
    mod[40].B_field = 9.18e-6; 
    mod[41].B_field = 9.07e-6; 
    mod[42].B_field = 8.96e-6; 
    mod[43].B_field = 8.85e-6; 
    mod[44].B_field = 8.74e-6; 
    mod[45].B_field = 8.63e-6; 
    mod[46].B_field = 8.52e-6; 
    mod[47].B_field = 8.41e-6; 
    mod[48].B_field = 8.3e-6; 
    mod[49].B_field = 8.19e-6; 
    mod[50].B_field = 8.07e-6; 
    mod[51].B_field = 7.96e-6; 

    
    z_min = 100.0; 
    for (counter = 1; counter <= 51; counter++)
    {
        z_diff = fabs(z - mod[counter].z);
        
       
        if (z_diff < z_min)
        {
            z_min = z_diff;
            counter_min = counter;
            
        }
    }


    magnetic_field_tube = mod[counter_min].B_field;
        

    return (magnetic_field_tube);

}

/***********************************************************************/
    
double velocity_field_south (double z)

{

    double z_diff, z_min, magnetic_flux, magnetic_field_tube;
    long counter, counter_min;


    mod[1].velocity =  5.93e9;
    mod[2].velocity =  5.67e9;
    mod[3].velocity =  5.29e9;
    mod[4].velocity =  4.97e9;
    mod[5].velocity =  4.73e9;
    mod[6].velocity =  4.67e9;
    mod[7].velocity =  4.76e9;
    mod[8].velocity =  4.85e9;
    mod[9].velocity =  4.88e9;
    mod[10].velocity = 4.92e9;
    mod[11].velocity = 4.90e9;
    mod[12].velocity = 4.82e9;
    mod[13].velocity = 4.77e9;
    mod[14].velocity = 4.75e9;
    mod[15].velocity = 4.75e9;
    mod[16].velocity = 4.74e9;
    mod[17].velocity = 4.69e9; 
    mod[18].velocity = 4.58e9; 
    mod[19].velocity = 4.43e9; 
    mod[20].velocity = 4.27e9; 
    mod[21].velocity = 4.13e9; 
    mod[22].velocity = 3.99e9; 
    mod[23].velocity = 3.87e9; 
    mod[24].velocity = 3.77e9; 
    mod[25].velocity = 3.69e9; 
    mod[26].velocity = 3.62e9; 
    mod[27].velocity = 3.54e9; 
    mod[28].velocity = 3.48e9; 
    mod[29].velocity = 3.41e9; 
    mod[30].velocity = 3.33e9; 
    mod[31].velocity = 3.21e9; 
    mod[32].velocity = 3.10e9; 
    mod[33].velocity = 2.96e9;
    mod[34].velocity = 2.83e9;
    mod[35].velocity = 2.68e9;
    mod[36].velocity = 2.55e9;
    mod[37].velocity = 2.41e9;
    mod[38].velocity = 2.28e9;
    mod[39].velocity = 2.15e9; 
    mod[40].velocity = 2.04e9; 
    mod[41].velocity = 1.92e9; 
    mod[42].velocity = 1.82e9; 
    mod[43].velocity = 1.71e9; 
    mod[44].velocity = 1.62e9; 
    mod[45].velocity = 1.54e9; 
    mod[46].velocity = 1.46e9; 
    mod[47].velocity = 1.39e9; 
    mod[48].velocity = 1.33e9;
    mod[49].velocity = 1.27e9; 
    mod[50].velocity = 1.21e9; 
    mod[51].velocity = 1.17e9; 

    
    z_min = 100.0; 
    for (counter = 1; counter <= 53; counter++)
    {
        z_diff = fabs(z - mod[counter].z);
        
       
        if (z_diff < z_min)
        {
            z_min = z_diff;
            counter_min = counter;
            
        }
    }


    magnetic_field_tube = mod[counter_min].B_field;
        

    return (magnetic_field_tube);
    
}


/***********************************************************************/


void set_interpolate_values_north (double z, int ii)

{

    double z_diff, z_min, radius_jet;
    int counter, counter_min;
    

    /* mod[1].z =  103.00; */
    /* mod[2].z =  120.19; */
    /* mod[3].z =  134.97; */
    /* mod[4].z =  149.07; */
    /* mod[5].z =  163.17; */
    /* mod[6].z =  179.32; */
    /* mod[7].z =  195.14; */
    /* mod[8].z =  208.20; */
    /* mod[9].z =  220.58; */
    /* mod[10].z = 232.61; */
    /* mod[11].z = 245.68; */
    /* mod[12].z = 260.12; */
    /* mod[13].z = 272.15; */
    /* mod[14].z = 286.25; */
    /* mod[15].z = 300.34; */
    /* mod[16].z = 315.47; */
    /* mod[17].z = 329.56; */
    /* mod[18].z = 342.63; */
    /* mod[19].z = 353.63; */
    /* mod[20].z = 367.73; */
    /* mod[21].z = 382.51; */
    /* mod[22].z = 395.92; */
    /* mod[23].z = 408.29; */
    /* mod[24].z = 420.67; */
    /* mod[25].z = 433.74; */
    /* mod[26].z = 445.08; */
    /* mod[27].z = 455.74; */
    /* mod[28].z = 466.74; */
    /* mod[29].z = 470.18; */
    /* mod[30].z = 488.40; */
    /* mod[31].z = 499.75; */
    /* mod[32].z = 512.81; */
    /* mod[33].z = 531.03; */
    /* mod[34].z = 551.66; */
    /* mod[35].z = 570.91; */
    /* mod[36].z = 592.57; */
 
    z_min = 100.0;
    for (counter = 1; counter <= 75; counter++)
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

    /* if (ii==35) */
    /*     printf("ii=%i, z=%g, diff1=%g, diff2=%g, diff3=%g\n", mod[counter_min].ii, mod[counter_min].z, fabs(cr[ii][1].z / parsec / 1.e3 -  mod[counter_min].z), fabs(cr[ii-1][1].z / parsec / 1.e3 -  mod[counter_min].z) , fabs(cr[ii+1][1].z / parsec / 1.e3 -  mod[counter_min].z) );   */
    
    

    

}

/***********************************************************************/

double interpolated_value (double value, double value_low, double value_high, int ii, int ii_mod)
{

    double interp;
    
    

    if ( mod[ii].z >= cr[ii_mod][1].z / parsec / 1.e3 )
        interp = value + (value_high  - value) * ( mod[ii].z  - cr[ii_mod][1].z / parsec / 1.e3 ) / (cr[ii_mod+1][1].z / parsec / 1.e3 - cr[ii_mod][1].z / parsec / 1.e3) ;
    else
        interp = value + (value  - value_low) * ( mod[ii].z  - cr[ii_mod][1].z / parsec / 1.e3 ) / (cr[ii_mod][1].z / parsec / 1.e3 - cr[ii_mod-1][1].z / parsec / 1.e3 );


    /* if (ii_mod == 39) */
    /*     printf("z=%g value=%g, value_low = %g, value_high = %g, interp=%g\n", mod[ii].z, value, value_low, value_high, interp) */;
    
       

    return (interp);
    

    
}

/***********************************************************************/


void set_interpolate_values_south (double z, int ii)

{

    double z_diff, z_min, radius_jet;
    int counter, counter_min;
    

    /* mod[1].z =  116.00;  */
    /* mod[2].z =  123.56;  */
    /* mod[3].z =  131.81;  */
    /* mod[4].z =  141.10;  */
    /* mod[5].z =  152.10;  */
    /* mod[6].z =  163.10;  */
    /* mod[7].z =  173.76;  */
    /* mod[8].z =  186.14;  */
    /* mod[9].z =  197.82;  */
    /* mod[10].z = 210.20;  */
    /* mod[11].z = 223.61;  */
    /* mod[12].z = 237.71;  */
    /* mod[13].z = 254.55;  */
    /* mod[14].z = 272.09;  */
    /* mod[15].z = 290.31;  */
    /* mod[16].z = 309.56;  */
    /* mod[17].z = 324.34;  */
    /* mod[18].z = 352.53;  */
    /* mod[19].z = 364.22;  */
    /* mod[20].z = 382.79;  */
    /* mod[21].z = 399.64;  */
    /* mod[22].z = 416.14;  */
    /* mod[23].z = 434.36;  */
    /* mod[24].z = 448.45;  */
    /* mod[25].z = 462.55;  */
    /* mod[26].z = 478.02;  */
    /* mod[27].z = 494.52;  */
    /* mod[28].z = 509.99;  */
    /* mod[29].z = 524.43;  */
    /* mod[30].z = 540.94;  */
    /* mod[31].z = 562.25;  */
    /* mod[32].z = 583.22;  */
    /* mod[33].z = 604.20;  */
    /* mod[34].z = 626.20;  */
    /* mod[35].z = 648.55;  */
    /* mod[36].z = 669.17;  */
    /* mod[37].z = 691.52;  */
    /* mod[38].z = 711.81;  */
    /* mod[39].z = 729.34;  */
    /* mod[40].z = 743.78;  */
    /* mod[41].z = 758.22;  */
    /* mod[42].z = 772.66;  */
    /* mod[43].z = 794.66;  */
 
    z_min = 100.0;
    for (counter = 1; counter <= 94; counter++)
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

    /* if (ii==50) */
    /*     printf("ii=%i, z=%g, diff1=%g, diff2=%g, diff3=%g\n", mod[counter_min].ii, mod[counter_min].z, fabs(cr[ii][1].z / parsec / 1.e3 -  mod[counter_min].z), fabs(cr[ii-1][1].z / parsec / 1.e3 -  mod[counter_min].z) , fabs(cr[ii+1][1].z / parsec / 1.e3 -  mod[counter_min].z) ); */
  
}

void read_intensity_model (void)
{
  FILE *myfile;
  double myvariable;
  long ii;
  long jj;
  long kk;
  long height, width;
  double value[100];

  kk=0;
  height = 76;
  width = 1;
//  printf("hier=%li\n", kk);

  myfile=fopen("int_interp2.dat", "r");
  kk=0;
  

  for(ii = 0; ii < height; ii++)
  {
    for (jj = 0 ; jj < width; jj++)
    {
      fscanf(myfile,"%lf",&myvariable);
      mod[kk].intensity2 = myvariable;
      kk = kk + 1;
//      printf("k=%li\n", kk);
      
    }
  }

  fclose(myfile);
}

void read_magnetic_field_model (void)
{
  FILE *myfile;
  double myvariable;
  long ii;
  double value[100];


  myfile=fopen("b2.dat", "r");

  for(ii = 0; ii <= grid_size+1; ii++)
  {
      fscanf(myfile,"%lf",&myvariable);
      B_field[ii] = myvariable;
//      printf("ii=%li B=%g\n", ii, B_field[ii]);
  }

  fclose(myfile);
}
