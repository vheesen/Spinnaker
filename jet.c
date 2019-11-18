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
            
    mod[0].z =  0.00  ;
    mod[1].z =  2.62  ; 
    mod[2].z =  5.28  ; 
    mod[3].z =  7.94  ; 
    mod[4].z =  10.95 ; 
    mod[5].z =  13.48 ; 
    mod[6].z =  15.84 ; 
    mod[7].z =  18.63 ; 
    mod[8].z =  21.12 ; 
    mod[9].z =  23.65 ; 
    mod[10].z = 26.35 ; 
    mod[11].z = 29.19 ; 
    mod[12].z = 31.72 ; 
    mod[13].z = 34.29 ; 
    mod[14].z = 37.17 ; 
    mod[15].z = 39.27 ; 
    mod[16].z = 41.27 ; 
    mod[17].z = 42.97 ; 
    mod[18].z = 45.59 ; 
    mod[19].z = 48.08 ; 
    mod[20].z = 50.04 ; 
    mod[21].z = 52.66 ; 
    mod[22].z = 55.28 ; 
    mod[23].z = 58.07 ; 
    mod[24].z = 61.60 ; 
    mod[25].z = 66.14 ; 
    mod[26].z = 70.11 ; 
    mod[27].z = 73.78 ; 
    mod[28].z = 77.96 ; 
    mod[29].z = 83.20 ; 
    mod[30].z = 90.57 ; 
    mod[31].z = 96.99 ; 
    mod[32].z = 102.35; 
    mod[33].z = 108.11; 
    mod[34].z = 112.69; 
    mod[35].z = 117.27; 
    mod[36].z = 122.29; 
    mod[37].z = 128.09; 
    mod[38].z = 133.55; 
    mod[39].z = 143.76; 
    mod[40].z = 151.13; 
    mod[41].z = 161.51; 
    mod[42].z = 170.55; 
    mod[43].z = 180.10; 
    mod[44].z = 191.75; 
    mod[45].z = 207.59; 
    mod[46].z = 219.37; 
    mod[47].z = 228.66; 
    mod[48].z = 238.26; 
    mod[49].z = 248.29; 
    mod[50].z = 259.59; 
    mod[51].z = 274.91; 
    mod[52].z = 290.22; 
    mod[53].z = 304.75; 
    mod[54].z = 326.69; 
    mod[55].z = 345.50; 
    mod[56].z = 361.51; 
    mod[57].z = 380.23; 
    mod[58].z = 395.58; 
    mod[59].z = 411.73; 
    mod[60].z = 430.70; 
    mod[61].z = 451.38; 
    mod[62].z = 467.83; 
    mod[63].z = 489.12; 
    mod[64].z = 506.57; 
    mod[65].z = 524.20; 
    mod[66].z = 548.89; 
    mod[67].z = 572.58; 
    mod[68].z = 593.40; 
    mod[69].z = 617.57; 
    mod[70].z = 641.13; 
    mod[71].z = 664.25; 
    mod[72].z = 686.85; 
    mod[73].z = 712.28; 
    mod[74].z = 736.98; 
    mod[75].z = 780.30; 


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
    mod[0].radius =  3.37 ;
    mod[1].radius =  3.30 ;
    mod[2].radius =  3.35 ;
    mod[3].radius =  3.88 ;
    mod[4].radius =  3.94 ;
    mod[5].radius =  4.19 ;
    mod[6].radius =  4.61 ;
    mod[7].radius =  4.40 ;
    mod[8].radius =  4.38 ;
    mod[9].radius =  4.50 ;
    mod[10].radius = 4.73 ;
    mod[11].radius = 5.26 ;
    mod[12].radius = 5.59 ;
    mod[13].radius = 5.47 ;
    mod[14].radius = 4.76 ;
    mod[15].radius = 4.38 ;
    mod[16].radius = 4.33 ;
    mod[17].radius = 4.01 ;
    mod[18].radius = 4.16 ;
    mod[19].radius = 4.16 ;
    mod[20].radius = 4.66 ;
    mod[21].radius = 4.76 ;
    mod[22].radius = 4.68 ;
    mod[23].radius = 3.85 ;
    mod[24].radius = 4.54 ;
    mod[25].radius = 5.35 ;
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
    mod[51].radius = 29.05; 
    mod[52].radius = 34.38; 
    mod[53].radius = 37.82; 
    mod[54].radius = 42.29; 
    mod[55].radius = 38.33; 
    mod[56].radius = 36.96; 
    mod[57].radius = 36.10; 
    mod[58].radius = 37.82; 
    mod[59].radius = 42.98; 
    mod[60].radius = 45.21; 
    mod[61].radius = 44.69; 
    mod[62].radius = 46.41; 
    mod[63].radius = 45.90; 
    mod[64].radius = 46.41; 
    mod[65].radius = 44.01; 
    mod[66].radius = 42.12; 
    mod[67].radius = 34.38; 
    mod[68].radius = 26.82; 
    mod[69].radius = 22.86; 
    mod[70].radius = 19.77; 
    mod[71].radius = 19.94; 
    mod[72].radius = 22.35; 
    mod[73].radius = 21.49; 
    mod[74].radius = 27.50; 
    mod[75].radius = 23.89; 


    mod[0].kappa =  5.61E-04;
    mod[1].kappa =  6.34E-04;
    mod[2].kappa =  6.54E-04;
    mod[3].kappa =  5.95E-04;
    mod[4].kappa =  6.01E-04;
    mod[5].kappa =  5.79E-04;
    mod[6].kappa =  5.03E-04;
    mod[7].kappa =  5.19E-04;
    mod[8].kappa =  4.80E-04;
    mod[9].kappa =  4.27E-04;
    mod[10].kappa = 4.02E-04;
    mod[11].kappa = 3.39E-04;
    mod[12].kappa = 2.93E-04;
    mod[13].kappa = 2.76E-04;
    mod[14].kappa = 2.81E-04;
    mod[15].kappa = 2.80E-04;
    mod[16].kappa = 2.67E-04;
    mod[17].kappa = 2.61E-04;
    mod[18].kappa = 2.22E-04;
    mod[19].kappa = 2.14E-04;
    mod[20].kappa = 1.67E-04;
    mod[21].kappa = 1.51E-04;
    mod[22].kappa = 1.45E-04;
    mod[23].kappa = 1.59E-04;
    mod[24].kappa = 1.50E-04;
    mod[25].kappa = 1.32E-04;
    mod[26].kappa = 1.31E-04;
    mod[27].kappa = 1.32E-04; 
    mod[28].kappa = 6.86E-05; 
    mod[29].kappa = 5.27E-05; 
    mod[30].kappa = 4.23E-05; 
    mod[31].kappa = 4.15E-05; 
    mod[32].kappa = 4.28E-05; 
    mod[33].kappa = 4.07E-05;
    mod[34].kappa = 3.69E-05;
    mod[35].kappa = 3.02E-05;
    mod[36].kappa = 2.88E-05;
    mod[37].kappa = 2.90E-05; 
    mod[38].kappa = 2.65E-05; 
    mod[39].kappa = 1.34E-05; 
    mod[40].kappa = 1.12E-05; 
    mod[41].kappa = 9.21E-06; 
    mod[42].kappa = 8.12E-06; 
    mod[43].kappa = 7.64E-06;
    mod[44].kappa = 7.61E-06; 
    mod[45].kappa = 6.84E-06; 
    mod[46].kappa = 4.66E-06; 
    mod[47].kappa = 4.34E-06; 
    mod[48].kappa = 3.80E-06; 
    mod[49].kappa = 3.49E-06; 
    mod[50].kappa = 3.48E-06; 
    mod[51].kappa = 2.62E-06; 
    mod[52].kappa = 1.93E-06; 
    mod[53].kappa = 1.60E-06; 
    mod[54].kappa = 1.27E-06; 
    mod[55].kappa = 1.37E-06; 
    mod[56].kappa = 1.42E-06; 
    mod[57].kappa = 1.23E-06; 
    mod[58].kappa = 9.97E-07; 
    mod[59].kappa = 7.63E-07; 
    mod[60].kappa = 6.00E-07; 
    mod[61].kappa = 6.42E-07; 
    mod[62].kappa = 5.71E-07; 
    mod[63].kappa = 5.55E-07; 
    mod[64].kappa = 4.87E-07; 
    mod[65].kappa = 4.55E-07; 
    mod[66].kappa = 4.32E-07; 
    mod[67].kappa = 4.74E-07; 
    mod[68].kappa = 5.81E-07; 
    mod[69].kappa = 7.03E-07; 
    mod[70].kappa = 6.58E-07; 
    mod[71].kappa = 4.84E-07; 
    mod[72].kappa = 3.90E-07; 
    mod[73].kappa = 4.05E-07; 
    mod[74].kappa = 2.78E-07; 
    mod[75].kappa = 3.16E-07;

/* 52-145 MHz spectral index */
    mod[0].alpha =  -0.42;
    mod[1].alpha =  -0.40;
    mod[2].alpha =  -0.40;
    mod[3].alpha =  -0.39;
    mod[4].alpha =  -0.39;
    mod[5].alpha =  -0.39;
    mod[6].alpha =  -0.40;
    mod[7].alpha =  -0.40;
    mod[8].alpha =  -0.42;
    mod[9].alpha =  -0.42;
    mod[10].alpha = -0.43;
    mod[11].alpha = -0.44;
    mod[12].alpha = -0.48;
    mod[13].alpha = -0.49;
    mod[14].alpha = -0.54;
    mod[15].alpha = -0.57;
    mod[16].alpha = -0.58;
    mod[17].alpha = -0.60;
    mod[18].alpha = -0.64;
    mod[19].alpha = -0.64;
    mod[20].alpha = -0.66;
    mod[21].alpha = -0.66;
    mod[22].alpha = -0.66;
    mod[23].alpha = -0.66;
    mod[24].alpha = -0.62;
    mod[25].alpha = -0.58;
    mod[26].alpha = -0.55;
    mod[27].alpha = -0.53; 
    mod[28].alpha = -0.55; 
    mod[29].alpha = -0.56; 
    mod[30].alpha = -0.53; 
    mod[31].alpha = -0.49; 
    mod[32].alpha = -0.46; 
    mod[33].alpha = -0.44;
    mod[34].alpha = -0.41;
    mod[35].alpha = -0.40;
    mod[36].alpha = -0.41;
    mod[37].alpha = -0.42; 
    mod[38].alpha = -0.43; 
    mod[39].alpha = -0.67; 
    mod[40].alpha = -0.73; 
    mod[41].alpha = -0.78; 
    mod[42].alpha = -0.80; 
    mod[43].alpha = -0.65;
    mod[44].alpha = -0.41; 
    mod[45].alpha = -0.32; 
    mod[46].alpha = -0.33; 
    mod[47].alpha = -0.38; 
    mod[48].alpha = -0.44; 
    mod[49].alpha = -0.50; 
    mod[50].alpha = -0.60; 
    mod[51].alpha = -0.70; 
    mod[52].alpha = -0.73; 
    mod[53].alpha = -0.81; 
    mod[54].alpha = -0.86; 
    mod[55].alpha = -0.82; 
    mod[56].alpha = -0.82; 
    mod[57].alpha = -0.92; 
    mod[58].alpha = -0.93; 
    mod[59].alpha = -0.86; 
    mod[60].alpha = -0.94; 
    mod[61].alpha = -0.91; 
    mod[62].alpha = -0.88; 
    mod[63].alpha = -0.88; 
    mod[64].alpha = -0.96; 
    mod[65].alpha = -1.02; 
    mod[66].alpha = -1.14; 
    mod[67].alpha = -1.20; 
    mod[68].alpha = -1.12; 
    mod[69].alpha = -0.91; 
    mod[70].alpha = -0.98; 
    mod[71].alpha = -1.39; 
    mod[72].alpha = -1.30; 
    mod[73].alpha = -1.30; 
    mod[74].alpha = -1.29; 
    mod[75].alpha = -1.33;

    counter_max = 75;

    return (counter_max);
    
    

}
/***********************************************************************/
int set_radius_south ( void )

{
    int counter_max;
    
    mod[0].z =  2.79  ;
    mod[1].z =  5.37  ; 
    mod[2].z =  8.33  ; 
    mod[3].z =  11.13 ; 
    mod[4].z =  13.92 ; 
    mod[5].z =  16.62 ; 
    mod[6].z =  19.59 ; 
    mod[7].z =  21.95 ; 
    mod[8].z =  24.78 ; 
    mod[9].z =  27.44 ; 
    mod[10].z = 30.71 ; 
    mod[11].z = 33.77 ; 
    mod[12].z = 37.04 ; 
    mod[13].z = 40.36 ; 
    mod[14].z = 43.45 ; 
    mod[15].z = 46.38 ; 
    mod[16].z = 48.91 ; 
    mod[17].z = 52.22 ; 
    mod[18].z = 55.32 ; 
    mod[19].z = 58.33 ; 
    mod[20].z = 61.91 ; 
    mod[21].z = 65.97 ; 
    mod[22].z = 69.68 ; 
    mod[23].z = 72.12 ; 
    mod[24].z = 74.61 ; 
    mod[25].z = 77.05 ; 
    mod[26].z = 79.84 ; 
    mod[27].z = 82.24 ; 
    mod[28].z = 84.60 ; 
    mod[29].z = 87.00 ; 
    mod[30].z = 89.40 ; 
    mod[31].z = 95.94 ; 
    mod[32].z = 102.66; 
    mod[33].z = 107.94; 
    mod[34].z = 113.61; 
    mod[35].z = 120.85; 
    mod[36].z = 127.92; 
    mod[37].z = 133.46; 
    mod[38].z = 139.44; 
    mod[39].z = 144.80; 
    mod[40].z = 152.53; 
    mod[41].z = 159.20; 
    mod[42].z = 166.36; 
    mod[43].z = 173.21; 
    mod[44].z = 183.81; 
    mod[45].z = 196.94; 
    mod[46].z = 213.39; 
    mod[47].z = 234.24; 
    mod[48].z = 251.21; 
    mod[49].z = 265.31; 
    mod[50].z = 283.28; 
    mod[51].z = 309.59; 
    mod[52].z = 330.97; 
    mod[53].z = 344.76; 
    mod[54].z = 360.16; 
    mod[55].z = 375.69; 
    mod[56].z = 396.19; 
    mod[57].z = 422.81; 
    mod[58].z = 446.50; 
    mod[59].z = 465.56; 
    mod[60].z = 483.49; 
    mod[61].z = 506.05; 
    mod[62].z = 525.38; 
    mod[63].z = 542.74; 
    mod[64].z = 560.98; 
    mod[65].z = 579.65; 
    mod[66].z = 597.54; 
    mod[67].z = 618.18; 
    mod[68].z = 639.25; 
    mod[69].z = 654.48; 
    mod[70].z = 671.54; 
    mod[71].z = 690.56; 
    mod[72].z = 709.93; 
    mod[73].z = 727.51; 
    mod[74].z = 744.79; 
    mod[75].z = 769.09; 
    mod[76].z = 789.73; 
    mod[77].z = 818.48; 
    mod[78].z = 845.05; 
    mod[79].z = 897.84; 


    mod[0].radius =  3.32 ;
    mod[1].radius =  3.58 ;
    mod[2].radius =  4.16 ;
    mod[3].radius =  4.68 ;
    mod[4].radius =  4.49 ;
    mod[5].radius =  4.69 ;
    mod[6].radius =  4.86 ;
    mod[7].radius =  5.81 ;
    mod[8].radius =  6.00 ;
    mod[9].radius =  6.41 ;
    mod[10].radius = 6.57 ;
    mod[11].radius = 6.19 ;
    mod[12].radius = 5.90 ;
    mod[13].radius = 5.62 ;
    mod[14].radius = 5.35 ;
    mod[15].radius = 5.62 ;
    mod[16].radius = 5.76 ;
    mod[17].radius = 5.84 ;
    mod[18].radius = 5.90 ;
    mod[19].radius = 5.93 ;
    mod[20].radius = 5.59 ;
    mod[21].radius = 6.12 ;
    mod[22].radius = 6.81 ;
    mod[23].radius = 8.25 ;
    mod[24].radius = 9.32 ;
    mod[25].radius = 10.83;
    mod[26].radius = 11.78;
    mod[27].radius = 12.15; 
    mod[28].radius = 12.27; 
    mod[29].radius = 10.83; 
    mod[30].radius = 11.17; 
    mod[31].radius = 7.74 ; 
    mod[32].radius = 8.25 ; 
    mod[33].radius = 8.25 ;
    mod[34].radius = 10.49;
    mod[35].radius = 15.90;
    mod[36].radius = 14.94;
    mod[37].radius = 14.16; 
    mod[38].radius = 13.08; 
    mod[39].radius = 12.31; 
    mod[40].radius = 11.91; 
    mod[41].radius = 10.98; 
    mod[42].radius = 10.54; 
    mod[43].radius = 10.49;
    mod[44].radius = 17.83; 
    mod[45].radius = 24.19; 
    mod[46].radius = 25.91; 
    mod[47].radius = 31.27; 
    mod[48].radius = 33.07; 
    mod[49].radius = 38.99; 
    mod[50].radius = 53.24; 
    mod[51].radius = 51.21; 
    mod[52].radius = 56.01; 
    mod[53].radius = 60.65; 
    mod[54].radius = 60.66; 
    mod[55].radius = 56.06; 
    mod[56].radius = 47.26; 
    mod[57].radius = 42.30; 
    mod[58].radius = 40.40; 
    mod[59].radius = 37.73; 
    mod[60].radius = 36.99; 
    mod[61].radius = 36.44; 
    mod[62].radius = 33.62; 
    mod[63].radius = 33.93; 
    mod[64].radius = 32.61; 
    mod[65].radius = 31.30; 
    mod[66].radius = 32.39; 
    mod[67].radius = 34.81; 
    mod[68].radius = 37.82; 
    mod[69].radius = 37.58; 
    mod[70].radius = 36.08; 
    mod[71].radius = 34.23; 
    mod[72].radius = 31.70; 
    mod[73].radius = 29.36; 
    mod[74].radius = 27.16; 
    mod[75].radius = 27.56; 
    mod[76].radius = 25.75; 
    mod[77].radius = 25.22; 
    mod[78].radius = 24.96; 
    mod[79].radius = 26.85; 

    mod[0].kappa =  5.44E-04;
    mod[1].kappa =  4.61E-04;
    mod[2].kappa =  3.70E-04;
    mod[3].kappa =  3.12E-04;
    mod[4].kappa =  3.11E-04;
    mod[5].kappa =  2.87E-04;
    mod[6].kappa =  2.64E-04;
    mod[7].kappa =  2.07E-04;
    mod[8].kappa =  2.00E-04;
    mod[9].kappa =  1.81E-04;
    mod[10].kappa = 1.68E-04;
    mod[11].kappa = 1.80E-04;
    mod[12].kappa = 1.85E-04;
    mod[13].kappa = 1.93E-04;
    mod[14].kappa = 2.07E-04;
    mod[15].kappa = 2.00E-04;
    mod[16].kappa = 1.98E-04;
    mod[17].kappa = 1.96E-04;
    mod[18].kappa = 1.99E-04;
    mod[19].kappa = 2.11E-04;
    mod[20].kappa = 2.29E-04;
    mod[21].kappa = 2.16E-04;
    mod[22].kappa = 1.92E-04;
    mod[23].kappa = 1.53E-04;
    mod[24].kappa = 1.37E-04;
    mod[25].kappa = 1.16E-04;
    mod[26].kappa = 1.02E-04;
    mod[27].kappa = 9.93E-05; 
    mod[28].kappa = 9.41E-05; 
    mod[29].kappa = 9.85E-05; 
    mod[30].kappa = 8.70E-05; 
    mod[31].kappa = 1.15E-04; 
    mod[32].kappa = 8.80E-05; 
    mod[33].kappa = 7.03E-05;
    mod[34].kappa = 4.35E-05;
    mod[35].kappa = 2.26E-05;
    mod[36].kappa = 2.10E-05;
    mod[37].kappa = 2.10E-05; 
    mod[38].kappa = 2.15E-05; 
    mod[39].kappa = 2.14E-05; 
    mod[40].kappa = 1.90E-05; 
    mod[41].kappa = 1.73E-05; 
    mod[42].kappa = 1.48E-05; 
    mod[43].kappa = 1.16E-05;
    mod[44].kappa = 4.26E-06; 
    mod[45].kappa = 2.23E-06; 
    mod[46].kappa = 1.70E-06; 
    mod[47].kappa = 1.33E-06; 
    mod[48].kappa = 1.28E-06; 
    mod[49].kappa = 1.10E-06; 
    mod[50].kappa = 8.68E-07; 
    mod[51].kappa = 1.15E-06; 
    mod[52].kappa = 1.13E-06; 
    mod[53].kappa = 1.08E-06; 
    mod[54].kappa = 1.01E-06; 
    mod[55].kappa = 9.46E-07; 
    mod[56].kappa = 1.02E-06; 
    mod[57].kappa = 9.42E-07; 
    mod[58].kappa = 8.57E-07; 
    mod[59].kappa = 9.58E-07; 
    mod[60].kappa = 9.89E-07; 
    mod[61].kappa = 9.44E-07; 
    mod[62].kappa = 8.87E-07; 
    mod[63].kappa = 7.53E-07; 
    mod[64].kappa = 7.91E-07; 
    mod[65].kappa = 8.67E-07; 
    mod[66].kappa = 7.64E-07; 
    mod[67].kappa = 6.69E-07; 
    mod[68].kappa = 5.86E-07; 
    mod[69].kappa = 5.73E-07; 
    mod[70].kappa = 5.39E-07; 
    mod[71].kappa = 5.35E-07; 
    mod[72].kappa = 5.80E-07; 
    mod[73].kappa = 6.31E-07; 
    mod[74].kappa = 7.04E-07; 
    mod[75].kappa = 5.79E-07;
    mod[76].kappa = 4.52E-07; 
    mod[77].kappa = 3.98E-07; 
    mod[78].kappa = 3.30E-07; 
    mod[79].kappa = 2.65E-07; 

/* 52-145 MHz spectral index */
    mod[0].alpha =  -0.43;
    mod[1].alpha =  -0.45;
    mod[2].alpha =  -0.47;
    mod[3].alpha =  -0.46;
    mod[4].alpha =  -0.46;
    mod[5].alpha =  -0.50;
    mod[6].alpha =  -0.51;
    mod[7].alpha =  -0.51;
    mod[8].alpha =  -0.51;
    mod[9].alpha =  -0.52;
    mod[10].alpha = -0.53;
    mod[11].alpha = -0.55;
    mod[12].alpha = -0.58;
    mod[13].alpha = -0.60;
    mod[14].alpha = -0.61;
    mod[15].alpha = -0.61;
    mod[16].alpha = -0.61;
    mod[17].alpha = -0.61;
    mod[18].alpha = -0.61;
    mod[19].alpha = -0.60;
    mod[20].alpha = -0.59;
    mod[21].alpha = -0.55;
    mod[22].alpha = -0.52;
    mod[23].alpha = -0.49;
    mod[24].alpha = -0.46;
    mod[25].alpha = -0.44;
    mod[26].alpha = -0.42;
    mod[27].alpha = -0.41; 
    mod[28].alpha = -0.40; 
    mod[29].alpha = -0.40; 
    mod[30].alpha = -0.42; 
    mod[31].alpha = -0.47; 
    mod[32].alpha = -0.50; 
    mod[33].alpha = -0.53;
    mod[34].alpha = -0.57;
    mod[35].alpha = -0.64;
    mod[36].alpha = -0.67;
    mod[37].alpha = -0.67; 
    mod[38].alpha = -0.64; 
    mod[39].alpha = -0.61; 
    mod[40].alpha = -0.58; 
    mod[41].alpha = -0.55; 
    mod[42].alpha = -0.55; 
    mod[43].alpha = -0.59;
    mod[44].alpha = -0.63; 
    mod[45].alpha = -0.38; 
    mod[46].alpha = -0.29; 
    mod[47].alpha = -0.58; 
    mod[48].alpha = -0.65; 
    mod[49].alpha = -0.72; 
    mod[50].alpha = -0.84; 
    mod[51].alpha = -0.74; 
    mod[52].alpha = -0.63; 
    mod[53].alpha = -0.56; 
    mod[54].alpha = -0.63; 
    mod[55].alpha = -0.65; 
    mod[56].alpha = -0.62; 
    mod[57].alpha = -0.86; 
    mod[58].alpha = -0.90; 
    mod[59].alpha = -0.85; 
    mod[60].alpha = -0.80; 
    mod[61].alpha = -0.83; 
    mod[62].alpha = -0.90; 
    mod[63].alpha = -0.98; 
    mod[64].alpha = -1.01; 
    mod[65].alpha = -0.90; 
    mod[66].alpha = -0.99; 
    mod[67].alpha = -1.06; 
    mod[68].alpha = -0.96; 
    mod[69].alpha = -0.91; 
    mod[70].alpha = -0.95; 
    mod[71].alpha = -1.03; 
    mod[72].alpha = -1.01; 
    mod[73].alpha = -0.87; 
    mod[74].alpha = -0.85; 
    mod[75].alpha = -1.02; 
    mod[76].alpha = -1.13; 
    mod[77].alpha = -1.30; 
    mod[78].alpha = -1.45; 
    mod[79].alpha = -1.28; 

    counter_max = 79;

    return (counter_max);
    
}

/***********************************************************************/

double radius (double z)

{

    double z_diff, z_min, radius_jet;
    int counter, counter_min;

    if ( model == 1 )
    {
        

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




    }

    else

    {

        if ( velocity_field == 3)
            radius_jet = R0 / kpc * pow(1.0 + pow(z / h_V, beta), 0.5);
        else
            radius_jet = R0 / kpc;


    }
    
        

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

    double z_diff, z_min, magnetic_field_strength, epsilon, kappa, alpha;
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
        kappa= mod[counter_min].kappa + (mod[counter_min+1].kappa - mod[counter_min].kappa) * ( z  - mod[counter_min].z ) / (mod[counter_min+1].z - mod[counter_min].z) ;
    else
        kappa = mod[counter_min].kappa + (mod[counter_min].kappa - mod[counter_min-1].kappa) * ( z - mod[counter_min].z ) / ( mod[counter_min].z - mod[counter_min-1].z );

    if ( z >= mod[counter_min].z )
        epsilon = mod[counter_min].epsilon + (mod[counter_min+1].epsilon - mod[counter_min].epsilon) * ( z  - mod[counter_min].z ) / (mod[counter_min+1].z - mod[counter_min].z) ;
    else
        epsilon = mod[counter_min].epsilon + (mod[counter_min].epsilon - mod[counter_min-1].epsilon) * ( z - mod[counter_min].z ) / ( mod[counter_min].z - mod[counter_min-1].z );

    if ( z >= mod[counter_min].z )
        alpha = mod[counter_min].alpha + (mod[counter_min+1].alpha - mod[counter_min].alpha) * ( z  - mod[counter_min].z ) / (mod[counter_min+1].z - mod[counter_min].z) ;
    else
        alpha = mod[counter_min].alpha + (mod[counter_min].alpha - mod[counter_min-1].alpha) * ( z - mod[counter_min].z ) / ( mod[counter_min].z - mod[counter_min-1].z );

    if (update_model == 1)
    {
//        printf("B=%g, intensity =%g, intensity2=%g, alpha=%g\n", B_field[i], intensity, intensity2, alpha);
        
        magnetic_field_strength =  B_field[i] * pow(kappa / (xi * epsilon), 1./(1.-alpha));
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
    for (counter = 0; counter <= number_of_data_points; counter++)
    {
        z_diff = fabs(z - mod[counter].z);
       
        if (z_diff < z_min)
        {
            z_min = z_diff;
            counter_min = counter;
        }
    }

    if (ii == 0)
    {
        mod[counter_min].ii = ii;
    }

    else
    {
        
        if ( ( fabs(cr[ii][0].z / kpc -  mod[counter_min].z) <= fabs(cr[ii-1][0].z / kpc -  mod[counter_min].z) ) && ( fabs(cr[ii][0].z / kpc -  mod[counter_min].z) <= fabs(cr[ii+1][0].z / kpc -  mod[counter_min].z) ) )
            mod[counter_min].ii = ii;
    }

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

  myfile=fopen("epsilon_interp2.dat", "r");
  kk=0;

  for(ii = 0; ii <= number_of_data_points; ii++)
  {
      fscanf(myfile,"%lf",&myvariable);
      mod[kk].epsilon = myvariable;
      kk = kk + 1;
//      printf("k=%i intensity=%g\n", kk, mod[kk].epsilon);
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
