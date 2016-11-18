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
    mod[0].radius =  3.37;
    mod[1].radius =  3.30;
    mod[2].radius =  3.35;
    mod[3].radius =  3.88;
    mod[4].radius =  3.94;
    mod[5].radius =  4.19;
    mod[6].radius =  4.61;
    mod[7].radius =  4.40;
    mod[8].radius =  4.38;
    mod[9].radius =  4.50;
    mod[10].radius = 4.73;
    mod[11].radius = 5.26;
    mod[12].radius = 5.59;
    mod[13].radius = 5.47;
    mod[14].radius = 4.76;
    mod[15].radius = 4.38;
    mod[16].radius = 4.33;
    mod[17].radius = 4.01;
    mod[18].radius = 4.16;
    mod[19].radius = 4.16;
    mod[20].radius = 4.66;
    mod[21].radius = 4.76;
    mod[22].radius = 4.68;
    mod[23].radius = 3.85;
    mod[24].radius = 4.54;
    mod[25].radius = 5.35;
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
    mod[64].alpha = -1.30; 
    mod[65].alpha = -1.40; 
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

    counter_max = 75;

    return (counter_max);
    
    

}
/***********************************************************************/
int set_radius_south ( void )

{
    int counter_max;
    
    mod[0].z =  2.20  ;
    mod[1].z =  4.23  ; 
    mod[2].z =  6.57  ; 
    mod[3].z =  8.77  ; 
    mod[4].z =  10.97 ; 
    mod[5].z =  13.10 ; 
    mod[6].z =  15.44 ; 
    mod[7].z =  17.29 ; 
    mod[8].z =  19.53 ; 
    mod[9].z =  21.63 ; 
    mod[10].z = 24.20 ; 
    mod[11].z = 26.61 ; 
    mod[12].z = 29.19 ; 
    mod[13].z = 31.80 ; 
    mod[14].z = 34.24 ; 
    mod[15].z = 36.55 ; 
    mod[16].z = 38.54 ; 
    mod[17].z = 41.15 ; 
    mod[18].z = 43.59 ; 
    mod[19].z = 45.97 ; 
    mod[20].z = 48.79 ; 
    mod[21].z = 51.98 ; 
    mod[22].z = 54.90 ; 
    mod[23].z = 56.83 ; 
    mod[24].z = 58.79 ; 
    mod[25].z = 60.72 ; 
    mod[26].z = 62.92 ; 
    mod[27].z = 64.81 ; 
    mod[28].z = 66.66 ; 
    mod[29].z = 71.82 ; 
    mod[30].z = 77.11 ; 
    mod[31].z = 81.27 ; 
    mod[32].z = 85.74 ; 
    mod[33].z = 91.59 ; 
    mod[34].z = 96.78 ; 
    mod[35].z = 103.07; 
    mod[36].z = 108.61; 
    mod[37].z = 114.31; 
    mod[38].z = 119.88; 
    mod[39].z = 124.25; 
    mod[40].z = 128.96; 
    mod[41].z = 133.19; 
    mod[42].z = 139.27; 
    mod[43].z = 144.53; 
    mod[44].z = 150.17; 
    mod[45].z = 155.57; 
    mod[46].z = 163.92; 
    mod[47].z = 174.27; 
    mod[48].z = 187.23; 
    mod[49].z = 203.67; 
    mod[50].z = 217.04; 
    mod[51].z = 228.15; 
    mod[52].z = 242.31; 
    mod[53].z = 263.04; 
    mod[54].z = 279.89; 
    mod[55].z = 290.75; 
    mod[56].z = 302.89; 
    mod[57].z = 315.13; 
    mod[58].z = 331.29; 
    mod[59].z = 352.26; 
    mod[60].z = 370.93; 
    mod[61].z = 385.95; 
    mod[62].z = 400.08; 
    mod[63].z = 417.85; 
    mod[64].z = 433.08; 
    mod[65].z = 446.77; 
    mod[66].z = 461.14; 
    mod[67].z = 475.85; 
    mod[68].z = 489.95; 
    mod[69].z = 506.21; 
    mod[70].z = 522.82; 
    mod[71].z = 534.82; 
    mod[72].z = 548.26; 
    mod[73].z = 563.25; 
    mod[74].z = 578.51; 
    mod[75].z = 592.37; 
    mod[76].z = 605.98; 
    mod[77].z = 625.13; 
    mod[78].z = 641.39; 
    mod[79].z = 664.05; 
    mod[80].z = 684.99; 
    mod[81].z = 726.59; 


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
    mod[29].radius = 19.79; 
    mod[30].radius = 20.01; 
    mod[31].radius = 18.81; 
    mod[32].radius = 19.17; 
    mod[33].radius = 16.04;
    mod[34].radius = 15.99;
    mod[35].radius = 16.31;
    mod[36].radius = 16.45;
    mod[37].radius = 15.90; 
    mod[38].radius = 14.94; 
    mod[39].radius = 14.16; 
    mod[40].radius = 13.08; 
    mod[41].radius = 12.31; 
    mod[42].radius = 11.91; 
    mod[43].radius = 10.98;
    mod[44].radius = 10.54; 
    mod[45].radius = 10.49; 
    mod[46].radius = 17.83; 
    mod[47].radius = 24.19; 
    mod[48].radius = 25.91; 
    mod[49].radius = 31.27; 
    mod[50].radius = 33.07; 
    mod[51].radius = 38.99; 
    mod[52].radius = 53.24; 
    mod[53].radius = 51.21; 
    mod[54].radius = 56.01; 
    mod[55].radius = 60.65; 
    mod[56].radius = 60.66; 
    mod[57].radius = 56.06; 
    mod[58].radius = 47.26; 
    mod[59].radius = 42.30; 
    mod[60].radius = 40.40; 
    mod[61].radius = 37.73; 
    mod[62].radius = 36.99; 
    mod[63].radius = 36.44; 
    mod[64].radius = 33.62; 
    mod[65].radius = 33.93; 
    mod[66].radius = 32.61; 
    mod[67].radius = 31.30; 
    mod[68].radius = 32.39; 
    mod[69].radius = 34.81; 
    mod[70].radius = 37.82; 
    mod[71].radius = 37.58; 
    mod[72].radius = 36.08; 
    mod[73].radius = 34.23; 
    mod[74].radius = 31.70; 
    mod[75].radius = 29.36; 
    mod[76].radius = 27.16; 
    mod[77].radius = 27.56; 
    mod[78].radius = 25.75; 
    mod[79].radius = 25.22; 
    mod[80].radius = 24.96; 
    mod[81].radius = 26.85; 


    mod[0].intensity =  1.99E+00;
    mod[1].intensity =  1.82E+00;
    mod[2].intensity =  1.70E+00;
    mod[3].intensity =  1.61E+00;
    mod[4].intensity =  1.54E+00;
    mod[5].intensity =  1.49E+00;
    mod[6].intensity =  1.42E+00;
    mod[7].intensity =  1.33E+00;
    mod[8].intensity =  1.33E+00;
    mod[9].intensity =  1.28E+00;
    mod[10].intensity = 1.22E+00;
    mod[11].intensity = 1.23E+00;
    mod[12].intensity = 1.20E+00;
    mod[13].intensity = 1.20E+00;
    mod[14].intensity = 1.23E+00;
    mod[15].intensity = 1.24E+00;
    mod[16].intensity = 1.26E+00;
    mod[17].intensity = 1.26E+00;
    mod[18].intensity = 1.30E+00;
    mod[19].intensity = 1.38E+00;
    mod[20].intensity = 1.41E+00;
    mod[21].intensity = 1.46E+00;
    mod[22].intensity = 1.44E+00;
    mod[23].intensity = 1.40E+00;
    mod[24].intensity = 1.41E+00;
    mod[25].intensity = 1.39E+00;
    mod[26].intensity = 1.33E+00;
    mod[27].intensity = 1.33E+00; 
    mod[28].intensity = 1.28E+00; 
    mod[29].intensity = 1.13E+00; 
    mod[30].intensity = 8.78E-01; 
    mod[31].intensity = 7.07E-01; 
    mod[32].intensity = 6.25E-01; 
    mod[33].intensity = 5.90E-01;
    mod[34].intensity = 5.78E-01;
    mod[35].intensity = 5.16E-01;
    mod[36].intensity = 4.28E-01;
    mod[37].intensity = 3.97E-01; 
    mod[38].intensity = 3.47E-01; 
    mod[39].intensity = 3.29E-01; 
    mod[40].intensity = 3.11E-01; 
    mod[41].intensity = 2.91E-01; 
    mod[42].intensity = 2.49E-01; 
    mod[43].intensity = 2.10E-01;
    mod[44].intensity = 1.72E-01; 
    mod[45].intensity = 1.35E-01; 
    mod[46].intensity = 8.39E-02; 
    mod[47].intensity = 5.97E-02; 
    mod[48].intensity = 4.87E-02; 
    mod[49].intensity = 4.58E-02; 
    mod[50].intensity = 4.70E-02; 
    mod[51].intensity = 4.74E-02; 
    mod[52].intensity = 5.11E-02; 
    mod[53].intensity = 6.52E-02; 
    mod[54].intensity = 7.00E-02; 
    mod[55].intensity = 7.24E-02; 
    mod[56].intensity = 6.75E-02; 
    mod[57].intensity = 5.86E-02; 
    mod[58].intensity = 5.31E-02; 
    mod[59].intensity = 4.40E-02; 
    mod[60].intensity = 3.82E-02; 
    mod[61].intensity = 4.00E-02; 
    mod[62].intensity = 4.04E-02; 
    mod[63].intensity = 3.80E-02; 
    mod[64].intensity = 3.30E-02; 
    mod[65].intensity = 2.82E-02; 
    mod[66].intensity = 2.85E-02; 
    mod[67].intensity = 3.00E-02; 
    mod[68].intensity = 2.73E-02; 
    mod[69].intensity = 2.57E-02; 
    mod[70].intensity = 2.45E-02; 
    mod[71].intensity = 2.38E-02; 
    mod[72].intensity = 2.15E-02; 
    mod[73].intensity = 2.03E-02; 
    mod[74].intensity = 2.03E-02; 
    mod[75].intensity = 2.05E-02;
    mod[76].intensity = 2.11E-02; 
    mod[77].intensity = 1.76E-02; 
    mod[78].intensity = 1.29E-02; 
    mod[79].intensity = 1.11E-02; 
    mod[80].intensity = 9.09E-03; 
    mod[81].intensity = 7.85E-03;

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
    mod[12].alpha = -0.58;
    mod[13].alpha = -0.60;
    mod[14].alpha = -0.61;
    mod[15].alpha = -0.62;
    mod[16].alpha = -0.62;
    mod[17].alpha = -0.62;
    mod[18].alpha = -0.62;
    mod[19].alpha = -0.60;
    mod[20].alpha = -0.59;
    mod[21].alpha = -0.56;
    mod[22].alpha = -0.55;
    mod[23].alpha = -0.55;
    mod[24].alpha = -0.55;
    mod[25].alpha = -0.55;
    mod[26].alpha = -0.55;
    mod[27].alpha = -0.55; 
    mod[28].alpha = -0.55; 
    mod[29].alpha = -0.54; 
    mod[30].alpha = -0.56; 
    mod[31].alpha = -0.54; 
    mod[32].alpha = -0.52; 
    mod[33].alpha = -0.52;
    mod[34].alpha = -0.52;
    mod[35].alpha = -0.55;
    mod[36].alpha = -0.59;
    mod[37].alpha = -0.64; 
    mod[38].alpha = -0.67; 
    mod[39].alpha = -0.67; 
    mod[40].alpha = -0.65; 
    mod[41].alpha = -0.62; 
    mod[42].alpha = -0.58; 
    mod[43].alpha = -0.55;
    mod[44].alpha = -0.55; 
    mod[45].alpha = -0.59; 
    mod[46].alpha = -0.63; 
    mod[47].alpha = -0.38; 
    mod[48].alpha = -0.29; 
    mod[49].alpha = -0.59; 
    mod[50].alpha = -0.65; 
    mod[51].alpha = -0.73; 
    mod[52].alpha = -0.84; 
    mod[53].alpha = -0.74; 
    mod[54].alpha = -0.64; 
    mod[55].alpha = -0.57; 
    mod[56].alpha = -0.63; 
    mod[57].alpha = -0.66; 
    mod[58].alpha = -0.63; 
    mod[59].alpha = -0.87; 
    mod[60].alpha = -0.91; 
    mod[61].alpha = -0.86; 
    mod[62].alpha = -0.81; 
    mod[63].alpha = -0.83; 
    mod[64].alpha = -0.91; 
    mod[65].alpha = -0.98; 
    mod[66].alpha = -1.01; 
    mod[67].alpha = -0.90; 
    mod[68].alpha = -1.00; 
    mod[69].alpha = -1.07; 
    mod[70].alpha = -0.96; 
    mod[71].alpha = -0.91; 
    mod[72].alpha = -0.96; 
    mod[73].alpha = -1.03; 
    mod[74].alpha = -1.01; 
    mod[75].alpha = -0.87; 
    mod[76].alpha = -0.85; 
    mod[77].alpha = -1.03; 
    mod[78].alpha = -1.13; 
    mod[79].alpha = -1.31; 
    mod[80].alpha = -1.46; 
    mod[81].alpha = -1.29;

    counter_max = 81;

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


  myfile=fopen("b2.dat", "r");

  for(ii = 0; ii <= grid_size+1; ii++)
  {
      fscanf(myfile,"%lf",&myvariable);
      B_field[ii] = myvariable;
//      printf("ii=%li B=%g\n", ii, B_field[ii]);
  }

  fclose(myfile);
}
/*EOF************************************************************************/
