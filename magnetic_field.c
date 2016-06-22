#include "dec.h"
#include "cosmic_rays.h"
#include "runge_kutta.h"
#include "read_parameters.h"
#include "synchrotron.h"
#include "magnetic_field.h"

/***********************************************************************/


double magnetic_field_north (double z)

{

    double z_diff, z_min, magnetic_flux, magnetic_field_tube;
    long counter, counter_min;

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
    mod[52].z = 117.; 
    mod[53].z = 119.;
 
    mod[1].B_field =  15.09e-6;
    mod[2].B_field =  14.44e-6;
    mod[3].B_field =  14.02e-6;
    mod[4].B_field =  13.74e-6;
    mod[5].B_field =  13.58e-6;
    mod[6].B_field =  13.48e-6;
    mod[7].B_field =  13.35e-6;
    mod[8].B_field =  13.19e-6;
    mod[9].B_field =  13.03e-6;
    mod[10].B_field = 12.87e-6;
    mod[11].B_field = 12.71e-6;
    mod[12].B_field = 12.54e-6;
    mod[13].B_field = 12.38e-6;
    mod[14].B_field = 12.25e-6;
    mod[15].B_field = 12.14e-6;
    mod[16].B_field = 12.09e-6;
    mod[17].B_field = 12.08e-6; 
    mod[18].B_field = 12.09e-6; 
    mod[19].B_field = 12.1e-6; 
    mod[20].B_field = 12.09e-6; 
    mod[21].B_field = 12.06e-6; 
    mod[22].B_field = 11.97e-6; 
    mod[23].B_field = 11.84e-6; 
    mod[24].B_field = 11.66e-6; 
    mod[25].B_field = 11.45e-6; 
    mod[26].B_field = 11.23e-6; 
    mod[27].B_field = 11.02e-6; 
    mod[28].B_field = 10.82e-6; 
    mod[29].B_field = 10.64e-6; 
    mod[30].B_field = 10.5e-6; 
    mod[31].B_field = 10.37e-6; 
    mod[32].B_field = 10.26e-6; 
    mod[33].B_field = 10.17e-6;
    mod[34].B_field = 10.08e-6;
    mod[35].B_field = 10e-6;
    mod[36].B_field = 9.91e-6;
    mod[37].B_field = 9.82e-6;
    mod[38].B_field = 9.73e-6;
    mod[39].B_field = 9.63e-6; 
    mod[40].B_field = 9.52e-6; 
    mod[41].B_field = 9.41e-6; 
    mod[42].B_field = 9.3e-6; 
    mod[43].B_field = 9.19e-6; 
    mod[44].B_field = 9.08e-6; 
    mod[45].B_field = 8.98e-6; 
    mod[46].B_field = 8.88e-6; 
    mod[47].B_field = 8.78e-6; 
    mod[48].B_field = 8.69e-6; 
    mod[49].B_field = 8.6e-6; 
    mod[50].B_field = 8.52e-6; 
    mod[51].B_field = 8.43e-6; 
    mod[52].B_field = 8.35e-6; 
    mod[53].B_field = 8.27e-6;


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



double velocity_field_north (double z)

{

    double z_diff, z_min;
    long counter, counter_min;

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
    mod[52].z = 117.; 
    mod[53].z = 119.;
 
    mod[1].velocity =  5.80e9;
    mod[2].velocity =  5.54e9;
    mod[3].velocity =  5.27e9;
    mod[4].velocity =  5.04e9;
    mod[5].velocity =  4.90e9;
    mod[6].velocity =  4.89e9;
    mod[7].velocity =  4.97e9;
    mod[8].velocity =  5.07e9;
    mod[9].velocity =  5.16e9;
    mod[10].velocity = 5.18e9;
    mod[11].velocity = 5.16e9;
    mod[12].velocity = 5.17e9;
    mod[13].velocity = 5.18e9;
    mod[14].velocity = 5.25e9;
    mod[15].velocity = 5.31e9;
    mod[16].velocity = 5.36e9;
    mod[17].velocity = 5.36e9; 
    mod[18].velocity = 5.33e9; 
    mod[19].velocity = 5.26e9;
    mod[20].velocity = 5.11e9; 
    mod[21].velocity = 4.88e9; 
    mod[22].velocity = 4.58e9; 
    mod[23].velocity = 4.23e9; 
    mod[24].velocity = 3.88e9; 
    mod[25].velocity = 3.57e9; 
    mod[26].velocity = 3.31e9; 
    mod[27].velocity = 3.09e9; 
    mod[28].velocity = 2.96e9; 
    mod[29].velocity = 2.89e9; 
    mod[30].velocity = 2.87e9;
    mod[31].velocity = 2.88e9; 
    mod[32].velocity = 2.91e9; 
    mod[33].velocity = 2.94e9;
    mod[34].velocity = 2.98e9;
    mod[35].velocity = 2.99e9;
    mod[36].velocity = 3.00e9;
    mod[37].velocity = 3.00e9;
    mod[38].velocity = 2.99e9;
    mod[39].velocity = 2.98e9; 
    mod[40].velocity = 2.96e9; 
    mod[41].velocity = 2.95e9; 
    mod[42].velocity = 2.94e9;
    mod[43].velocity = 2.92e9; 
    mod[44].velocity = 2.91e9; 
    mod[45].velocity = 2.90e9; 
    mod[46].velocity = 2.89e9; 
    mod[47].velocity = 2.87e9; 
    mod[48].velocity = 2.85e9; 
    mod[49].velocity = 2.84e9;
    mod[50].velocity = 2.82e9; 
    mod[51].velocity = 2.81e9; 
    mod[52].velocity = 2.79e9; 
    mod[53].velocity = 2.76e9;
   

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

//    if (j==1)
//        printf("z=%g, V=%g\n", z, mod[counter_min].velocity);
    


    return (mod[counter_min].velocity);

    


}



/***********************************************************************/


double magnetic_field_south (double z)

{

    double z_diff, z_min, magnetic_flux, magnetic_field_tube;
    long counter, counter_min;

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
    mod[52].z = 117.; 
    mod[53].z = 119.;
 
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
    mod[52].B_field = 7.85e-6; 
    mod[53].B_field = 7.73e-6;


    
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
    
double velocity_field_south (double z)

{

    double z_diff, z_min, magnetic_flux, magnetic_field_tube;
    long counter, counter_min;

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
    mod[52].z = 117.; 
    mod[53].z = 119.;
 
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
    mod[52].velocity = 1.12e9; 
    mod[53].velocity = 1.08e9;


    
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
