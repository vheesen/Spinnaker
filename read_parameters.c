#include "dec.h"
#include "cosmic_rays.h"
#include "runge_kutta.h"
#include "read_parameters.h"
#include "synchrotron.h"
#include "magnetic_field.h"
#include "adiabatic.h"

/*************************************************************************/
float get_int_parameter(struct int_parameter param, FILE *f)
{

    
    char line[100];
    char string_value[100];
    int string_equal = *"=";
    int string_garden = *"#";
    int ii = 0;
    char *start;
    char *stop;


   
/*Check for comment lines starting with '#'*/
    stop = line;
    while (stop - line == 0)
    {
	fgets(line, sizeof(line), f);
	start = strchr(line, string_equal);
	stop = strchr(line, string_garden);
	if (ii > 10)
	{
	    printf ("Error while reading comment lines in parameter file. Stop.\n");
	    exit (0);//something wrong
	    
	}
	ii++;
    }
    strncpy(string_value, start + 1, stop - start);
    strncpy(param.string_file, line, start - line - 1);
    param.value = atof (string_value);

/* #if DEBUG == 2 */
/*     printf("string_file = %s\n", param.string_file); */
/*     printf("param.value = %g\n", param.value); */
/*     printf("string start = %s stop = %s\n", start, stop); */
/*     printf("string stop - start %d\n", stop - start); */
/* //    printf ("c_light = %s", value); */
/*     exit (0); */
/* #endif */
//    }

#if DEBUG == 1
    printf("%s = %d\n", param.string_search, param.value);
#endif

/*     if (strstr(param.string_search, param.string_file) == NULL) */
/*     { */
/* 	printf("Error while reading parameter %s.\n",  */
/* 	       param.string_search); */
/* 	printf("Tried to read %s in parameter file.\n", param.string_file); */
/* 	printf("Searched for %s\n", param.string_search); */
/* 	printf("Check, whether the string in the parameter file  */
/* is correct!\n"); */
/* 	exit(0); */
/*     } */

//    fgets(param.string_file, sizeof(param.string_file), f);
//    sscanf(param.string_file, "%g", &param.value);



    return param.value;


}
/*************************************************************************/
float get_float_parameter(struct float_parameter param, FILE *f)
{
    char line[100];
    char string_value[100];
    int string_equal = *"=";
    int string_garden = *"#";
    int ii = 0;
    char *start;
    char *stop;
    
/*Check for comment lines starting with '#'*/
    stop = line;
    while (stop - line == 0)
    {
	fgets(line, sizeof(line), f);
	start = strchr(line, string_equal);
	stop = strchr(line, string_garden);
	if (ii > 10)
	{
	    printf ("Error while reading comment lines in parameter file. Stop.\n");
	    exit (0);//something wrong
	    
	}
	ii++;
    }
  
    strncpy(string_value, start + 1, stop - start);
    strncpy(param.string_file, line, start - line);
    param.value = atof (string_value);

/* #if DEBUG == 1 */
/*     printf("ii = %d", ii); */
/*     printf("string_file = %s\n", param.string_file); */
/*     printf("param.value = %g\n", param.value); */
/*     printf("string start = %s stop = %s\n", start, stop); */
/*     printf("string stop - start %d\n", stop - start); */
/*     printf("param.string_search = %s \n", param.string_search); */
/*     printf("param.string_file = %s \n", param.string_file); */
/*     exit(0); */
/* //    printf ("c_light = %s", value); */
    
/* #endif */
//    }
/*     if (strstr(param.string_search, line) == NULL) */
/*     { */
/* 	printf("Error while reading parameter %s.\n",  */
/* 	       param.string_search); */
/* 	printf("Tried to read %s in parameter file.\n", param.string_file); */
/* 	printf("Searched for %s\n", param.string_search); */
/* 	printf("output string search = %d\n", */
/* 	       strstr(param.string_search, param.string_file)); */
/* 	printf("Check, whether the string in the parameter file  */
/* is correct!\n"); */
/* 	exit(0); */
/*     } */

//    fgets(param.string_file, sizeof(param.string_file), f);
//    sscanf(param.string_file, "%g", &param.value);

#if DEBUG == 1
    printf("%s = %g\n", param.string_search, param.value);
#endif


    return param.value;


}
/**************************************************************************/
/* Read parameters from parameter file */
void read_parameters(void)
{
    FILE *f;

/* This stuff is no needed. However, it doesn't work without on my machine!*/
/*-------------------------------------------------------------------------*/
    char line[100];
    char dummy[100];
/*-------------------------------------------------------------------------*/
    struct int_parameter parami;
    struct float_parameter paramf;

    f=fopen("./parameters", "r");
    if (f == NULL)
    {
	printf("Cannot open parameter file. Stop.");
	exit(0);
    }


/*--------------------------------------------------------------------------*/
    strcpy(paramf.string_search, "c_light");
    c_light = get_float_parameter(paramf, f); 
/*--------------------------------------------------------------------------*/
    strcpy(paramf.string_search, "parsec");
    parsec = get_float_parameter(paramf, f); 
/*--------------------------------------------------------------------------*/
    strcpy(paramf.string_search, "pi");
    pi = get_float_parameter(paramf, f); 
/*--------------------------------------------------------------------------*/
    strcpy(paramf.string_search, "sigma_t");
    sigma_t = get_float_parameter(paramf, f); 
/*--------------------------------------------------------------------------*/
    strcpy(paramf.string_search, "m_electron");
    m_electron = get_float_parameter(paramf, f); 
/*--------------------------------------------------------------------------*/
    strcpy(paramf.string_search, "e_elem");
    e_elem = get_float_parameter(paramf, f); 
/*--------------------------------------------------------------------------*/
    strcpy(parami.string_search, "grid_size");
    grid_size = get_int_parameter (parami, f);
/*--------------------------------------------------------------------------*/
    strcpy(parami.string_search, "grid_delta");
    grid_delta = get_int_parameter (parami, f);
/*--------------------------------------------------------------------------*/
    strcpy(paramf.string_search, "z_halo");
    z_halo_parsec = get_float_parameter(paramf, f);
/*--------------------------------------------------------------------------*/
    strcpy(parami.string_search, "mode_0");
    mode_0 = get_int_parameter (parami, f);
/*--------------------------------------------------------------------------*/
    strcpy(parami.string_search, "mode_1");
    mode_1 = get_int_parameter (parami, f);
/*--------------------------------------------------------------------------*/
    strcpy(parami.string_search, "mode_2");
    mode_2 = get_int_parameter (parami, f);
/*--------------------------------------------------------------------------*/
    strcpy(paramf.string_search, "D0");
    D0 = get_float_parameter(paramf, f);  
/*--------------------------------------------------------------------------*/
    strcpy(paramf.string_search, "D1");
    D1 = get_float_parameter(paramf, f);  
/*--------------------------------------------------------------------------*/
    strcpy(paramf.string_search, "D2");
    D2 = get_float_parameter(paramf, f);  
/*--------------------------------------------------------------------------*/
    strcpy(paramf.string_search, "mu_diff");
    mu_diff = get_float_parameter(paramf, f); 
/*--------------------------------------------------------------------------*/
    strcpy(paramf.string_search, "v0");
    v0 = get_float_parameter(paramf, f); 
/*--------------------------------------------------------------------------*/
    strcpy(paramf.string_search, "v1");
    v1 = get_float_parameter(paramf, f); 
/*--------------------------------------------------------------------------*/
    strcpy(paramf.string_search, "v2");
    v2 = get_float_parameter(paramf, f); 
/*--------------------------------------------------------------------------*/
    strcpy(paramf.string_search, "nu_low");
    nu_low = get_float_parameter(paramf, f); 
/*--------------------------------------------------------------------------*/
    strcpy(paramf.string_search, "nu_high");
    nu_high = get_float_parameter(paramf, f); 
/*--------------------------------------------------------------------------*/
    strcpy(paramf.string_search, "nu_channel");
    nu_channel = get_float_parameter(paramf, f); 
/*--------------------------------------------------------------------------*/
    strcpy(paramf.string_search, "gamma_in");
    gamma_in = get_float_parameter(paramf, f);
/*--------------------------------------------------------------------------*/
    strcpy(paramf.string_search, "B0");
    B0 = get_float_parameter(paramf, f); 
/*--------------------------------------------------------------------------*/
    strcpy(paramf.string_search, "rad_field");
    rad_field = get_float_parameter(paramf, f); 
/*--------------------------------------------------------------------------*/
    strcpy(paramf.string_search, "h_B0");
    h_B0 = get_float_parameter(paramf, f); 
/*--------------------------------------------------------------------------*/
    strcpy(paramf.string_search, "z0");
    z0 = get_float_parameter(paramf, f); 
/*--------------------------------------------------------------------------*/
    strcpy(paramf.string_search, "h_B1");
    h_B1 = get_float_parameter(paramf, f); 
/*--------------------------------------------------------------------------*/
    strcpy(paramf.string_search, "z1");
    z1 = get_float_parameter(paramf, f); 
/*--------------------------------------------------------------------------*/
    strcpy(paramf.string_search, "h_B2");
    h_B2 = get_float_parameter(paramf, f); 
/*--------------------------------------------------------------------------*/
    strcpy(parami.string_search, "power_law");
    power_law = get_int_parameter (parami, f);
/*--------------------------------------------------------------------------*/
    strcpy(paramf.string_search, "DF0");
    DF0 = get_float_parameter(paramf, f); 
/*--------------------------------------------------------------------------*/
    strcpy(paramf.string_search, "beta0");
    beta0 = get_float_parameter(paramf, f); 
/*--------------------------------------------------------------------------*/
    strcpy(paramf.string_search, "beta1");
    beta1 = get_float_parameter(paramf, f); 
/*--------------------------------------------------------------------------*/
    strcpy(paramf.string_search, "beta2");
    beta2 = get_float_parameter(paramf, f); 
/*--------------------------------------------------------------------------*/
    strcpy(paramf.string_search, "z_red");
    z_red = get_float_parameter(paramf, f); 
/*--------------------------------------------------------------------------*/
    strcpy(parami.string_search, "adiabatic_losses");
    adiabatic_losses = get_int_parameter(parami, f); 
/*--------------------------------------------------------------------------*/
    strcpy(paramf.string_search, "t_ad");
    t_ad = get_float_parameter(paramf, f); 
/*--------------------------------------------------------------------------*/
    strcpy(parami.string_search, "model");
    model = get_int_parameter(parami, f); 
/*--------------------------------------------------------------------------*/
    strcpy(parami.string_search, "model_north");
    model_north = get_int_parameter(parami, f); 
/*--------------------------------------------------------------------------*/
    strcpy(parami.string_search, "galaxy_mode");
    galaxy_mode = get_int_parameter(parami, f); 
/*--------------------------------------------------------------------------*/
    strcpy(paramf.string_search, "B1");
    B1 = get_float_parameter(paramf, f); 
/*--------------------------------------------------------------------------*/
    strcpy(paramf.string_search, "B2");
    B2 = get_float_parameter(paramf, f); 
/****************************************************************************/

//    printf("magnetic flux = %g\n", magnetic_flux);
    
    if (nu_channel < 2)
    {
	printf("You need at least to separate frequencies to compute\n");
	printf("the spectral index: nu_channel has to be >= 2! Stop.\n");
	exit (0); 
    }

    if (nu_high <= nu_low)
	{
	    printf ("Error: nu_high <= nu_low. Stop.\n");
	    exit(0);
	}

    z_halo = z_halo_parsec * parsec;

    printf("z0=%g\n", z0);

//Convert adiabatic loss time scale from Myr to s    
    t_ad = 3.15e13 * t_ad;
    

    fclose(f);

}
