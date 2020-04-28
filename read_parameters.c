#include "dec.h"
#include "cosmic_rays.h"
#include "runge_kutta.h"
#include "read_parameters.h"
#include "synchrotron.h"
#include "jet.h"

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

    printf("Input parameter file = '%s'\n", parameter_file_name);

/* The old version with a static parameter file name    
    f=fopen("./parameters", "r"); */
    f=fopen(parameter_file_name, "r");
    if (f == NULL)
    {
	printf("Cannot open parameter file. Stop.\n");
	exit(0);
    }

/* Include some fundamental constants */
    c_light = 2.99792458e10;//[cm s^-1]
    kpc = 3.085677582e21;//[cm]
    pi = 3.141565492;
    sigma_t = 6.6524616e-25;//[cm^2]
    m_electron = 9.1093897e-28;//[g]
    e_elem = 1.60217733e-19;

/*--------------------------------------------------------------------------*/
/*Setup of the 2-dimensional grid*/
/*--------------------------------------------------------------------------*/
    strcpy(parami.string_search, "grid_size");
    grid_size = get_int_parameter (parami, f);
/*--------------------------------------------------------------------------*/
    strcpy(paramf.string_search, "nu_channel");
    nu_channel = get_float_parameter(paramf, f); 
/*--------------------------------------------------------------------------*/
    strcpy(parami.string_search, "grid_delta");
    grid_delta = get_int_parameter (parami, f);
/*--------------------------------------------------------------------------*/
    strcpy(paramf.string_search, "z_halo");
    z_halo_kpc = get_float_parameter(paramf, f);
/*--------------------------------------------------------------------------*/
    strcpy(parami.string_search, "first_data_point_at_0kpc");
    first_data_point_at_0kpc = get_int_parameter (parami, f);
/*--------------------------------------------------------------------------*/
    strcpy(parami.string_search, "normalize_intensities");
    normalize_intensities = get_int_parameter (parami, f);
/*--------------------------------------------------------------------------*/
/*Output */
/*--------------------------------------------------------------------------*/
    strcpy(paramf.string_search, "nu_1");
    nu_1 = get_float_parameter(paramf, f); 
/*--------------------------------------------------------------------------*/
    strcpy(paramf.string_search, "nu_2");
    nu_2 = get_float_parameter(paramf, f); 
/*--------------------------------------------------------------------------*/
    strcpy(paramf.string_search, "nu_3");
    nu_3 = get_float_parameter(paramf, f); 
/*--------------------------------------------------------------------------*/
    strcpy(paramf.string_search, "nu_4");
    nu_4 = get_float_parameter(paramf, f); 
/*--------------------------------------------------------------------------*/
    strcpy(parami.string_search, "mode");
    mode = get_int_parameter (parami, f);
/*--------------------------------------------------------------------------*/
    strcpy(parami.string_search, "epsilon");
    epsilon = get_int_parameter (parami, f);
/*--------------------------------------------------------------------------*/
    strcpy(paramf.string_search, "FWHM_effective_beam");
    FWHM_effective_beam = get_float_parameter (paramf, f);
/*--------------------------------------------------------------------------*/
/*Setup of the advection and diffusion model*/
/*--------------------------------------------------------------------------*/
    strcpy(paramf.string_search, "gamma_in");
    gamma_in = get_float_parameter(paramf, f);
/*--------------------------------------------------------------------------*/
    strcpy(paramf.string_search, "rad_field");
    rad_field = get_float_parameter(paramf, f); 
/*--------------------------------------------------------------------------*/
    strcpy(paramf.string_search, "V0");
    V0 = get_float_parameter(paramf, f); 
/*--------------------------------------------------------------------------*/
    strcpy(parami.string_search, "velocity_field");
    velocity_field = get_int_parameter (parami, f);
/*--------------------------------------------------------------------------*/
    strcpy(paramf.string_search, "h_V");
    h_V = get_float_parameter(paramf, f); 
/*--------------------------------------------------------------------------*/
    strcpy(parami.string_search, "adiabatic_losses");
    adiabatic_losses = get_int_parameter(parami, f); 
/*--------------------------------------------------------------------------*/
    strcpy(paramf.string_search, "D0");
    D0 = get_float_parameter(paramf, f);  
/*--------------------------------------------------------------------------*/
    strcpy(paramf.string_search, "mu_diff");
    mu_diff = get_float_parameter(paramf, f); 
/*--------------------------------------------------------------------------*/
    strcpy(paramf.string_search, "V_rot");
    V_rot = get_float_parameter(paramf, f); 
/*--------------------------------------------------------------------------*/
/*Magnetic field setup*/
/*--------------------------------------------------------------------------*/
    strcpy(parami.string_search, "galaxy_mode");
    galaxy_mode = get_int_parameter(parami, f); 
/*--------------------------------------------------------------------------*/
    strcpy(paramf.string_search, "z1");
    z1 = get_float_parameter(paramf, f); 
/*--------------------------------------------------------------------------*/
    strcpy(paramf.string_search, "B0");
    B0 = get_float_parameter(paramf, f); 
/*--------------------------------------------------------------------------*/
    strcpy(paramf.string_search, "B1");
    B1 = get_float_parameter(paramf, f); 
/*--------------------------------------------------------------------------*/
    strcpy(paramf.string_search, "h_B1");
    h_B1 = get_float_parameter(paramf, f); 
/*--------------------------------------------------------------------------*/
    strcpy(paramf.string_search, "h_B2");
    h_B2 = get_float_parameter(paramf, f); 
/*--------------------------------------------------------------------------*/
     strcpy(paramf.string_search, "beta");
    beta = get_float_parameter(paramf, f); 
/*--------------------------------------------------------------------------*/
    strcpy(paramf.string_search, "R0");
    R0 = get_float_parameter(paramf, f); 
/*--------------------------------------------------------------------------*/
/*#Use a magnetic field model (needs edit of the source files)*/
/*--------------------------------------------------------------------------*/
    strcpy(parami.string_search, "model");
    model = get_int_parameter(parami, f); 
/*--------------------------------------------------------------------------*/
    strcpy(parami.string_search, "initialze_model");
    initialize_model = get_int_parameter(parami, f); 
/*--------------------------------------------------------------------------*/
    strcpy(parami.string_search, "model_north");
    model_north = get_int_parameter(parami, f); 
/*--------------------------------------------------------------------------*/
    strcpy(parami.string_search, "update_model");
    update_model = get_int_parameter(parami, f); 
/*--------------------------------------------------------------------------*/
    strcpy(paramf.string_search, "xi");
    xi = get_float_parameter(paramf, f); 
/*--------------------------------------------------------------------------*/
/****************************************************************************/
    z_halo = z_halo_kpc * kpc;
    
    if (nu_channel < 2)
    {
	printf("You need at least to separate frequencies to compute\n");
	printf("the spectral index: nu_channel has to be >= 2! Stop.\n");
	exit (0); 
    }

    if (nu_2 <= nu_1)
	{
	    printf ("Error: nu_2 <= nu_1. Stop.\n");
	    exit(0);
	}

/* For diffusion the velocity field should be set to zero
 so that the CRe number density is calculated correctly */
    if (mode != 1)
        velocity_field = 0;
    
    
    if (model == 1 || initialize_model == 1)
    {
        
        if (model_north == 1)
            number_of_data_points = set_radius_north ();
        else
            number_of_data_points = set_radius_south ();
        if (update_model == 1)
            read_intensity_model ();
        if  (mod[number_of_data_points].z > z_halo_kpc)
        {
            printf("Halo size is smaller than model size. Stop.\n");
            printf("Halo_size = %g kpc, Model size = %g kpc\n",  z_halo_kpc, mod[number_of_data_points].z);
            exit(0);
        }
    }

    grid_size = (int) 2.0 * grid_size;

    R0 = R0 * kpc;
    h_grav = 4.0 * R0 / kpc;
    
    fclose(f);

}
