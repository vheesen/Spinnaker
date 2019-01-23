#include "dec.h"
#include "cosmic_rays.h"
#include "runge_kutta.h"
#include "read_parameters.h"
#include "read_parameters2.h"
#include "synchrotron.h"
#include "jet.h"

/*************************************************************************/
/* Read parameters from command line input, given by main as char array */

//reads values between the searched variable and the # character
char regit(char input[], char *buf, char var[])
{
	char *pch = strstr(input,var);
	char temp[100]={'\0'};
	strncpy(temp,pch,50);
	char *chr=strchr(temp,'#');
	int len = strlen(var);
	strncpy(buf,temp+len,chr-temp-len);
}

//returns a float variable
double getfloat(char input[], char *regstr)
{
	char buf[50] = {'\0'};
	double var_d=0; 
	regit(input,buf,regstr);
	sscanf(buf,"%lg",&var_d);
	return var_d;
}

//returns an int variable
int getint(char input[], char *regstr)
{
	char buf[50] = {'\0'};
	int var_i=0;
	regit(input,buf,regstr);
	sscanf(buf,"%d",&var_i);
	return var_i;
}

void read_parameters2(char input[])
{

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
	grid_size = getint(input,"grid_size=");
/*--------------------------------------------------------------------------*/
	nu_channel = getfloat(input,"nu_channel=");
/*--------------------------------------------------------------------------*/
	grid_delta = getint(input,"grid_delta=");
/*--------------------------------------------------------------------------*/
	z_halo_kpc = getfloat(input,"z_halo=");
/*--------------------------------------------------------------------------*/
	first_data_point_at_0kpc = getint(input,"first_data_point_at_0kpc=");
/*--------------------------------------------------------------------------*/
	normalize_intensities = getint(input,"normalize_intensities=");
/*--------------------------------------------------------------------------*/
/*Output */
/*--------------------------------------------------------------------------*/
	nu_1 = getfloat(input,"nu_1="); 
/*--------------------------------------------------------------------------*/
	nu_2 = getfloat(input,"nu_2="); 
/*--------------------------------------------------------------------------*/
	nu_3 = getfloat(input,"nu_3="); 
/*--------------------------------------------------------------------------*/
	nu_4= getfloat(input,"nu_4="); 
/*--------------------------------------------------------------------------*/
	mode = getint(input,"mode=");
/*--------------------------------------------------------------------------*/
	epsilon = getint(input,"epsilon=");
/*--------------------------------------------------------------------------*/
	FWHM_effective_beam= getfloat(input,"FWHM_effective_beam="); 
/*--------------------------------------------------------------------------*/
/*Setup of the advection and diffusion model*/
/*--------------------------------------------------------------------------*/
	gamma_in= getfloat(input,"gamma_in="); 
/*--------------------------------------------------------------------------*/
	rad_field= getfloat(input,"rad_field="); 
/*--------------------------------------------------------------------------*/
	V0= getfloat(input,"V0="); 
/*--------------------------------------------------------------------------*/
	velocity_field = getint(input,"velocity_field=");
/*--------------------------------------------------------------------------*/
	h_V= getfloat(input,"h_V="); 
/*--------------------------------------------------------------------------*/
	adiabatic_losses = getint(input,"adiabatic_losses="); 
/*--------------------------------------------------------------------------*/
	D0= getfloat(input,"D0=");
/*--------------------------------------------------------------------------*/
	mu_diff= getfloat(input,"mu_diff=");
/*--------------------------------------------------------------------------*/
/*Magnetic field setup*/
/*--------------------------------------------------------------------------*/
	galaxy_mode = getint(input,"galaxy_mode="); 
/*--------------------------------------------------------------------------*/
	z1= getfloat(input,"z1=");
/*--------------------------------------------------------------------------*/
 	B0= getfloat(input,"B0=");
/*--------------------------------------------------------------------------*/
	B1= getfloat(input,"B1=");
/*--------------------------------------------------------------------------*/
	h_B1= getfloat(input,"h_B1=");
/*--------------------------------------------------------------------------*/
	h_B2= getfloat(input,"h_B2=");
/*--------------------------------------------------------------------------*/
/*#Use a magnetic field model (needs edit of the source files)*/
/*--------------------------------------------------------------------------*/
	model = getint(input,"model=");
/*--------------------------------------------------------------------------*/
	initialize_model = getint(input,"initialize_model="); 
/*--------------------------------------------------------------------------*/
	model_north = getint(input,"model_north="); 
/*--------------------------------------------------------------------------*/
	update_model = getint(input,"update_model=");  
/*--------------------------------------------------------------------------*/
	xi= getfloat(input,"xi="); 
/*--------------------------------------------------------------------------*/
	beta= getfloat(input,"beta="); 
/*--------------------------------------------------------------------------*/
	R0= getfloat(input,"R0="); 
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
            
 //   fclose(f);

}
