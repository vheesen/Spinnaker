#include<stdio.h>
#include<string.h>
#include<stdlib.h>
#include<math.h>
#include<unistd.h>

#define DEBUG 0

#ifndef DEC_H
#define DEC_H

/*Global integer variables*/
int i, j;

/*Definitions for the parameters read from file*/
char parameter_file_name[25];
int grid_size, nu_channel, grid_delta, first_data_point_at_0kpc;
int epsilon, normalize_intensities;
double z_halo, FWHM_effective_beam;
double nu_1, nu_2, nu_3, nu_4;
int mode;
double gamma_in, rad_field, V0;
int velocity_field;
double h_V;
int adiabatic_losses;
double D0, mu_diff;
int galaxy_mode;
double z1, B0, B1, h_B1, h_B2;
int model, initialize_model, model_north, update_model;
double xi;
double beta, R0;

/* Output file */
double intensity_nu1[402], intensity_nu2[402], intensity_nu3[402], intensity_nu4[402];


/*Physical constants used*/
double pi, c_light, sigma_t, m_electron, e_elem, kpc;


/* Definitions for the cosmic ray propagation code */
double delta_z, z_halo_kpc, delta_nu, delta_nu_factor;
double u_CMB, b;
double B_field[402], u_B[402];
double v_z[402], t_ad[402], t_ad_r[402], t_ad_v[402];
double i_syn_spec[7][402];


/* Modifications for the jet model*/
double factor_model;
int model, model_north, update_model, number_of_data_points;


/** Definitions for the Runge-Kutta**/
struct ODE
{
    double h;
    double x, y1, y2, y3;
    double dy1_dx, dy2_dx, dy3_dx;
    double k1_y1, k1_y2, k1_y3;
    double k2_y1, k2_y2, k2_y3;
    double k3_y1, k3_y2, k3_y3;
    double k4_y1, k4_y2, k4_y3;
};
struct ODE rk;

/** Definitions for the two-dimentionsal grid**/
struct grid_1d
{
    double z, delta_z;
    double E;
    double y, N;
    double nu;
    double t_adv, gamma, dN_dE, gamma_rk, dN_dE_rk;
    double k1_N, k2_N, k3_N, k4_N;
    int integrate_true;
};

struct grid_1d cr[402][402];


/* Definitions for parameter file */ 
struct float_parameter 
{
    char string_search[100];
    double value;
    char string_file[100];
};

struct int_parameter
{
    char string_search[100];
    int value;
    char string_file[100];
};

/* Definitions for calculation of synchrotron intensities */ 
struct syn_emission
{
    double f, g;
    double x;
};

struct syn_emission syn[60];

/* Definitions for the jet model */ 
struct model
{
    double z;
    double B_field;
    double velocity;
    double radius;
    double epsilon, kappa, alpha;
    int ii;
};

struct model mod[86], mod2[60];


#endif
