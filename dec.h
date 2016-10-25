#include<stdio.h>
#include<string.h>
#include<stdlib.h>
#include<math.h>

#define DEBUG 0

#ifndef DEC_H
#define DEC_H

/*global integer variables*/
long i, j;
long choice_rk, mode, mode_0, mode_1, mode_2, power_law;
long adiabatic_losses, model, model_north, galaxy_mode;


double pi, c_light, sigma_t, m_electron, e_elem, gamma_in;
double rad_field;//=U_rad/U_B
double B0, B1, B2;//in units of mikro Gauss
double z_red;//Redshift
double t_ad;


/* Definitions for the Cosmic Ray propagation code */
long nu_channel;
double nu_gyro;
double u_b;
double delta_z, z_halo, z_halo_parsec, parsec, delta_nu, delta_nu_factor;
double v0, v1, v2, b;
double diff_0, D0, D1, D2, mu_diff;
double nu_low, nu_high;
double h_B0, h_B1, h_B2, z0, z1, DF0, beta0, beta1, beta2;//Magnetic field setup
double bz[402][402];
double u_CMB;
double B_field[402], u_B[402], v_z[402];
double i_syn_spec[7][402];

/** Definitions for the Runge-Kutta**/
long grid_size, grid_delta;

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

struct grid_1d
{
    double z, delta_z;//in parsec
    double E;//electron energy in units of GeV
    double y, N;
    double nu;//in units of MHz
//    double alpha; //computed spectral index alpha(nu_low, nu_high)
//    double gamma;//computed energy index gamma(nu_low, nu_high)
    double Q;//source term for cosmic ray acceleration
    double alpha, gamma, dN_dE, dgamma_dE, gamma_rk, dN_dE_rk;
    double k1_N, k2_N, k3_N, k4_N;
};

struct grid_1d cr[402][402];


/* parameter file */ 

struct float_parameter 
{
    char string_search[100];
    double value;
    char string_file[100];
};

struct int_parameter
{
    char string_search[100];
    long value;
    char string_file[100];
};

struct syn_emission
{
    double f, g;
    double x;//nu/nu_crit
};

struct syn_emission syn[60];

struct model
{
    double z;
    double B_field;
    double velocity;
    double radius;
    int ii;
    
};

struct model mod[86];


#endif
