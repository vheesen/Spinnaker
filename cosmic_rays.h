double dN_dz (double z, double N, double E, double gamma, double dN_dE);

double d2N_dz2 (double z, double N, double y, double E, 
		double gamma, double dN_dE);

void gamma_cr (void);

void dN_dE (void);

void adiabatic (void);

struct grid_1d setup_initial_grid (void);

void output_file (int i_max);

void output_stdout (int i_max);

double interpolate_frequency (double value, double value_low, double value_high, double nu, int jj);
