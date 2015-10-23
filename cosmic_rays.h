double dN_dz (double z, double N, double E, double gamma, double dN_dE);

double d2N_dz2 (double z, double N, double y, double E, 
		double gamma, double dN_dE);

void gamma_cr (void);

void dN_dE (void);

struct grid_1d setup_initial_grid (void);

void output_file (long i_max);


