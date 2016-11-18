int set_radius_north ( void );
int set_radius_south ( void );
double radius (double z);
double dr_dz (double z);
double magnetic_field (double z);
void set_interpolate_values (double z, int ii);
double interpolated_value (double value, double value_low, double value_high, int ii, int ii_mod);
void read_intensity_model (void);
void read_magnetic_field_model (void);

