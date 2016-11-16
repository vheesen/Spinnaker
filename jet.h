void set_radius_north ( void );
double radius_north (double z);
double dr_dz_north (double z);
void set_radius_south ( void );
double radius_south (double z);
double dr_dz_south (double z);
double magnetic_field_north (double z);
double velocity_field_north (double z);
double magnetic_field_south (double z);
double velocity_field_south (double z);
void set_interpolate_values_north (double z, int ii);
double interpolated_value (double value, double value_low, double value_high, int ii, int ii_mod);
void set_interpolate_values_south (double z, int ii);
void read_intensity_model (void);
void read_magnetic_field_model (void);

