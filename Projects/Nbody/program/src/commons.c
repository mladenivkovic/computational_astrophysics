//====================================
// Define common global variables
// Add one useless function
//====================================


#include <stdio.h>
#include <stdarg.h>


//---------------
// VARIABLES
//---------------

//runtime vars
int verbose = 0; //verbose = false by default

int npart = 0;
int ngaspart = 0;
int nstarpart = 0;



// units
double scale_m = 0;
double scale_t = 0;
double scale_l = 0;

// a_SI = scale_a * a_CODE
// Units of choice:
// set G_CODE = 1
// set M_tot_CODE = 1 => scale_m = M_tot_SI
// set R_max_CODE = 1 => scale_l = R_max_SI
// => scale_t = [ scale_l**3 / (G * scale_m) ]^(1/2) 



//---------------
//ARRAYS
//---------------

double *m = 0;
double *x = 0;
double *y = 0;
double *z = 0;
double *vx = 0;
double *vy = 0;
double *vz = 0;
double *softening_data = 0;
double *potential_data = 0;


double *r = 0;  // distance from origin
double *fx = 0; // force in x direction
double *fy = 0; // force in y direction
double *fz = 0; // force in z direction
double *phi = 0;// potential




//---------------
// FUNCITONS
//---------------

