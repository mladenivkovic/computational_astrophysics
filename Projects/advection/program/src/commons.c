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
int verbose = 0;          // verbose = false by default
int density_profile = -1; 
int method = -1;
double t_end = 0;
int use_2d = 0;
double u = 1;
double v = 0;




// global params
int nx = 0;
int ny = 0;
double dx = 0;
double dy = 0;
double t = 0;
double dt = 0;
double t_out[7] = {1, 2, 5, 10, 20, 50, 100};
// double t_out[11] = {0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1};
// double t_out[17] = {0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1, 2, 5, 10, 20, 50, 100};
// double t_out[16] = {3.5, 3.6, 3.7, 3.8, 3.9, 4.0, 4.1, 4.2, 4.3, 4.4, 4.5, 4.6, 4.7, 4.8, 4.9, 5};
int t_out_step = 0;
double courant_factor = 1;



//---------------
//ARRAYS
//---------------

double *rho = 0;
double *rho_old = 0;

double **rho2d = 0;
double **rho2d_old = 0;

