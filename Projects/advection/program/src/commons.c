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
double v = 1;




// global params
int nx = 0;
double dx = 0;
double t = 0;
double dt = 0;
double t_out[7] = {1, 2, 5, 10, 20, 50, 100};
int t_out_step = 0;
double courant_factor = 1;


//---------------
//ARRAYS
//---------------

double *u = 0;
double *u_old = 0;

