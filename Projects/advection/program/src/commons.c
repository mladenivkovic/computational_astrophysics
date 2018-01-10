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
int density_profile = 0; 
double t_end = 0;
double v = 1;




// global params
int nx = 0;
double dx = 0;
double t = 0;
double dt = 0;


//---------------
//ARRAYS
//---------------

double *u = 0;

