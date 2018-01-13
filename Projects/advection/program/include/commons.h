//====================================
// Define common global variables
// Add one useless function
//====================================


//---------------
// VARIABLES
//---------------



//runtime vars
extern int verbose;         //verbose = false by default
extern int density_profile; // which density profile to use 
                            // 0 : step function
                            // 1 : linear step
                            // 2 : gauss
extern int method;          // which method to use
                            // 0 : piecewise constant
                            // 1 : 
extern double t_end;        // do multiplole calculation



// global params
extern int nx;          // number of cells
extern double dx;       // cell width
extern double t;        // current time
extern double dt;       // time step
extern double t_out[];  // next output time 
extern int t_out_step;  // step for output time
extern double v;        // global velocity
extern double courant_factor;


//---------------
//ARRAYS
//---------------

extern double *u;
extern double *u_old;


