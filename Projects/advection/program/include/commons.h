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
extern double t_end;        // do multiplole calculation



// global params
extern int nx;        // number of cells
extern double dx;     // cell width
extern double t;      // current time
extern double dt;     // time step
extern double t_out;  // next output time 
extern double t_out_step; // step for output time
extern double v;      // global velocity


//---------------
//ARRAYS
//---------------

extern double *u;
extern double *u_old;


