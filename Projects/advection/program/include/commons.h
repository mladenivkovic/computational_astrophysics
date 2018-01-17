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
                            // 1 : piecewise linear
                            // 2 : piecewise linear with minmod slope limiter
                            // 3 : piecewise linear with VanLeer slope limiter

extern double t_end;        // time when to end
extern int use_2d;          // use 2d instead of 1d advection

extern double u;        // global velocity in x direction
extern double v;        // global velocity in y direction


// global params
extern int nx;          // number of cells
extern int ny;          // number of cells
extern double dx;       // cell width
extern double dy;       // cell height
extern double t;        // current time
extern double dt;       // time step
extern double t_out[];  // next output time 
extern int t_out_step;  // step for output time
extern double courant_factor; // factor for courant condition; Must be <= 1


//---------------
//ARRAYS
//---------------

extern double *rho;
extern double *rho_old;

extern double **rho2d;
extern double **rho2d_old;


