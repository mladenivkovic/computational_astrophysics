//====================================
// Define common global variables
// Add one useless function
//====================================


//---------------
// VARIABLES
//---------------



//runtime vars
extern int verbose;         //verbose = false by default
extern double f_softening;  // softening factor 
extern int direct_force;    // do direct force calculation
extern int calc_potential;  // do potential calculation



// global params
extern int npart;
extern int ngaspart;
extern int nstarpart;

extern double softening;




// units
extern double scale_m;
extern double scale_t;
extern double scale_l;

// a_SI = scale_a * a_CODE
// Units of choice:
// set G_CODE = 1
// set M_tot_CODE = 1 => scale_m = M_tot_SI
// set R_max_CODE = 1 => scale_l = R_max_SI
// => scale_t = [ scale_l**3 / (G * scale_m) ]^(1/2) 



//---------------
// Arrays
//---------------

extern double *m;
extern double *x;
extern double *y;
extern double *z;
extern double *vx;
extern double *vy;
extern double *vz;
extern double *softening_data;
extern double *potential_data;


extern double *r;  // distance from origin
extern double *fx; // force in x direction
extern double *fy; // force in y direction
extern double *fz; // force in z direction
extern double *phi;// potential


//---------------
// constants
//---------------


extern const double pi;
extern const double G;
