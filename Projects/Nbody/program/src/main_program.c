//=================================
// Contains main program.
//=================================



#include <stdio.h>         
#include <stdlib.h>
#include <math.h>
#include <omp.h>
#include <time.h>
#include "commons.h"
#include "io.h"
#include "direct_forces.h"
#include "multipole.h"



//=====================================
// functions defined below
//=====================================
void initialise(int argc, char *argv[]);
void set_units();
void get_scales();
void check_root();







//=====================================
int main(int argc, char *argv[])    
//=====================================
{

  if (verbose)
  {
    printf("Started program.\n");
  }
 

  //-----------------------
  // Initialise program
  //-----------------------
  initialise(argc, argv);


  // set up timing
  clock_t start, end;
  double cpu_time_used;

  //-----------------------
  // get direct force calc
  //-----------------------
  if (direct_force) {
    if (verbose) {printf("Started direct force calculation.\n");}

    start = clock();
    get_direct_force();
    end = clock();

    write_output(1);

    if (multipole){
      // clean up forces
      for (int i = 0; i < npart; i++){
        fx[i] = 0;
        fy[i] = 0;
        fz[i] = 0;
      }
    }
  }








  //--------------------------
  // Get multipole force calc
  //--------------------------
  if (multipole) {
    if (verbose) {printf("Started multipole force calculation.\n");}
    printf("Softening for direct force parts: %g\n", softening);
    build_tree();
    
    start = clock();
    // Calculate multipole stuff
#pragma omp parallel
    {
      int nthreads = omp_get_num_threads();
#pragma omp master 
      { printf("Entered parallel region with %d threads.\nCalculating multipoles.\n", nthreads); }

#pragma omp for
      for (int i = 0; i<8; i++){
        get_multipole(i);
      }

#pragma omp master 
      { 
        printf("Calculating multipole forces.\n");
        printf("Using multipole order %d, bucket=%d, theta_max=%g\n", multipole_order, ncellpartmax, theta_max); 
      }

#pragma omp for
      // calculate actual forces
      for (int i = 0; i<8; i++){
        calculate_multipole_forces(i);
      }
    }

    end = clock();


    // write_cellparticles();

    write_output(2);
    // check_root();
  }





  //------------------
  // write run info
  //------------------
  cpu_time_used = (double) (end - start) /  CLOCKS_PER_SEC;
  write_info(cpu_time_used);



  if (verbose) { printf("I'm finished!\n"); }
  return(0);

}










//=======================================
void initialise(int argc, char *argv[])
//=======================================
{
  
  //----------------------
  // Set everything up.
  //----------------------


  //read parameters
  readparams(argc, argv);

  //read data
  readdata(argv);

  //transform read data to units
  set_units();

  // get softening
  get_softening();

}












//==================
void set_units()
//==================
{

  //-------------------------------------------------------
  // Sets units of read data to chosen units
  //--------------------------------------------------------
  

  //------------------------
  // first compute scales
  //------------------------

  get_scales();

  if (verbose) {
    printf("Found scales: scale_l = %g, scale_m = %g, scale_t = %g\n", scale_l, scale_m, scale_t  );
  }



  //-----------------------------------
  // translate values to code units
  //-----------------------------------

  if (verbose) { printf("Setting units\n"); }

  double fmass = 1.0/scale_m;
  double flen = 1.0/scale_l;
  double fvel = scale_t/scale_l;
  
  double *phys_arrays[8] = {m, x, y, z, r, vx, vy, vz};
  double factor[8] = {fmass, flen, flen, flen, flen, fvel, fvel, fvel};

  for (int array = 0; array < 8; array++){
    for (int part = 0; part < npart; part++){
      phys_arrays[array][part] *= factor[array];
    }
  }

  double mtot = 0;
  for (int part = 0; part < npart; part++){
    mtot += m[part];
  }
  printf("Total mass after translation %g\n", mtot);


}












//=========================
void get_scales()
//=========================
{

  // ==================================================
  // Computes the scales for code units.
  //
  // a_SI = scale_a * a_CODE
  // Units of choice:
  // set G_CODE = 1
  // set M_tot_CODE = 1 => scale_m = M_tot_SI
  // set R_max_CODE = 1 => scale_l = R_max_SI
  // => scale_t = [ scale_l**3 / (G * scale_m) ]^(1/2) 
  // ==================================================
  

  if (verbose) { printf("Computing unit scales\n"); }

  //-----------------
  // Get mass scale
  //-----------------
  
  double mtot = 0;


  for (int i = 0; i<npart; i++)
  {
    mtot += m[i];
  }

  scale_m = mtot;






  //--------------------
  // Get length scale
  //--------------------
  
  double rmax=0;


  if (scale_cube){
    for (int i = 0; i<npart; i++)
    {
      r[i] = sqrt( pow(x[i], 2) + pow(y[i], 2) + pow(z[i], 2));
      
      if (fabs(x[i])>rmax)
      {
        rmax = fabs(x[i]);
      }
    }
  }
  else {
    for (int i = 0; i<npart; i++)
    {

      r[i] = sqrt( pow(x[i], 2) + pow(y[i], 2) + pow(z[i], 2));
      
      if (r[i]>rmax)
      {
        rmax = r[i];
      }

    }
  }

  scale_l = rmax;

  // set boxlen
  if (scale_cube){
    boxlen = 2;
  }
  else{
    boxlen = 2;
  }




  //--------------------
  // Get time scale
  //--------------------

  scale_t = sqrt(pow(scale_l,3) / (G * scale_m));

}







//===================
void check_root()
//===================
{

  //---------------------------------------------
  // Prints the properties of the root nodes.
  //---------------------------------------------

  printf("\n");
  printf("CHECKING ROOT VALUES\n");

  double totmass = 0;

  for (int i = 0; i<8; i++){
    node * n = cells[i];
    printf("\n");
    printf("=============================\n");
    printf("ROOT %d of level %d\n", n->cellindex, n->level);
    printf("=============================\n");
    printf("has parent: %d\nChildren: ", n->parent);
    for (int j = 0; j<8; j++){
      printf("%3d ", n->child[j]);
    }
    printf("\nCell center: ");
    for (int j = 0; j<3; j++){
      printf("%7g ", n->center[j]);
    }
    printf("\nCenter of mass: ");
    for (int j = 0; j<3; j++){
      printf("%7g ", n->centre_of_mass[j]);
    }
    printf("\nGroup size: %7g np: %5d mass: %7g\n", n->groupsize, n->np, n->mass);
    printf("Multipole vector: ");
    for (int j = 0; j<3; j++){
      printf("%7g ", n->multip_vector[j]);
    }
    printf("\nMultipole matrix:\n");
    for (int j = 0; j<3; j++){
      for (int k = 0; k<3; k++){
        printf("%7g ", n->multip_matrix[j][k]);
      }
      printf("\n");
    }
  printf("\n");


  totmass += n->mass;
  }




  printf("Root total mass: %7g\n\n", totmass);

}





