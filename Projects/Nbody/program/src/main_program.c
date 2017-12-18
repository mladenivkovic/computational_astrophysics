//=================================
// Contains main program.
//=================================



#include <stdio.h>         
#include <stdlib.h>
#include <math.h>
#include "commons.h"
#include "io.h"
#include "direct_forces.h"



//=====================================
// functions defined below
//=====================================
void initialise(int argc, char *argv[]);
void set_units();
void get_scales();







//=====================================
int main(int argc, char *argv[])    
//=====================================
{

  if (verbose)
  {
    printf("Started program.\n");
  }
 


  // Initialise program
  initialise(argc, argv);


  // get direct force calc
  if (direct_force) {
    get_direct_force();
  }
    output_direct_force();



  // write run info
  write_info();



  if (verbose)
  {
    printf("I'm finished!\n");
  }
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

  for (int i = 0; i<npart; i++)
  {

    r[i] = sqrt( pow(x[i], 2) + pow(y[i], 2) + pow(z[i], 2));
    
    if (r[i]>rmax)
    {
      rmax = r[i];
    }

  }

  scale_l = rmax;






  //--------------------
  // Get time scale
  //--------------------

  scale_t = sqrt(pow(scale_l,3) / (G * scale_m));

}











