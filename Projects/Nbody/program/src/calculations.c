//=========================================
// Includes the functions for calculations
//=========================================

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "commons.h"


//=============================
// Functions not in header
//=============================
void get_scales();




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

  r = (double *)calloc(npart, sizeof(double));

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

 
  double G = 6.674083e-11; //gravitational constant

  scale_t = sqrt(pow(scale_l,3) / (G * scale_m));

}
