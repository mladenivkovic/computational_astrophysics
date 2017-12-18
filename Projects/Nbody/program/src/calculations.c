//=========================================
// Includes the functions for calculations
//=========================================

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <omp.h>
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










//=========================
void get_direct_force()
//=========================
{

  //---------------------------------------------
  // Compute the direct force for each particle
  //---------------------------------------------

  if (verbose) { printf("Computing direct forces.\n"); }

  // get softening
  void get_softening();
  get_softening();


  if (verbose) { printf("Found softening: %g\n", softening); }


  double force_fact;
  double rsq;
  double rij;




#pragma omp parallel
  {

    int nthreads = omp_get_num_threads();
#pragma omp master
    {
      if (verbose){
        printf("Started parallel region. Nthreads: %d\n", nthreads);
      }
    }

#pragma omp for     
    for (int i = 0; i < npart; i++)
    {
      for (int j = 0; j < npart; j++)
      {
        if (!(i == j))
        {
          rsq = pow(x[i]-x[j], 2) + pow(y[i]-y[j], 2) + pow(z[i]-z[j], 2);
          rij = sqrt(rsq);
          force_fact = - m[i]*m[j] / sqrt(rsq + pow(softening, 2));
          fx[i] += force_fact * (x[i]-x[j])/rij;
          fy[i] += force_fact * (y[i]-y[j])/rij;
          fz[i] += force_fact * (z[i]-z[j])/rij;
        }
      }
    }

  }
}









//========================================
void get_softening()
//========================================
{
  // ------------------------------------------------
  // Computes the mean interparticle distance and
  // the softening.
  // ------------------------------------------------



  // first get mean interparticle distance
  double rmax = 0;
  
  for (int i = 0; i<npart; i++)
  {
    if (r[i] > rmax){
      rmax = r[i];
    }
  }
 
  double full_volume = 4.0/3.0*pi*pow(rmax,3);
  double particle_volume = full_volume / npart;
  double mean_interparticle_distance = pow(3*particle_volume/(4 * pi), 1/3.);
  softening = f_softening * mean_interparticle_distance;
  

}


