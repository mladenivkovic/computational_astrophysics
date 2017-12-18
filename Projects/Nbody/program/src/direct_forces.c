//=========================================
// Includes the functions for calculations
//=========================================

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <omp.h>
#include "commons.h"







//=========================
void get_direct_force()
//=========================
{

  //---------------------------------------------
  // Compute the direct force for each particle
  //---------------------------------------------


  // get softening
  void get_softening();
  get_softening();


  if (verbose) { printf("Found softening: %g\n", softening); }
  if (verbose) { printf("Computing direct forces.\n"); }

  double force_fact;
  double rsq;




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
          force_fact = - m[i]*m[j] / pow(rsq + pow(softening, 2), 1.5);
          fx[i] += force_fact * (x[i]-x[j]);
          fy[i] += force_fact * (y[i]-y[j]);
          fz[i] += force_fact * (z[i]-z[j]);
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


