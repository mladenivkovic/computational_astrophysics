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




#pragma omp parallel 
  {
    double force_fact;
    double rsq;


#pragma omp master
    {
      if (verbose){
        printf("Started parallel region.\n");
      }
    }

#pragma omp for 
    for (int i = 0; i < npart; i++)
    {
      for (int j = i+1; j < npart; j++)
      {
        double dx = x[i]-x[j];
        double dy = y[i]-y[j];
        double dz = z[i]-z[j];

        rsq = pow(dx, 2) + pow(dy, 2) + pow(dz, 2);
        force_fact = - m[i]*m[j] / pow(rsq + pow(softening, 2), 1.5);
#pragma omp critical
        {
          fx[i] += force_fact * dx;
          fy[i] += force_fact * dy;
          fz[i] += force_fact * dz;
          fx[j] += force_fact * (-dx);
          fy[j] += force_fact * (-dy);
          fz[j] += force_fact * (-dz);
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


