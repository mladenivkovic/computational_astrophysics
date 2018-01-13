#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include "commons.h"


double get_flux(int i, double *some_u);



//=============================
void get_timestep()
//=============================
{
  //--------------------------------------------
  // Compute the next timestep according to CFL
  //--------------------------------------------

  dt = courant_factor*fabs(dx / v);

  // check if you're jumping over output time step
  if ( t + dt >= t_out[t_out_step] ){
    dt = t_out[t_out_step] - t;
  }

  // if (verbose) {
  //   printf("dt =%6g, t= %6g\n", dt, t+dt);
  // }
}










//================
void advect()
//================
{
  //-----------------------------------------
  // Do the integration of the advection eqn
  //-----------------------------------------
  

  //----------------
  // Preparation 
  //----------------
  double *u_inter = malloc((nx+2)*sizeof(double));

#pragma omp parallel
  {
    // store old values
#pragma omp for
    for (int i = 0; i<nx+2; i++){
      u_old[i] = u[i];
    }


    
    // Loop over cells
#pragma omp for
    for (int i = 1; i<nx+1; i++){
      // u_inter[i] = u_old[1]+get_flux(i, u_old)/dx*dt;
      // u[i] = 0.5*(u_old[i] + dt/dx*get_flux(i, u_inter));
      u[i] = u_old[i] + dt/dx*get_flux(i, u_old);
    }

  }


  //implement periodic boundary condition
  u[0] = u[nx];
  u[nx+1] = u[1];


  free(u_inter);

}





double get_flux(int i, double *some_u)
{
  double f = 0;

  if (v < 0){
    f = v*(some_u[i]-some_u[i+1]);
  }
  else{
    f = v*(some_u[i-1]-some_u[i]);
  }

  return (f);
}


