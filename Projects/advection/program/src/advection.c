#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include "commons.h"


void get_interface_values(double *rho_left, double *rho_right, int i);
double get_flux_inter(int i);
double flux_final(int i, double *u_inter);



//=============================
void get_timestep()
//=============================
{
  //--------------------------------------------
  // Compute the next timestep according to CFL
  //--------------------------------------------

  dt = 0.1*fabs(dx / v);

  // check if you're jumping over output time step
  if ( t + dt > t_out ){
    dt = t_out - t;
  }

  if (verbose) {
    printf("dt =%g, t= %g\n", dt, t+dt);
  }
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
  double *flux_inter = malloc((nx+2)*sizeof(double));

  // store old values
  for (int i = 0; i<nx+2; i++){
    u_old[i] = u[i];
  }

  // compute intermediate fluxes
  for (int i = 1; i<=nx; i++){
    flux_inter[i] = get_flux_inter(i);
  }

  // apply periodicity
  flux_inter[0] = flux_inter[nx];
  flux_inter[nx+1] = flux_inter[1];



  double rho_left, rho_right;
  double rho_left_star, rho_right_star;
  
  // Loop over cells
  for (int i = 1; i<=nx; i++){
    // get_interface_values(&rho_left, &rho_right, i);
    // rho_left_star = rho_left + v * dt/2 *(rho_left - rho_right) / dx;
    // rho_right_star = rho_right + v * dt/2 *(rho_left - rho_right) / dx;
    // u[i] = u_old[i] - v * dt / dx * (rho_right_star - rho_left_star);
    // u_inter[i] = u[i] + (dt/2 * flux_inter[i]/dx);
    u_inter[i] = u_old[i] + 0.5*dt/dx * get_flux_inter(i);
    u[i] = u_old[i] + dt/dx*flux_final(i, u_inter);
  }

  //implement periodic boundary condition
  u[0] = u[nx];
  u[nx+1] = u[1];


  free(u_inter);
}





//========================================================================
void get_interface_values(double *rho_left, double *rho_right, int i)
//========================================================================
{
  //--------------------------
  // Get the interface values
  //--------------------------
  
  if (v < 0){
    *rho_left = u_old[i];
    *rho_right = u_old[i+1];
  }
  else{
    *rho_left = u_old[i-1];
    *rho_right = u_old[i];
  }

  return;
}



double get_flux_inter(int i)
{
  double f = 0;

  if (v < 0){
    f = v*(u_old[i]-u_old[i+1]);
  }
  else{
    f = v*(u_old[i-1]-u_old[i]);
  }

  return (f);
}


double flux_final(int i, double *u_inter)
{
  double f = 0;

  if (v < 0){
    f = v*(u_inter[i]-u_inter[i+1]);
  }
  else{
    f = v*(u_inter[i-1]-u_inter[i]);
  }

  return (f);
}
