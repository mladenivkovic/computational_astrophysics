#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include "commons.h"
#include "io.h"


void set_boundary();
void get_pwlin_flux(int i, double *flux1, double *flux2);
double get_slope(int i);
double get_minmod_slope(int i);
double VanLeer_limiter1(int i);
double VanLeer_limiter2(int i);







//=======================================
void initialise1d()
//=======================================
{
  //----------------------
  // Set everything up.
  //----------------------



  //-------------------------
  // initialize variables
  //-------------------------

  if (nx == 0){
    printf("Got nx = 0, can't work with that. quitting.\n");
    exit(1);
  }

#pragma omp single
  { rho = malloc((nx+4)*sizeof(double)); }
#pragma omp single
  { rho_old = malloc((nx+4)*sizeof(double)); }
#pragma omp single
  { rho_inter = malloc((nx+4)*sizeof(double)); }
#pragma omp single
  { 
    dx = 1.0/((double) nx); 
    u = 1;
  }



  //------------------------------
  // Initialise density profile
  //------------------------------

  if (density_profile == 0){
    //--------------------------------------------
    printf("Using step density profile.\n");
    //--------------------------------------------
#pragma omp for
    for (int i = 0; i<nx; i++){
      if (i*dx <= 0.3){
        rho[i+2] = 1;
      }
      else if (i*dx <= 0.6){
        rho[i+2] = 2;
      }
      else{
        rho[i+2] = 1;
      }
    }
  }
  else if (density_profile == 1){
    //--------------------------------------------
    printf("Using linear step density profile.\n");
    //--------------------------------------------
#pragma omp for
    for (int i = 0; i<nx; i++){
      if (i*dx <= 0.3){
        rho[i+2] = 1;
      }
      else if (i*dx <= 0.6){
        rho[i+2] = 2+1.7*(i*dx-0.3);
      }
      else{
        rho[i+2] = 1;
      }
    }
  }
  else if (density_profile == 2){
    //--------------------------------------------
    printf("Using gauss density profile.\n");
    //--------------------------------------------
#pragma omp for
    for (int i = 0; i<nx; i++){
      rho[i+2] = 1 + exp(-pow((i*dx - 0.5), 2)/0.1);
    }
  }
  else {
#pragma omp single
    {
      printf("Not recognized density_profile = %d\n", density_profile);
    }
  }


  set_boundary();

}














//=============================
void get_timestep1d()
//=============================
{
  //--------------------------------------------
  // Compute the next timestep according to CFL
  //--------------------------------------------

  dt = courant_factor*fabs(dx / u);

  // check if you're jumping over output time step
  if ( t + dt >= t_out[t_out_step] ){
    dt = t_out[t_out_step] - t;
  }

  // if (verbose) {
  //   printf("dt =%6g, t= %6g\n", dt, t+dt);
  // }
}










//================
void advect1d()
//================
{
  //-----------------------------------------
  // Do the integration of the advection eqn
  //-----------------------------------------
  

  //----------------
  // Preparation 
  //----------------

  // store old values
#pragma omp for
  for (int i = 0; i<nx+4; i++){
    rho_old[i] = rho[i];
  }


  //---------------------- 
  // Loop over cells
  //---------------------- 
  
  double c = dt/dx * u;

  if (method == 0){
    //--------------------------------
    // piecewise constant method
    //--------------------------------
#pragma omp for
    for (int i = 2; i<nx+2; i++){
      rho[i] = rho_old[i] - c*(rho_old[i]-rho_old[i-1]);
    }
  }


  else if (method == 1){
    //--------------------------------
    // piecewise linear method
    //--------------------------------
    double slope_l, slope_r;
  
    if (u >= 0){
      
#pragma omp for
      for (int i = 2; i<nx+2; i++){
        slope_l = get_slope(i-1);
        slope_r = get_slope(i);
        rho[i] = rho_old[i] - c *(rho_old[i] - rho_old[i-1]) - 0.5*c*(slope_r - slope_l)*(dx - u * dt);
      }
    }
    else {
      printf("Can't handle u<0. Aborting\n");
      exit(123);
    }
  }



  else if (method == 2){
    //-----------------------------------------------------
    // piecewise linear method with minmod slope limiter
    //-----------------------------------------------------

    double slope_l, slope_r;
   
    if (u >= 0) {
#pragma omp for
      for (int i = 2; i<nx+2; i++){
        slope_l = get_minmod_slope(i-1);
        slope_r = get_minmod_slope(i);
        rho[i] = rho_old[i] - c *(rho_old[i] - rho_old[i-1]) - 0.5*c*(slope_r - slope_l)*(dx - u * dt);
      }
    }
    else{
      printf("Can't handle u<0. Aborting\n");
      exit(123);
    }
  }




  else if (method == 3){
    //-----------------------------------------------------
    // piecewise linear method with VanLeer flux limiter
    //-----------------------------------------------------

    double fluxleft, fluxright; 

    if (u > 0) {
#pragma omp for
      for (int i = 2; i<nx+2; i++){
        fluxleft = u * rho_old[i-1] + 0.5*u*(1 - c) * VanLeer_limiter1(i) * (rho_old[i]-rho_old[i-1]);
        fluxright = u * rho_old[i] + 0.5*u*(1 - c) * VanLeer_limiter1(i+1) * (rho_old[i+1]-rho_old[i]);
        rho[i] = rho_old[i] + c * (fluxleft - fluxright);
      }
    }
    else{
#pragma omp for
      for (int i = 2; i<nx+2; i++){
        fluxleft = u * rho_old[i] - 0.5*u*(1 + c) * VanLeer_limiter2(i) * (rho_old[i]-rho_old[i-1]);
        fluxright = u * rho_old[i+1] - 0.5*u*(1 + c) * VanLeer_limiter2(i+1) * (rho_old[i+1]-rho_old[i]);
        rho[i] = rho_old[i] + c * (fluxleft - fluxright);
      }
    }
  }




  set_boundary();

}






//======================
void set_boundary()
//======================
{
  //------------------------------------
  // Set periodic boundary conditions
  //------------------------------------

#pragma omp single
  {
    rho[0] = rho[nx];
    rho[1] = rho[nx+1];
    rho[nx+2] = rho[2];
    rho[nx+3] = rho[3];
  }

  return;

}













//========================================================
void get_pwlin_flux(int i, double *flux1, double *flux2)
//========================================================
{
  //-----------------------------------
  // Get the flux for piecewise linear 
  // method.
  //-----------------------------------
  if (u < 0){
    *flux1 = (rho_old[i-1] + 3*rho_old[i] - 5*rho_old[i+1] + rho_old[i+2]);
    *flux2 = (rho_old[i-1] - rho_old[i] - rho_old[i+1] + rho_old[i+2]);
  }
  else{
    *flux1 = (rho_old[i+1] + 3*rho_old[i] - 5*rho_old[i-1] + rho_old[i-2]);
    *flux2 = (rho_old[i+1] - rho_old[i] - rho_old[i-1] + rho_old[i-2]);
  }
}



//===========================
double get_slope(int i)
//===========================
{
  //-------------------------------------------
  // Get centered slope for piecewise linear
  // method (Fromm’s method)
  //-------------------------------------------
  
  double slope = rho_old[i+1]-rho_old[i-1];
  slope = slope/(2 * dx);
  return (slope);
}






//================================
double get_minmod_slope(int i)
//================================
{
  
  //---------------------------------
  // Computes the slope for the 
  // minmod slope limiter
  //---------------------------------
  
  double a = rho_old[i] - rho_old[i-1];
  double b = rho_old[i+1] - rho_old[i];

  a = a/dx;
  b = b/dx;

  if (a * b > 0){
    if (fabs(a) < fabs(b)){
      return (a);
    } 
    else{
      return (b);
    }
  }
  else{
    return (0.0);
  }
}








//================================
double VanLeer_limiter1(int i)
//================================
{
  //---------------------------------------
  // computes the Van Leer flux limiter
  // for u >= 0
  //---------------------------------------



  if (rho_old[i]-rho_old[i-1] != 0){
    double r = (rho_old[i-1] - rho_old[i-2])/(rho_old[i]-rho_old[i-1]);
    double absr = fabs(r);
    double limiter = (r + absr)/(1 + absr);

    return (limiter);

  }
  else{
    return (0.0);
  }

}



//================================
double VanLeer_limiter2(int i)
//================================
{
  //---------------------------------------
  // computes the Van Leer flux limiter
  // for u < 0
  //---------------------------------------

  if (rho_old[i]-rho_old[i-1] != 0){
    double r = (rho_old[i+1] - rho_old[i])/(rho_old[i]-rho_old[i-1]);
    double absr = fabs(r);
    double limiter = (r + absr)/(1 + absr);
    return (limiter);
  }
  else{
    return (0.0);
  }

}
