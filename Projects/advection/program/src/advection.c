#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include "commons.h"


double get_pwconst_flux(int i);
void get_pwlin_flux(int i, double *flux1, double *flux2);
double get_slope(int i);
double get_minmod_slope(int i);
double VanLeer_limiter1(int i);
double VanLeer_limiter2(int i);


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

#pragma omp parallel
  {
    // store old values
#pragma omp for
    for (int i = 0; i<nx+4; i++){
      u_old[i] = u[i];
    }


    //---------------------- 
    // Loop over cells
    //---------------------- 
    
    if (method == 0){
      //--------------------------------
      // piecewise constant method
      //--------------------------------
      double a = dt/dx * v;
#pragma omp for
      for (int i = 2; i<nx+2; i++){
        u[i] = u_old[i] + a*get_pwconst_flux(i);
      }
    }


    else if (method == 1){
      //--------------------------------
      // piecewise linear method
      //--------------------------------
//       double a1 = v * dt / (4 * dx);
//       double a2 = a1 * a1*4;
//       double f1 = 0, f2 = 0;
// #pragma omp for
//       for (int i = 2; i<nx+2; i++){
//         get_pwlin_flux(i, &f1, &f2);
//         u[i] = u_old[i] - a1*f1 + a2*f2;
      double c = v * dt / dx;
      double slope_l, slope_r;
     
      if (v > 0) {
#pragma omp for
        for (int i = 2; i<nx+2; i++){
          slope_l = get_slope(i-1);
          slope_r = get_slope(i);
          u[i] = u_old[i] - c *(u_old[i] - u_old[i-1]) - 0.5*c*(slope_r - slope_l)*(dx - v * dt);
        }
      }
      else{
#pragma omp for
        for (int i = 2; i<nx+2; i++){
          slope_l = get_slope(i);
          slope_r = get_slope(i+1);
          u[i] = u_old[i] - c *(u_old[i+1] - u_old[i]) - 0.5*c*(slope_r - slope_l)*(dx - v * dt);
        }
      }
    }



    else if (method == 2){
      //-----------------------------------------------------
      // piecewise linear method with minmod slope limiter
      //-----------------------------------------------------

      double c = v * dt / dx;
      double slope_l, slope_r;
     
      if (v > 0) {
#pragma omp for
        for (int i = 2; i<nx+2; i++){
          slope_l = get_minmod_slope(i-1);
          slope_r = get_minmod_slope(i);
          u[i] = u_old[i] - c *(u_old[i] - u_old[i-1]) - 0.5*c*(slope_r - slope_l)*(dx - v * dt);
        }
      }
      else{
#pragma omp for
        for (int i = 2; i<nx+2; i++){
          slope_l = get_minmod_slope(i);
          slope_r = get_minmod_slope(i+1);
          u[i] = u_old[i] - c *(u_old[i+1] - u_old[i]) - 0.5*c*(slope_r - slope_l)*(dx - v * dt);
        }
      }
    }




    else if (method == 3){
      //-----------------------------------------------------
      // piecewise linear method with VanLeer flux limiter
      //-----------------------------------------------------

      double fluxleft, fluxright; 
      double c = dt / dx * v;

      if (v > 0) {
#pragma omp for
        for (int i = 2; i<nx+2; i++){
          fluxleft = v * u_old[i-1] + 0.5*v*(1 - c) * VanLeer_limiter1(i) * (u_old[i]-u_old[i-1]);
          fluxright = v * u_old[i] + 0.5*v*(1 - c) * VanLeer_limiter1(i+1) * (u_old[i+1]-u_old[i]);
          u[i] = u_old[i] + c * (fluxleft - fluxright);
        }
      }
      else{
#pragma omp for
        for (int i = 2; i<nx+2; i++){
          fluxleft = v * u_old[i] - 0.5*v*(1 + c) * VanLeer_limiter2(i) * (u_old[i]-u_old[i-1]);
          fluxright = v * u_old[i+1] - 0.5*v*(1 + c) * VanLeer_limiter2(i+1) * (u_old[i+1]-u_old[i]);
          u[i] = u_old[i] + c * (fluxleft - fluxright);
        }
      }
    }
  } // end parallel region




  //implement periodic boundary condition
  u[0] = u[nx];
  u[1] = u[nx+1];
  u[nx+2] = u[2];
  u[nx+3] = u[3];



}




//=========================================
double get_pwconst_flux(int i)
//=========================================
{
  //---------------------------------------
  // computes the flux for the piecewise
  // constant method
  //---------------------------------------

  double f = 0;

  if (v < 0){
    f = (u_old[i]-u_old[i+1]);
  }
  else{
    f = (u_old[i-1]-u_old[i]);
  }

  return (f);
}







//========================================================
void get_pwlin_flux(int i, double *flux1, double *flux2)
//========================================================
{
  if (v < 0){
    *flux1 = (u_old[i-1] + 3*u_old[i] - 5*u_old[i+1] + u_old[i+2]);
    *flux2 = (u_old[i-1] - u_old[i] - u_old[i+1] + u_old[i+2]);
  }
  else{
    *flux1 = (u_old[i+1] + 3*u_old[i] - 5*u_old[i-1] + u_old[i-2]);
    *flux2 = (u_old[i+1] - u_old[i] - u_old[i-1] + u_old[i-2]);
  }
}



//===========================
double get_slope(int i){
//===========================
  double slope = u_old[i+1]-u_old[i-1];
  slope = slope/(2 * dx);
  return (slope);
}






//================================
double get_minmod_slope(int i){
//================================
  
  //---------------------------------
  // Computes the slope for the 
  // minmod slope limiter
  //---------------------------------
  
  double a = u_old[i] - u_old[i-1];
  double b = u_old[i+1] - u_old[i];

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
  // for v >= 0
  //---------------------------------------



  if (u_old[i]-u_old[i-1] != 0){
    double r = (u_old[i-1] - u_old[i-2])/(u_old[i]-u_old[i-1]);
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
  // for v < 0
  //---------------------------------------

  if (u_old[i]-u_old[i-1] != 0){
    double r = (u_old[i+1] - u_old[i])/(u_old[i]-u_old[i-1]);
    double absr = fabs(r);
    double limiter = (r + absr)/(1 + absr);
    return (limiter);
  }
  else{
    return (0.0);
  }

}
