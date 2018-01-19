#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include "commons.h"
#include "io.h"

void set_boundaries();

void piecewise_constant_advection();

void piecewise_linear_advection();
double slope_x(int i, int j);
double slope_y(int i, int j);

void minmod_advection();
double minmod_slope_x(int i, int j);
double minmod_slope_y(int i, int j);

void VanLeer_advection();
double VanLeer_limiter_x(int i, int j);
double VanLeer_limiter_y(int i, int j);







//=======================================
void initialise2d()
//=======================================
{
  //----------------------
  // Set everything up.
  //----------------------



  //-------------------------
  // initialize variables
  //-------------------------

  if (nx == 0 || ny == 0){
    printf("Got nx = %d, ny = %d, can't work with that. quitting.\n", nx, ny);
    exit(1);
  }

#pragma omp single
  {
    dx = 1.0/((double) nx);
    dy = 1.0/((double) ny);
  }
#pragma omp single
  { rho2d = malloc((nx+4)*sizeof(double *)); }
#pragma omp single
  { rho2d_old = malloc((nx+4)*sizeof(double *)); }

#pragma omp barrier

#pragma omp for
  for ( int i = 0; i<nx+4; i++ ){
    rho2d[i] = calloc((ny+4), sizeof(double));
    rho2d_old[i] = malloc((ny+4)*sizeof(double));
  }





  //------------------------------
  // Initialise density profile
  //------------------------------

  if (density_profile == 0){
#pragma omp master
    {
      //--------------------------------------------
      printf("Using step density profile.\n");
      //--------------------------------------------
    }
#pragma omp for
    for (int i = 0; i<nx+4; i++){
      for (int j = 0; j<ny+4; j++){

        if (i*dx <= 0.3){
          rho2d[i][j] += 1;
        }
        else if (i*dx <= 0.6){
          rho2d[i][j] += 2;
        }
        else{
          rho2d[i][j] += 1;
        }

        if (j*dy <= 0.3){
          rho2d[i][j] += 1;
        }
        else if (j*dy <= 0.6){
          rho2d[i][j] += 2;
        }
        else{
          rho2d[i][j] += 1;
        }
      }
    }
  }





  else if (density_profile == 1){
#pragma omp master
    {
      //--------------------------------------------
      printf("Using linear step density profile.\n");
      //--------------------------------------------
    }
#pragma omp for
    for (int i = 0; i<(nx+4); i++){
      for (int j = 0; j<(ny+4); j++){

        if (i*dx <= 0.3){
          rho2d[i][j] += 1;
        }
        else if (i*dx <= 0.6){
          rho2d[i][j] += 2+1.7*(i*dx-0.3);
        }
        else{
          rho2d[i][j] += 1;
        }

        if (j*dy <= 0.3){
          rho2d[i][j] += 1;
        }
        else if (j*dy <= 0.6){
          rho2d[i][j] += 2+1.7*(j*dy-0.3);
        }
        else{
          rho2d[i][j] += 1;
        }
      }
    }
  }



  else if (density_profile == 2){
#pragma omp master
    {
      //--------------------------------------------
      printf("Using gauss density profile.\n");
      //--------------------------------------------
    }
#pragma omp for
    for (int i = 0; i<(nx+2); i++){
      for (int j = 0; j<(nx+4); j++){
        rho2d[i][j] = 2 + exp(-pow((i*dx - 0.5), 2)/0.1) + exp(-pow((j*dy - 0.5),2)/0.1);
      }
    }
  }


  else {
    printf("Not recognized density_profile = %d\n", density_profile);
  }



  set_boundaries();

}














//=============================
void get_timestep2d()
//=============================
{
  //--------------------------------------------
  // Compute the next timestep according to CFL
  //--------------------------------------------


  dt = courant_factor/(u/dx + v/dy);
  // dt = 0.5*courant_factor *dx/ (sqrt(u*u + v*v));

  // check if you're jumping over output time step
  if ( t + dt >= t_out[t_out_step] ){
    if ( fabs(t_out_step -t -dt) <= 1e-4 ){
      // if difference is too small, cheat a bit
      t = t_out[t_out_step];
    }
    else{
      // reduce dt so that you get exact timestep for output
      dt = t_out[t_out_step] - t;
    }
  }


  // if (verbose) {
  //   printf("dt =%6g, t= %6g, dx= %6g, dy= %6g\n", dt, t+dt, dx, dy);
  // }
}










//================
void advect2d()
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
    for (int j = 0; j<ny+4; j++){
      rho2d_old[i][j] = rho2d[i][j];
    }
  }




  
  if (u >= 0) {
    if (v >= 0){
      
      if (method == 0){
        
        //--------------------------------
        // piecewise constant method
        //--------------------------------
        piecewise_constant_advection();
      }

      else if (method == 1){
        
        //--------------------------------
        // piecewise linear method
        //--------------------------------
        piecewise_linear_advection();
      }


      else if (method == 2){

        //-----------------------------------------------------
        // piecewise linear method with minmod slope limiter
        //-----------------------------------------------------
        
        minmod_advection();

      }

      else if (method == 3){
      
        //-----------------------------------------------------
        // piecewise linear method with VanLeer flux limiter
        //-----------------------------------------------------
        VanLeer_advection();
      }

    }
    else {
      printf("Can't handle v<0. Aborting.\n");
      exit(2);
    }
  }
  else {
    printf("Can't handle u<0. Aborting.\n");
    exit(2);
  }



  set_boundaries();

}





//===========================
void set_boundaries()
//===========================
{
  //------------------------------
  // Set periodic boundaries
  //------------------------------
  
#pragma omp for
  for (int i = 0; i<nx+4; i++){
    // for all x, copy boundary y
    rho2d[i][0] = rho2d[i][ny];
    rho2d[i][1] = rho2d[i][ny+1];
    rho2d[i][ny+2] = rho2d[i][2];
    rho2d[i][ny+3] = rho2d[i][3];
  }

#pragma omp for
  for (int j = 0; j<ny+4; j++){
    // for all y, copy boundary x
    rho2d[0][j] = rho2d[nx][j];
    rho2d[1][j] = rho2d[nx+1][j];
    rho2d[nx+2][j] = rho2d[2][j];
    rho2d[nx+3][j] = rho2d[3][j];
  }


}












//========================================
void piecewise_constant_advection()
//========================================
{
  //--------------------------------
  // piecewise constant method
  //--------------------------------

  double udtdx = u*dt/dx;
  double vdtdy = v*dt/dy;

#pragma omp for
  for (int i = 1; i<nx+2; i++){
    for (int j = 1; j<ny+2; j++){
      rho2d[i][j] = rho2d_old[i][j] + 
        udtdx*(rho2d_old[i-1][j]-rho2d_old[i][j]) +
        vdtdy*(rho2d_old[i][j-1]-rho2d_old[i][j]);
    }
  }
}









//========================================
void piecewise_linear_advection()
//========================================
{
  //--------------------------------
  // piecewise linear method
  //--------------------------------

  double udtdx = u*dt/dx;
  double vdtdy = v*dt/dy;

#pragma omp for
  for (int i = 2; i<nx+2; i++){
    for (int j = 2; j<ny+2; j++){
      rho2d[i][j] = rho2d_old[i][j] - 
        udtdx *(rho2d_old[i][j] - rho2d_old[i-1][j]) - 0.5*udtdx*(slope_x(i, j) - slope_x(i-1, j))*(dx - u * dt) -
        vdtdy *(rho2d_old[i][j] - rho2d_old[i][j-1]) - 0.5*vdtdy*(slope_y(i, j) - slope_y(i, j-1))*(dy - v * dt);
    }
  }
  
}



//====================================
double slope_x(int i, int j)
//====================================
{
  double slope = rho2d_old[i+1][j] - rho2d_old[i][j];
  slope = slope/dx;
  return (slope);
}




//====================================
double slope_y(int i, int j)
//====================================
{
  double slope = rho2d_old[i][j+1] - rho2d_old[i][j];
  slope = slope/dy;
  return (slope);
}









//=============================
void minmod_advection()
//=============================
{

  //----------------------------------------------------
  // piecewise linear method with minmod slope limiter
  //----------------------------------------------------

  double udtdx = u*dt/dx;
  double vdtdy = v*dt/dy;

#pragma omp for
  for (int i = 2; i<nx+2; i++){
    for (int j = 2; j<ny+2; j++){
      rho2d[i][j] = rho2d_old[i][j] - 
        udtdx *(rho2d_old[i][j] - rho2d_old[i-1][j]) - 0.5*udtdx*(minmod_slope_x(i, j) - minmod_slope_x(i-1, j))*(dx - u * dt) -
        vdtdy *(rho2d_old[i][j] - rho2d_old[i][j-1]) - 0.5*vdtdy*(minmod_slope_y(i, j) - minmod_slope_y(i, j-1))*(dy - v * dt);
    }
  }

}








//========================================================
double minmod_slope_x(int i, int j)
//========================================================
{  
  //---------------------------------
  // Computes the slope for the 
  // minmod slope limiter
  //---------------------------------
  
  double a = rho2d_old[i][j] - rho2d_old[i-1][j];
  double b = rho2d_old[i+1][j] - rho2d_old[i][j];

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




//========================================================
double minmod_slope_y(int i, int j)
//========================================================
{   
  //---------------------------------
  // Computes the slope for the 
  // minmod slope limiter
  //---------------------------------
  
  double a = rho2d_old[i][j] - rho2d_old[i][j-1];
  double b = rho2d_old[i][j+1] - rho2d_old[i][j];

  a = a/dy;
  b = b/dy;

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







//========================
void VanLeer_advection()
//========================
{

  double fluxleftx, fluxrightx, fluxlefty, fluxrighty;
  
  double udtdx = u * dt/dx;
  double vdtdy = v * dt/dy;

#pragma omp for
  for (int i = 2; i<nx+2; i++){
    for (int j = 2; j < ny+2; j++){
      fluxleftx = u * rho2d_old[i-1][j] + 
        0.5*u*(1 - udtdx) * VanLeer_limiter_x(i, j) * (rho2d_old[i][j]-rho2d_old[i-1][j]);
      fluxrightx = u * rho2d_old[i][j] + 
        0.5*u*(1 - udtdx) * VanLeer_limiter_x(i+1, j) * (rho2d_old[i+1][j]-rho2d_old[i][j]);
      fluxlefty = v * rho2d_old[i][j-1] + 
        0.5*v*(1 - vdtdy) * VanLeer_limiter_y(i, j) * (rho2d_old[i][j]-rho2d_old[i][j-1]);
      fluxrighty = v * rho2d_old[i][j] + 
        0.5*v*(1 - vdtdy) * VanLeer_limiter_y(i, j+1) * (rho2d_old[i][j+1]-rho2d_old[i][j]);
      rho2d[i][j] = rho2d_old[i][j] + udtdx * (fluxleftx - fluxrightx) + vdtdy * (fluxlefty - fluxrighty);
    }
  }


}






//======================================
double VanLeer_limiter_x(int i, int j)
//======================================
{
  //---------------------------------------
  // computes the Van Leer flux limiter
  // for u >= 0
  //---------------------------------------



  if (rho2d_old[i][j]-rho2d_old[i-1][j] != 0){
    double r = (rho2d_old[i-1][j] - rho2d_old[i-2][j])/(rho2d_old[i][j]-rho2d_old[i-1][j]);
    double absr = fabs(r);
    double limiter = (r + absr)/(1 + absr);

    return (limiter);
  }
  else{
    return (0.0);
  }
}





//======================================
double VanLeer_limiter_y(int i, int j)
//======================================
{
  //---------------------------------------
  // computes the Van Leer flux limiter
  // for u >= 0
  //---------------------------------------



  if (rho2d_old[i][j]-rho2d_old[i][j-1] != 0){
    double r = (rho2d_old[i][j-1] - rho2d_old[i][j-2])/(rho2d_old[i][j]-rho2d_old[i][j-1]);
    double absr = fabs(r);
    double limiter = (r + absr)/(1 + absr);

    return (limiter);

  }
  else{
    return (0.0);
  }
}




