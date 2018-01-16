#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include "commons.h"
#include "io.h"

void set_boundaries();
void piecewise_constant_advection();
void piecewise_linear_advection();
void minmod_advection();
double get_minmod_slope_x(int i, int j, double **rho_arr);
double get_minmod_slope_y(int i, int j, double **rho_arr);
double VanLeer_limiter1_x(int i, int j);
double VanLeer_limiter1_y(int i, int j);
double VanLeer_limiter2_x(int i, int j);
double VanLeer_limiter2_y(int i, int j);







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
    dx = 1.0/((double) nx + 1);
    dy = 1.0/((double) ny + 1);
    rho2d = malloc((nx+4)*sizeof(double *));
    rho2d_old = malloc((nx+4)*sizeof(double *));
    rho2d_inter = malloc((nx+4)*sizeof(double *));
  }

#pragma omp barrier

#pragma omp for
  for ( int i = 0; i<nx+4; i++ ){
    rho2d[i] = calloc((ny+4), sizeof(double));
    rho2d_old[i] = malloc((ny+4)*sizeof(double));
    rho2d_inter[i] = malloc((ny+4)*sizeof(double));
  }


  // u = 1;
  // v = 1;



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
      for (int j = 0; j<(ny + 4); j++){

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


  dt = courant_factor/(fabs(u)/dx + fabs(v)/dy);
  // dt = 0.5*courant_factor *dx/ (sqrt(u*u + v*v));

  // check if you're jumping over output time step
  if ( t + dt >= t_out[t_out_step] ){
    dt = t_out[t_out_step] - t;
  }

  stepcounter += 1;

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



  double dtdx = dt/dx;
  double dtdy = dt/dy;
  double slope_lx, slope_rx, slope_ly, slope_ry;





  
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





//   else if (method == 2){
//     //-----------------------------------------------------
//     // piecewise linear method with minmod slope limiter
//     //-----------------------------------------------------
//
//
//     if (u >= 0) {
//       if (v >= 0) {
// #pragma omp for
//       for (int i = 2; i<nx+2; i++){
//         for (int j = 2; j<ny+2; j++){
//           slope_lx = get_minmod_slope_x(i-1, j);
//           slope_rx = get_minmod_slope_x(i, j);
//           slope_ly = get_minmod_slope_y(i, j-1);
//           slope_ry = get_minmod_slope_y(i, j);
//           rho2d[i][j] = rho2d_old[i][j] -
//             u*dtdx *(rho2d_old[i][j] - rho2d_old[i-1][j]) - 0.5*u*dtdx*(slope_rx - slope_lx)*(dx - u * dt) -
//             v*dtdy *(rho2d_old[i][j] - rho2d_old[i][j-1]) - 0.5*v*dtdy*(slope_ry - slope_ly)*(dy - v * dt);
//         }
//       }
//     }
//       else {
//         printf("Can't handle v<0. Aborting.\n");
//         exit(2);
//       }
//     }
//     else {
//       printf("Can't handle u<0. Aborting.\n");
//       exit(2);
//     }
//   }
//
//
//
//
//   else if (method == 3){
//     //-----------------------------------------------------
//     // piecewise linear method with VanLeer flux limiter
//     //-----------------------------------------------------
//
//     double fluxleft, fluxright;
//     double c = dt / dx * u;
//
//     if (u > 0) {
// #pragma omp for
//       for (int i = 2; i<nx+2; i++){
//         fluxleft = u * rho2d_old[i][i-1] + 0.5*u*(1 - c) * VanLeer_limiter1(i) * (rho2d_old[i][i]-rho2d_old[i][i-1]);
//         fluxright = u * rho2d_old[i][i] + 0.5*u*(1 - c) * VanLeer_limiter1(i+1) * (rho2d_old[i][i+1]-rho2d_old[i][i]);
//         rho2d[i][i] = rho2d_old[i][i] + c * (fluxleft - fluxright);
//       }
//     }
//     else{
// #pragma omp for
//       for (int i = 2; i<nx+2; i++){
//         fluxleft = u * rho2d_old[i][i] - 0.5*u*(1 + c) * VanLeer_limiter2(i) * (rho2d_old[i][i]-rho2d_old[i][i-1]);
//         fluxright = u * rho2d_old[i][i+1] - 0.5*u*(1 + c) * VanLeer_limiter2(i+1) * (rho2d_old[i][i+1]-rho2d_old[i][i]);
//         rho2d[i][i] = rho2d_old[i][i] + c * (fluxleft - fluxright);
//       }
//     }
//   }
//

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

  double dtdx = dt/dx;
  double dtdy = dt/dy;

  if (stepcounter % 2 == 0){
    //-----------------
    // sweep x first
    //-----------------

#pragma omp for
    for (int i = 1; i<nx+3; i++){
      for (int j = 0; j<ny+4; j++){
        //sweep x
        rho2d_inter[i][j] = rho2d_old[i][j] + 
          dtdx*u*(rho2d_old[i-1][j]-rho2d_old[i][j]);
      }
    }

#pragma omp for
    for (int i = 0; i<nx+4; i++){
      for (int j = 1; j<ny+3; j++){
        //sweep y
        rho2d[i][j] = rho2d_inter[i][j] + 
          dtdy*v*(rho2d_inter[i][j-1]-rho2d_inter[i][j]);
      }
    }
  }


  else{
    //-----------------
    // sweep y first
    //-----------------
    
#pragma omp for
    for (int i = 0; i<nx+4; i++){
      for (int j = 1; j<ny+3; j++){
        //sweep y
        rho2d_inter[i][j] = rho2d_old[i][j] + 
          dtdy*v*(rho2d_old[i][j-1]-rho2d_old[i][j]);
      }
    }

#pragma omp for
    for (int i = 1; i<nx+3; i++){
      for (int j = 0; j<ny+4; j++){
        //sweep x
        rho2d[i][j] = rho2d_inter[i][j] + 
          dtdx*u*(rho2d_inter[i-1][j]-rho2d_inter[i][j]);
      }
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

  double dtdx = dt/dx;
  double dtdy = dt/dy;
  double rho_xl, rho_xr, rho_yl, rho_yr;


  if (stepcounter % 2 == 0){
    //-----------------
    // sweep x first
    //-----------------

#pragma omp for
    for (int i = 2; i<nx+2; i++){
      for (int j = 0; j<ny+4; j++){
        // sweep x
        rho_xl = 0.5 * (rho2d_old[i][j] + rho2d_old[i-1][j]) - 
          0.5*u*dtdx*(rho2d_old[i][j] - rho2d_old[i-1][j]);
        rho_xr = 0.5 * (rho2d_old[i+1][j] + rho2d_old[i][j]) - 
          0.5*u*dtdx*(rho2d_old[i+1][j] - rho2d_old[i][j]);
        
        rho2d_inter[i][j] = rho2d_old[i][j] - u*dtdx*(rho_xr - rho_xl);
      }
    }

#pragma omp for
    for (int i = 0; i<nx+4; i++){
      for (int j = 2; j<ny+2; j++){
        // sweep y
        rho_yl = 0.5 * (rho2d_inter[i][j] + rho2d_inter[i][j-1]) - 
          0.5*v*dtdy*(rho2d_inter[i][j] - rho2d_inter[i][j-1]);
        rho_yr = 0.5 * (rho2d_inter[i][j+1] + rho2d_inter[i][j]) - 
          0.5*v*dtdy*(rho2d_inter[i][j+1] - rho2d_inter[i][j]);
        
        rho2d[i][j] = rho2d_inter[i][j] - v*dtdx*(rho_yr - rho_yl);
      }
    }
  }
  

  else {
    //-----------------
    // sweep y first
    //-----------------

#pragma omp for
    for (int i = 0; i<nx+4; i++){
      for (int j = 2; j<ny+2; j++){
        // sweep y
        rho_yl = 0.5 * (rho2d_old[i][j] + rho2d_old[i][j-1]) - 
          0.5*v*dtdy*(rho2d_old[i][j] - rho2d_old[i][j-1]);
        rho_yr = 0.5 * (rho2d_old[i][j+1] + rho2d_old[i][j]) - 
          0.5*v*dtdy*(rho2d_old[i][j+1] - rho2d_old[i][j]);
        
        rho2d_inter[i][j] = rho2d_old[i][j] - v*dtdx*(rho_yr - rho_yl);
      }
    }

#pragma omp for
    for (int i = 2; i<nx+2; i++){
      for (int j = 0; j<ny+4; j++){
        // sweep x
        rho_xl = 0.5 * (rho2d_inter[i][j] + rho2d_inter[i-1][j]) - 
          0.5*u*dtdx*(rho2d_inter[i][j] - rho2d_inter[i-1][j]);
        rho_xr = 0.5 * (rho2d_inter[i+1][j] + rho2d_inter[i][j]) - 
          0.5*u*dtdx*(rho2d_inter[i+1][j] - rho2d_inter[i][j]);
        
        rho2d[i][j] = rho2d_inter[i][j] - u*dtdx*(rho_xr - rho_xl);
      }
    }
  }
}








//=============================
void minmod_advection()
//=============================
{

  //----------------------------------------------------
  // piecewise linear method with minmod slope limiter
  //----------------------------------------------------

  double dtdx = dt/dx;
  double dtdy = dt/dy;
  double rho_xl, rho_xr, rho_yl, rho_yr;


  if (stepcounter % 2 == 0){
    //-----------------
    // sweep x first
    //-----------------

#pragma omp for
    for (int i = 2; i<nx+2; i++){
      for (int j = 0; j<ny+4; j++){
        // sweep x
        rho_xl = 0.5 * (rho2d_old[i][j] + rho2d_old[i-1][j]) - 
          0.5*u*dtdx*get_minmod_slope_x(i-1, j, rho2d_old);
        rho_xr = 0.5 * (rho2d_old[i+1][j] + rho2d_old[i][j]) - 
          0.5*u*dtdx*get_minmod_slope_x(i, j, rho2d_old);
        
        rho2d_inter[i][j] = rho2d_old[i][j] - u*dtdx*(rho_xr - rho_xl);
      }
    }

#pragma omp for
    for (int i = 0; i<nx+4; i++){
      for (int j = 2; j<ny+2; j++){
        // sweep y
        rho_yl = 0.5 * (rho2d_inter[i][j] + rho2d_inter[i][j-1]) - 
          0.5*v*dtdy*get_minmod_slope_y(i, j-1, rho2d_inter);
        rho_yr = 0.5 * (rho2d_inter[i][j+1] + rho2d_inter[i][j]) - 
          0.5*v*dtdy*get_minmod_slope_y(i, j, rho2d_inter);
        
        rho2d[i][j] = rho2d_inter[i][j] - v*dtdy*(rho_yr - rho_yl);
      }
    }
  }
  

  else {
    //-----------------
    // sweep y first
    //-----------------

#pragma omp for
    for (int i = 0; i<nx+4; i++){
      for (int j = 2; j<ny+2; j++){
        // sweep y
        rho_yl = 0.5 * (rho2d_old[i][j] + rho2d_old[i][j-1]) - 
          0.5*v*dtdy*get_minmod_slope_y(i, j-1, rho2d_old);
        rho_yr = 0.5 * (rho2d_old[i][j+1] + rho2d_old[i][j]) - 
          0.5*v*dtdy*get_minmod_slope_y(i, j, rho2d_old);
        
        rho2d_inter[i][j] = rho2d_old[i][j] - v*dtdy*(rho_yr - rho_yl);
      }
    }

#pragma omp for
    for (int i = 2; i<nx+2; i++){
      for (int j = 0; j<ny+4; j++){
        // sweep x
        rho_xl = 0.5 * (rho2d_inter[i][j] + rho2d_inter[i-1][j]) - 
          0.5*u*dtdx*get_minmod_slope_x(i-1, j, rho2d_inter);
        rho_xr = 0.5 * (rho2d_inter[i+1][j] + rho2d_inter[i][j]) - 
          0.5*u*dtdx*get_minmod_slope_x(i, j, rho2d_inter);
        
        rho2d[i][j] = rho2d_inter[i][j] - u*dtdx*(rho_xr - rho_xl);
      }
    }
  }
}












//========================================================
double get_minmod_slope_x(int i, int j, double **rho_arr)
//========================================================
{  
  //---------------------------------
  // Computes the slope for the 
  // minmod slope limiter
  //---------------------------------
  
  double a = rho_arr[i][j] - rho_arr[i-1][j];
  double b = rho_arr[i+1][j] - rho_arr[i][j];

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
double get_minmod_slope_y(int i, int j, double **rho_arr)
//========================================================
{  
  //---------------------------------
  // Computes the slope for the 
  // minmod slope limiter
  //---------------------------------
  
  double a = rho_arr[i][j] - rho_arr[i][j-1];
  double b = rho_arr[i][j+1] - rho_arr[i][j];

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








//======================================
double VanLeer_limiter1_x(int i, int j)
//======================================
{
  //---------------------------------------
  // computes the Van Leer flux limiter
  // for u >= 0
  //---------------------------------------



  if (rho2d_old[i][i]-rho2d_old[i][i-1] != 0){
    double r = (rho2d_old[i][i-1] - rho2d_old[i][i-2])/(rho2d_old[i][i]-rho2d_old[i][i-1]);
    double absr = fabs(r);
    double limiter = (r + absr)/(1 + absr);

    return (limiter);

  }
  else{
    return (0.0);
  }
}





//======================================
double VanLeer_limiter1_y(int i, int j)
//======================================
{
  //---------------------------------------
  // computes the Van Leer flux limiter
  // for u >= 0
  //---------------------------------------



  if (rho2d_old[i][i]-rho2d_old[i][i-1] != 0){
    double r = (rho2d_old[i][i-1] - rho2d_old[i][i-2])/(rho2d_old[i][i]-rho2d_old[i][i-1]);
    double absr = fabs(r);
    double limiter = (r + absr)/(1 + absr);

    return (limiter);

  }
  else{
    return (0.0);
  }
}




//=======================================
double VanLeer_limiter2_x(int i, int j)
//=======================================
{
  //---------------------------------------
  // computes the Van Leer flux limiter
  // for u < 0
  //---------------------------------------

  if (rho2d_old[i][i]-rho2d_old[i][i-1] != 0){
    double r = (rho2d_old[i][i+1] - rho2d_old[i][i])/(rho2d_old[i][i]-rho2d_old[i][i-1]);
    double absr = fabs(r);
    double limiter = (r + absr)/(1 + absr);
    return (limiter);
  }
  else{
    return (0.0);
  }
}


//=======================================
double VanLeer_limiter2_y(int i, int j)
//=======================================
{
  //---------------------------------------
  // computes the Van Leer flux limiter
  // for u < 0
  //---------------------------------------

  if (rho2d_old[i][i]-rho2d_old[i][i-1] != 0){
    double r = (rho2d_old[i][i+1] - rho2d_old[i][i])/(rho2d_old[i][i]-rho2d_old[i][i-1]);
    double absr = fabs(r);
    double limiter = (r + absr)/(1 + absr);
    return (limiter);
  }
  else{
    return (0.0);
  }
}
