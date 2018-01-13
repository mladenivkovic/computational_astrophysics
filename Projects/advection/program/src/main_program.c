//=================================
// Contains main program.
//=================================



#include <stdio.h>         
#include <stdlib.h>
#include <math.h>
#include <omp.h>
#include "commons.h"
#include "io.h"
#include "advection.h"



//=====================================
// functions defined below
//=====================================
void initialise(int argc, char *argv[]);







//=====================================
int main(int argc, char *argv[])    
//=====================================
{

  printf("Started advection program.\n");
 

  //-----------------------
  // Initialise program
  //-----------------------
  initialise(argc, argv);
  if (verbose) {printf("Courant factor: %g\n", courant_factor);}


  //-----------------------
  // Main advection loop
  //-----------------------
  while (t < t_end){
    
    // compute next timestep
    get_timestep();

    // integrate density advection
    advect();

    // move timestep along, write output if necessary
    t += dt;
    if (t == t_out[t_out_step]){
      write_output(0);
      t_out_step += 1;
    }

    // t = t_end; // stop at first step
  }




  if (verbose) { printf("\nI'm finished!\n"); }
  return(0);

}










//=======================================
void initialise(int argc, char *argv[])
//=======================================
{
  //----------------------
  // Set everything up.
  //----------------------



  //--------------------
  //read parameters
  //--------------------
  readparams(argc, argv);


  //-------------------------
  // initialize variables
  //-------------------------

  if (nx == 0){
    printf("Got nx = 0, can't work with that. quitting.\n");
    exit(1);
  }

  u = malloc((nx+2)*sizeof(double));
  u_old = malloc((nx+2)*sizeof(double));
  dx = 1.0/((double) nx);
  v = 1;



  //------------------------------
  // Initialise density profile
  //------------------------------

  if (density_profile == 0){
    //--------------------------------------------
    printf("Using step density profile.\n");
    //--------------------------------------------
    for (int i = 0; i<(nx+2); i++){
      if (i*dx <= 0.3){
        u[i] = 1;
      }
      else if (i*dx <= 0.6){
        u[i] = 2;
      }
      else{
        u[i] = 1;
      }
    }
  }
  else if (density_profile == 1){
    //--------------------------------------------
    printf("Using linear step density profile.\n");
    //--------------------------------------------
    for (int i = 0; i<(nx+2); i++){
      if (i*dx <= 0.3){
        u[i] = 1;
      }
      else if (i*dx <= 0.6){
        u[i] = 2+1.7*(i*dx-0.3);
      }
      else{
        u[i] = 1;
      }
    }
  }
  else if (density_profile == 2){
    //--------------------------------------------
    printf("Using gauss density profile.\n");
    //--------------------------------------------
    for (int i = 0; i<(nx+2); i++){
      u[i] = 1 + exp(-pow((i*dx - 0.5), 2)/0.1);
    }
  }
  else {
    printf("Not recognized density_profile = %d\n", density_profile);
  }

}

