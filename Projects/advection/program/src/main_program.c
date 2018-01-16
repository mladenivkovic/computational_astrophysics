//=================================
// Contains main program.
//=================================



#include <stdio.h>         
#include <stdlib.h>
#include <math.h>
#include <omp.h>
#include "commons.h"
#include "io.h"
#include "advection1d.h"
#include "advection2d.h"






//=====================================
// functions defined below
//=====================================
void run1d();
void run2d();






//=====================================
int main(int argc, char *argv[])    
//=====================================
{

  printf("Started advection program.\n");
 
  //--------------------
  //read parameters
  //--------------------
  readparams(argc, argv);


  if (use_2d){
    run2d();
  }
  else{
    run1d();
  }



  if (verbose) { printf("\nI'm finished!\n"); }
  return(0);

}





//==========================
void run1d()
//==========================
{
 
  
#pragma omp parallel shared (rho, rho_old, rho_inter, dt, dx, nx, courant_factor)
  {

    //-----------------------
    // Initialise program
    //-----------------------
    
    initialise1d();
#pragma omp master
    {
      if (verbose) {printf("Courant factor: %g\n", courant_factor);}
    }


    //-----------------------
    // Main advection loop
    //-----------------------
    while (t < t_end){
      
      // compute next timestep
#pragma omp single
      {
        get_timestep1d();
      }
#pragma omp barrier

      // integrate density advection
      advect1d();

      // move timestep along, write output if necessary
#pragma omp single
      {
        t += dt;
        if (t == t_out[t_out_step]){
          write_output();
          t_out_step += 1;
        }
      }
    }
  } // end parallel
}






//==========================
void run2d()
//==========================
{
  

#pragma omp parallel
  {

    //-----------------------
    // Initialise program
    //-----------------------
    initialise2d();
#pragma omp single
    {
      if (verbose) {printf("Courant factor: %g\n", courant_factor);}
    }

    


    //-----------------------
    // Main advection loop
    //-----------------------
    while (t < t_end){
#pragma omp single
      { 
        // compute next timestep
        get_timestep2d(); 
      }
#pragma omp barrier

      // integrate density advection
      advect2d();

#pragma omp single
      {
        // move timestep along, write output if necessary
        t += dt;
        if (t == t_out[t_out_step]){
          write_output();
          t_out_step += 1;
        }
      }
    }
  } // end parallel region
}




