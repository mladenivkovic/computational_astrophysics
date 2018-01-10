//=================================
// Contains main program.
//=================================



#include <stdio.h>         
#include <stdlib.h>
#include <math.h>
#include <omp.h>
#include "commons.h"
#include "io.h"



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



  write_output(0);

  if (verbose) { printf("I'm finished!\n"); }
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
  dx = 1.0/((double) nx);



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
}

