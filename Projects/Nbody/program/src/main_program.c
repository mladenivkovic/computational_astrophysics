//=================================
// Contains main program.
//=================================



#include <stdio.h>         
#include <stdlib.h>
#include <math.h>
#include "commons.h"
#include "io.h"
#include "calculations.h"



//=====================================
// functions defined below
//=====================================
void initialise(int argc, char *argv[]);








//=====================================
int main(int argc, char *argv[])    
//=====================================
{

  if (verbose)
  {
    printf("Started program.\n");
  }
  
  initialise(argc, argv);


  if (verbose)
  {
    printf("I'm finished!\n");
  }
  return(0);
}










//=======================================
void initialise(int argc, char *argv[])
//=======================================
{
  
  //----------------------
  // Set everything up.
  //----------------------


  //read parameters
  readparams(argc, argv);

  //read data
  readdata(argv);

  //transform read data to units
  set_units();

}




