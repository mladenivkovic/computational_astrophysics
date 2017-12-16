//=================================
// Contains main program.
//=================================



#include <stdio.h>         
#include <stdlib.h>
#include "commons.h"
#include "io.h"



//=====================================
// functions defined below
//=====================================
void initialise(int argc, char *argv[]);








//=====================================
int main(int argc, char * argv[])    
//=====================================
{

  if (verbose);
  {
    printf("Started program.\n");
  }
  
  initialise(argc, argv);


  if (verbose);
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


  //read params
  readparams(argc, argv);

  //read data
  readdata(argv);

  //transform read data to units
  void set_units();
  set_units();

}





//==================
void set_units()
//==================
{

  //--------------------------------------------
  // Sets units of read data to chosen units
  //--------------------------------------------

}

