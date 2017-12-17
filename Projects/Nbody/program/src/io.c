//================================
// INPUT AND OUTPUT 
//================================



#include <stdio.h>         
#include <string.h>     
#include <stdlib.h>    
#include "commons.h"


//======================================
void readparams(int argc, char *argv[]) 
//======================================
{
 
  //------------------------------------------
  // reading parameters from parameters file
  //------------------------------------------



  if (verbose)
  {
    printf("Reading parameters.\n");
  }




  if (argc!=3)
  {
    // check if called correctly
    printf("ERROR: Usage ./nbody nbody-params ../files/data.ascii\n");
    exit(1);
  }
  else
  {

    //open file
    FILE *params = fopen(argv[1], "r");

    // check if file exists
    if (params == NULL) { 
      printf("Error: file '%s' not found.\n", argv[1]);
      exit(1);
    }

    char varname[80] ;
    char varvalue[80] ;
    char tempbuff[80] ;

    


    while (fgets(tempbuff,80,params))
    // fgets(str_buff, n,filepointer) :
    // gets n characters from file in filepointer and stores them
    // in str_buff.
    // returns 0 if EoF is reached.
    
    {
      sscanf(tempbuff, "%20s = %56[^\n]\n", varname, varvalue);
      // reads formatted input from a string, writes it in
      // the variables given after the format string.
      // The format used is <string> separator <=> <string> ends with <;>
    

      if (strcmp(varname,"verbose")==0) {
        verbose = atoi(varvalue);
      // atoi/atof: convert string to integer/float
      // from stdlib.h
      } 
      else if (strcmp(varname, "//")==0) {
        // ignore comments
        continue;
      }
      else if (strcmp(varname, "/*")==0) {
        // ignore comments
        continue;
      }
      else{
        printf("Unrecongized parameter : \"%s\"\n", varname);
      }
    }

    fclose(params);


  }


}











//==========================
void readdata(char *argv[])
//==========================
{

  //-----------------------------------------------------------
  // Reading in the data, allocating data arrays and filling 
  // them up.
  //-----------------------------------------------------------


  if (verbose)
  {
    printf("Reading in data.\n");
  }



 
  
  //------------------
  // Read in header
  //------------------


  //open file
  FILE *data = fopen(argv[2], "r");


  // check if file exists
  if (data == NULL) { 
    printf("Error: file '%s' not found.\n", argv[2]);
    exit(1);
    }

  int info;
  info = fscanf(data, "%d %d %d\n", &npart, &ngaspart, &nstarpart);
  


  //--------------------------
  // Initialize arrays
  //--------------------------
  
 //TODO: remove cast? 
  m = calloc(npart, sizeof(double));
  x = calloc(npart, sizeof(double));
  y = calloc(npart, sizeof(double));
  z = calloc(npart, sizeof(double));
  vx = calloc(npart, sizeof(double));
  vy = calloc(npart, sizeof(double));
  vz = calloc(npart, sizeof(double));
  softening_data = calloc(npart, sizeof(double));
  potential_data = calloc(npart, sizeof(double));



  // --------------------------
  // Read in the actual data
  // --------------------------


  double *data_arr[9] = {m, x, y, z, vx, vy, vz, softening_data, potential_data};


  for (int value=0; value<9; value++){
    for (int i = 0; i < npart; i++){
      info = fscanf(data, "%lf\n", &data_arr[value][i]);
    }
  }
 


}
