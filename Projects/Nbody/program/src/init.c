/*
 * INITIALISATION AND SET UP
 */


#include <stdio.h>      /* input, output    */
#include <string.h>     /* string manipulation */
#include <stdlib.h>     /* some other string manipulation, alloc */
#include "commons.h"
#include "init.h"


void readparams(int argc, char *argv[]) // to pass on param file as cmd line arg
{
  /* 
   * reading parameters from parameters file
   */


  if (argc!=2)
  {
    // check if called correctly
    printf("ERROR: Usage ./my_program params.txt\n");
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
      sscanf(tempbuff, "%20s : %56[^;];", varname, varvalue);
      // reads formatted input from a string, writes it in
      // the variables given after the format string.
      // The format used is <string> separator <:> <string> ends with <;>


      if (strcmp(varname,"Nx")==0) {
        Nx = atoi(varvalue);
      // atoi/atof: convert string to integer/float
      // from stdlib.h
      } 
      else if (strcmp(varname,"Ny")==0) {
        Ny = atoi(varvalue);
      }
      else if (strcmp(varname,"Nz")==0) {
        Nz = atoi(varvalue);
      }
      else if (strcmp(varname,"dx")==0) {
        dx = atof(varvalue);
      }
      else if (strcmp(varname,"dy")==0) {
        dy = atof(varvalue);
      }
      else if (strcmp(varname,"dz")==0) {
        dz = atof(varvalue);
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





void initialise(int argc, char *argv[]){
  /*
   * Set everything up.
   */


  //read params
  readparams(argc, argv);

  // Nx = Ny = Nz = 1 initially
  grid = calloc(Nx*Ny*Nz, sizeof(int)); 
  


  // fill array with random shit to test

  void fillarray(int *arr);
  fillarray(grid);

}









void fillarray(int *arr){
  // fill array with random shit to test

  int i, j, k;
  int index;


  //fill up array and print it
  for (k = 0; k<Nz; k++){
    printf("\n\n\n z = %d \n\n", k);
    for (j = 0; j<Ny; j++){
      for (i = 0; i<Nx; i++){
        index=getiind(i,j,k);
        arr[index] = i + 10*j + 100*k;
        printf(" %5d ", arr[getiind(i,j,k)]);
      }
    printf("\n");
    }
  }

}


