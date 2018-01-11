//================================
// INPUT AND OUTPUT 
//================================



#include <stdio.h>         
#include <string.h>     
#include <stdlib.h>    
#include <math.h>
#include "commons.h"


//======================================
void readparams(int argc, char *argv[]) 
//======================================
{
 
  //------------------------------------------
  // reading parameters from parameters file
  //------------------------------------------



  printf("Reading parameters.\n");




  if (argc!=2)
  {
    // check if called correctly
    printf("ERROR: Usage ./advection advection-params\n");
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
      else if (strcmp(varname, "nx") == 0){
        nx = atoi(varvalue);
      }
      else if (strcmp(varname, "density_profile")==0){
        density_profile = atoi(varvalue);
      }
      else if (strcmp(varname, "t_end")==0){
        t_end = atof(varvalue);
      }
      else if (strcmp(varname, "t_out_step")==0){
        t_out_step = atof(varvalue);
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











//==================================
void write_output(int output_case)
//==================================
{

  //----------------------------------------
  // Writes direct forces results to file.
  // Cases:
  //----------------------------------------



  //---------------------
  // get filename
  //---------------------

  char filename[80] = "output_"; 

  if (density_profile == 0) {


    strcat(filename, "step_function-");
    // char softening_str[10];
    // sprintf(softening_str, "%.4g", f_softening);
    // strcat(filename, softening_str);
  }
  // else if (output_case == 2){
  //
  //
  //   // get filename
  //   strcat(filename, "multipole_");
  //   char multipole_str[10];
  //   char thetamax_str[10];
  //   sprintf(multipole_str, "%1d", multipole_order);
  //   sprintf(thetamax_str, "%3g", theta_max);
  //   strcat(filename, multipole_str);
  //   strcat(filename, "-");
  //   strcat(filename, thetamax_str);
  //   strcat(filename, ".dat");
  // }
  else{
    printf("Something went wrong with output. Got case=%d, which I dont recognise.\n", output_case);
    return;
  }
  
  char time_str[8] = "";
  sprintf(time_str, "%06.2f", t);
  strcat(filename, time_str);
  strcat(filename, ".dat");



  //-----------------------
  // write output to file
  //-----------------------
  
  FILE *outfilep = fopen(filename, "w");
  
  fprintf(outfilep, "%15.7g\n", t);
  for (int i = 0; i<nx; i++){
    fprintf(outfilep, "%15.7g\n", u[i+1] );
  }

  fclose(outfilep);

}












