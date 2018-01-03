//================================
// INPUT AND OUTPUT 
//================================



#include <stdio.h>         
#include <string.h>     
#include <stdlib.h>    
#include <math.h>
#include "commons.h"
#include "multipole.h"


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
      else if (strcmp(varname, "f_softening") == 0){
        f_softening = atof(varvalue);
      }
      else if (strcmp(varname, "direct_force")==0){
        direct_force = atoi(varvalue);
      }
      else if (strcmp(varname, "ncellmax")==0){
        ncellmax = atoi(varvalue);
      }
      else if (strcmp(varname, "multipole")==0){
        multipole = atoi(varvalue);
      }
      else if (strcmp(varname, "ncellpartmax")==0){
        ncellpartmax = atoi(varvalue);
      }
      else if (strcmp(varname, "scale_cube")==0){
        scale_cube = atoi(varvalue);
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
  
  m = calloc(npart, sizeof(double));
  x = calloc(npart, sizeof(double));
  y = calloc(npart, sizeof(double));
  z = calloc(npart, sizeof(double));
  vx = calloc(npart, sizeof(double));
  vy = calloc(npart, sizeof(double));
  vz = calloc(npart, sizeof(double));
  softening_data = calloc(npart, sizeof(double));
  potential_data = calloc(npart, sizeof(double));

  fx = calloc(npart, sizeof(double));
  fy = calloc(npart, sizeof(double));
  fz = calloc(npart, sizeof(double));
  r = calloc(npart, sizeof(double));
  phi = calloc(npart, sizeof(double));


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





//============================
void output_direct_force()
//============================
{

  //-------------------------------------
  // Writes direct forces results to file.
  //-------------------------------------
 
  if(verbose){ printf("Writing direct forces to file.\n"); }

  // get filename
  char filename[80] = "output_direct_force_"; 
  char softening_str[10];
  sprintf(softening_str, "%.4g", f_softening);
  strcat(filename, softening_str);
  strcat(filename, ".dat");




  // write to file
  FILE *outfilep = fopen(filename, "w");

  fprintf(outfilep, "%15s   %15s   %15s   %15s   %15s   %15s   %15s   %15s   %15s\n", "x ", "y ", "z ", "r ", "m ", "fx", "fy", "fz", "ftot");
  for (int i = 0; i<npart; i++){
    fprintf(outfilep, "%15g   %15g   %15g   %15g   %15g   %15g   %15g   %15g   %15g\n", x[i], y[i], z[i], r[i], m[i], fx[i], fy[i], fz[i], sqrt(pow(fx[i],2) + pow(fy[i], 2) + pow(fz[i], 2) ) );
  }

  fclose(outfilep);

}









//============================
void output_multipole()
//============================
{

  //-------------------------------------
  // Writes direct forces results to file.
  //-------------------------------------
 
  if(verbose){ printf("Writing multipole forces to file.\n"); }


  // get filename
  char filename[80] = "output_multipole_"; 
  char multipole_str[10];
  sprintf(multipole_str, "%.4g", 0.0);
  strcat(filename, multipole_str);
  strcat(filename, ".dat");




  // write to file
  FILE *outfilep = fopen(filename, "w");

  fprintf(outfilep, "%15s   %15s   %15s   %15s   %15s   %15s   %15s   %15s   %15s\n", "x ", "y ", "z ", "r ", "m ", "fx", "fy", "fz", "ftot");
  for (int i = 0; i<npart; i++){
    fprintf(outfilep, "%15g   %15g   %15g   %15g   %15g   %15g   %15g   %15g   %15g\n", x[i], y[i], z[i], r[i], m[i], fx[i], fy[i], fz[i], sqrt(pow(fx[i],2) + pow(fy[i], 2) + pow(fz[i], 2) ) );
  }

  fclose(outfilep);

}







//==============================
void write_info()
//==============================
{

  //-------------------------------------
  // Writes run info to file.
  //-------------------------------------
 
  if(verbose){ printf("Writing info to file.\n"); }
  // get filename
  char filename[80] = "info_"; 

  if (multipole) {
    char multipole_str[10];
    sprintf(multipole_str, "%.4g", 0.0);
    strcat(filename, "multipole_");
    strcat(filename, multipole_str);
  }
  else {
    char softening_str[10];
    sprintf(softening_str, "%.4g", f_softening);
    strcat(filename, "direct_"); 
    strcat(filename, softening_str);
  }

  strcat(filename, ".txt");


  // write to file
  FILE *outfilep = fopen(filename, "w");

  fprintf(outfilep, "%15s \t %15s \t %15s \t %15s\n", "scale_m", "scale_l", "scale_t", "softening");
  fprintf(outfilep, "%15g \t %15g \t %15g \t %15g\n", scale_m, scale_l, scale_t, softening);

  fclose(outfilep);


}













//============================
void write_cellparticles()
//============================
{

  //-----------------------------------------------------
  // Writes particle-cell division to file.
  //-----------------------------------------------------
 
  if(verbose){ printf("Writing cell particle info to file.\n"); }



  // write to file
  FILE *outfilep = fopen("cellparticles.dat", "w");


  fprintf(outfilep, "%15s   %15s   %15s   %15s  \n", "x ", "y ", "z ", "cell ");
  for (int p = 0; p < npart; p++){
    fprintf(outfilep, "%15g   %15g   %15g   %15d \n", x[p], y[p], z[p], partcell[p] );
  }

  fclose(outfilep);

}


