/*
 * Define common global variables
 * Add one useless function
 */


#include <stdio.h>
#include <stdarg.h>


// VARIABLES

int    Nx = 1, Ny = 1, Nz = 1;
double dx = 1, dy = 1, dz = 1;

const double pi = 3.14159;
const double c = 2.998E8;


int verbose = 0; //verbose = false by default



//ARRAYS
int *grid;




// FUNCITONS

int getiind(int x, ...){

  va_list ap;
  int temp;
  int y = 0, z = 0;

  va_start(ap,x);

  temp = va_arg(ap, int);
  if (temp){
    y = temp;
  }

  temp = va_arg(ap, int);
  if (temp){
    z = temp;
  }

  va_end(ap);

#if NDIM == 2
  z = 0;
#endif
#if NDIM == 1
  y = 0;
#endif

    int result = x + Ny * y + Nz * z;
    return (result);
  }




