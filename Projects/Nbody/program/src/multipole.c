//=====================
// Multipole Method
//=====================

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "commons.h"
#include "multipole.h"






//===============
// Global Vars
//===============

node * cells = 0;
unsigned int lastcell = 8;
int * partcell = 0;

int * thislevel_refine = 0;
int * nextlevel_refine = 0;

int nthislevel = 0;
int nnextlevel = 0;

int to_refine = 1;





//================
//Functions
//================

void build_root();
void refine(node parentnode);









//=====================
void build_tree(){
//=====================


  //-------------
  // PREPARATION
  //-------------
 

  if (verbose) { printf("Building tree.\n"); }

  // check whether ncellmax is defined
  if (ncellmax == 0) { ncellmax = 20*8*npart; }

  // allocate cells array
  cells = calloc(ncellmax, sizeof(node));
  
  // set boxlen = 2*rmax
  boxlen = 2 ;




  //----------------
  // BUILD ROOT
  //----------------

  build_root();




  //------------------
  // Start refinement
  //------------------


  int levelcounter = 0;
  while (to_refine)
  {
    levelcounter += 1;
    to_refine = 0;

    for (int n = 0; n < nthislevel; n++){

// printf("Going for cell %d with %d particles at level %d\n", n, cells[n].np, levelcounter-1);
     
      refine(cells[thislevel_refine[n]]);
    }

    
printf("Finished level %d with nthis=%d and nnext=%d\n", levelcounter-1, nthislevel, nnextlevel);
printf("nthis:\n");
for (int i = 0; i<nthislevel; i++){
  printf("%d ", thislevel_refine[i]);
}
printf("\nnnext:\n");
for (int i = 0; i<nnextlevel; i++){
  printf("%d ", nextlevel_refine[i]);
}
printf("\ntotal cells:%d\n\n",lastcell);


    // check whether another step is necessary
    if (nnextlevel > 0){
      
      // repeat on
      to_refine = 1;

      // copy next level cells into this level cell array
      thislevel_refine = calloc(pow(8,levelcounter), sizeof(int));
      for (int n = 0; n<nnextlevel; n++){
        thislevel_refine[n]=nextlevel_refine[n];
      }
      nthislevel = nnextlevel;
      
      // initialize new next level cell array
      nextlevel_refine = calloc(pow(8,levelcounter + 1), sizeof(int));

      // reset number of next level particles
      nnextlevel = 0;
    }


  }




  printf("Finished refinement. Max level: %d, ncells: %d\n", max_refinement_level, lastcell);
  

}










//======================
void build_root()
//======================

{


  if (verbose) { printf("Building root.\n"); }





  //-----------------
  // Preparation
  //-----------------
 
  // initialize cell particles are in
  partcell = calloc(npart, sizeof(int));


  double c = boxlen/2;
  
  // centres for level one:
  double cs[8][3] = { {-c, -c, -c}, {c, -c, -c}, {-c, c, -c}, {c, c, -c}, {-c, -c, c}, {c, -c, c}, {-c, c, c}, {c, c, c} };


  // refinement steps array initialisation
  thislevel_refine = calloc(8, sizeof(int));
  nextlevel_refine = calloc(64, sizeof(int));



  
  //------------------
  // Initialise root
  //------------------
  
  for (int i =0; i<8; i++) {

    int nochild[] = {-1, -1, -1, -1, -1, -1, -1, -1};

    double thiscenter[] = {cs[i][0], cs[i][1], cs[i][2]};
    int * parray = calloc(npart, sizeof(int));

    node root = {-1, i, 1, nochild, thiscenter , 0, parray } ;


    cells[i] = root;

    // printf("TESTING %d\n", i);
    // printf("parent %d\nchildren ", cells[i].parent);
    // for (int j = 0; j<8; j++) { printf("%d ", cells[i].child[j]); }
    // printf("\nlevel %d\ncenter ", cells[i].level);
    // for (int j = 0; j<3; j++) { printf("%g ", cells[i].center[j]); }
    // printf("\n\n");
  
  }





  //-----------------------------
  // DIVIDE PARTICLES INTO ROOT
  //-----------------------------
  
  // assume center is at (0,0,0)

  int i = 0;
  int j = 0;
  int k = 0;
  int ind;


  for (int p = 0; p < npart; p++){
    // find in which root a particle belongs
    i = 0;
    j = 0;
    k = 0;
    if (x[p] > 0){ i = 1; }
    if (y[p] > 0){ j = 1; }
    if (z[p] > 0){ k = 1; }

    ind = i + 2*j + 4*k;

    partcell[p] = ind;

    // add particle to root linked list
    cells[ind].particles[ cells[ind].np ] = p;

    // raise particle counter for root cell
    cells[ind].np += 1;

  }


  // find out which root cells need to be refined

  for (int i = 0; i<8; i++){
    if (cells[i].np > ncellpartmax){
// printf("roottest %d %d %d\n", i, cells[i].np, cells[i].cellindex);
      thislevel_refine[nthislevel] = i;
      nthislevel += 1;
    }
  }


// printf("Finished building root with nthis=%d and nnext=%d\n", nthislevel,  nnextlevel);
// printf("\nnthis:\n");
// for (int i = 0; i<nthislevel; i++){
//   printf("%d ", thislevel_refine[i]);
// }


}













//===========================
void refine(node parent)
//===========================
{
  

  if (verbose) {printf("Refining cell %d of level %d\n", parent.cellindex, parent.level);}



  //---------------------
  // Preparation
  //---------------------

  node * newchildren = calloc(8, sizeof(node));

  double xp = parent.center[0];
  double yp = parent.center[1];
  double zp = parent.center[2];
  int npar = parent.np;
  int parind = parent.cellindex;


  // get child level, update max ref level so far
  int childlevel = parent.level + 1;

  if (max_refinement_level < childlevel) { 
    max_refinement_level = childlevel; 
  }

  
  // centres for child level:
  double c = boxlen/pow(2, childlevel);

printf("boxlen: %g c: %g\n", boxlen, c);
  double cs[8][3] = { {xp-c, yp-c, zp-c}, {xp+c, yp-c, zp-c}, {xp-c, yp+c, zp-c}, {xp+c, yp+c, zp-c}, {xp-c, yp-c, zp+c}, {xp+c, yp+-c, zp+c}, {xp-c, yp+c, zp+c}, {xp+c, yp+c, zp+c} };

  



  //----------------------
  // Initialise children
  //----------------------
  
  for (int i =0; i<8; i++) {

    int nochild[] = {-1, -1, -1, -1, -1, -1, -1, -1};
    double thiscenter[] = {cs[i][0], cs[i][1], cs[i][2]};
    int * parray = calloc(npar, sizeof(int));

    node thischild = {parind, -1, childlevel, nochild, thiscenter, 0, parray } ;


    newchildren[i] = thischild;

  }




  //-----------------------------------
  // DIVIDE PARTICLES AMONGST CHILDREN
  //-----------------------------------
  

  int i = 0;
  int j = 0;
  int k = 0;
  int ind;


  node *thischild;
  for (int p = 0; p < npar; p++){

    i = 0;
    j = 0;
    k = 0;
    if (x[p] > xp){ i = 1; }
    if (y[p] > yp){ j = 1; }
    if (z[p] > zp){ k = 1; }
printf("x %7g xp %7g y %7g yp %7g z %7g zp %7g c %7g\n", x[p], xp, y[p], yp, z[p], zp, c);
    ind = i + 2*j + 4*k;

    thischild = &newchildren[ind];
    (*thischild).particles[(*thischild).np] = p;
    (*thischild).np += 1;

  }




  
 
  
  //------------------------------
  // ADD CHILDREN TO CELL ARRAY
  //------------------------------
  

  for (int child = 0; child<8; child++){
   
    // if it has any particle, add to list
    if (newchildren[child].np > 0){

printf("Found new child, local index %d, global index %d, ptcl %d\n", child, lastcell, newchildren[child].np);
      // give child a unique index
      int thischild_ind = lastcell;
      lastcell += 1;
      newchildren[child].cellindex = thischild_ind;

      // copy child to cell array
      cells[thischild_ind] = newchildren[child];

      // tell parent what index its child has
      parent.child[child] = thischild_ind;

      // tell particles what cell they belong to
      for (int p = 0; p < cells[thischild_ind].np; p++){
        partcell[p] = thischild_ind;
      }


      if ( cells[thischild_ind].np > ncellpartmax ){
        // this child needs to be refined
        nextlevel_refine[nnextlevel] = thischild_ind;
        nnextlevel += 1;
      }


    }
  }


  





}
