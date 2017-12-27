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

node ** cells = 0;
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
void refine(node * parentnode);









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
  cells = malloc(ncellmax * sizeof(node*));
  
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

  if (verbose) {printf("Refining level %3d -> %3d: %5d cells to refine, ", levelcounter, levelcounter+1, nthislevel); }

    for (int n = 0; n < nthislevel; n++){
      refine(cells[thislevel_refine[n]]);
    }

  if (verbose) { printf("%5d children need refinement\n", nnextlevel); }


    // check whether another step is necessary
    if (nnextlevel > 0){

      // repeat on
      to_refine = 1;

      // copy next level cells into this level cell array
      thislevel_refine = realloc(thislevel_refine, nnextlevel*sizeof(int));
      for (int n = 0; n<nnextlevel; n++){
        thislevel_refine[n]=nextlevel_refine[n];
      }

      nthislevel = nnextlevel;

      // initialize new next level cell array
      nextlevel_refine = realloc(nextlevel_refine,nthislevel*8*sizeof(int));

      // reset number of next level particles
      nnextlevel = 0;
    }

  }



  // cleanup
  free(thislevel_refine);
  free(nextlevel_refine);


  printf("Finished refinement. Max level: %d, ncells: %d\n", max_refinement_level, lastcell);



// write to file
FILE *outfilep = fopen("cellcentres.dat", "w");



fprintf(outfilep, "%15s   %15s   %15s   %15s  \n", "x ", "y ", "z ", "cell ");
for (unsigned int i = 0; i<lastcell; i++){
  node *cell = cells[i];
  double xc = cell->center[0];
  double yc = cell->center[1];
  double zc = cell->center[2];

  fprintf(outfilep, "%15g   %15g   %15g   %15d \n", xc, yc, zc, i );
}


fclose(outfilep);







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



  double c = boxlen/4;
  
  // centres for level one:
  double cs[8][3] = { {-c, -c, -c}, {c, -c, -c}, {-c, c, -c}, {c, c, -c}, {-c, -c, c}, {c, -c, c}, {-c, c, c}, {c, c, c} };


  // refinement steps array initialisation
  thislevel_refine = calloc(8, sizeof(int));
  nextlevel_refine = calloc(64, sizeof(int));



  
  //------------------
  // Initialise root
  //------------------
  
  for (int i =0; i<8; i++) {

    int * nochild = malloc(8*sizeof(int));
    for (int j = 0; j<8; j++) { nochild[j] = -1; }

    double *thiscenter = malloc(3*sizeof(double));
    for (int j = 0; j<3; j++) { thiscenter[j] = cs[i][j]; }

    int *parray = malloc(npart * sizeof(int));

    node *root = malloc(sizeof(node));

    root->parent = -1;
    root->cellindex = i;
    root->level = 1;
    root->child = nochild;
    root->center = thiscenter;
    root->np = 0;
    root->particles = parray;


    cells[i] = root;

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
    cells[ind]->particles[ cells[ind]->np ] = p;

    // raise particle counter for root cell
    cells[ind]->np += 1;

  }


  // find out which root cells need to be refined

  for (int i = 0; i<8; i++){
  if (cells[i]->np > ncellpartmax){
      thislevel_refine[nthislevel] = i;
      nthislevel += 1;
    }
  }



}













//===========================
void refine(node * parent)
//===========================
{
  

  // if (verbose) {printf("Refining cell %d of level %d with npart %d\n", parent->cellindex, parent->level, parent->np);}



  //---------------------
  // Preparation
  //---------------------

  node ** newchildren = malloc(8 * sizeof(node*));

  double xp = parent->center[0];
  double yp = parent->center[1];
  double zp = parent->center[2];
  int npar = parent->np;
  int parind = parent->cellindex;


  // get child level, update max ref level so far
  int childlevel = parent->level + 1;

  if (max_refinement_level < childlevel) { 
    max_refinement_level = childlevel; 
  }

  
  // centres for child level:
  double c = boxlen/pow(2, childlevel+1);


  double cs[8][3] = { {xp-c, yp-c, zp-c}, {xp+c, yp-c, zp-c}, {xp-c, yp+c, zp-c}, {xp+c, yp+c, zp-c}, {xp-c, yp-c, zp+c}, {xp+c, yp+-c, zp+c}, {xp-c, yp+c, zp+c}, {xp+c, yp+c, zp+c} };

  



  //----------------------
  // Initialise children
  //----------------------
  
  for (int i =0; i<8; i++) {



    int * nochild = malloc(8*sizeof(int));
    for (int j = 0; j<8; j++) { nochild[j] = -1; }

    double *thiscenter = malloc(3*sizeof(double));
    for (int j = 0; j<3; j++) { thiscenter[j] = cs[i][j]; }

    int *parray = malloc(npart * sizeof(int));

    node *thischild = malloc(sizeof(node));


    thischild->parent = parind;
    thischild->cellindex = -1;
    thischild->level = childlevel;
    thischild->child = nochild;
    thischild->center = thiscenter;
    thischild->np = 0;
    thischild->particles = parray;

    newchildren[i] = thischild;
  }




  //-----------------------------------
  // DIVIDE PARTICLES AMONGST CHILDREN
  //-----------------------------------
  

  int i = 0;
  int j = 0;
  int k = 0;
  int ind;
  int thispart;


  node *thischild;
  for (int p = 0; p < npar; p++){

    thispart = parent->particles[p];
    i = 0;
    j = 0;
    k = 0;
    if (x[thispart] > xp){ i = 1; }
    if (y[thispart] > yp){ j = 1; }
    if (z[thispart] > zp){ k = 1; }
    ind = i + 2*j + 4*k;

    thischild = newchildren[ind];
    thischild->particles[thischild->np] = thispart;
    thischild->np += 1;

  }



  
 
  
  //------------------------------
  // ADD CHILDREN TO CELL ARRAY
  //------------------------------
  

  for (int child = 0; child<8; child++){
   
    // if it has any particle, add to list
    if (newchildren[child]->np > 0){

      // give child a unique index
      newchildren[child]->cellindex = lastcell;

      // copy child to cell array
      cells[lastcell] = newchildren[child];

      // tell parent what index its child has
      parent->child[child] = lastcell;

      // tell particles what cell they belong to
      // do it here, when cell has obtained an index
      for (int p = 0; p < cells[lastcell]->np; p++){
        partcell[cells[lastcell]->particles[p]] = lastcell;
      }


      if ( cells[lastcell]->np > ncellpartmax ){
        // this child needs to be refined
        nextlevel_refine[nnextlevel] = lastcell;
        nnextlevel += 1;
      }

      lastcell += 1;

    }
  }




  //-----------------
  // CLEANUP
  //-----------------
  
  free(newchildren);
  




}
