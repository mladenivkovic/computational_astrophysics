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

int multipole_order = 2;




//================
//Functions
//================

void build_root();
void refine(node * parentnode);
void calc_direct(node *target, node *source);
void walk_tree(node * target, node * source);
double theta(double y, node * source);
void calc_multipole(double d[3], node * target, node * source);
int is_ancestor(node* target, node* source);




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

  // for (int i = 0; i<nthislevel; i++){
  //   printf("Child%5d has %5d particles\n", thislevel_refine[i], cells[thislevel_refine[i]]->np);
  // }
  // printf("\n");
    if (verbose) {printf("Refining level %3d -> %3d: %5d cells to refine, ", levelcounter, levelcounter+1, nthislevel); }

    for (int n = 0; n < nthislevel; n++){
      refine(cells[thislevel_refine[n]]);
    }

    if (verbose) { 
      printf("%5d children need refinement. Ncells = %d\n", nnextlevel, lastcell); 
    }



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

  } //while






  //--------------
  // cleanup
  //--------------
  
  free(thislevel_refine);
  free(nextlevel_refine);


  printf("Finished refinement. Max level: %d, ncells: %d\n", max_refinement_level, lastcell);



  // // write to file
  // FILE *outfilep = fopen("cellcentres.dat", "w");
  //
  //
  //
  // fprintf(outfilep, "%15s   %15s   %15s   %15s  \n", "x ", "y ", "z ", "cell ");
  // for (unsigned int i = 0; i<lastcell; i++){
  //   node *cell = cells[i];
  //   double xc = cell->center[0];
  //   double yc = cell->center[1];
  //   double zc = cell->center[2];
  //
  //   fprintf(outfilep, "%15g   %15g   %15g   %15d \n", xc, yc, zc, i );
  // }
  // fclose(outfilep);
  //
  //
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



  // assuming simulation is centered around origin
  double c = boxlen/4;
  
  // centres for level one:
  double cs[8][3] = { 
    {-c, -c, -c}, 
    {c, -c, -c}, 
    {-c, c, -c}, 
    {c, c, -c}, 
    {-c, -c, c}, 
    {c, -c, c}, 
    {-c, c, c}, 
    {c, c, c} 
  };


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

    // generate empty matrix
    double ** matrix = malloc(3*sizeof(double*));
    for (int j = 0; j<3; j++){
      double * column = calloc(3 , sizeof(double));
      matrix[j] = column;
    }



    node *root = malloc(sizeof(node));

    root->parent = -1;
    root->cellindex = i;
    root->level = 1;
    root->isleaf = 1;
    root->child = nochild;
    root->center = thiscenter;
    root->np = 0;
    root->particles = malloc(npart * sizeof(int));
    root->centre_of_mass = calloc(3, sizeof(double));
    root->groupsize = 0; 
    root->mass = 0;
    root->multip_vector = calloc(3, sizeof(double));
    root->multip_matrix = matrix;
    root->multip_sq = 0;


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
    if (x[p] >= 0){ i = 1; }
    if (y[p] >= 0){ j = 1; }
    if (z[p] >= 0){ k = 1; }

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
      cells[i]->isleaf = 0;
      thislevel_refine[nthislevel] = i;
      nthislevel += 1;
    }
  }


}













//===========================
void refine(node * parent)
//===========================
{

  //==============================================
  // Refine a parent node into up to 8 children,
  // that contain at least 1 particle
  //==============================================
  
  // if (verbose) {printf("Refining cell %d of level %d\n", parent->cellindex, parent->level);}


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
  double c = boxlen/pow(2.0, (float) (childlevel+1));

  double cs[8][3] = { 
    {xp-c, yp-c, zp-c}, 
    {xp+c, yp-c, zp-c}, 
    {xp-c, yp+c, zp-c}, 
    {xp+c, yp+c, zp-c}, 
    {xp-c, yp-c, zp+c}, 
    {xp+c, yp-c, zp+c}, 
    {xp-c, yp+c, zp+c}, 
    {xp+c, yp+c, zp+c} 
  };

  



  //----------------------
  // Initialise children
  //----------------------
  
  for (int i =0; i<8; i++) {



    int * nochild = malloc(8*sizeof(int));
    for (int j = 0; j<8; j++) { nochild[j] = -1; }

    double *thiscenter = malloc(3*sizeof(double));
    for (int j = 0; j<3; j++) { thiscenter[j] = cs[i][j]; }


    // generate empty matrix
    double ** matrix = malloc(3*sizeof(double*));
    for (int j = 0; j<3; j++){
      double * column = calloc(3 , sizeof(double));
      matrix[j] = column;
    }



    node *thischild = malloc(sizeof(node));

    thischild->parent = parind;
    thischild->cellindex = -1;
    thischild->level = childlevel;
    thischild->isleaf = 1;
    thischild->child = nochild;
    thischild->center = thiscenter;
    thischild->np = 0;
    thischild->particles = malloc(npart * sizeof(int));
    thischild->centre_of_mass = calloc(3, sizeof(double));
    thischild->groupsize = 0;
    thischild->mass = 0;
    thischild->multip_vector = calloc(3, sizeof(double));
    thischild->multip_matrix = matrix;
    thischild->multip_sq = 0;


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
    if (x[thispart] >= xp){ i = 1; }
    if (y[thispart] >= yp){ j = 1; }
    if (z[thispart] >= zp){ k = 1; }
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
      int newind;
#pragma omp critical
      {
        newind = lastcell;
        lastcell += 1;
      }
      newchildren[child]->cellindex = newind;

      // copy child to cell array
      cells[newind] = newchildren[child];

      // tell parent what index its child has
      parent->child[child] = newind;

      // tell particles what cell they belong to
      // do it here, when cell has obtained an index
      for (int p = 0; p < cells[newind]->np; p++){
        partcell[cells[newind]->particles[p]] = newind;
      }


      if ( cells[newind]->np > ncellpartmax ){
        // this child needs to be refined
        cells[newind]->isleaf = 0;
        nextlevel_refine[nnextlevel] = newind;
        nnextlevel += 1;
      }


    }
  }




  //-----------------
  // CLEANUP
  //-----------------
  
  free(newchildren);
  

}





















//=====================================
void get_multipole(int index)
//=====================================
{

  //==================================================
  // This function calculates all necessary parts for
  // the multipoles.
  //==================================================
  
  node * thiscell = cells[index];


  // loop over children recursively

  int isleaf = thiscell->isleaf;
  int thispart = 0;
  double com[] = {0,0,0};

  if (!isleaf){
    for (int i = 0; i < 8; i++){
      if (thiscell->child[i] > 0){
        get_multipole(thiscell->child[i]);
      }
    }
  }
  else{
    // explicitly compute centre of mass for leaf cells
    for (int p = 0; p < thiscell->np; p++){
      thispart = thiscell->particles[p];
      com[0]+=m[thispart]*x[thispart];
      com[1]+=m[thispart]*y[thispart];
      com[2]+=m[thispart]*z[thispart];
#pragma omp atomic
      thiscell->mass += m[thispart];
    }

#pragma omp critical
    {
     for (int i = 0; i<3; i++){
        thiscell->centre_of_mass[i] = com[i];
      }
    }
  }

  
  // pass centre of mass data to parent, if parent exists
  if (thiscell->parent >= 0) {
    node * parent = cells[thiscell->parent];
#pragma omp atomic
    parent->mass += thiscell->mass;
    for (int i = 0; i<3; i++){
#pragma omp atomic
      parent->centre_of_mass[i] += thiscell->centre_of_mass[i]; 
    }
  }
  
  


  if (thiscell->np > 0){

    //--------------------------------
    // Get Centre of Mass
    //--------------------------------
    for (int i = 0; i<3; i++){
#pragma omp atomic
      thiscell->centre_of_mass[i] /= thiscell->mass;
    }




    //---------------------------------------------
    // calculate ( s - x_i ) vector related stuff
    //---------------------------------------------
    
    double vec[3] = {0, 0, 0};
   
    for (int p = 0; p < thiscell->np; p++){
      int pind = thiscell->particles[p];

      vec[0] = thiscell->centre_of_mass[0] - x[pind];
      vec[1] = thiscell->centre_of_mass[1] - y[pind];
      vec[2] = thiscell->centre_of_mass[2] - z[pind];

      //store vector
      for (int i = 0; i<3; i++){
#pragma omp atomic
        thiscell->multip_vector[i] += vec[i]; 
      }
      
      //get and store square
      double v =  pow(vec[0], 2) + pow(vec[1], 2) + pow(vec[2], 2);
#pragma omp atomic
      thiscell->multip_sq += v;

      // get group size
      v = sqrt(v);
#pragma omp critical
      if (v > thiscell->groupsize){
        thiscell->groupsize = v;
      }
      
      //get and store matrix
      // TODO: doublecheck
#pragma omp critical
      {
        for (int i = 0; i<3; i++){
          for (int j = 0; j<3; j++){
            thiscell->multip_matrix[i][j] += vec[i]*vec[j];
          }
        }
      }

  // if (thiscell->cellindex == 8){
  // printf("Got vec: %g %g %g\n", vec[0], vec[1], vec[2]);
  // printf("Multipole matrix:\n");
  // double check = 0;
  // for (int j = 0; j<3; j++){
  //   for (int k = 0; k<3; k++){
  //     printf("%7g ", thiscell->multip_matrix[j][k]);
  //   }
  //   printf("\n");
  //   check += thiscell->multip_matrix[j][j];
  // }
  //   printf("mp vector square: %g, check: %g\n\n", thiscell->multip_sq, check);
  // }

    }
  }
}












//=========================================
void calculate_multipole_forces(int index)
//=========================================
{
  //=====================================
  // Descend until you find the leaves.
  // Then calculate force.
  //=====================================


  node * thiscell = cells[index];


  //-------------------------------------
  // loop over children recursively
  // until you reach leaves
  //-------------------------------------

  int isleaf = thiscell->isleaf;

  
  if (!isleaf){
    for (int i = 0; i < 8; i++){
      // if there is a child: descend into child
      if (thiscell->child[i] > 0){
        calculate_multipole_forces(thiscell->child[i]);
      }
    }
  }
  else{
    //------------------------------------
    // calculate forces if this is a leaf
    // walk the tree, starting with root
    //------------------------------------

    for (int root = 0; root < 8; root++){
      if (cells[root]->np > 0){
        walk_tree(thiscell, cells[root]);
      }
    }
  }
}








//==========================================
void walk_tree(node * target, node * source)
//==========================================
{

  //============================================
  // Walk the tree and calculate the forces.
  //============================================
 



  if (is_ancestor(target, source)){
    // if the source contains the target
    if (target->cellindex == source->cellindex){
      // if source is target: calculate direct if necessary, then return
      if (target->np > 1){
        calc_direct(target, target);
      }
      return;
    }
    else{ // move along.
      for (int c = 0; c < 8; c++){
        if (source->child[c]>0){
          walk_tree( target, cells[source->child[c]] );
        }
      }
      return;
    }
    

  }

  
  // Otherwise: do your thing
  double absd, absd_sq;
  double d[3] = {0, 0, 0};

  // get distance
  absd_sq = 0;
  for (int i = 0; i<3; i++){
    d[i] = target->centre_of_mass[i] - source->centre_of_mass[i];
    absd_sq += pow( d[i], 2.0 );
  }
  
  if (absd_sq > 0){

    absd = sqrt(absd_sq);

    if (theta(absd, source) <= theta_max){ 
      // if multipole approx condition satisfied
      calc_multipole(d, target, source); 
    }
    else {
      // check whether source is a leaf cell
      if (!(source->isleaf)){
        for (int c = 0; c < 8; c++){
          if (source->child[c] > 0){
            // if it isn't leaf, try children for multipole approach
            walk_tree( target, cells[source->child[c]]);
          }
        }
      }
      else {
        // if it was a leaf cell, do direct force calculation
        calc_direct(target, source);
      }
    }

  }
  else{
    printf("!!!!!!!!!!!! I SHOULDN'T BE HERE!!!!!!!!!!!!\n");
  }
  

  return;

}







//============================================
double theta(double y, node * source)
//============================================
{
  //===================================================
  // Calculate approximate angle of source wrt target 
  //===================================================
  
  double theta = source->groupsize / y;
  return (theta);
}










//===================================================================
void calc_multipole(double d[3], node * target, node * source)
//===================================================================

{
 
  //==================================================================
  // This function does the actual force calculation of the multipole
  //==================================================================
  

  double absd, absd_sq=0, absd_cube, absd_five, absd_seven;

  for (int i = 0; i<3; i++){
    absd_sq += pow(d[i], 2.0);
  }
  absd = sqrt(absd_sq);
  absd_cube = pow(absd, 3.0);
  absd_five = pow(absd, 5.0);
  absd_seven = pow(absd, 7.0);



  double force[3] = {0, 0, 0};

  
  //-------------------
  // calc monopole 
  //-------------------
  for (int j = 0; j<3; j++){
    force[j] -= d[j]/absd_cube;
  }


  // if (multipole_order > 0){
  //   //-------------------
  //   // calc dipole
  //   //-------------------
  //   double scalar_product = 0;
  //   for (int k = 0; k<3; k++){
  //     scalar_product += d[k] * source->multip_vector[k];
  //   }
  //   for (int j = 0; j <3; j++){
  //     force[j] += source->multip_vector[j]/absd_cube - 3.0* scalar_product * d[j]/ absd_five;
  //   }
  // }


  if (multipole_order > 1){
    //-------------------
    // Calc quadrupole
    //-------------------
    double sum1[3] = {0, 0, 0}; 
    double sum2 = 0;

    for (int j = 0; j<3; j++){
      for (int k = 0; k<3; k++){
        sum1[j] += 3*(source->multip_matrix[k][j])*d[j];
      }
    }

    for (int k = 0; k<3; k++){
      sum2 += sum1[k]*d[k];
    }

    for (int j = 0; j<3; j++){
      force[j] += 
        (sum1[j] - (source->multip_sq)*d[j])/absd_five -
          2.5 * (sum2 - absd_sq*(source->multip_sq))* d[j] / absd_seven;
    }
  }




  
  //--------------------------------------
  // apply calculated force to particles
  //--------------------------------------
  int pind;
#pragma omp critical
  {
    for (int p = 0; p<target->np; p++){
      pind = target->particles[p];
      fx[pind] += force[0] * m[pind] * source->mass;
      fy[pind] += force[1] * m[pind] * source->mass;
      fz[pind] += force[2] * m[pind] * source->mass;
    }
  }


}









//===========================================================
void calc_direct(node *target, node *source)
//===========================================================
{
  //===================================================
  // Calculate the direct forces between particles
  // of list1 with length len1 and particles of list2
  // with length len2
  // works for same-cell same-cell interactions and
  // different cells interacting (e.g. neighbours)
  //===================================================


  int len_tar = target->np;
  int len_src = source->np;

  int *plist_tar = target->particles;
  int *plist_src = source->particles;

  double force_fact, rsq, dx, dy, dz;
  int p_i, p_j;
  for (int i = 0; i < len_tar; i++){
    for (int j = 0; j < len_src; j++){
      p_i = plist_tar[i];
      p_j = plist_src[j];

      if (p_i != p_j) {
        dx = x[p_i]-x[p_j];
        dy = y[p_i]-y[p_j];
        dz = z[p_i]-z[p_j];

        rsq = dx*dx + dy*dy + dz*dz;

        force_fact = - m[p_i]*m[p_j] / pow(rsq + softening*softening, 1.5);
#pragma omp critical
        {
          fx[p_i] += force_fact * (x[p_i]-x[p_j]);
          fy[p_i] += force_fact * (y[p_i]-y[p_j]);
          fz[p_i] += force_fact * (z[p_i]-z[p_j]);
        }
      }
    }
  }
}














//=================================================
int is_ancestor(node* target, node* source)
//=================================================
{
  //=======================================
  // checks if source is ancestor of target
  //=======================================
  
  int tlev = target->level;
  int slev = source->level;
  int srcind = source->cellindex;
  int ancind = target->cellindex; //initialise ancestor index as first parent



  for (int i = slev; i < tlev; i++){
    ancind = cells[ancind]->parent;
  }

  if (srcind == ancind){
    return 1;
  }
  else{
    return 0;
  }
}
