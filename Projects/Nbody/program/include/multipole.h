//=====================
// Multipole Method
//=====================




//=======================
// Structs
//=======================



//------------------
// tree node
//------------------

typedef struct {

  int parent;
  int cellindex;
  int level;
  int * child; // contains 8 children
  double * center; // x, y, z of cell centre

  int np; // number of particles
  int * particles; // list of particle indices

  // multipole expansion stuff
  double * centre_of_mass;
  double diagonal;
  double mass;
  double * multip_vector;
  double ** multip_matrix;
  double multip_sq;


} node;




//=======================
// Global vars and arrays
//=======================
extern unsigned int lastcell; // index of first free cell
extern node ** cells;           // array of all cells of all levels
extern int * partcell;        // most refined cell a particle belongs to

extern int * thislevel_refine; // cells to refine on this level
extern int * nextlevel_refine; // cells to refine on next level
extern int nthislevel;
extern int nnextlevel;
extern int to_refine; // whether there is a next level to refine
extern int multipole_order; // up to which multipole order to use






//======================
// Functions
//======================

extern void build_tree();
extern void get_multipole(int index);
extern void calculate_multipole_forces(int index);
