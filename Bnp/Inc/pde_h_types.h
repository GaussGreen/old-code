/* ========================================================================
   FILENAME:     utpdetypes.h
   PURPOSE:     all the types used when playing with PDE's
   ======================================================================== */
#ifndef PDE_H_TYPES
#define PDE_H_TYPES

/* Which method used for the PDE */
typedef enum {

  EXPLICIT,
  THETASCHEME_CENTRED,
  THETASCHEME_DECENTRED,
  SHUPACK
} SrtPdeSolvingScheme;

typedef enum { EXPONENTIAL, REGULAR, GAUSSIAN_CENTRED } SrtPdeMeshType;

typedef enum { DIRICHLET = 0, NEUMANN = 1 } SrtPdeBoundaryCond;

typedef enum { LINEAR = 0, NONLINEAR = 1 } SrtPdePDEType;
#endif
/*-------------------------------------------------------------------*/
/*--- End of File ---*/
