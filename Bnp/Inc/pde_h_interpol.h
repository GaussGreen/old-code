/* ========================================================================
   FILENAME:     utpdeinterpol.h

  PURPOSE:     interpol the solution of the PDE on the new mesh
   ======================================================================== */

#ifndef PDE_H_INTERPOL
#define PDE_H_INTERPOL

#include "pde_h_struct.h"
#include "utallhdr.h"

/* ======================================================================== */

Err pde_interpol(SrtPdeObject *pde, double M[2][2], double M1[2][2]);

#endif
/*-------------------------------------------------------------------*/
/*--- End of File ---*/
