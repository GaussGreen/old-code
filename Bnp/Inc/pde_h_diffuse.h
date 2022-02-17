/* ========================================================================
   FILENAME:     utpdediffuse.h

  PURPOSE:     solves the PDE for one time step (diffuse it)
   ======================================================================== */

#ifndef PDE_H_DIFFUSE
#define PDE_H_DIFFUSE

#include "pde_h_struct.h"
#include "utallhdr.h"

/* ======================================================================== */

Err pde_onestep_diffuse(SrtPdeObject *pde, int i);

void swap_pde_new_old(double ******U, double ******V);

Err solve_in_X(SrtPdeObject *pde);

Err solve_in_Y(SrtPdeObject *pde);

#endif
/*-------------------------------------------------------------------*/
/*--- End of File ---*/
