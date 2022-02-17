/* ========================================================================
   FILENAME:     utpdegrid.h

  PURPOSE:     Allocation and Desallocation functions  for the work space
   ======================================================================== */

#ifndef PDE_H_GRID
#define PDE_H_GRID

#include "pde_h_struct.h"
#include "utallhdr.h"

/* ----------------------------------------------------------------------------------------------------------------------------------------
 */
/* allocation of workspace for pde solving */

Err pde_new_grid(SrtPdeObject *pde, long num_points[3]);

Err pde_delete_grid(SrtPdeGridObject *grid);

#endif
/*-------------------------------------------------------------------*/
/*--- End of File ---*/
