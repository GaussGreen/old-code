/* ========================================================================
   FILENAME:     utpdesetup.h

  PURPOSE:     Allocation and setup functions
   ======================================================================== */

#ifndef PDE_H_SETUP
#define PDE_H_SETUP

#include "pde_h_struct.h"
#include "utallhdr.h"

/*-------------------------------------------------------------------*/

Err pde_set_up(
    int                 dimension,
    int                 nb_payoffs,
    int                 nb_second_members,
    SrtPdeMeshType      mesh_type[3],
    SrtPdeBoundaryCond  bound_cond[3],
    SrtPdeSolvingScheme scheme,
    double              theta,
    long                num_points[3],
    void*               info,
    PdeSetUpModelFunc   the_pde_set_up,
    SrtPdeObject*       pde);

Err pde_clear(SrtPdeObject* pde);

Err pde_attach_info_model(SrtPdeObject* pde, void* info);

#endif
/*-------------------------------------------------------------------*/
/*--- End of File ---*/
