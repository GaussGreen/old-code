/* ========================================================================
   FILENAME:     utpdesetup.c

  PURPOSE:     Allocation and setup functions
   ======================================================================== */

#include "pde_h_grid.h"
#include "pde_h_setup.h"

/* For Solving PDE up to three space dimensions */
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
    SrtPdeObject*       pde)
{
    Err  err = NULL;
    long taillealloc;

    /* Associate inputs with Structures fields */
    pde->dimension = dimension;

    pde->nb_payoffs = nb_payoffs;

    pde->nb_second_members = nb_second_members;

    taillealloc = sizeof(SrtPdeGridObject);
    pde->grid   = (SrtPdeGridObject*)malloc(taillealloc);
    err         = pde_new_grid(pde, num_points);

    pde->mesh_type_x = mesh_type[0];
    pde->mesh_type_y = mesh_type[1];
    pde->mesh_type_z = mesh_type[2];

    pde->bound_cond_x = bound_cond[0];
    pde->bound_cond_y = bound_cond[1];
    pde->bound_cond_z = bound_cond[2];

    pde->scheme = scheme;
    pde->theta  = theta;

    pde->info = info;

    /* updating info pointer function (for time, etc)
    pde->pde_attach_info_model = &pde_attach_info_model;  */

    /* the procedure of change-of-variables and localization (space boundaries setting)
    pde->pde_set_up = the_pde_set_up;
    err = (*(pde->pde_set_up))(pde);	*/

    the_pde_set_up(pde);

    return err;
} /* END Err pde_set_up(...) */

/* ----------------------------------------------------------------------------------------------------------------------------------------
 */
/* ----------------------------------------------------------------------------------------------------------------------------------------
 */

Err pde_clear(SrtPdeObject* pde)
{
    Err err = NULL;
    err     = pde_delete_grid(pde->grid);
    free(pde->grid);

    return err;
} /* END Err pde_clear(...) */

/* ----------------------------------------------------------------------------------------------------------------------------------------
 */
/* updating of info pointer   */

Err pde_attach_info_model(SrtPdeObject* pde, void* info)
{
    Err err   = NULL;
    pde->info = info;

    return err;
}

/*-------------------------------------------------------------------*/
/*--- End of File ---*/
