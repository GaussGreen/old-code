/* ========================================================================
   FILENAME:     pde_h_discretise.h

  PURPOSE:     all the functions to discretise a pde
   ======================================================================== */

#ifndef PDE_H_DISCRETISE
#define PDE_H_DISCRETISE

#define CUTPDEGRIDATSTDEV 5.0

#include "utallhdr.h"
#include "pde_h_struct.h"

/* ----------------------------------------------------------------------------------------------------------------------------------------------------- */

Err pde_discretise_space(
				SrtPdeObject        *pde);


Err pde_discretise_space_dim1(
				SrtPdeMeshType         mesh_typ,
				int                         mx,
				double                   min,
				double                   max,
				double                   *vec,
				double                   *inc); 


Err pde_discretise_operator(
						SrtPdeObject	*pde,
						double dt);


 #endif
 /*-------------------------------------------------------------------*/
/*--- End of File ---*/

