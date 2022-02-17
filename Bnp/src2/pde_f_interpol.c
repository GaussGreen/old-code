/* ========================================================================
   FILENAME:     utpdeinterpol.c

  PURPOSE:     2D interpolation

  AUTHOR:        Eric Fournie

  CREATION:     22 - january - 1998

  MODIFICATION:     22 - january - 1998
   ======================================================================== */

#include "pde_h_interpol.h"
#include "pde_h_diffuse.h"
#include "srt_h_all.h"

/* ======================================================================== */

Err pde_interpol(
						SrtPdeObject	*pde,
						double   M[2][2], double  M1[2][2])
{
	Err err = NULL;
	SrtPdeGridObject*  grid = pde->grid;
	SrtStpPtr 			              cur;				                           /* Current time step */
	SrtTwoFacPdeStepInfo 	*cur_pde_info, *next_pde_info;    /* PdeInfo at current and next step */
    int  mxx = grid->num_x, myy = grid->num_y, nbcoeff = grid->nb_coeffs;

	int  i = 0, j = 0, k = 0, m = 0, l = 0;
	int  ip, jp, imin, imax, jmin, jmax;
	double xp, yp, ximin, ximax, yjmin, yjmax, xk, yk;
	double thetax, thetay, ux, ux1, uxy;


	cur = (SrtStpPtr) pde->info; 
	cur_pde_info = (SrtTwoFacPdeStepInfo*) (cur->trinf);
   	next_pde_info = (SrtTwoFacPdeStepInfo*) (cur->next->trinf);


	for (i = 0 ; i <= mxx + 1; i++)
	{
		for (j = 0 ; j <= myy + 1; j++)
		{
				xk =  M1[0][0]*grid->vec_x[i]	+  M1[0][1]*grid->vec_y[j];
				yk =  M1[1][0]*grid->vec_x[i]	+  M1[1][1]*grid->vec_y[j];
				xp =  M[0][0]*xk	+  M[1][0]*yk;
				yp =  M[0][1]*xk	+  M[1][1]*yk;  

				/* --------------- */
				/* BASIC DICHOTOMIE TO FIND INDEX OF THE INTERVAL CONTAINING xp, yp */
				imin =  0;
				imax = mxx +1;
				while ((imax - imin) > 1)
				{
					ximax =grid->vec_x[imax]; 
					ximin =grid->vec_x[imin]; 
					k = (imax + imin)/2;
					xk = grid->vec_x[k];
					if (xp > xk) imin = k;
					else imax = k;
				}
				ip = imin;

				jmin =  0;
				jmax = myy +1;
				while ((jmax - jmin) > 1)
				{
					yjmax =grid->vec_y[jmax]; 
					yjmin =grid->vec_y[jmin]; 
					k = (jmax + jmin)/2;
					yk = grid->vec_y[k];
					if (yp > yk) jmin = k;
					else jmax = k;
				}
				jp = jmin;

				/* --------------- */

				if (xp >= grid->vec_x[mxx + 1])		ip = mxx;
				if (xp <= grid->vec_x[0])					ip = 0;
 				if (yp >= grid->vec_y[mxx + 1])		jp = myy;
				if (yp <= grid->vec_y[0])					jp = 0;

				/* --------------- */

				thetax = (xp - grid->vec_x[ip])/(grid->vec_x[ip+1] - grid->vec_x[ip]);
				ux1 = thetax * grid->new_value[0][ip+1][jp+1][0][0] + (1.0 - thetax) * grid->new_value[0][ip][jp+1][0][0];
				ux = thetax * grid->new_value[0][ip+1][jp][0][0] + (1.0 - thetax) * grid->new_value[0][ip+1][jp][0][0];

				thetay = (yp - grid->vec_y[jp])/(grid->vec_y[jp+1] - grid->vec_y[jp]);
				uxy = thetay * ux1 + (1.0 - thetay) * ux;

				grid->old_value[0][i][j][0][0] = uxy; 
		}
	}

	/*---- swap the two grids old and new ----*/

	swap_pde_new_old(&(grid->new_value), &(grid->old_value));

    return err;
} /* END Err  pde_onestep_diffuse(...) */

/*-------------------------------------------------------------------*/
/*--- End of File ---*/
