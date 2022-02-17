/* ========================================================================
   FILENAME:     pde_f_discretise.c

  PURPOSE:     all the functions to discretise a pde
   ======================================================================== */
/* ----------------------------------------------------------------------------------------------------------------------------------------------------- */
/* Define a discrete mesh in space, for each coordinate being within min/max point according to PDE */
#include "srt_h_all.h"
#include "pde_h_discretise.h"

Err pde_discretise_space(
				SrtPdeObject        *pde)
{
	  Err err = NULL;

	  SrtPdeGridObject  *grid   = pde->grid;

	  if (pde->dimension >=1)
	         pde_discretise_space_dim1(pde->mesh_type_x, grid->num_x, 
			                        pde->min_x, pde->max_x, grid->vec_x, grid->inc_x);
	  if (pde->dimension >=2)
	         pde_discretise_space_dim1(pde->mesh_type_y, grid->num_y, 
			                        pde->min_y, pde->max_y, grid->vec_y, grid->inc_y);
	  if (pde->dimension >=3)
	         pde_discretise_space_dim1(pde->mesh_type_z, grid->num_z, 
			                        pde->min_z, pde->max_z, grid->vec_z, grid->inc_z);
	  return err;
}
	
Err pde_discretise_space_dim1(
				SrtPdeMeshType    mesh_typ,
				int                         mx,
				double                   min,
				double                   max,
				double                   *vec,
				double                   *inc)
{
	Err err = NULL;
	int         i;
    double   dx, sigma;

    switch ( mesh_typ )
    {
        case REGULAR:
             dx = (max - min)/(mx + 1);
             for ( i = 0; i <= mx + 1; i++)
             {
                 vec[i] = min + i * dx;
             } 
             break;

        case GAUSSIAN_CENTRED:
			 sigma = 2.0*(max - min)/(2.0*CUTPDEGRIDATSTDEV);
             dx = 1.0/(mx + 1);
			 vec[0] = min; 
			 vec[mx + 1] = max;
             for ( i = 1; i <= mx; i++)
             {
                 vec[i] = sigma*inv_cumnorm_fast(i*dx);
             } 
             break;

        default:
             break;
    }
    
    for ( i = 0; i <= mx; i++)
    {
        inc[i] = vec[i+1] - vec[i];
    }

	return err;
}

/* ----------------------------------------------------------------------------------------------------------------------------------------------------- */
/* ----------------------------------------------------------------------------------------------------------------------------------------------------- */

Err pde_discretise_operator(
						SrtPdeObject	*pde,
						double dt)
{
	Err err = NULL;
    SrtPdeGridObject*    grid = pde->grid;
	double theta = pde->theta;
		
	int           i = 0, j = 0, k = 0, l = 0, m = 0;

    int           mxx = grid->num_x, 
		          myy = grid->num_y;

    double    dxd = 0.0, dx = 0.0, 								   
		          dyd = 0.0, dy = 0.0,
                  c0 = 0.0, cx = 0.0, cxx = 0.0, 
				  cy = 0.0, cyy = 0.0, cxy = 0.0,
				  aa, bb, cc;

    double     alpha = 0.0, 
		          sumd = 0.0, 
				  sumd2 = 0.0;

    double     sgn, g1, g2, /* lam , */ cobeta = 1000.0; 

	double point[3] = {0.0, 0.0, 0.0};

    /*------- Building the discretized operator in X ---------*/

    for (j = 0; j <= myy + 1 ; j++)
    {
        point[1] =  grid->vec_y[j]; 

        grid->BX[0][j][0] = 1.0;

        grid->CX[0][j][0] = - pde->bound_cond_x;

        for (i = 1 ; i <= mxx ; i++)
        {
            point[0] =  grid->vec_x[i];  

			cx = pde->coeff_X(pde->info, point, theta); 	
			cxx = pde->coeff_XX(pde->info, point, theta); 

			dxd = grid->inc_x[i-1];
			dx = grid->inc_x[i];

            alpha = (dx - dxd)/(dx + dxd);
	        sumd2 = (dx*dx + dxd*dxd)/2.0;

			if (pde->scheme == THETASCHEME_DECENTRED)
			{
				/*
				lam = 0.5*(1.0 + cobeta*cx/(1.0 + cobeta*fabs(cx))); 
				g1 = cx*(1.0-lam)/dxd;
				g2 = cx*lam/dx;	  
				
				aa = dt * (cxx * (1.0 + alpha)/sumd2 - g1);
				bb = dt *(- 2.0 * cxx/sumd2 + g1 - g2);
				cc = dt * (cxx * (1.0 - alpha)/sumd2 + g2);	
			
				sgn = (cx >0) ? 1.0 : -1.0;
				g1 = cx*0.5*(1.0 - sgn)/dxd;           UPWIND
				g2 = cx*0.5*(1.0 + sgn)/dx;	
				
				g1 = cx*0.5*(1.0 + sgn)/dxd;		  DOWNWIND
				g2 = cx*0.5*(1.0 - sgn)/dx;	 
				*/

				sgn = (cx >0) ? 1.0 : -1.0;
				g1 = cx*0.5*(1.0 - sgn)/dxd;
				g2 = cx*0.5*(1.0 + sgn)/dx;	 

				aa = dt * (cxx * (1.0 + alpha)/sumd2 - g1);
				bb = dt *(- 2.0 * cxx/sumd2 + g1 - g2);
				cc = dt * (cxx * (1.0 - alpha)/sumd2 + g2);			

			}
			else
			{
				sumd = dx + dxd;

				aa = dt * (cxx * (1.0 + alpha)/sumd2 - cx/sumd);
				bb = - 2.0 * dt * cxx/sumd2;
				cc = dt * (cxx * (1.0 - alpha)/sumd2 + cx/sumd);
			}

            /*----------*/

            grid->AX[i][j][0] = theta*aa;
            grid->BX[i][j][0] = 1.0 + theta*bb;
            grid->CX[i][j][0] = theta*cc;

            grid->AX1[i][j][0] = -(1.0 - theta)*aa;
            grid->BX1[i][j][0] = 1.0 - (1.0 - theta)*bb;
            grid->CX1[i][j][0] = -(1.0 - theta)*cc;

            /*----------*/
        }

        grid->AX[mxx+1][j][0] = - pde->bound_cond_x;

        grid->BX[mxx+1][j][0] = 1.0;
    }


    /*------- Building the discretized operator in Y ---------*/

    for (i = 0; i <= mxx +1; i++)
    {
        point[0] =  grid->vec_x[i]; 

        grid->BY[i][0][0] = 1.0;

        grid->CY[i][0][0] = - pde->bound_cond_y;

        for (j = 1 ; j <= myy ; j++)
        {
            point[1] =  grid->vec_y[j]; 

			c0 = pde->coeff_0(pde->info, point, theta);
			cy = pde->coeff_Y(pde->info, point, theta); 
			cyy = pde->coeff_YY(pde->info, point, theta); 

			dyd = grid->inc_y[j-1];
			dy = grid->inc_y[j];

            alpha = (dy - dyd)/(dx + dyd);
			sumd2 = (dy*dy + dyd*dyd)/2.0;

			if (pde->scheme == THETASCHEME_DECENTRED)
			{
				sgn = (cy > 0) ? 1.0 : -1.0;
				g1 = cy*0.5*(1.0 - sgn)/dyd;
				g2 = cy*0.5*(1.0 + sgn)/dy;

				aa = dt * (cyy * (1.0 + alpha)/sumd2 - g1);
				bb = dt *(c0 - 2.0 * cyy/sumd2 + g1 - g2);
				cc = dt * (cyy * (1.0 - alpha)/sumd2 + g2);	 

				/*
				am = 0.5*(1.0 + cobeta*cy/(1.0 + cobeta*fabs(cy))); 
				g1 = cy*(1.0-lam)/dyd;
				g2 = cy*lam/dy;	  
				
				aa = dt * (cyy * (1.0 + alpha)/sumd2 - g1);
				bb = dt *(c0 - 2.0 * cyy/sumd2 + g1 - g2);
				cc = dt * (cyy * (1.0 - alpha)/sumd2 + g2);
				*/
			}
			else
			{
				sumd = dy + dyd;

				aa =dt * (cyy * (1.0 + alpha)/sumd2 - cy/sumd);
				bb = dt * (c0 - 2.0 * cyy/sumd2);
				cc = dt * (cyy * (1.0 - alpha)/sumd2 + cy/sumd);
			}

            /*----------*/

            grid->AY[i][j][0] = theta*aa;
            grid->BY[i][j][0] = 1.0 + theta*bb;
            grid->CY[i][j][0] = theta*cc;

            grid->AY1[i][j][0] = -(1.0 - theta)*aa;
            grid->BY1[i][j][0] = 1.0 - (1.0 - theta)*bb;
            grid->CY1[i][j][0] = -(1.0 - theta)*cc;

           /*----------*/

        }

        grid->AY[i][myy+1][0] = - pde->bound_cond_y;

        grid->BY[i][myy+1][0] = 1.0;
    }


    /*------- Building the discretized operator in XY ---------*/

    for (i = 1 ; i <= mxx ; i++)
    {
        point[0] =  grid->vec_x[i];

        dxd = grid->inc_x[i-1];
        dx = grid->inc_x[i];

        for (j = 1 ; j <= myy ; j++)
        {
            point[1] =  grid->vec_y[j]; 

			dyd = grid->inc_y[j-1];
			dy = grid->inc_y[j];

			cxy = pde->coeff_XY(pde->info, point, theta);

            /*----------*/

            grid->EXY[i][j][0] = - dt * cxy/(2.0 * (dy + dyd) * (dx + dxd));

            /*----------*/

        }
    }

    /*---------- Decomposition of the discretized operator X and Y
	                           for implicit solving                                  ----------*/


    for (j = 0 ; j <= myy +1; j++)
    {
        grid->D1[0][j][0] = grid->BX[0][j][0];

        for (i = 1 ; i <= mxx + 1 ; i++)
        {
            grid->L1[i][j][0] = grid->AX[i][j][0] / grid->D1[i-1][j][0];

            grid->D1[i][j][0] = grid->BX[i][j][0] 
                       - grid->L1[i][j][0] * grid->CX[i-1][j][0];
        }
    }

    for (i = 0 ; i <= mxx + 1 ; i++)
    {
        grid->D2[i][0][0] = grid->BY[i][0][0];

        for (j = 1 ; j <= myy + 1 ; j++)
        {
            grid->L2[i][j][0] = grid->AY[i][j][0] / grid->D2[i][j-1][0];

            grid->D2[i][j][0] = grid->BY[i][j][0] 
                         - grid->L2[i][j][0] * grid->CY[i][j-1][0];
        }
    }

	return err;
}


/*-------------------------------------------------------------------*/
/*--- End of File ---*/

