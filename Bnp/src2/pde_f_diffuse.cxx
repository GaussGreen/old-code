/* ========================================================================
   FILENAME:     utpdediffuse.cxx

  PURPOSE:     solves the PDE for one time step (diffuse it)

  AUTHOR:        Eric Fournie

  CREATION:     4 - january - 1998

  MODIFICATION:     14 - january - 1998
   ======================================================================== */

#include "pde_h_diffuse.h"

/* ======================================================================== */

Err pde_onestep_diffuse(SrtPdeObject *pde, int i) {
  Err err = NULL;
  SrtPdeGridObject *grid = pde->grid;

  swap_pde_new_old(&(grid->new_value), &(grid->old_value));

  if (i == 1) {
    err = solve_in_X(pde);

    swap_pde_new_old(&(grid->new_value), &(grid->old_value));

    err = solve_in_Y(pde);
  } else {
    err = solve_in_Y(pde);

    swap_pde_new_old(&(grid->new_value), &(grid->old_value));

    err = solve_in_X(pde);
  }

  return err;

} /* END Err  pde_onestep_diffuse(...) */

/*------------------------------------------------------------------------------------------*/

Err solve_in_Y(SrtPdeObject *pde) {
  Err err = NULL;
  SrtPdeGridObject *grid = pde->grid;
  int mxx = grid->num_x, myy = grid->num_y, nbcoeff = grid->nb_coeffs;

  int i = 0, j = 0, k = 0, m = 0, l = 0;

  for (m = 0; m < grid->nb_payoffs; m++) {
    /*---- step 1: solving the system in Y ----*/

    for (l = 0; l < nbcoeff; l++) {
      /*------- computation of the explicit part of the theta scheme -------*/

      for (i = 1; i <= mxx; i++) {
        for (j = 1; j <= myy; j++) {
          grid->new_value[m][i][j][0][l] =
              grid->AX1[i][j - 1][0] * grid->old_value[m][i][j - 1][0][l] +
              grid->BX1[i][j][0] * grid->old_value[m][i][j][0][l] +
              grid->CX1[i][j + 1][0] * grid->old_value[m][i][j + 1][0][l];
        }
      }

      swap_pde_new_old(&(grid->new_value), &(grid->old_value));

      /*------- computation of the implicit part of the theta scheme -------*/
      /*------- we use a LU decomposition method for linear operator  -------*/

      for (i = 0; i <= mxx + 1; i++) {
        /*--------------*/

        grid->temp_value[m][i][0][0][l] =
            grid->old_value[m][i][0][0][l] -
            pde->bound_cond_y * grid->old_value[m][i][1][0][l];

        for (j = 1; j <= myy; j++) {
          if ((i == 0) || (i == mxx + 1)) {
            grid->temp_value[m][i][j][0][l] =
                grid->old_value[m][i][j][0][l] -
                grid->L2[i][j][0] * grid->temp_value[m][i][j - 1][0][l];
          } else {
            grid->temp_value[m][i][j][0][l] =
                grid->old_value[m][i][j][0][l] +
                grid->EXY[i][j][0] * (grid->old_value[m][i + 1][j + 1][0][l] -
                                      grid->old_value[m][i - 1][j + 1][0][l] -
                                      grid->old_value[m][i + 1][j - 1][0][l] +
                                      grid->old_value[m][i - 1][j - 1][0][l]) -
                grid->L2[i][j][0] * grid->temp_value[m][i][j - 1][0][l];
          }
        }

        grid->temp_value[m][i][myy + 1][0][l] =
            grid->old_value[m][i][myy + 1][0][l] -
            pde->bound_cond_y * grid->old_value[m][i][myy][0][l] -
            grid->L2[i][myy + 1][0] * grid->temp_value[m][i][myy][0][l];

        /*--------------*/

        grid->new_value[m][i][myy + 1][0][l] =
            grid->temp_value[m][i][myy + 1][0][l] / grid->D2[i][myy + 1][0];

        for (j = myy; j >= 0; j--) {
          grid->new_value[m][i][j][0][l] =
              (grid->temp_value[m][i][j][0][l] -
               grid->CY[i][j][0] * grid->new_value[m][i][j + 1][0][l]) /
              grid->D2[i][j][0];
        }

        /*--------------*/
      }
    }
  }

  return err;
} /* END Err  solve_in_Y(...) */

/*------------------------------------------------------------------------------------------*/

Err solve_in_X(SrtPdeObject *pde) {
  Err err = NULL;
  SrtPdeGridObject *grid = pde->grid;
  int mxx = grid->num_x, myy = grid->num_y, nbcoeff = grid->nb_coeffs;

  int i = 0, j = 0, k = 0, m = 0, l = 0;

  for (m = 0; m < grid->nb_payoffs; m++) {
    /*---- step  2: solving the systeme in X ----*/
    for (l = 0; l < nbcoeff; l++) {
      /*------- computation of the explicit part of the theta scheme -------*/

      for (j = 1; j <= myy; j++) {
        for (i = 1; i <= mxx; i++) {
          grid->new_value[m][i][j][0][l] =
              grid->AY1[i - 1][j][0] * grid->old_value[m][i - 1][j][0][l] +
              grid->BY1[i][j][0] * grid->old_value[m][i][j][0][l] +
              grid->CY1[i + 1][j][0] * grid->old_value[m][i + 1][j][0][l];
        }
      }

      swap_pde_new_old(&(grid->new_value), &(grid->old_value));

      /*------- computation of the implicit part of the theta scheme -------*/
      /*------- we use a LU decomposition method for linear operator  -------*/

      for (j = 0; j <= myy + 1; j++) {
        /*--------------*/

        grid->temp_value[m][0][j][0][l] =
            grid->old_value[m][0][j][0][l] -
            pde->bound_cond_x * grid->old_value[m][1][j][0][l];

        for (i = 1; i <= mxx; i++) {
          if ((j == 0) || (j == myy + 1)) {
            grid->temp_value[m][i][j][0][l] =
                grid->old_value[m][i][j][0][l] -
                grid->L1[i][j][0] * grid->temp_value[m][i - 1][j][0][l];
          } else {
            grid->temp_value[m][i][j][0][l] =
                grid->old_value[m][i][j][0][l] +
                grid->EXY[i][j][0] * (grid->old_value[m][i + 1][j + 1][0][l] -
                                      grid->old_value[m][i - 1][j + 1][0][l] -
                                      grid->old_value[m][i + 1][j - 1][0][l] +
                                      grid->old_value[m][i - 1][j - 1][0][l]) -
                grid->L1[i][j][0] * grid->temp_value[m][i - 1][j][0][l];
          }
        }

        grid->temp_value[m][mxx + 1][j][0][l] =
            grid->old_value[m][mxx + 1][j][0][l] -
            pde->bound_cond_x * grid->old_value[m][mxx][j][0][l] -
            grid->L1[mxx + 1][j][0] * grid->temp_value[m][mxx][j][0][l];

        /*--------------*/

        grid->new_value[m][mxx + 1][j][0][l] =
            grid->temp_value[m][mxx + 1][j][0][l] / grid->D1[mxx + 1][j][0];

        for (i = mxx; i >= 0; i--) {
          grid->new_value[m][i][j][0][l] =
              (grid->temp_value[m][i][j][0][l] -
               grid->CX[i][j][0] * grid->new_value[m][i + 1][j][0][l]) /
              grid->D1[i][j][0];
        }

        /*--------------*/
      }
    }
  }

  return err;
} /* END Err  solve_in_X(...) */

/*-------------------------------------------------------------------*/
/*-------------------------------------------------------------------*/

void swap_pde_new_old(double ******U, double ******V) {
  double *****W;

  W = *V;
  *V = *U;
  *U = W;
}

/*-------------------------------------------------------------------*/
/*--- End of File ---*/
