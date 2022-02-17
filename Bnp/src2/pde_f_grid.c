/* ========================================================================
   FILENAME:     utpdegrid.c

  PURPOSE:     Allocation and Desallocation functions  for the work space
   ======================================================================== */

#include "pde_h_grid.h"

#define MEMSET(p, c, sz) memset((char *)p, c, sz)

/* ----------------------------------------------------------------------------------------------------------------------------------------
 */
/* allocation of workspace for pde solving */

Err pde_new_grid(SrtPdeObject *pde, long num_points[3]) {
  Err err = NULL;
  SrtPdeGridObject *grid = pde->grid;

  int taillealloc;
  int mxx = 0, nxx = 1, myy = 0, nyy = 1, mzz = 0, nzz = 1;
  int i, j, k, m;

  if ((pde->dimension == 0) || (pde->nb_payoffs == 0)) {
    return serror("ill-defined pde");
  }

  grid->nb_payoffs = pde->nb_payoffs;
  grid->dimension = pde->dimension;
  grid->num_x = num_points[0];
  grid->num_y = num_points[1];
  grid->num_z = num_points[2];
  grid->nb_coeffs = pde->nb_second_members;

  nxx = grid->num_x + 2;
  mxx = grid->num_x + 1;

  nyy = grid->num_y + 2;

  nzz = grid->num_z + 2;

  taillealloc = nxx * sizeof(double);

  grid->vec_x = (double *)malloc(taillealloc);
  MEMSET(grid->vec_x, 0, taillealloc);
  grid->inc_x = (double *)malloc(taillealloc);
  MEMSET(grid->inc_x, 0, taillealloc);

  taillealloc = nxx * sizeof(double **);

  grid->AX = (double ***)malloc(taillealloc);
  grid->BX = (double ***)malloc(taillealloc);
  grid->CX = (double ***)malloc(taillealloc);
  grid->AX1 = (double ***)malloc(taillealloc);
  grid->BX1 = (double ***)malloc(taillealloc);
  grid->CX1 = (double ***)malloc(taillealloc);

  grid->L1 = (double ***)malloc(taillealloc);
  grid->D1 = (double ***)malloc(taillealloc);

  if (grid->dimension >= 2) {
    taillealloc = nyy * sizeof(double);

    grid->vec_y = (double *)malloc(taillealloc);
    MEMSET(grid->vec_y, 0, taillealloc);
    grid->inc_y = (double *)malloc(taillealloc);
    MEMSET(grid->inc_y, 0, taillealloc);

    taillealloc = nxx * sizeof(double **);

    grid->AY = (double ***)malloc(taillealloc);
    grid->BY = (double ***)malloc(taillealloc);
    grid->CY = (double ***)malloc(taillealloc);
    grid->AY1 = (double ***)malloc(taillealloc);
    grid->BY1 = (double ***)malloc(taillealloc);
    grid->CY1 = (double ***)malloc(taillealloc);

    grid->EXY = (double ***)malloc(taillealloc);
    grid->L2 = (double ***)malloc(taillealloc);
    grid->D2 = (double ***)malloc(taillealloc);
  }

  if (grid->dimension == 3) {
    taillealloc = nzz * sizeof(double);

    grid->vec_z = (double *)malloc(taillealloc);
    MEMSET(grid->vec_z, 0, taillealloc);
    grid->inc_z = (double *)malloc(taillealloc);
    MEMSET(grid->inc_z, 0, taillealloc);

    taillealloc = nxx * sizeof(double **);

    grid->AZ = (double ***)malloc(taillealloc);
    grid->BZ = (double ***)malloc(taillealloc);
    grid->CZ = (double ***)malloc(taillealloc);
    grid->AZ1 = (double ***)malloc(taillealloc);
    grid->BZ1 = (double ***)malloc(taillealloc);
    grid->CZ1 = (double ***)malloc(taillealloc);

    grid->EXZ = (double ***)malloc(taillealloc);
    grid->EYZ = (double ***)malloc(taillealloc);

    grid->L3 = (double ***)malloc(taillealloc);
    grid->D3 = (double ***)malloc(taillealloc);
  }

  i = 0;
  while (i <= mxx) {
    nyy = (grid->dimension >= 2) ? grid->num_y + 2 : 1;

    taillealloc = nyy * sizeof(double *);

    grid->AX[i] = (double **)malloc(taillealloc);
    grid->BX[i] = (double **)malloc(taillealloc);
    grid->CX[i] = (double **)malloc(taillealloc);
    grid->AX1[i] = (double **)malloc(taillealloc);
    grid->BX1[i] = (double **)malloc(taillealloc);
    grid->CX1[i] = (double **)malloc(taillealloc);

    grid->L1[i] = (double **)malloc(taillealloc);
    grid->D1[i] = (double **)malloc(taillealloc);

    if (grid->dimension >= 2) {
      grid->AY[i] = (double **)malloc(taillealloc);
      grid->BY[i] = (double **)malloc(taillealloc);
      grid->CY[i] = (double **)malloc(taillealloc);
      grid->AY1[i] = (double **)malloc(taillealloc);
      grid->BY1[i] = (double **)malloc(taillealloc);
      grid->CY1[i] = (double **)malloc(taillealloc);
      grid->EXY[i] = (double **)malloc(taillealloc);

      grid->L2[i] = (double **)malloc(taillealloc);
      grid->D2[i] = (double **)malloc(taillealloc);
    }

    if (grid->dimension == 3) {
      grid->AZ[i] = (double **)malloc(taillealloc);
      grid->BZ[i] = (double **)malloc(taillealloc);
      grid->CZ[i] = (double **)malloc(taillealloc);
      grid->AZ1[i] = (double **)malloc(taillealloc);
      grid->BZ1[i] = (double **)malloc(taillealloc);
      grid->CZ1[i] = (double **)malloc(taillealloc);
      grid->EXZ[i] = (double **)malloc(taillealloc);
      grid->EYZ[i] = (double **)malloc(taillealloc);

      grid->L3[i] = (double **)malloc(taillealloc);
      grid->D3[i] = (double **)malloc(taillealloc);
    }

    myy = (grid->dimension >= 2) ? grid->num_y + 1 : 0;
    j = 0;
    while (j <= myy) {
      nzz = (grid->dimension >= 3) ? grid->num_z + 2 : 1;

      taillealloc = nzz * sizeof(double);

      grid->AX[i][j] = (double *)malloc(taillealloc);
      MEMSET(grid->AX[i][j], 0, taillealloc);
      grid->BX[i][j] = (double *)malloc(taillealloc);
      MEMSET(grid->BX[i][j], 0, taillealloc);
      grid->CX[i][j] = (double *)malloc(taillealloc);
      MEMSET(grid->CX[i][j], 0, taillealloc);
      grid->AX1[i][j] = (double *)malloc(taillealloc);
      MEMSET(grid->AX1[i][j], 0, taillealloc);
      grid->BX1[i][j] = (double *)malloc(taillealloc);
      MEMSET(grid->BX1[i][j], 0, taillealloc);
      grid->CX1[i][j] = (double *)malloc(taillealloc);
      MEMSET(grid->CX1[i][j], 0, taillealloc);

      grid->L1[i][j] = (double *)malloc(taillealloc);
      MEMSET(grid->L1[i][j], 0, taillealloc);
      grid->D1[i][j] = (double *)malloc(taillealloc);
      MEMSET(grid->D1[i][j], 0, taillealloc);

      if (grid->dimension >= 2) {
        taillealloc = nzz * sizeof(double);

        grid->AY[i][j] = (double *)malloc(taillealloc);
        MEMSET(grid->AY[i][j], 0, taillealloc);
        grid->BY[i][j] = (double *)malloc(taillealloc);
        MEMSET(grid->BY[i][j], 0, taillealloc);
        grid->CY[i][j] = (double *)malloc(taillealloc);
        MEMSET(grid->CY[i][j], 0, taillealloc);
        grid->AY1[i][j] = (double *)malloc(taillealloc);
        MEMSET(grid->AY1[i][j], 0, taillealloc);
        grid->BY1[i][j] = (double *)malloc(taillealloc);
        MEMSET(grid->BY1[i][j], 0, taillealloc);
        grid->CY1[i][j] = (double *)malloc(taillealloc);
        MEMSET(grid->CY1[i][j], 0, taillealloc);
        grid->EXY[i][j] = (double *)malloc(taillealloc);
        MEMSET(grid->EXY[i][j], 0, taillealloc);

        grid->L2[i][j] = (double *)malloc(taillealloc);
        MEMSET(grid->L2[i][j], 0, taillealloc);
        grid->D2[i][j] = (double *)malloc(taillealloc);
        MEMSET(grid->D2[i][j], 0, taillealloc);
      }

      if (grid->dimension == 3) {
        taillealloc = nzz * sizeof(double);

        grid->AZ[i][j] = (double *)malloc(taillealloc);
        MEMSET(grid->AZ[i][j], 0, taillealloc);
        grid->BZ[i][j] = (double *)malloc(taillealloc);
        MEMSET(grid->BZ[i][j], 0, taillealloc);
        grid->CZ[i][j] = (double *)malloc(taillealloc);
        MEMSET(grid->CZ[i][j], 0, taillealloc);
        grid->AZ1[i][j] = (double *)malloc(taillealloc);
        MEMSET(grid->AZ1[i][j], 0, taillealloc);
        grid->BZ1[i][j] = (double *)malloc(taillealloc);
        MEMSET(grid->BZ1[i][j], 0, taillealloc);
        grid->CZ1[i][j] = (double *)malloc(taillealloc);
        MEMSET(grid->CZ1[i][j], 0, taillealloc);
        grid->EXZ[i][j] = (double *)malloc(taillealloc);
        MEMSET(grid->EXZ[i][j], 0, taillealloc);
        grid->EYZ[i][j] = (double *)malloc(taillealloc);
        MEMSET(grid->EYZ[i][j], 0, taillealloc);

        grid->L3[i][j] = (double *)malloc(taillealloc);
        MEMSET(grid->L3[i][j], 0, taillealloc);
        grid->D3[i][j] = (double *)malloc(taillealloc);
        MEMSET(grid->D3[i][j], 0, taillealloc);
      }

      j++;
    }

    i++;
  }

  taillealloc = grid->nb_payoffs * sizeof(double ****);

  grid->new_value = (double *****)malloc(taillealloc);
  grid->old_value = (double *****)malloc(taillealloc);

  grid->temp_value = (double *****)malloc(taillealloc);

  m = 0;
  while (m < grid->nb_payoffs) {
    nxx = grid->num_x + 2;
    taillealloc = nxx * sizeof(double ***);

    grid->new_value[m] = (double ****)malloc(taillealloc);
    grid->old_value[m] = (double ****)malloc(taillealloc);

    grid->temp_value[m] = (double ****)malloc(taillealloc);

    mxx = grid->num_x + 1;
    i = 0;
    while (i <= mxx) {
      nyy = (grid->dimension >= 2) ? grid->num_y + 2 : 1;
      taillealloc = nyy * sizeof(double **);

      grid->new_value[m][i] = (double ***)malloc(taillealloc);
      ;
      grid->old_value[m][i] = (double ***)malloc(taillealloc);

      grid->temp_value[m][i] = (double ***)malloc(taillealloc);

      myy = (grid->dimension >= 2) ? grid->num_y + 1 : 0;
      j = 0;
      while (j <= myy) {
        nzz = (grid->dimension >= 3) ? grid->num_z + 2 : 1;
        taillealloc = nzz * sizeof(double *);

        grid->new_value[m][i][j] = (double **)malloc(taillealloc);
        grid->old_value[m][i][j] = (double **)malloc(taillealloc);
        ;

        grid->temp_value[m][i][j] = (double **)malloc(taillealloc);

        mzz = (grid->dimension >= 3) ? grid->num_z + 1 : 0;
        k = 0;
        while (k <= mzz) {
          taillealloc = grid->nb_coeffs * sizeof(double);

          grid->new_value[m][i][j][k] = (double *)malloc(taillealloc);
          MEMSET(grid->new_value[m][i][j][k], 0, taillealloc);
          grid->old_value[m][i][j][k] = (double *)malloc(taillealloc);
          MEMSET(grid->old_value[m][i][j][k], 0, taillealloc);

          grid->temp_value[m][i][j][k] = (double *)malloc(taillealloc);
          MEMSET(grid->temp_value[m][i][j][k], 0, taillealloc);

          k++;
        }

        j++;
      }

      i++;
    }

    m++;
  }

  return err;

} /* END Err pde_new_grid(...) */

/* ----------------------------------------------------------------------------------------------------------------------------------------
 */
/* allocation of workspace for pde solving */

Err pde_delete_grid(SrtPdeGridObject *grid) {
  Err err = NULL;
  int mxx = 0, myy = 0, mzz = 0;
  int i, j, m, k;

  mxx = grid->num_x + 1;
  i = 0;
  while (i <= mxx) {
    myy = (grid->dimension >= 2) ? grid->num_y + 1 : 0;
    j = 0;
    while (j <= myy) {

      free(grid->L1[i][j]);
      free(grid->D1[i][j]);

      free(grid->AX[i][j]);
      free(grid->BX[i][j]);
      free(grid->CX[i][j]);
      free(grid->AX1[i][j]);
      free(grid->BX1[i][j]);
      free(grid->CX1[i][j]);

      if (grid->dimension >= 2) {
        free(grid->AY[i][j]);
        free(grid->BY[i][j]);
        free(grid->CY[i][j]);
        free(grid->AY1[i][j]);
        free(grid->BY1[i][j]);
        free(grid->CY1[i][j]);
        free(grid->EXY[i][j]);

        free(grid->L2[i][j]);
        free(grid->D2[i][j]);
      }

      if (grid->dimension == 3) {
        free(grid->AZ[i][j]);
        free(grid->BZ[i][j]);
        free(grid->CZ[i][j]);
        free(grid->AZ1[i][j]);
        free(grid->BZ1[i][j]);
        free(grid->CZ1[i][j]);
        free(grid->EXZ[i][j]);
        free(grid->EYZ[i][j]);

        free(grid->L3[i][j]);
        free(grid->D3[i][j]);
      }

      j++;
    }

    free(grid->L1[i]);
    free(grid->D1[i]);

    free(grid->AX[i]);
    free(grid->BX[i]);
    free(grid->CX[i]);
    free(grid->AX1[i]);
    free(grid->BX1[i]);
    free(grid->CX1[i]);

    if (grid->dimension >= 2) {
      free(grid->AY[i]);
      free(grid->BY[i]);
      free(grid->CY[i]);
      free(grid->AY1[i]);
      free(grid->BY1[i]);
      free(grid->CY1[i]);
      free(grid->EXY[i]);

      free(grid->L2[i]);
      free(grid->D2[i]);
    }

    if (grid->dimension == 3) {
      free(grid->AZ[i]);
      free(grid->BZ[i]);
      free(grid->CZ[i]);
      free(grid->AZ1[i]);
      free(grid->BZ1[i]);
      free(grid->CZ1[i]);
      free(grid->EXZ[i]);
      free(grid->EYZ[i]);

      free(grid->L3[i]);
      free(grid->D3[i]);
    }

    i++;
  }

  free(grid->L1);
  free(grid->D1);

  free(grid->vec_x);
  free(grid->inc_x);

  free(grid->AX);
  free(grid->BX);
  free(grid->CX);
  free(grid->AX1);
  free(grid->BX1);
  free(grid->CX1);

  if (grid->dimension >= 2) {
    free(grid->vec_y);
    free(grid->inc_y);

    free(grid->AY);
    free(grid->BY);
    free(grid->CY);
    free(grid->AY1);
    free(grid->BY1);
    free(grid->CY1);
    free(grid->EXY);

    free(grid->L2);
    free(grid->D2);
  }

  if (grid->dimension == 3) {
    free(grid->vec_z);
    free(grid->inc_z);

    free(grid->AZ);
    free(grid->BZ);
    free(grid->CZ);
    free(grid->AZ1);
    free(grid->BZ1);
    free(grid->CZ1);
    free(grid->EXZ);
    free(grid->EYZ);

    free(grid->L3);
    free(grid->D3);
  }

  m = 0;
  while (m < grid->nb_payoffs) {

    mxx = grid->num_x + 1;
    i = 0;
    while (i <= mxx) {

      myy = (grid->dimension >= 2) ? grid->num_y + 1 : 0;
      j = 0;
      while (j <= myy) {

        mzz = (grid->dimension >= 3) ? grid->num_z + 1 : 0;
        k = 0;
        while (k <= mzz) {
          free(grid->new_value[m][i][j][k]);
          free(grid->old_value[m][i][j][k]);

          free(grid->temp_value[m][i][j][k]);

          k++;
        }

        free(grid->new_value[m][i][j]);
        free(grid->old_value[m][i][j]);

        free(grid->temp_value[m][i][j]);

        j++;
      }

      free(grid->new_value[m][i]);
      free(grid->old_value[m][i]);

      free(grid->temp_value[m][i]);

      i++;
    }

    free(grid->new_value[m]);
    free(grid->old_value[m]);

    free(grid->temp_value[m]);

    m++;
  }

  free(grid->new_value);
  free(grid->old_value);

  free(grid->temp_value);

  return err;
} /* END Err pde_delete_grid(...) */

/*-------------------------------------------------------------------*/
/*--- End of File ---*/
