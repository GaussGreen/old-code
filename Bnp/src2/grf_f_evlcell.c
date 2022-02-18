#include "grf_h_all.h"

#define CELL(col, row) cell[row][col]

static int qcmp(const void* vp, const void* vq)
{
    const double* p    = vp;
    const double* q    = vq;
    double        diff = *p - *q;

    return ((diff >= 0.0) ? ((diff > 0.0) ? -1 : 0) : +1);
}

/* ------------------------------------------------------- */

/* Function called when evaluating the COMLL */

void grf_f_evlcllfct(
    double**  cell,
    int       c,
    int       r1,
    int       r2,
    int       ind, /* used for ROWSORT, COLSORT  */
    COMLLType fun_type,
    double*   answer)
{
    double  x;
    int     i;
    double* tmp_row;
    double* tmp_col;

    switch (fun_type)
    {
        /* Maximum column entry */
    case COMLL_C_COLMAX:

        x = CELL(c, r1);

        for (i = r1 + 1; i <= r2; i++)
            x = DMAX(x, CELL(c, i));

        break;

    /* Minimum column entry */
    case COMLL_C_COLMIN:

        x = CELL(c, r1);

        for (i = r1 + 1; i <= r2; i++)
            x = DMIN(x, CELL(c, i));

        break;

        /* Sum of column entries */
    case COMLL_C_COLSUM:

        x = CELL(c, r1);

        for (i = r1 + 1; i <= r2; i++)
            x += CELL(c, i);

        break;

        /* Average of column entries */
    case COMLL_C_COLAVG:

        x = CELL(c, r1);

        for (i = r1 + 1; i <= r2; i++)
            x += CELL(c, i);

        x /= (r2 - r1 + 1.0);
        break;

        /* Sort the values of a column and return the requested ind th largest value */
    case COMLL_C_COLSORT:

        /* Make a temporary copy of the current cell values of column # c  */
        tmp_col = dvector(0, r2 - r1);
        for (i = 0; i <= (r2 - r1); i++)
            tmp_col[i] = CELL(c, r1 + i);

        /* Sort the values of this column (for rows from c1 to c2 ) */
        qsort(tmp_col, r2 - r1 + 1, sizeof(double), qcmp);

        /* Gets the ith index in the tableau */
        x = tmp_col[(int)ind];

        free_dvector(tmp_col, 0, r2 - r1);

        break;

        /* Maximum column entry */
    case COMLL_C_ROWMAX:

        x = CELL(r1, c);

        for (i = r1 + 1; i <= r2; i++)
            x = DMAX(x, CELL(i, c));

        break;

        /* and its index */
    case COMLL_C_ROWMAXIDX:                         /* added by C. Godart for Regis Benichou */
        tmp_row  = (double*)malloc(sizeof(double)); /* it's there, so use it */
        *tmp_row = CELL(r1, c);
        x        = c;

        for (i = r1 + 1; i <= r2; i++)
        {
            if (*tmp_row < CELL(i, c))
            {
                x        = i;
                *tmp_row = CELL(i, c);
            }
        }
        free(tmp_row);
        break;

    /* Minimum column entry */
    case COMLL_C_ROWMIN:

        x = CELL(r1, c);

        for (i = r1 + 1; i <= r2; i++)
            x = DMIN(x, CELL(i, c));

        break;
        /* and its index */
    case COMLL_C_ROWMINIDX: /* added by C. Godart for Regis Benichou */
        tmp_row  = (double*)malloc(sizeof(double));
        *tmp_row = CELL(r1, c);
        x        = r1;

        for (i = r1 + 1; i <= r2; i++)
        {
            if (*tmp_row > CELL(i, c))
            {
                x        = (double)i;
                *tmp_row = CELL(i, c);
            }
        }
        free(tmp_row);
        break;

        /* Sum of rows entries */
    case COMLL_C_ROWSUM:

        x = CELL(r1, c);

        for (i = r1 + 1; i <= r2; i++)
            x += CELL(i, c);

        break;

        /* Average of column entries */
    case COMLL_C_ROWAVG:

        x = CELL(r1, c);

        for (i = r1 + 1; i <= r2; i++)
            x += CELL(i, c);

        x /= (r2 - r1 + 1.0);
        break;

        /* Sort the values of a row and return the requested ind th largest value*/
    case COMLL_C_ROWSORT:

        /* Here, row and column inputs are really the column and row inputs.
           So c is really r, and (r1,r2) are really (c1,c2) */

        /* Make a temporary copy of the current cell values of row # c  */
        tmp_row = dvector(0, r2 - r1);
        for (i = 0; i <= (r2 - r1); i++)
            tmp_row[i] = CELL(r1 + i, c);

        /* Sort the values of this row (for cells from r1 to r2 ) */
        qsort(tmp_row, r2 - r1 + 1, sizeof(double), qcmp);

        /* Gets the ith index in the tableau */
        x = tmp_row[ind];

        free_dvector(tmp_row, 0, r2 - r1);

        break;

        /* Default case: do nothing */
    default:

        x = 0.0;

        break;
    }

    /* The result to return */
    *answer = x;
}

/* --------------------------------------------------------------------------------------- */
#define PVREF(c1) pv_copy[c1] /* Gets a previously computed PV */

/*  Approximate y = f(x) by linear interpolation on a row ( only for COMLL_C_ROWINTERP) */
void grf_f_evlcllfct_interp(
    double*   x_values,
    double**  cell,
    int       r,
    int       c1,
    int       c2,
    double    x,
    COMLLType type,
    double*   answer)
{
    int     i;
    double* tmp_row;
    double  y;
    int     n;

    switch (type)
    {
    case COMLL_C_ROWINTERP:

        /* Make a temporary copy of the current cell values of row # c  */
        tmp_row = dvector(0, c2 - c1);
        for (i = 0; i <= (c2 - c1); i++)
            tmp_row[i] = CELL(c1 + i, r);

        /* Sort the values of this row (for cells from r1 to r2 ) */
        qsort(tmp_row, c2 - c1 + 1, sizeof(double), qcmp);

        /* Interpolate y=f(x) over the values of this row using the x_values input */
        interp(x_values, tmp_row, c2 - c1 + 1, x, 0, &y);

        /* Free the temporary array */
        free_dvector(tmp_row, 0, c2 - c1);
        break;

    case COMLL_X_PVINTERP:
        /* Make a temporary copy of the current cell values of row # c  */
        n       = (int)(*answer);
        tmp_row = dvector(0, n - 1);
        for (i = 0; i < n; i++)
            tmp_row[i] = CELL(c1 + i * c2, r);
        //				tmp_row[i] = PVREF(c1+i*c2);

        /* Sort the values of this row (for cells from r1 to r2 ) */
        /*	qsort(tmp_row,c2-c1+1,sizeof(double),qcmp); */

        /* Interpolate y=f(x) over the values of this row using the x_values input */
        interp(x_values, tmp_row, n, x, 0, &y);

        /* Free the temporary array */
        free_dvector(tmp_row, 0, n - 1);
        break;

    default:

        y = 0.0;

        break;
    }

    /* The value to return */
    *answer = y;
}

/* ------------------------------------------------------------------------------------- */

#undef CELL