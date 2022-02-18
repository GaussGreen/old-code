/******************************************************************************\
*                Copyright (c) 1995 PARIBAS Capital Markets Group              *
********************************************************************************
*
*   SYSTEM          :   GRF
*
*   MODULE NAME     :   grf_f_types
*
*   PURPOSE         :   Memory handling of all GRFN data types
*
*   AUTHOR          :
*
*   DATE            :
*
*   VERSION         :
*
*   DESCRIPTION     :
*
*   FUNCTIONS USED  :
*
*   PARAMETERS      :
*
*   RETURNS         :
*
*
********************************************************************************
*                           Amendment History                                  *
********************************************************************************
*
*   AMENDED BY      :
*
*   DATE            :
*
*   VERSION         :
*
*   REASON          :
*
*   REQUEST/BUG NO  :
*
*   DESCRIPTION
*
\******************************************************************************/

/******************************************************************************\
*                               Import Include Files                           *
\******************************************************************************/

#include "grf_h_all.h"

/* Checks there is at least one cell in the tableau, and that the dates are increasing strictly */
static Err grf_validate_args(Date* dates, long m, long n)

{
    int i;

    /* Check there is at least one cell */
    if (m < 1 || n < 1)
        return GRERR_BAD_EVDIM;

    /* Checks the dates are increasing */
    for (i = 0; i < m - 1; i++)
    {
        if (dates[i] >= dates[i + 1])
            return serror("GRFN: Event dates must be strictly increasing.");
    }
    return NULL;
}

/* ------------------------------------------------------------------------------------- */

/* Check from the bottom of the tableau that there are no empty rows; if YES: remove them */
static Err grf_check_and_rem_null_rows(long* nrows, long ncols, GrfnCell** sprdsht)
{
    long i, j;
    int  null_row;
    Err  err = NULL;

    /* Check for NULL rows and remove them */
    for (i = *nrows - 1; i > 0; i--)
    {
        null_row = 1;
        for (j = 0; j < ncols; j++)
        {
            null_row *= (sprdsht[i][j].type == GRFNBCELL);
        }
        if (null_row)
            (*nrows)--;
        else
            break;
    }

    return NULL;
}

/* ------------------------------------------------------------------------------------- */

/* -------------------------------------------------------------------------------------
 FUNCNAME        :GrfnCellmatrix
  AUTHOR          :allocate an array of grfn cells.
  if strlen is greater than zero, we allocate a character
  array of len strlen for the sval field of each GrfnCell.
  otherwise the sval field is set to NULL.  Everything else is set to
  zero.  If the sval field is allocated, the SRT_Boolean field str_alloced
  will be set to SRT_YES; otherwise it will be set to SRT_NO.
  Returns NULL if it fails.
  DESCRIPTION     :
  MODIFIES        :
  CALL            :

  ------------------------------------------------------------------------------------- */

GrfnCell** GrfnCellmatrix(long nrow, long ncol, long slen)
{
    int        i, j;
    GrfnCell **a, *b;

    a = (GrfnCell**)srt_calloc(nrow, sizeof(GrfnCell*));
    /* check for failure */
    if (!a)
        return NULL;

    b = (GrfnCell*)srt_calloc(nrow * ncol, sizeof(GrfnCell));

    /* check for failure */
    if (!b)
    {
        srt_free(a);
        return NULL;
    }

    for (i = 0; i < nrow; i++)
    {
        a[i] = b + i * ncol;
    }

    /* allocate strings */
    if (slen > 0)
    {
        for (i = 0; i < nrow; i++)
        {
            for (j = 0; j < ncol; j++)
            {
                a[i][j].sval = (String)srt_calloc(1, (slen + 1) * sizeof(char));
                /* check for failure */
                if (!a[i][j].sval)
                {
                    srt_free(b);
                    grfn_free_GrfnCellmatrix(a, nrow, ncol);
                    return NULL;
                }
                a[i][j].str_alloced = SRT_YES;
            }
        }
    }
    return a;
}

/* -------------------------------------------------------------------------------------------- */

/*  -------------------------------------------------------------------------
  FUNCNAME        :grfn_free_GrfnCellmatrix
  AUTHOR          :E.Auld
  DESCRIPTION     :free a matrix of grfn cells.  Note that each GrfnCell has a
boolean field called str_alloced.  If the sval field of a GrfnCell is not NULL,
then this function will look at the str_alloced field to determined whether or
no the sval field was allocated (in which case it will be freed).
  MODIFIES        :
  CALL            :

  ---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------
  AMENDMENTS      :
  Reference       :
  Author          :
  Date            :
  Description     :
-----------------------------------------------------------------------------*/

void grfn_free_GrfnCellmatrix(GrfnCell** m, long nrow, long ncol)
{
    int i, j;

    for (i = 0; i < nrow; i++)
    {
        for (j = 0; j < ncol; j++)
        {
            if (m[i][j].sval && m[i][j].str_alloced == SRT_YES)
                srt_free(m[i][j].sval);
        }
    }
    srt_free(m[0]);
    srt_free(m);
}

/* ---------------------------------------------------------------------- */

void grfn_free_GrfnEvent(GrfnEvent* event)
{
    int i;

    if (event == NULL)
        return;

    /* Free the discount factor (and related) arrays for all underlyings */
    for (i = 0; i < MAXUNDERLYING; i++)
    {
        if (event[0].df[i] != NULL)
            srt_free(event[0].df[i]);
        if (event[0].dft[i] != NULL)
            srt_free(event[0].dft[i]);
        if (event[0].dfd[i] != NULL)
            srt_free(event[0].dfd[i]);
        if (event[0].yp[i] != NULL)
            srt_free(event[0].yp[i]);
    }

    if (event[0].icom != NULL)
        comll_free_list(event[0].icom);

    srt_free(event);
}

/* ----------------------------------------------------------------------- */

/* -----------------------------------------------------------------------

  FUNCNAME        :grfn_free_inGrfnDeal
  AUTHOR          :E.Auld
  DESCRIPTION     :FREE everything inside the GrfnDeal.  Do not free
  the GrfnDeal itself.  If the sptr field of the GrfnDeal is not NULL
  and the stp->e fields within are not empty, they are assumed to be
  pointers to GrfnEvents and freed as such.
  MODIFIES        :
  CALL            :

   ----------------------------------------------------------------------- */

void grfn_free_inGrfnDeal(GrfnDeal* gd)
{
    if (gd->event_dates)
        srt_free(gd->event_dates);

    if (gd->auxlen)
        free_lngvector(gd->auxlen, 0, gd->auxwidth - 1);
    if (gd->aux)
    {
        srt_free(gd->aux[0]);
        srt_free(gd->aux);
    }

    if (gd->rowstatus)
        srt_free(gd->rowstatus);

    if (gd->cells)
        free_matrix(gd->cells, 0, gd->sslength - 1, 0, gd->sswidth - 1);

    if (gd->gcells)
    {
        grfn_free_GrfnCellmatrix(gd->gcells, gd->sslength, gd->sswidth);
        gd->gcells = NULL;
    }

    if (gd->pay_dates)
        srt_free(gd->pay_dates);
    if (gd->event_pay_dates)
        srt_free(gd->event_pay_dates);
    if (gd->pay_amounts)
        srt_free(gd->pay_amounts);

   /* if (gd->outfileptr)
        fclose(gd->outfileptr);*/

    memset(gd, 0, sizeof(GrfnDeal));
}

/* -------------------------------------------------------------------------- */

/* --------------------------------------------------------------------------
  FUNCNAME        :grfn_check_tableau_and_dates
  DESCRIPTION:      checks the tableau dates, null rows and at least one cell
                    if NULL rows are detected: removes them
  --------------------------------------------------------------------------- */

Err grfn_check_tableau_and_dates(long* nrows, long ncols, GrfnCell** sprdsht, Date* eventdates)
{
    Err err = NULL;

    /* Check for and removes NULL rows (and corrects the number of rows ) */
    err = grf_check_and_rem_null_rows(nrows, ncols, sprdsht);
    if (err)
        return err;

    /*  Check that the tableau is not empty, check the dates are increasing */
    err = grf_validate_args(eventdates, *nrows, ncols);
    if (err)
        return err;

    return NULL;
}

/* -------------------------------------------------------------------------------------------- */

/* -----------------------------------------------------------------
  FUNCNAME        :grfn_copy_to_GrfnDeal
  AUTHOR          :E.Auld
  DESCRIPTION     :make copies of all the inputs and place these copies
  inside the GrfnDeal.  Copies are made in order to make
  the GrfnDeal autonomous.  If the bottom of the GRFN Tableau is
  padded with blank rows these will be discarded.
  FIX: this does not as yet do anything with the grng field.
  MODIFIES        :
  CALL            :

  ---------------------------------------------------------------------*/

Err grfn_copy_to_GrfnDeal(
    GrfnDeal*  gd,
    long       nrows,
    long       ncols,
    GrfnCell** sprdsht,
    Date*      eventdates,
    long       numgrng,
    GrfnRng*   grng,
    long       auxwidth,
    long*      auxlen,
    double**   aux)
{
    int  i, j;
    long auxtoalloc = 0;

    /* INITIALISATION: set  the GrfnDeal to NULL  */
    memset(gd, 0, sizeof(GrfnDeal));

    /* Starts populating the GrfnDeal with the tableau dimensions, dates,... */
    gd->sswidth  = ncols;
    gd->sslength = nrows;

    gd->event_dates = (Date*)srt_malloc(nrows * sizeof(Date));
    memcpy(gd->event_dates, eventdates, nrows * sizeof(Date));

    gd->rowstatus = (GrfnCellStatus*)srt_calloc(nrows, sizeof(GrfnCellStatus));

    /* Copies the tableau (as strings) into the GrfnDeal */
    gd->gcells = GrfnCellmatrix(nrows, ncols, 0);
    memcpy(&gd->gcells[0][0], &sprdsht[0][0], nrows * ncols * sizeof(GrfnCell));
    for (i = 0; i < nrows; i++)
    {
        for (j = 0; j < ncols; j++)
        {
            if (gd->gcells[i][j].sval != NULL && gd->gcells[i][j].str_alloced == SRT_YES)
            {
                gd->gcells[i][j].sval = new_string(gd->gcells[i][j].sval);
            }
        }
    }

    /* Allocates space for the GrfnCell values in the GrfnDeal (values of the cells once evaluated)
     */
    gd->cells = dmatrix(0, nrows - 1, 0, ncols - 1);
    memset(&gd->cells[0][0], 0, nrows * ncols * sizeof(double));

    /* Allocate memory for and set up the auxiliary ranges into the GrfnDeal */
    if (aux && auxlen && auxwidth > 0)
    {
        gd->auxwidth = auxwidth;
        gd->auxlen   = lngvector(0, auxwidth - 1);
        memcpy(gd->auxlen, auxlen, auxwidth * sizeof(long));
        for (i = 0; i < auxwidth; i++)
        {
            auxtoalloc += auxlen[i];
        }
        gd->aux    = (double**)srt_calloc(auxwidth, sizeof(double*));
        gd->aux[0] = (double*)srt_calloc(auxtoalloc, sizeof(double));
        for (i = 1; i < auxwidth; i++)
        {
            gd->aux[i] = gd->aux[i - 1] + auxlen[i - 1];
        }
        for (i = 0; i < auxwidth; i++)
        {
            memcpy(gd->aux[i], aux[i], auxlen[i] * sizeof(double));
        }
    }

    return NULL;
}

/* -------------------------------------------------------------------------------------------- */

/* Extend tableau when all events are in the past */

Err grfn_extend_tableau(GrfnCell*** tab, long* nrow, long ncol)
{
    /*..Local variables..*/
    long       i, j, m, n;
    int        s_len;
    GrfnCell** tableau; /* Extended tableau */

    /*..Initialise local variables..*/

    s_len = 0;
    m     = *nrow;
    n     = ncol;

    /* Allocate memory for a new EXTENDED tableau */
    tableau = GrfnCellmatrix(*nrow + 1, ncol, 0);

    /* Copy old tableau to new one, row by row, cell by cell */
    for (i = 0; i < m; i++)
    {
        for (j = 0; j < n; j++)
        {
            tableau[i][j].type = (*tab)[i][j].type;
            GrfnCopyStatus(tableau[i][j].status, (*tab)[i][j].status);
            tableau[i][j].dval = (*tab)[i][j].dval;
            tableau[i][j].sval = NULL;

            s_len = 0;
            if ((*tab)[i][j].type != GRFNBCELL)
            {
                /* If the old tableau has allocated for strings,
                       allocate space in new tableau and copy string
                        */
                if ((*tab)[i][j].sval != NULL && (*tab)[i][j].str_alloced == SRT_YES)
                {
                    s_len = strlen((*tab)[i][j].sval);
                }

                if (s_len > 0)
                {
                    tableau[i][j].sval = (String)srt_calloc(s_len + 1, sizeof(char));

                    /* Check allocation status */
                    if (!tableau[i][j].sval)
                    {
                        grfn_free_GrfnCellmatrix(tableau, m, n);
                        tableau = NULL;
                        return serror("Memory Allocation failure in extend_tableau");
                    }
                    strcpy(tableau[i][j].sval, (*tab)[i][j].sval);
                    tableau[i][j].str_alloced = SRT_YES;

                } /*end if */
            }     /* end if */
        }         /*end for */
    }             /* end for */

    /* Now the copy is done: frees the old tableau and points to the new one */
    grfn_free_GrfnCellmatrix(*tab, *nrow, ncol);
    (*tab) = tableau;

    /* Add the extra cell in the last row of the tableau: 0.0 */
    tableau[m][n - 1].type = GRFNDCELL;
    tableau[m][n - 1].dval = 0.0;

    /* Corrects the number of rows */
    *nrow = m + 1;

    return NULL;

} /* END Err grfn_extend_tableau(...) */

/*---------------------------------------------------------------------------*/
