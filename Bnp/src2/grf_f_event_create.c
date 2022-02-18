/******************************************************************************\
*                Copyright (c) 1995 PARIBAS Capital Markets Group              *
********************************************************************************
*
*   SYSTEM          :   GRF
*
*   MODULE NAME     :   GRF_F_GET_EVENT
*
*   PURPOSE         :
*
*   AUTHOR          :
*
*   DATE            :
*
*   VERSION         :
*
*   DESCRIPTION     :   Implementation of the GRFN Language
*
*   FUNCTIONS USED  :   YYPARSE
*
*   PARAMETERS      :
*
*   RETURNS         :
*
*
********************************************************************************



/******************************************************************************\
*                               Import Include Files                           *
\******************************************************************************/

#include "grf_h_all.h"
#include "math.h"

/******************************************************************************\
*                               Global Variables                               *
\******************************************************************************/

/* Global variables used in evaluation of strings in grfn */

EXTERND String grfn_input_string; /*  Source of input to lex_yy() */

EXTERND SRT_Boolean grfn_lang_cur_amcell; /*  Is set to SRT_YES if parsed cell
                                          contains AM.
                                          Due to YACCs limited look-ahead a
                                          pattern of the form AM(...)
                                          will NOT work; In this case (...)
                                          will be evaluated prior to AM
                                          being taken into account. */

/******************************************************************************\
*                           Macro Functions                                    *
\******************************************************************************/

#define set_lattice_flg                                                                           \
    {                                                                                             \
        lattice_flg =                                                                             \
            (GrfnIsStatus(GrfnCGCell(gd).status, GRFNCSPVREF) ? COMLL_T_SET : COMLL_T_INCREMENT); \
    }

/******************************************************************************\
*                           Private Function Definitions                       *
\******************************************************************************/

/******************************************************************************\
*
*   FUNCTION        :   join_adjacent_cells
*
*   PARAMETERS      :   COMLL_PTR   c1
*                       COMLL_PTR   c2
*                       COMLLType   type
*                       int         i, j
*                       double      **cells
*
*   RETURNS         :   COMLL_PTR   c2
*
*   CALLED BY       :   create_event
*
*   DESCRIPTION     :   Joins together the COMLLs for two adjacent cells in
*                       the same row of a GRFN Tableau.

                        c1 = Command linked list for cell[0:j-1,i]
                        c2 =           ...           cell[j,i]

                        Note: [j,i] means column j, row i.

                        The two are joined together as follows:

                        c1 --> POP --> c2 --> SET/INCREMENT --> ASSIGN.

                        COMLL_POP  => Remove the value of c1 from the stack
                                      during evaluation.

                        COMLL_T_...=> Return the result of the evaluation to
                                      the calling tree (if we have a tree).

                    NB: The obscure rule of the GRFN Language says that if
                        the PV[...] function occurs in a cell within ANY of
                                                of the tableau (with the same column reference or
NOT) then the value being discounted through the tree is REPLACED by the value at that cell after
evaluation (like American style option to represent an exercise, or a choice depending on an otpion
value) otherwise it is INCREMENTED by that amount (like for the sum of deterministic cash flows, or
the final discounted sum of all last column cells) The "treeset" will be the value of the "type"
input to this function; it will not be inserted if it is not one of COMLL_T_SET and
COMLL_T_INCREMENT.

                        COMLL_ASSIGN => Assign saved value of the expression in
                                        cell[j,i] i.e. c2.

    SPECIAL CASES   :   If c2 is an atom of type real, skip COMLL_ASSIGN and
                        place its value directly into **cells.

                        If c2 is NULL or c1 and c2 are both NULL do nothing.

                        If c1 is NULL, then only do:
                        c2-->COMLL_T_SET/INCREMENT-->COMLL_ASSIGN.




\******************************************************************************/

static COMLL_PTR join_adjacent_cells(
    COMLL_PTR c1, COMLL_PTR c2, COMLLType type, int i, int j, double** cells)

{
    /*..Local variables..*/
    COMLL_PTR temp = NULL, temp1 = NULL;

    /*..Main Body..*/
    /* Both cells are empty */
    if (c1 == NULL && c2 == NULL)
        return NULL;

    /* One cell is empty    */
    if (c2 == NULL)
        return (c1);

    /* Beginning of c2      */
    temp = comll_gotobot(c2);

    /*  If c2 is NOT a "single" node of type real */
    if (!comll_atom(c2, COMLL_REAL))
    {
        temp             = comll_insert_after(temp);
        temp->type       = COMLL_ASSIGN;
        temp->gdcells    = &cells[i][j];
        temp->evalstatus = temp->prev->evalstatus;
    }
    else
    {
        cells[i][j] = c2->dval;
    }

    if (type == COMLL_T_SET || type == COMLL_T_INCREMENT)
    {
        temp             = comll_insert_after(temp);
        temp->type       = type;
        temp->ivec[0]    = j;
        temp->evalstatus = temp->prev->evalstatus;
    }

    /* Always ends the first cell COMLL by a POP comll to remove all value from the stack */
    if (c1 != NULL)
    {
        temp             = comll_gotobot(c1);
        temp             = comll_insert_after(temp);
        temp->type       = COMLL_POP;
        temp->evalstatus = temp->prev->evalstatus;
        temp1            = comll_gototop(c2);
        c2               = comll_join(c1, c2);
        temp             = c2;
        while (temp->nextcol)
            temp = temp->nextcol;
        temp->nextcol = temp1;
    }

    return (c2);
}

/* -------------------------------------------------------------------------------

   FUNCTION        :   parse_cell

   PARAMETERS      :   GrfnDeal    *gd
                       COMLL_PTR   *icomptr

   CALLED BY       :   create_event

   CALLS           :   yyparse

   DESCRIPTION     :   Parse a particular cell (GrfnCell) in the GRFN tableau.

                       For cells of type GRFNDCELL (cell containing a real
                       number) a single COMLL of type COMLL_REAL is created.

                       Cells of type GRFNBCELL are treated as being of
                       type GRFNDCELL (real) with value zero.

                       For cells of type GRFNSCELL the global variable
                       GRFN_INPUT_STRING is set to the cell's SVAL field,
                       and YYPARSE is called.  If YYPARSE sets the global
                       variable GRFN_LANG_CUR_AMCELL to be SRT_YES (meaning
                       this cell contains the "AM" qualifier), this fact is
                       recorded in the STATUS fields of the GrfnDeal for
                       the cell, row and the entire deal.

   Note            :   Cell types are set in the calling application i.e. for
                       swap_tools in SRT_F_GRF2020UTIL.C etc.

  ------------------------------------------------------------------------------- */

static Err parse_cell(GrfnDeal* gd, COMLL_PTR* icomptr)
{
    /*..Local variables..*/
    Err       error = NULL;
    COMLL_PTR icom  = NULL;
    COMLL_PTR bot   = NULL;

    /* Set the following GLOBAL variable to SRT_NO:
           Gets reset by the parser if it finds an AM  */
    grfn_lang_cur_amcell = SRT_NO;

    /*..Main Body..*/
    switch (GrfnCGCell(gd).type)
    {
        /* Cell contains nothing */
    case GRFNBCELL:

        GrfnCGCell(gd).dval = 0.0;
        GrfnCGCell(gd).type = GRFNDCELL;

    /* Cell contains a double */
    case GRFNDCELL:

        icom       = comll_insert_after(NULL);
        icom->type = COMLL_REAL;
        icom->dval = GrfnCGCell(gd).dval;
        break;

    /* Cell contains a string */
    case GRFNSCELL:

        /* Copy the cell from the GRFN tableau into the input string taken by YACC */
        grfn_input_string = GrfnCGCell(gd).sval;

        /* Call to YACC parser (This in turn calls LEX) */
        error = yyparse(&icom, gd);

        /* If parser found PAY in the cell, make sure nothing ele is hanging around outside it in
         * the cell */
        if (GrfnIsStatus(GrfnCGCell(gd).status, GRFNCSPAYREF))
        {
            bot = comll_gotobot(icom);

            /* If the last element is not of the PAY.. type, there is an operation outside PAY */
            if ((bot->type != COMLL_F_PAYFROMPAST) && (bot->type != COMLL_F_PAY) &&
                (bot->type != COMLL_F_PAYINPAST))
                return serror("[GRFN Error] %s in c[%d,%d]", GRERR_PAY_ALONE, gd->J, gd->I);
        }

        /* If parser found AM then set status of GRFN cell accordingly */
        if (grfn_lang_cur_amcell == SRT_YES)
        {
            GrfnAddStatus(GrfnCGCell(gd).status, GRFNCSAMERICAN);
        }

        /* Aggregate the current cell status into the row and full deal status */
        GrfnAggregateStatus(gd->rowstatus[gd->I], GrfnCGCell(gd).status);
        GrfnAggregateStatus(gd->sumstatus, GrfnCGCell(gd).status);

    } /* END switch loop on cell type */

    /* Returns the COMLL created by the LEX/YACC parser */
    *icomptr = icom;

    return error;

} /* END parse_cell */

/* ---------------------------------------------------------------------------------------- */

/* ----------------------------------------------------------------------------------------

    FUNCTION        :   create_event

    DESCRIPTION     :   Creates a GrfnEvent for a particular time.

    CALLED BY       :   grfn_attach_event

    DESCRIPTION     :   For each event date do the following:

                        o   If all the GrfnCells in the current row are of
                            type GRFNBCELL (blank), set *new_event_ptr = NULL

                        o   If gd->american==SRT_NO, we are actually at the
                            current row (and not between it and the next one)
                            Then for each cell do the following:

                                * Create a comll        (parse_cell)
                                * Join adjacent comlls  (join_adjacent_cells)

                        o   If gd->american==SRT_YES, we are generating an
                            event between two actual event dates. ONLY call
                            yyparse for cells that actually contain AM.

                            Since we are heading forwards in time, we must
                            have evaluated this row before when gd->american
                            SRT_NO. Hence we can look at the status fields
                            corresponding to this row within the GrfnDeal in
                            order to determine which cells we wish to reparse.

                        o   Set up all the df dates for this event
                                                    (grfn_set_df_dates_in_event )

                                                o   This function is called for PAST and FUTURE
  events, with the restriction that for PAST events, the gd->american flag is set to SRT_NO

     Modification    :   Eric Fournie
                                                Modification de la liste comll grfn qui va
  maintenant contenir un  champ permettant de se deplacer par colonne lors de l'evaluation.

  ---------------------------------------------------------------------------------------- */

static Err create_event(GrfnEvent** new_event_ptr, GrfnDeal* gd)
{
    /*..Local variables..*/
    COMLL_PTR  icom = NULL, temp = NULL;
    COMLLType  lattice_flg = COMLL_NULL;
    GrfnEvent* new_event;
    Err        err;
    int        null_row = 1;
    int        und_index;
    String     s;

    /* If all the cells in row are blank, then return NULL */
    for (gd->J = 0; gd->J < gd->sswidth; gd->J++)
        null_row *= (GrfnCGCell(gd).type == GRFNBCELL);

    if (null_row == 1)
    {
        *new_event_ptr = NULL;
        return (NULL);
    }

    /* Allocate space for the new event and check allocation status */
    new_event = (GrfnEvent*)srt_calloc(1, sizeof(GrfnEvent));

    if (new_event == NULL)
        return serror("Allocation failure in create_event");

    /* If it is not an American event */
    if (gd->american == SRT_NO)
    {
        /* Parse and links all cells in the current row (i.e. column 0 to n) */
        for (gd->J = 0; gd->J < gd->sswidth; gd->J++)
        {
            /* Sets the status of the current cell to its default values */
            GrfnResetStatus(GrfnCGCell(gd).status);

            /* Parse the currenct cell, using Lex & Yacc */
            err = parse_cell(gd, &icom);
            if (err)
                return err;

            /* add infos for ammc */
            temp = icom;
            if (!(GrfnCGCell(gd).status[GRFNCSPASTREF]) &&
                (GrfnCGCell(gd).status[GRFNCSPVREF] || GrfnCGCell(gd).status[GRFNCSFUTUREREF]))
                while (icom)
                {
                    icom->evalstatus = BWDCELL;
                    icom             = icom->next;
                }
            else
                while (icom)
                {
                    icom->evalstatus = FWDCELL;
                    icom             = icom->next;
                }
            icom = temp;

            /* Aggregate the current cell status in the event status */
            GrfnAggregateStatus(new_event->status, GrfnCGCell(gd).status);

            /* Sets the correspondiong flag if GRFNCSPVREF : COMLL_T_INCREMENT or COMLL_T_SET */
            set_lattice_flg;

            /* Adds the current cell COMLL to the already existing one of the event (for previous
             * cells) */
            new_event->icom =
                join_adjacent_cells(new_event->icom, icom, lattice_flg, gd->I, gd->J, gd->cells);

        } /* END loop on columns of event row */

    } /* END gd->american == SRT_NO */

    /* If it is an American event (this will not be called for historical events) */
    else
    {
        /* Parse all cells in the current row (i.e. column 0 to n) */

        for (gd->J = 0; gd->J < gd->sswidth; gd->J++)
        {
            /* Check if there is a call to AM in the cell ( it has already been analysed ) */
            if (GrfnIsStatus(GrfnCGCell(gd).status, GRFNCSAMERICAN))
            {
                /* 	Reparse the cell, to duplicate the event on each discretisation date */
                err = parse_cell(gd, &icom);
                if (err)
                    return err;

                /* add infos for ammc */
                temp = icom;
                if (!(GrfnCGCell(gd).status[GRFNCSPASTREF]) &&
                    (GrfnCGCell(gd).status[GRFNCSPVREF] || GrfnCGCell(gd).status[GRFNCSFUTUREREF]))
                    while (icom)
                    {
                        icom->evalstatus = BWDCELL;
                        icom             = icom->next;
                    }
                else
                    while (icom)
                    {
                        icom->evalstatus = FWDCELL;
                        icom             = icom->next;
                    }
                icom = temp;
                /* Set the flag for the COMLL_T: SET or INCREMENT if GRFNCSPVREF */
                set_lattice_flg;

                /* Joints the current cell to the previous ones, taking this flag into account */
                new_event->icom = join_adjacent_cells(
                    new_event->icom, icom, lattice_flg, gd->I, gd->J, gd->cells);
            }
            else
                /* Make sure that there is a NULL cash-flow for empty cell if in the CashFlow (last)
                   column */
                if (gd->J == gd->sswidth - 1)
            {
                /* Create a new node */
                icom = comll_insert_after(NULL);

                /* Add a POP (the previous value from the stack) - As we are in
           this loop we are guaranteed that there will be something on the stack */
                icom->type       = COMLL_POP;
                icom->evalstatus = BWDCELL;
                icom             = comll_insert_after(icom);
                icom->dval       = 0.0;
                icom->type       = COMLL_REAL;
                icom->evalstatus = BWDCELL;
                /* Joins this POP -> ( REAL = 0 ) element to the existing COMLL */
                new_event->icom = comll_join(new_event->icom, icom);
            }

        } /* END loop on all cells of current row */

    } /* END gd->american == SRT_YES */

    /* If we have a Historical American event then do not initialise event */
    /* OVE: AND WHAT IF THE CALL TO AM IS NOT IN THE LAST PARSED CELL ? */
    /* && GrfnIsStatus(gd->rowstatus[gd->I], GRFNCSAMERICAN); */
    if (gd->is_history == SRT_YES && grfn_lang_cur_amcell == SRT_YES)
    {
        *new_event_ptr = new_event;
        return NULL;
    }

    /* Set up and aggregates the dates (for events and df) in newly created event for the domestic
     * underlying first (index = 0) */
    und_index = 0;
    err       = grfn_set_df_dates_in_event(new_event, gd, und_index);
    if (err)
    {
        srt_free(new_event);
        return err;
    }

    /* Loop on all following underlyings referenced in the tableau */
    for (und_index = 1; (s = get_grf_und_name(und_index)) != NULL; und_index++)
    {
        /* Frees the memory allocated by get_grf_und_name(...) */
        if (s)
            srt_free(s);

        /* Set up and aggregate df dates for that event for that underlying */
        err = grfn_set_df_dates_in_event(new_event, gd, und_index);
        if (err)
        {
            srt_free(new_event);
            return err;
        }
    }

    /* The newly created event to be returned */
    *new_event_ptr = new_event;

    return err;

} /* END create_event */

/* --------------------------------------------------------------------------- */

/* -----------------------------------------------------------------------------------

    FUNCTION        :   grfn_create_future_event

    DESCRIPTION     :   When dealing with future events,

                        GrfnDeal->I corresponds to the index of the previous event date
                                                ( == previous row of the tableau )
                                                It is incremented only if the step date is on an
   EventDate
   ----------------------------------------------------------------------------------- */

Err grfn_create_future_event(
    GrfnDeal* gd, GrfnEvent** event, double step_ddate, double next_step_ddate)
{
    Err    error = NULL;
    double x;

    /* If the current step date is the same as the next event date in the Grfn tableau */
    if (step_ddate == gd->event_dates[gd->I + 1])
    {
        /* Move on to the next event in the tableau (the one that will be treated) */
        gd->I++;

        /* This is not an American (repeated) event */
        gd->american = SRT_NO;
        gd->prvdt    = gd->nowdt;
        gd->nowdt    = step_ddate;

        if (next_step_ddate != 0)
        {
            gd->nxtstpdt = next_step_ddate;
            gd->nxtevdt  = gd->event_dates[gd->I + 1];
        }

        /* Create an event for the ith event date */
        error = create_event(event, gd);
    }
    else
        /* If there is a call to AM in the previous row (already parsed) (if there is a previous
           row) */
        if ((gd->I >= 0) && (GrfnIsStatus(gd->rowstatus[gd->I], GRFNCSAMERICAN)))
    {
        /* Create AMERICAN event between the (i-1)th and ith event date */
        gd->prvdt    = gd->nowdt;
        gd->nowdt    = step_ddate;
        gd->american = SRT_YES;

        /* If the date is in the future */
        if (step_ddate > 0)
        {
            gd->amdt   = 0.0;
            gd->amdays = (long)step_ddate - gd->event_dates[gd->I];
        }
        /* If the step date is in the past */
        else
        {
            x          = (gd->nowdt - gd->event_dates[gd->I]);
            gd->amdays = DTOL(x);
            gd->amdt   = DMAX(x - (double)gd->amdays, 0.0);
        }

        if (next_step_ddate)
        {
            gd->nxtstpdt = next_step_ddate;
            gd->nxtevdt  = gd->event_dates[gd->I + 1];
        }

        /* ReParse the Ith event stored in the GrfnDeal, and returns an evaluable one */
        error = create_event(event, gd);

        /* Resets the American flag to NO */
        gd->american = SRT_NO;
    }
    /* This is just a discretisation date : no event in the tableau, no AM before : nothing */
    else
    {
        event = NULL;
    }

    return error;

} /* END Err grfn_create_future_event(...) */

/* ----------------------------------------------------------------------------------- */

/* -----------------------------------------------------------------------------------

    FUNCTION        :   grfn_create_historical_event

    DESCRIPTION     :   When dealing with historical events, the creation ( and
                            then the evaluation ) has to be done in order to set
                                                variables that might have an impact on future events
                                                (auxiliary variables: va, vb...)

                        For historical events, the treatment of AM has to be
                                                specific: only for the event dates in the tableau so
   far

                                                The GrfnDeal->ssparam and GrfnDeal->cells will be
                                                populated accordingly !

                        When fixings are encountered, they are going to be evaluated
                                                using the HistoricalFixings functions...

                        GrfnDeal->I corresponds to the index of the previous event date
                                                and  tableau row associated with the past date we
   are looking at
   ----------------------------------------------------------------------------------- */

Err grfn_create_historical_event(
    GrfnDeal* gd, GrfnEvent** event, double step_date, double next_step_date)
{
    Err error = NULL;

    /* Move on to the next event date (in the tableau) (i.e. the one to treat) */
    if (step_date == gd->event_dates[gd->I + 1])
    {
        /* The next row corresponds to the current step: go to this row  */
        gd->I++;
    }
    else
    {
        /* The current step is on a discretised broken date ( during today's day with e.o.d. flag
         * on!): no event to attach */
        event = NULL;

        return NULL;
    }

    /* Set the previous dates and now */
    gd->prvdt = gd->nowdt;
    gd->nowdt = gd->event_dates[gd->I];

    /* Set the next step dates and next event dates  */
    if (gd->I + 1 < gd->first_unkn_index)
    {
        /* If the next event date is before the first unknown date */
        gd->nxtstpdt = gd->nxtevdt = gd->event_dates[gd->I + 1];
    }
    else
    {
        /* If the current event date is today and the tableau is not finished */
        if ((gd->event_dates[gd->I] == gd->today) && (gd->I + 1 < gd->sslength))
        {
            gd->nxtstpdt = next_step_date;
            gd->nxtevdt  = gd->event_dates[gd->I + 1];
        }
        else
        /* If the current event date is before today or the tableau is finished */
        {
            gd->nxtevdt  = gd->today;
            gd->nxtstpdt = gd->today;
        }
    }

    /* By default, assume no Repeated (American) events in the past */
    gd->american = SRT_NO;

    /* Create the COMLL from the tableau and attaches it to the corresponding event */
    error = create_event(event, gd);
    if (error)
        return error;

    /* Return an error or a success message */
    return error;

} /* END grfn_treat_historical_event(...) */

/* ---------------------------------------------------------------------------- */

/* -----------------------------------------------------------------------------------

    FUNCTION        :   grfn_attach_past_to_future

    DESCRIPTION     :   After evaluating historical events, GrfnDeal->ssparams are
                            populated with values (va, vb...) implied by past events.

                        In order to keep track of these values that might impact on
                                                future cash flows, this function copies the ssparams
   produced by the evaluation of historical events into the GrfnDeal->hist_ssparams.

                                                In order to force the use of these historical
   values, a new RESET_STATE node will be atttached in the COMLL of today's GrfnEvent (passed as an
   input). If today's event does not exist (or is empty), it will be created with just this node in
   the COMLL

   ----------------------------------------------------------------------------------- */
Err grfn_attach_past_info_to_future(GrfnDeal* gd, GrfnEvent** today_event)
{
    Err err = NULL;

    /* If there are dates in the past, and variables (va, vb...)  have been referenced */
    if ((gd->first_unkn_index > 0) && GrfnIsStatus(gd->sumstatus, GRFNCSVARREF))
    {
        /* Copies the SsParams produced by the past evaluation into the historical ones for future
         * use */
        memcpy(&gd->hist_ssparam, &gd->ssparam, sizeof(GrfnSSParam));

        /* If today's event is empty, create a ZERO value one (for the POP() at the end of
         * eval_event) */
        if (!(*today_event))
        {
            (*today_event) = (GrfnEvent*)srt_calloc(1, sizeof(GrfnEvent));
            if ((*today_event) == NULL)
                return serror(" Allocation failure in grfn_attach_past_info_to_today ");

            /* Attaches a ZERO value REAL COMLL to the new event */
            (*today_event)->icom       = comll_insert_after(NULL);
            (*today_event)->icom->type = COMLL_REAL;
            (*today_event)->icom->dval = 0.0;

        } /* end if */

        /* Add a node in the COMLL of today's event to force the copy of hist_ssparam */
        (*today_event)->icom       = comll_insert_before((*today_event)->icom);
        (*today_event)->icom->type = COMLL_RESET_STATE;

    } /* END if */

    return NULL;

} /* END Err grfn_attach_past_to_future(...) */

/* ---------------------------------------------------------------------------- */
