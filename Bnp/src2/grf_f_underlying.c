/******************************************************************************\
*                Copyright (c) 1995 PARIBAS Capital Markets Group              *
********************************************************************************

    SYSTEM          :   GRF

    MODULE NAME     :   GRF_F_UNDERLYING

    PURPOSE         :   Retrieval of underlying information from GRFN tableau

    AUTHOR          :   M NADJM

    DATE            :

    VERSION         :

    DESCRIPTION     :

    FUNCTIONS USED  :

    EXT FUNCTIONS   :


********************************************************************************
*                           Amendment History                                  *
********************************************************************************

    AMENDED BY      :

    DATE            :

    VERSION         :

    REASON          :

    REQUEST/BUG NO  :

    DESCRIPTION     :

*******************************************************************************/

/*******************************************************************************
 *                           Import Include Files                               *
 *******************************************************************************/

#include "grf_h_lang.h"
#include "grf_h_lex.h"
#include "grf_h_types.h"
#include "utallhdr.h"

/*******************************************************************************
 *                           External variables                                 *
 *******************************************************************************/

EXTERND String grfn_input_string; /*  Source of input to lex_yy() */

/*******************************************************************************
 *                           Private Function Definitions                       *
 *******************************************************************************/

/* Underlying names found in the GRFN tableau */

#define GRFN_UND_LIST_SIZE 64 /* Maximum number of underlyings allowed in the list */

static String _grfn_und_list[GRFN_UND_LIST_SIZE];

static Err Check_Brackets(String s, int i, int j)

{
    String temp;  /* temporary storage for input string  */
    int    count; /* counter for left and right brackets */

    temp  = s;
    count = 0;

    /* Look for left/right brackets in string and increment/decrement
       counter accordingly
    */

    while (*temp != '\0')
    {
        if (*temp == '(')
            count++;
        else if (*temp == ')')
            count--;

        temp++;

    } /* end while */

    /* Since we are simply matching left and right brackets the counter i,
       must be zero if all brackets match, > 0 if more right brackets than
       left and < 0 otherswise
    */

    if (count < 0)
        return (serror("Error: Too many ) in c[%d,%d]", j, i));
    else if (count > 0)
        return (serror("Error: Too many ( in c[%d,%d]", j, i));
    else
        return (NULL);
}

static void reset_tableau(void)
{
    memset(_grfn_und_list, 0, GRFN_UND_LIST_SIZE * sizeof(String));
}

/*******************************************************************************

    FUNCTION    :   get_grf_und_index

    DESCRIPTION :   Given a string s, will return the positional index of s in
                    the list of underlying names if found, -1 otherwise.

*******************************************************************************/

int get_grf_und_index(String s)
{
    int i;

    i = 0;

    /* While not at the end of list */
    while (_grfn_und_list[i] != NULL)
    {
        /* If underlying name found then return index, else increment i */
        if (strcmp(_grfn_und_list[i], s) == 0)
            return (i);
        i++;
    }

    /* We have reached the end of the list (i.e. NULL) so return -1 */

    return (-1); /* Invalid index signifying error */

} /* end und_index */

/*******************************************************************************

    FUNCTION    :   get_grf_und_name

    DESCRIPTION :   Given an integer i, returns the string in the static array
                    _grfn_und_list[i] (if valid).

                    Returns NULL in case of error.

*******************************************************************************/

String get_grf_und_name(int index)
{
    String s; /* Local variable for returning result */
    int    s_len;

    /* Check range of index */

    if (index < 0 || index >= GRFN_UND_LIST_SIZE)
        return (NULL);

    /* If we have a valid name, then get its length,
       allocate space for string and copy it */

    if (_grfn_und_list[index] != NULL)
    {
        s_len = strlen(_grfn_und_list[index]);

        /* Allocate space for string - gets freed in calling routine */
        s = (String)srt_calloc(s_len + 1, sizeof(char));

        /* Check allocation */
        if (s == NULL)
        {
            return (NULL);
        }
        else
        {
            strcpy(s, _grfn_und_list[index]);
            return (s);
        }

    } /*end if */

    return (NULL);

} /* end get_grf_und_name */

/*******************************************************************************
 *                           Public Function Definitions                        *
 *******************************************************************************/

/*******************************************************************************

    FUNCTION    :   grfn_store_und_name

    DESCRIPTION :   Saves the name of the underlying curve found in the
                    GRFN tableau in the static string list _grfn_und_list.
                                        Names are stored in the same order as that found in the
                                        tableau, except domestic that is stored as _grfn_und_list[0]

*******************************************************************************/

int grfn_store_und_name(String s, String domestic)
{
    /*..Local variables..*/
    int i, s_len, d_len;

    /*..Initialise local variables..*/
    s_len = strlen(s);
    d_len = strlen(domestic);

    /*..Main Body..*/

    /* Store name of domestic underlying in list */

    if (_grfn_und_list[0] == NULL)
    {
        /* Allocate space for new string */
        _grfn_und_list[0] = (String)srt_calloc(d_len + 1, sizeof(char));

        /* Check allocation */
        if (_grfn_und_list[0] == NULL)
            return (-1);

        /* Copy domestic underlying name as first element of the list */
        strcpy(_grfn_und_list[0], domestic);
    }

    i = 0;

    /* Find the first available slot in the list to store new name */
    while ((_grfn_und_list[i] != NULL) && (i < GRFN_UND_LIST_SIZE))
    {
        /* If name is already in the list then return index */
        if (strcmp(_grfn_und_list[i], s) == 0)
            return (i);
        i++;
    }

    if (i >= GRFN_UND_LIST_SIZE)
        return (-1);

    /* Allocate space for new string */
    _grfn_und_list[i] = (String)srt_calloc(s_len + 1, sizeof(char));

    /* Check allocation */
    if (_grfn_und_list[i] == NULL)
        return (-1);

    /* Copy string to list */
    strcpy(_grfn_und_list[i], s);

    return (i);

} /* end grfn_store_und_name */

/*******************************************************************************

    FUNCTION        :   grfn_free_und_name

    PURPOSE         :   free the list of underlying names specified in the
                        GRFN tableau.

    CALLED BY       :   srt_f_grfn

    DESCRIPTION     :   This function frees the list of underlying used
                                                by GRFN to price a product.

*******************************************************************************/

void grfn_free_und_name()
{
    int i;

    for (i = 0; i < GRFN_UND_LIST_SIZE && _grfn_und_list[i] != NULL; i++)
        srt_free(_grfn_und_list[i]);
}

/*******************************************************************************

    FUNCTION        :   grfn_list_deal_underlyings

    PURPOSE         :   Returns a list of underlying names specified in the
                        GRFN tableau.

    CALLED BY       :   srt_f_grfn

    DESCRIPTION     :   This function (indirectly) creates a (static) list of
                        underlying names found in the GRFN tableau by parsing
                        all the cells row by row. The static list is created
                        by the call to grfn_interp_fin_func. This function
                        will know (according to the value of gd->pass) whether
                        or not a list needs to be created, i.e. if gd->pass = 1
                        then we are just parsing and extracting information
                        from the tableau and NOT creating events.

*******************************************************************************/

Err grfn_list_deal_underlyings(char** list, GrfnDeal* gd)

{
    String    s;
    int       i, s_len;
    COMLL_PTR temp;
    Err       error;

    /* Set gd->pass to 1 => We are ONLY parsing */
    gd->pass = 1;

    /* Reset list of underlying names i.e. delete previous names */
    reset_tableau();

    /* Parse all cells in the tableau row by row */

    for (gd->I = 0; gd->I < gd->sslength; gd->I++)
    {
        /* is used by IFFIX */
        gd->nowdt = gd->event_dates[gd->I];
        for (gd->J = 0; gd->J < gd->sswidth; gd->J++)
        {
            /* If we have a cell containing a string ... */
            if (GrfnCGCell(gd).type == GRFNSCELL)
            {
                /* Get string taken from the GRFN tableau */
                grfn_input_string = GrfnCGCell(gd).sval;

                error = Check_Brackets(GrfnCGCell(gd).sval, gd->I, gd->J);
                if (error)
                    return error;

                /* Call YACC parser (This in turn calls LEX) */
                error = yyparse(&temp, gd);
                if (error)
                    return error;

                comll_free_list(temp);

            } /* end if */

        } /* end for */

    } /* end for */

    /*We start getting the underlying names at index 1 NOT 0 ( the domestic one is in position 0)*/
    i = 1;

    while (1)
    {
        /* Get name of underlying : Free after use */
        s = get_grf_und_name(i);

        if (s == NULL)
        {
            /* We have reached the end of the list of names : make it a blank one (not to loose
             * address) */
            strcpy(list[i], "");
            break;
        }
        else
        {
            /* Determine length of string */
            s_len = strlen(s);

            if (list[i] != NULL)
            {
                /* Copy name  list */
                strcpy(list[i], s);

                /* Free s */
                srt_free(s);

                i++;
            }

        } /* end if */

    } /* end while */

    /* Set gd->pass to 2 => We will be creating events */

    gd->pass = 2;

    return (NULL);

} /* end grfn_list_underlying */
