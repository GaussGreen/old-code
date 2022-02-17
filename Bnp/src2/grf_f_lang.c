/******************************************************************************\
*                Copyright (c) 1995 PARIBAS Capital Markets Group              *
********************************************************************************
                                                                              
    SYSTEM          :   GRF
                                                                              
    MODULE NAME     :   GRF_F_LANG                                        

    PURPOSE         :   GRFN Language Syntax Checking

    AUTHOR          :

    DATE            :

    VERSION         :


    DESCRIPTION     :   For functions of the form f(a,b,c,...) we have already
                        evaluated the arguments and joined them together. Hence
                        COMLL_PTR argument to these functions is a list:
                                        "c"-->"b"-->"a".  

                        For other functions we simply check the number of
                        arguments and then append the function name: 
                                        "c"-->"b"-->"a"-->"f";
  
                        In other cases we will evaluate or partially evaluate 
                        the arguments, replacing what was there before with a 
                        single "f(a,b,c)". Back in YACC, the pattern matching 
                        code contains a pointer to the TOP of this list, which
                        we can NOT give up. Therefore in the second case we
                        will free this part of the list: 
                                            -->"b"-->"a"  
                        and replace "c" with "f(a,b,c)". 

                        Other important modules called from this one include:

                        o grf_f_lngsymtab.c (GRFN language symbol table)
  
                        o grf_f_lang_finance.c (Calls to financial functions)
  
                        o grf_f_lang_pvrng.c (grfn function pvrng())


    FUNCTIONS USED  :   grfn_interp_name        (Variable names)
                        grfn_interp_func        (Math/Underlying)
                        grfn_interp_cell_func   (Cell functions)
                        
    PROBLEMS        :   Error checking may not respond well to character 
                        strings used in the wrong place, e.g. "annual" + 2 
                        should return an error BUT it will not!

                        Memory allocation in comll_insert_after/before 
                        are not checked for failure.

                        What is needed is the following:
                        1) A GrfnSymbol should contain a more formal function
                           prototype for the types of its arguments.
                        2) Need a single function to check the arguments to 
                           any function, taking as inputs a prototype and a 
                           COMLL representing the arguments.
  
   
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
*   DESCRIPTION     :
*
\******************************************************************************/
 


/******************************************************************************\
*                           Import Include Files                               *
\******************************************************************************/

 
#include "grf_h_all.h"
#include "math.h"
#include "OPFNCTNS.H>
#include "RainbowOpt.h"

/******************************************************************************\
*                           Private Function Definitions                       *
\******************************************************************************/


static int n_end_r_args(COMLL_PTR args)

{
    int n = 0;
    COMLL_PTR last = comll_gotobot(args);

    while (last->prev !=NULL && last->type==COMLL_REAL) {
        n++;
        last = last->prev;
    }

    n += (last->type==COMLL_REAL);

    return n;

}

/* ------------------------------------------------------------- */
/* Gets the nargs first arguments in a COMLL (from the bottom up) */
static double *get_det_args(int nargs, COMLL_PTR args)

{
    int i;
    double *arg_vec = NULL;
    arg_vec = srt_calloc(nargs,sizeof(double));
    args = comll_gotobot(args);

    for (i=0; i < nargs; i++)
	{
        arg_vec[i] = args->dval;
        args = args->prev;
    }

    return arg_vec;

} 

/* ------------------------------------------------------------- */

static Err check_string_arg(COMLL_PTR args, COMLLType type)
{
    Err             error;
    SrtCallPutType callputtype;

    args = comll_gototop(args);

	if ((type==COMLL_A_BLKSCHLS) || (type==COMLL_A_BLKSCHLSNRM) || (type==COMLL_A_BSSABR) 
		|| (type==COMLL_A_SPREADNUM))
    {
		if ((args->type != COMLL_STRING)
		   || (error = interp_call_put(args->sval,&callputtype))) 
		{
			return serror("%s: %s",GRERR_BAD_CALL_PUT,args->sval);
		}
    
		args->type = COMLL_REAL;
	
		if (callputtype == SRT_CALL) {
			args->dval = 0.0;
		}
		else {
			args->dval = 1.0;
		}   
    }

	return NULL;
}


/******************************************************************************\
 
    FUNCTION        :   check_cell_coords
    DESCRIPTION     :   Check cell co-ordinates col, row1, row2 specified
                        in cell functions. All cells in the implied range 
                        must be already computed either on the same row and
                        to the left of, or above, the current cell.

\******************************************************************************/

                        
static Err check_cell_coords    (   
				GrfnDeal    *gd, 
				int         *detargs,
				String      name        )
{
long    c , r1, r2;
long    r , c1, c2, ind;

    if (  (strcmp(name,"ROWSORT") == 0) 
		|| (strcmp(name,"ROWINTERP") == 0) )
	{ 
	/* Gets the parameters of the ROWSORT function */ 
		r   = detargs[0];
		c1  = detargs[1];
		c2  = detargs[2];
		ind = detargs[3];

	/* If the row to sort is before the current event: reference to past */
		if (r < gd->I)
		{
			GrfnAddStatus(GrfnCGCell(gd).status, GRFNCSPASTREF);
		}
		
	/* If arguments are within the normal bounds, nothing to do */
		if (0 <= r && r <= gd->I && c1 >= 0 && c2 >= c1 && c2 < gd->J 
			&& ind <= (c2-c1) && ind >= 0)
			return NULL;
	
	/* If there is a reference to the future, mention it */
		if ( ( r > gd->I) ||
			(( r == gd->I ) && ( c2 >= gd->J ) ) )
		return serror("%s: %s[%d,%d,%d,%d]",
				GRERR_ILLEGAL_REF,name,r,c1,c2, ind);	

	/* If the arguments are out of the tableau range, mention it */	
		if (r < 0 || r > gd->sslength-1 || c1 < 0 || c1 > gd->sswidth-1
			|| c2 < 0 || c2 > gd->sswidth-1)
			return serror("%s: %s[%d,%d,%d,%d]",
				GRERR_REF_OUT_OF_RANGE,
				name,r,c1,c2,ind);

	/* If the second column reference is before the first one */
		if (c2 < c1)
			return serror("%s: %s[%d,%d,%d,%d]",
				GRERR_COL_OUT_OF_RANGE,
				name,r,c1,c2,ind);

	/* If the index required is not within the number of columns on which to sort, mention it */ 
		if ((ind < 0 ) || ( ind > (c2 - c1 )) )
			return  serror("%s: %s[%d,%d,%d,%d]",
				GRERR_INDEX_OUT_OF_RANGE,
				name,r,c1,c2,ind);
	}
	else if ( (strcmp(name,"ROWMIN") == 0) 
			  || (strcmp(name,"ROWMAX") == 0)  
			  || (strcmp(name,"ROWSUM") == 0)  
			  || (strcmp(name,"ROWAVG") == 0) 
			  || (strcmp(name,"ROWMAXIDX") == 0)
			  || (strcmp(name,"ROWMINIDX") == 0))
	{ 
	/* Gets the parameters of the ROWSORT function */ 
		r   = detargs[0];
		c1  = detargs[1];
		c2  = detargs[2];

	/* If the row to sort is before the current event: reference to past */
		if (r < gd->I)
		{
			GrfnAddStatus(GrfnCGCell(gd).status, GRFNCSPASTREF);
		}
		
	/* If arguments are within the normal bounds, nothing to do */
		if (0 <= r && r <= gd->I && c1 >= 0 && c2 >= c1 && c2 < gd->J )
			return NULL;
	
	/* If there is a reference to the future, mention it */
		if ( (r > gd->I) || (( r == gd->I ) && ( c2 >= gd->J )) )
		return serror("%s: %s[%d,%d,%d]",
				GRERR_ILLEGAL_REF,name,r,c1,c2);	

	/* If the arguments are out of the tableau range, mention it */	
		if (r < 0 || r > gd->sslength-1 || c1 < 0 || c1 > gd->sswidth-1
			|| c2 < 0 || c2 > gd->sswidth-1)
			return serror("%s: %s[%d,%d,%d]",
				GRERR_REF_OUT_OF_RANGE,
				name,r,c1,c2);

	/* If the second column reference is before the first one */
		if (c2 < c1)
			return serror("%s: %s[%d,%d,%d]",
				GRERR_COL_OUT_OF_RANGE,
				name,r,c1,c2);

	}
	else  /* perform check on COLXXX functions */
	{
	/* Gets the parameters of the COL... function: Column, RowStart, RowEnd */ 
		c   = detargs[0];
		r1  = detargs[1];
		r2  = detargs[2];

	/* If there is a reference to the past, adds it to the cell status */
		if (r1 < gd->I)
			GrfnAddStatus(GrfnCGCell(gd).status, GRFNCSPASTREF);

	/* If parameters are within bound, nothing to do */
		if ((-1 < c && c < gd->J && r1 > -1 && r2 >= r1 && r2 <= gd->I) ||
		   (gd->J <= c && c < gd->sswidth && r1 > -1 && r2 >= r1 && r2 < gd->I))
			return NULL;
		
	/* If parameters are out of bound, mention it */
		if (c < 0 || c > gd->sswidth-1 || r1 < 0 || r1 > gd->sslength-1
			|| r2 < 0 || r2 > gd->sslength-1)
			return serror("%s: %s[%d,%d,%d]",
				GRERR_REF_OUT_OF_RANGE,
				name,c,r1,r2);
	
	/* If the second row reference is before the first one */
		if (r2 < r1)
			return serror("%s: %s[%d,%d,%d,%d]",
				GRERR_ROW_OUT_OF_RANGE,
				name,c,r1,r2);

	/* If there is an illegel reference to the future, mention it */
		if ( ( r2 > gd->I) ||
			(( c == gd->I ) && ( r2 >= gd->J ) ) )
		return serror("%s: %s[%d,%d,%d]",
				GRERR_ILLEGAL_REF,
				name,c,r1,r2);
	}
	return NULL;
}


/******************************************************************************\
*                       Public Function Definitions                            *
\******************************************************************************/



/******************************************************************************\

    FUNCTION        :   grfn_interp_cell_func        

    DESCRIPTION     :   Check a cell function.

                        A cell function is one which takes square brackets
                        which take deterministic arguments

						Called by the LEX/YACC parser when encountering a 
						NAME[...] type of input string

						Arguments are passed BACKWARDS by the parser 
						(see the 'Args' section in Grfn_Yacc.y)
\******************************************************************************/


Err grfn_interp_cell_func   (   String      name, 
                                COMLL_PTR   args,
                                GrfnDeal    *gd     )

{

/*..Local variables..*/

Err         error;
GrfnSymbol  gs;
COMLL_PTR   lastarg;
int         i, c, r1, r2, ind, dateindex, rindex, cindex;
int         detargs[4];

/* Check function name is in GRFN symbol table. If so return 
   information relating to name (including COMLL type ) */
    if (error = grfn_interp_comll_name(name, &gs)) 
		return error;

/* Point to first argument (they are linked backwards by the parser) */
    lastarg = comll_gotobot(args);

/* Check number of arguments between parsed and required (by the GrfnSymbol) */
    if (args->nargs != gs.ndargs + gs.nrargs)
        return serror("%s %s %d, requires %d args",
					   name,GRERR_WRONG_NUM_ARGS,
					   args->nargs,gs.ndargs+gs.nrargs); 

/* Look at nametype and take appropriate action */
 
    switch(gs.nametype)
	{

        case CLLNCellFunc:

        /* Makes sure all the arguments are REAL and collects them by function order (back of COMLL) */
			for (i=0; i < gs.ndargs;i++)
			{
                test_real(i,name,lastarg);
                detargs[i] = DTOL(lastarg->dval);
                lastarg = lastarg->prev;
            }
		
		/* Checks nothing wrong has been input (illegal references,...) */
            if (error = check_cell_coords(gd,detargs,name)) 
				return error;

		/* C is the column index , except for ROWSORT/ROWINTERP where it is the row index */
            c  = detargs[0];

		/* R... is the row index, except for ROWSORT/ROWINTERP where it is the column index */
            r1 = detargs[1];
            r2 = detargs[2];

		/* The final argument of ROWSORT, COLSORT or ROWINTERP: the index wanted */
			if ( gs.type == COMLL_C_COLSORT || gs.type == COMLL_C_ROWSORT   || gs.type == COMLL_C_ROWINTERP  )
				ind = detargs[3];   
			else
				ind = 0;
		
		/* Sets the COMLL type */
	        args->type = gs.type;
                    
        /* Treats ROWINTERP separately first */
			if (gs.type == COMLL_C_ROWINTERP)
			{
			/* Checks the auxiliary range has the right length for interpolation (and that ind is not outside of range) */
				CheckAuxRange1(ind,gd);
				if (gd->auxlen[ind] != (r2-r1))
					serror("%s: lengths %d,%d", GRERR_INTERP_SAME_SIZE, gd->auxlen[ind], (r2 - r1) );

			/* For ROWINTERP, the value to interpolate (x) is not deterministic : need a TWO nodes COMLL */
				args->next->next->prev = NULL;
				comll_free_list(args->next->next);
				args->next->next = NULL;

			/* The _ROW_INTERP COMLL will be the second node in the comll */
				args->next->type = gs.type;

			/* Populate the COMLL node with the indexes: aux, row, col_s, col_e */
				memcpy(args->next->ivec, detargs, 3 * sizeof(int));

			/* The (complex) argument at the top of the list is x (args points on the last argument of "x") */
						
				return NULL;        

			}
			else
			{
			/* We need only one element in this function COMLL (storing everything in ivec) */
				args->next->prev = NULL;
				comll_free_list(args->next);
				args->next = NULL;

		   /* Populate the COMLL for this Cell function */
				memcpy(args->ivec, detargs, 4 * sizeof(int));

				return NULL;
			}
        case CLLNCellRef:
		/* Gets the Cell Rows and Columns */
            test_real(0,name,lastarg);
            test_real(1,name,lastarg->prev);
            cindex = DTOL(lastarg->dval);
            rindex = DTOL(lastarg->prev->dval);

        /* Checks the Rows and Columns Vs the Tableau size */
			if ( cindex < 0 || cindex > gd->sswidth-1 ||
                rindex < 0 || rindex > gd->sslength-1)
                return serror("%s C[%d,%d]",
                    GRERR_REF_OUT_OF_RANGE,cindex,rindex);

		/* Adds the Past Reference status if reference to a row before current */
			if (rindex < gd->I)
               GrfnAddStatus(GrfnCGCell(gd).status,GRFNCSPASTREF);
        
		/* Returns an error if illegal reference to future (no logical meaning) */
			if (rindex > gd->I)
               return serror("%s: C[%d,%d]",
				GRERR_ILLEGAL_REF,cindex, rindex);
	
		/* Builds the corresponding COMLL_CONST element refering to the gd->cells */
			lastarg = lastarg->prev;
            lastarg->next->prev = NULL;
            comll_free_list(lastarg->next);
            lastarg->next = NULL;
		    lastarg->type = COMLL_CONST;
            lastarg->gdcells = &gd->cells[rindex][cindex];
        
            return NULL;

        case CLLNAuxRef:

            test_real(0,name,lastarg);
            test_real(1,name,lastarg->prev);
            cindex = DTOL(lastarg->dval);
            rindex = DTOL(lastarg->prev->dval);
            if ( cindex < 0 || rindex < 0 ||
                cindex > gd->auxwidth-1 ||
                rindex > gd->auxlen[cindex]-1 )
                    return serror("%s A[%d,%d]",
                    GRERR_REF_OUT_OF_RANGE,cindex,rindex);
            lastarg = lastarg->prev;
            lastarg->next->prev = NULL;
            comll_free_list(lastarg->next);
            lastarg->next = NULL;
            lastarg->type = COMLL_REAL;
            lastarg->dval = gd->aux[cindex][rindex];
            return NULL;

        case CLLNDateRef:

            test_real(0,name,lastarg);
            dateindex = DTOL(lastarg->dval);

            if ( dateindex < 0 || 
              dateindex >= gd->sslength)
                    return serror("%s: D[%d]",
                    GRERR_BAD_DATE_REF,dateindex);
            lastarg->dval = (double)gd->event_dates[dateindex];
            lastarg->type = COMLL_REAL;
            return NULL;

        case CLLNPvRef:

/*            test_real(0,name,lastarg);*/
            cindex = DTOL(lastarg->dval);
            if (cindex < 0 || cindex > gd->sswidth-1)
                return serror("%s : PV[%d]",
                    GRERR_ILLEGAL_REF,cindex);
            lastarg->ivec[0] = cindex;
            lastarg->type = COMLL_T_REF;
        /* Sets both PV and Future Reference Status...(so far equivalent) */
			GrfnAddStatus(GrfnCGCell(gd).status,GRFNCSFUTUREREF);
			GrfnAddStatus(GrfnCGCell(gd).status,GRFNCSPVREF);

            return NULL;        
// IMPLEMENTATION OF PVR IN GRFN		

		case CLLNPvRRef:

			test_real(0,name,lastarg);
			cindex = DTOL(lastarg->dval);
            if (cindex < 0 || cindex > gd->sswidth-1)
                return serror("%s : PVR[%d]",
                    GRERR_ILLEGAL_REF,cindex);
			lastarg->ivec[0] = cindex;
            lastarg->type = COMLL_T_RISKY_REF;
        // Sets past Reference Status...(so far equivalent)
			GrfnAddStatus(GrfnCGCell(gd).status,GRFNCSFUTUREREF);
			GrfnAddStatus(GrfnCGCell(gd).status,GRFNCSPVREF);
			
			return NULL;

        case CLLNCumRef:

 			test_real(0,name,lastarg);
            cindex = DTOL(lastarg->dval);
            if (cindex < 0 || cindex > gd->sswidth-1)
                return serror("%s : PV[%d]",
                    GRERR_ILLEGAL_REF,cindex);
			lastarg->ivec[0] = cindex;
            lastarg->type = COMLL_T_CUMREF;
        /* Sets past Reference Status...(so far equivalent) */
			GrfnAddStatus(GrfnCGCell(gd).status,GRFNCSFUTUREREF);
			GrfnAddStatus(GrfnCGCell(gd).status,GRFNCSPVREF);
 
            return NULL;        

       case CLLNCurRef:

 			test_real(0,name,lastarg);
            cindex = DTOL(lastarg->dval);
            if (cindex < 0 || cindex > gd->sswidth-1)
                return serror("%s : PV[%d]",
                    GRERR_ILLEGAL_REF,cindex);
			lastarg->ivec[0] = cindex;
            lastarg->type = COMLL_T_CURREF;
        /* Sets PV Reference Status...(so far equivalent) */
			GrfnAddStatus(GrfnCGCell(gd).status,GRFNCSFUTUREREF);
			GrfnAddStatus(GrfnCGCell(gd).status,GRFNCSPVREF);
 
            return NULL;        

	   default:

            return serror("%s %s[]",GRERR_UNKNOWN_SYMBOL,name);

    } /* END switch gs.nametype */

} /* end grfn_interp_cell_func */
                   



/******************************************************************************\
*
*   FUNCTION        :   grfn_interp_func

    PARAMETERS      :   String      name
                        COMLL_PTR   args
                        GrfnDeal    *gd


    DESCRIPTION     :   Interprets f(...) where f is a reserved GRFN function
						Only the name of the function is passed, as well as the
						COMLL that corresponds only to the ARGUMENTS of the function
						No specific COMLL element is yet built for the function name
  
						The COMLL_PTR args returned ALWAYS points to the top of the 
						COMLL corresponding to this function (not to break the
						full COMLL construction in the YACC parsing)
						This top correponds to the LAST argument input to the 
						function (see the Grfn_Yacc.y file: reverse building)

                        
						It also performs error checking on its arguments, and 
                        modifies its arguments  appropriately.  If the result 
                        of the function can be computed now, this is done.
						

                        Example: If name is "max", args would have to describe
                        two  expressions, "a"-->"b".  For "a" and "b" both real 
                        numbers, we replace this with one node of type 
                        COMLL_REAL and value max(a,b).  
                        Otherwise, we would place another node at the end of 
                        the function of type COMLL_A_MAX: "a"-->"b"-->"max". 


    PROBLEM         :   In order to pre-compute results when arguments are 
                        known, there is duplication within this function of 
                        evaluation of all arithmetic functions with names, 
                        like max, exp, etc. This is confusing in that other
                        arithmetic functions are precomputed directly by YACC.
                        ( e.g. +, -, * ...). Also, if a new arithmetic function
                        with a name were added (sin(x) say), 3 modules rather
                        than 2 would need modification. (grf_f_eval_event.c 
                        and grf_f_lang.c to say how to evaluate it, and 
                        grf_f_lnfsymtab.c to declare it).

                        This could lead to possible inconsistencies in the 
                        evaluation of the "same" function in different places.

						Function called by the YACC parser when encountering a
						NAME(...) type of input

\******************************************************************************/



Err grfn_interp_func(String name, COMLL_PTR args, GrfnDeal *gd)

{

/*..Local variables..*/
    
Err         error;
GrfnSymbol  gs;                     /* Local copy of symbol table entries */
COMLL_PTR   lastarg, new;
long        detargs[4],l;
int         i, numenddetargs;       
int         flag;                   /* Set to 1 if we can calculate 
                                       arithmetic functions, 0 otherwise  */
int         index;                  /* Underlying index                   */
double      *arg_vec = NULL;        /* Temp storage for det arguments     */
double      x, f, k, v, t, d;
double		fwdx, fwdy, sigx, sigy, a, b, rho;
SrtCallPutType call_put;

/*..Main Body..*/


/* Check to see if name is a valid GRFN function - If so then
   return the relevant information contained in the symbol table
   e.g. symbol & COMLL types, number of (non)deterministic arguments */
	if ( error = grfn_interp_comll_name(name, &gs)) 
		return error;


/* Point to the FIRST argument of this function ( the last one put into the COMLL ) */
    lastarg = comll_gotobot(args);


/* Check number of arguments n:
       det_args + non_det_args <= n <= det_args + non_det_args + opt_args  */
	if ((args->nargs < gs.ndargs + gs.nrargs) ||
        (args->nargs > gs.ndargs + gs.nrargs  + gs.no_opt_args)) 
        return (serror("%s %s %d, requires %d args",
                name, GRERR_WRONG_NUM_ARGS,
                args->nargs, gs.ndargs+gs.nrargs)); 


    switch (gs.nametype)
	{

       /* Auxiliary Range Function */
           case CLLNAuxFunc:
		    
            if (gs.type==COMLL_X_INTERP)
			{      /* Interpolate Function */

                /* NB. Since we have interp(r1,r2,x), then r1 and r2 must be 
                   deterministic, but x may not be. Hence r1 and r2 will 
                   be the last two nodes of args; there may be any number of 
                   nodes above representing x */
        
                
                /* If all the arguments are real then store them in detargs */

                for (i=0; i < gs.ndargs; i++)
				{
                    test_real(i, name, lastarg);
                    detargs[i] = DTOL(lastarg->dval);
                    lastarg = lastarg->prev;
                }

                lastarg = lastarg->next;
                lastarg->next->prev = NULL;
                comll_free_list(lastarg->next);
                lastarg->next = NULL;
                lastarg->type = gs.type;

                /* The following are (row,col) in auxiliary ranges */

                lastarg->ivec[0] = detargs[0];
                lastarg->ivec[1] = detargs[1];

                /* Are the auxiliary ranges valid */

                CheckAuxRange2(detargs[0], detargs[1], gd);

                return NULL;
                break;

            }
            else if (gs.type==COMLL_X_PVINTERP)
			{      /* Interpolate Function */

                /* PVInterp(AuxIndex,StartCol,ColSkip,x).  The first three arguments are
                   deterministic, but x may not be. Hence the deterministic arguments will be in the last 
                   be the last three nodes of args; there may be any number of 
                   nodes above representing x */
        
                
                /* If all the arguments are real then store them in detargs */
                for (i=0; i < gs.ndargs; i++)
				{
                    test_real(i, name, lastarg);
                    detargs[i] = DTOL(lastarg->dval);
                    lastarg = lastarg->prev;
                }

                lastarg = lastarg->next;
                lastarg->next->prev = NULL;
                comll_free_list(lastarg->next);
                lastarg->next = NULL;
                lastarg->type = gs.type;

                /* Store the indices in the in local index vector */
                lastarg->ivec[0] = detargs[0]; /* auxilliary index */
//                lastarg->ivec[1] = detargs[1]; /* row */
                lastarg->ivec[1] = detargs[1]; /* first col */
                lastarg->ivec[2] = detargs[2]; /* skip col */

                /* Are the auxiliary ranges valid */
                CheckAuxRange1(detargs[0], gd);

/* Stuff from PV to store result
            lastarg->ivec[0] = cindex;
            lastarg->type = COMLL_T_REF;*/
			GrfnAddStatus(GrfnCGCell(gd).status,GRFNCSFUTUREREF);
			GrfnAddStatus(GrfnCGCell(gd).status,GRFNCSPVREF);



                return NULL;
                break;

            }



	/* Arithmetic functions */
	    case CLLNArithFunc:
			
            /* Is last argument for Black-Scholes a string? */
            
            if ( ((gs.type==COMLL_A_BLKSCHLS) || (gs.type==COMLL_A_BLKSCHLSNRM) || (gs.type==COMLL_A_BSSABR) 
				|| (gs.type==COMLL_A_SPREADNUM))
				 && (error=check_string_arg(args, gs.type))) 
                return error; 

            /* If the arguments are deterministic, then compute result */

            flag			= 0;
            arg_vec			= NULL;
            numenddetargs	= n_end_r_args(lastarg);

            if (numenddetargs == gs.ndargs + gs.nrargs) {

                /* Copy pre-determined args of function to the array arg_vec 
                   NB. This function allocates memory for arg_vec */

                arg_vec = get_det_args(numenddetargs, lastarg);

                if (arg_vec==NULL)
                    return (serror("** Allocation failure in get_det_args **"));

                /* IMPORTANT: If flag = 1, then it is implied that the value of 
                              the arithmetic function can be computed.

                              If flag = 0, then we can NOT evaluate the function
                              flag is NOT currently set to 0. However, we may
                              need this functionality in future for arithmetic
                              functions which may not be pre-computable. */

                /* Perform arithmetic operation according to function type */
                
                switch (gs.type) {
            
                case COMLL_A_POW:
                    flag = 1;
                    x = pow(arg_vec[0], arg_vec[1]);
                    break;               

                case COMLL_A_EXP:
                    flag = 1;
                    x = exp(arg_vec[0]);
                    break;

                case COMLL_A_LOG:
                    flag = 1;
                    x = log(arg_vec[0]);
                    break;

                case COMLL_A_MAX:
                    flag = 1;
                    x = DMAX(arg_vec[0], arg_vec[1]);
                    break;

                case COMLL_A_MIN:
                    flag = 1;
                    x = DMIN(arg_vec[0], arg_vec[1]);
                    break;

                case COMLL_A_ABS:
                    flag = 1;
                    x = fabs(arg_vec[0]);
                    break;

                case COMLL_A_FLOOR:
                    flag = 1;
                    x = arg_vec[0]/arg_vec[1];
                    l = DTOL(x);
                    x = (double)l*arg_vec[1];
                    break;

                case COMLL_A_CEILING:
                    flag = 1;
                    x = -arg_vec[0]/arg_vec[1];
                    l = DTOL(x);
                    x = (double)l* -arg_vec[1];
                    break;

                case COMLL_A_NORM:
                    flag = 1;
                    x = norm(arg_vec[0]);
                    break;

                case COMLL_A_IF:
                    flag = 1;
                    x = (arg_vec[0] > 0.5 ? arg_vec[1] : arg_vec[2]);
                    break;

                case COMLL_A_BLKSCHLS:
                    flag = 1;
                    f  = arg_vec[0];
                    k  = arg_vec[1];
                    v  = arg_vec[2];
                    t  = arg_vec[3];
                    v *= v * t;
                    if (v <=0)
					{
						if (arg_vec[4] < 0.5)
		                    x =  f - k;
			            else  
				            x = k - f;
						if ( x < 0 )
							x = 0;
					}
					else
					{
						v = sqrt(v) ;
						d  = (log(f/k) - 0.5 * v * v ) / v;

	                    if (arg_vec[4] < 0.5)
		                    x =  f*norm( d + v ) - k*norm( d);
			            else  
				            x = -f*norm( - d - v ) + k*norm(-d);
                    }
					break;
 
				case COMLL_A_BLKSCHLSNRM:
                    flag = 1;
                    f  = arg_vec[0];
                    k  = arg_vec[1];
                    v  = arg_vec[2];
                    t  = arg_vec[3];
                    v *= v * t;
					if ( v <= 0 )
					{
						if (arg_vec[4] < 0.5)
		                    x =  f - k;
			            else  
				            x = k - f;
						if ( x < 0 )
							x = 0;
					}
                    else
					{
						v = sqrt(v);
						d  = (f-k)/v;

						if (arg_vec[4] < 0.5)
							x = v*gauss(d) + (f-k)*norm(d);
						else  
							x = v*gauss(d) - (f-k)*norm(-d);
                    }
					break;

				case COMLL_A_SPREADNUM:
                    flag = 1;

                    fwdx = arg_vec[0];
					fwdy = arg_vec[1];
					sigx = arg_vec[2];
					sigy = arg_vec[3];
					rho  = arg_vec[4];
					a	 = arg_vec[5];
					b	 = arg_vec[6];
					k    = arg_vec[7];
					t	 = arg_vec[8];
					call_put = (int) arg_vec[9];

					error = OptSpreadNew(
									fwdx,
									a,
									sigx,
									fwdy,
									-b,
									sigy,
									k,
									t, 
									rho,
									call_put,
									&x);

					if (error)
						return error;
					
					break;
                
                } /* end switch */
          
            } /* end if */

            /* Free memory allocated for temporary array arg_vec */

            if (arg_vec) srt_free(arg_vec);

            /* If we have succeeded in evaluating the function, then
               update the COMLL (i.e. replace the function with the 
               computed result) and cut the link to the function and 
               its arguments. Otherwise append the argument to the COMLL .
            */

            if (flag) {
                args = comll_gototop(args);
                args->dval = x;
                args->type = COMLL_REAL;
                comll_cut(args);
                return NULL;
            }
            else {
                new = comll_insert_after(lastarg);
                new->type = gs.type;
                return NULL;
             } /* end if */      
          
            break;
        
	/* Financial functions */
        case CLLNFinFunc:

            return grfn_interp_fin_func(name, &gs, gd, lastarg);
            break;

	/* Utility functions   */
		case CLLNUtil:
			
            /* We currently have only one utility function specified
               viz. Histogram function. Consequently we know that
               the first argument is a string and the second is a positive
               real number.
            */
            
            /* Point to first argument */
            lastarg = comll_gotobot(args);

            /* If argument is not a string, then return an error */
            if (lastarg->type != COMLL_STRING)
                return serror("%s: %s",GRERR_BAD_TYPE, lastarg->sval);

            lastarg->type = COMLL_HIST;
            return NULL;
            
            break;        

	/* S(...): Equity/Forex Underlying */
        case CLLNSample:
            
            /* Check argument type - must be a string */

            if (lastarg->type != COMLL_STRING) 
                return serror("%s: %s", name, GRERR_WRONG_TYPE_OPT_ARG);
		
            if (gd->pass == 1) 
			{
			/* Store the underlying name in the grfn list of underlyings  */
                i = grfn_store_und_name(lastarg->sval, (String)gd->domestic_und);
				if (i < 0) 
					return serror("Too many underlyings used in Grfn: cannot list them");

            }
            else if (gd->pass == 2) 
			{
                /* Set index of underlying */
                index = get_grf_und_index(lastarg->sval);
                args = comll_gototop(args);
                args->type = gs.type;
                if (index >= 0) 
					args->ivec[0] = index;
                comll_cut(args);
                return NULL;
            }

            break;
// IMPLEMENTATION OF PVR IN GRFN
		case CLLNSourceAssign:

			new = comll_insert_after(lastarg);
            new->type = gs.type;
            new->ivec[0] = gd->J;
			return NULL;

			break;

        default:
            return serror("%s: %s()",GRERR_UNKNOWN_SYMBOL,name);
    
    } /* end switch */

    return NULL;
                             
} /* end grfn_interp_func */                                          




/******************************************************************************\

    FUNCTION        :   grfn_interp_name 

    CALLED BY       :   yyparse

    DESCRIPTION     :   Interprets a string like "CPN" or "I" or "SEMI".  
                        In this case, c points to a blank COMLL_PTR allocated 
                        from within the yacc code.

\******************************************************************************/


Err grfn_interp_name    (   COMLL_PTR   c, 
                            String      name, 
                            GrfnDeal    *gd     )

{
Err         error;
GrfnSymbol  gs;

/* Check to see if name is in the GRFN symbol table.
   If so then return information relating to the symbol. */
	error = grfn_interp_comll_name(name, &gs);
	if (error)
		return error;

/* Look at the nametype of the symbol and take appropriate action */
    switch (gs.nametype)
	{

        case CLLNVar:

            c->gdcells = &gd->ssparam[gs.ival];
            c->type = COMLL_VAR;
            GrfnAddStatus(GrfnCGCell(gd).status, GRFNCSVARREF);
            break;

        case CLLNI:         

            c->type = COMLL_REAL;
            c->dval = (double)gd->I;
            break;          

        case CLLNJ:

            c->type = COMLL_REAL;
            c->dval = (double)gd->J;
            break;

        case CLLNNow:

            c->type = COMLL_REAL;
            c->dval = (double)gd->nowdt;
            break;

        case CLLNNxt:

            c->type = COMLL_REAL;
            if (gd->I > gd->sslength-2) return serror("Invalid NXT keyword");

            if (grfn_lang_cur_amcell)
                c->dval = (double)gd->nxtstpdt;
            else 
                c->dval= (double)gd->nxtevdt;

            break;

        case CLLNPrv:

            if (gd->I <= 0) return serror("Invalid PRV keyword");
            c->type = COMLL_REAL;
            c->dval = (double)gd->prvdt;
            break;

        case CLLNSample:

            c->type = gs.type;
            c->ivec[0] = gs.ival;
            break;

        default: 

            return serror("%s %s",GRERR_UNKNOWN_SYMBOL,name); 

    } /* END switch */

    return NULL;

} /* END grfn_interp_name */

/* -------------------------------------------------------------------------------- */

/* Resets a status to its default values */
void GrfnResetStatus(GrfnCellStatus status)
{
int i;

	for (i =0;i <GRFNCSLASTSTATUS;i++)
		status[i] = 0;
}

/* -------------------------------------------------------------------------------------------- */
/* Aggregates into a full status (row, tableau...) a single status (for the cell) */
void GrfnAggregateStatus(GrfnCellStatus full_status,GrfnCellStatus cell_status)  
{
int i;

	for (i =0;i <GRFNCSLASTSTATUS;i++)
		full_status[i] = ( full_status[i] || cell_status[i] );
}

/* -------------------------------------------------------------------------------------------- */

/* Copy into a new status an old status  */
void GrfnCopyStatus(GrfnCellStatus new_status,GrfnCellStatus old_status)  
{
int i;

	for (i =0;i <GRFNCSLASTSTATUS;i++)
		new_status[i] = old_status[i];
}
/* -------------------------------------------------------------------------------------------- */