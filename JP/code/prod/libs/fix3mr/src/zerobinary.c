/****************************************************************************/
/*      Zero binary payoff													*/
/****************************************************************************/
/*      STATEVAR.c                                                          */
/****************************************************************************/




#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "fix123head.h"



/*****  ZeroBinary_t  ****************************************************/
/*
*   Calculates the zero binary leg price, no dev, just acretes notional.
*	Note that memory is allocated here for the state variable vector
*	pNotlsCurrent.  The responsibility for freeing this at the end rests
*	with the caller.
*/
int ZeroBinary_t(
	double				**ZeroBinary,		/* (O) Prices for all states						*/
	long				ResetFlagSt,		/* (I) Reset at this timestep flag					*/
	double				*Step,				/* (I) Step for ladder								*/
	double				DcfSt,				/* (I) Day count fract								*/
	char				CompSt,				/* (I) Simple/Compound pmt							*/
	int					NbStates,			/* (I) Nb of notional states						*/
	double				**pNotlsCurrent,	/* (I) Notionals at current zero binary reset date	*/
	double				MaxState,			/* (I) Upper notional state bound, this timestep	*/
	double				MinState,			/* (I) Lower notional state bound, this timestep	*/
	int					t,					/* (I) Current time point							*/
	TREE_DATA			*tree_data)			/* (I) Tree data structure							*/
{
    /* Local slice pointers */
    double				*ZeroBinaryL[MAXNBSTATES+1];
    double				*StepL;

	/* Local notional pointers */
	double				*NotlsTPlus1 = NULL;
	double				*NotlsCurrent = NULL;

    /* state-variable variables */
    double				D[MAXNBSTATES][3];  /* Precomputed quadratic coeffs    */
    int					s;                  /* State variable index            */

    /* payoff variables */
	double				NewNotional;		/* Subsequent notional for this state */
    double				ZBPayoff[MAXNBSTATES];   /* Zero binary payoff due to future steps */
    double				Accretion;
    int					isSimpleSt;

    /* Tree variables */
    int					Top1, Bottom1;      /* Tree limits (1rst dim) */
    int					*Top2, *Bottom2;    /* Tree limits (2nd dim)  */
    int					**Top3, **Bottom3;  /* Tree limits (3rd dim)  */

    int					i, j, k;            /* Node indices           */
	int					rc;					/* Result code			  */
    int					offset;             /* Node offset            */
    int					status = FAILURE;   /* Error status           */

	/* Find the upper and lower array bounds. */
    Top1    = tree_data->Top1[t];
    Top2    = tree_data->Top2[t];
    Top3    = tree_data->Top3[t];
    Bottom1 = tree_data->Bottom1[t];
    Bottom2 = tree_data->Bottom2[t];
    Bottom3 = tree_data->Bottom3[t];

	/* No good if we have too many states. */
    if (NbStates > MAXNBSTATES)
    {
        DR_Error("ZeroBinary_t: nb of state variables exceeds limit (200)!");
        goto RETURN;
    }

	/* What kind of compounding? */
    isSimpleSt = (CompSt == 'S');

	/* ZERO BINARY LEG */

	/****************************************************************/
	/* Reset date													*/
	/*																*/
    /* Work backwards through tree applying the rule below.			*/
	/*		ZB[s,i,t] = ZB[s',i,t+1]								*/
	/*		where N(s) = N(s') * (1 + Rate(t,i))					*/
	/****************************************************************/
    if (ResetFlagSt)
    {
        /* ------------------------------------------------ */
        /*   Prepare notionals                              */
        /* ------------------------------------------------ */

		/* Store a copy of notionals from t+1 because we'll need them. */
		NotlsTPlus1 = *pNotlsCurrent;

		/* Generate the notionals for this timestep */
		if ((rc = StateVar_Generate(
						pNotlsCurrent,		/* (I/O) Notionals at current zero binary reset date	*/
						NbStates,			/* (I) Nb of notional states							*/
						MaxState,			/* (I) Upper notional state bound, this timestep		*/
						MinState,			/* (I) Lower notional state bound, this timestep		*/
						'G'					/* (I) Geometric generation of notionals				*/
						)) != SUCCESS)
	    {
			DR_Error("ZeroBinary_t: nb of state variables exceeds limit (200)!");
			goto RETURN;
		}

		/* Dereference */
		NotlsCurrent = *pNotlsCurrent;

        /* ------------------------------------------------ */
        /*   Prepare prev states for interp if appropriate  */
        /* ------------------------------------------------ */

        if (!IS_EQUAL(NotlsTPlus1[0], NotlsTPlus1[NbStates-1]))
        {
            for (s = 0; s < NbStates-2; s++)
            {
                D[s][0] = 1. / ((NotlsTPlus1[s] - NotlsTPlus1[s+1])
                               *(NotlsTPlus1[s] - NotlsTPlus1[s+2]));
                D[s][1] = 1. / ((NotlsTPlus1[s+1] - NotlsTPlus1[s])
                               *(NotlsTPlus1[s+1] - NotlsTPlus1[s+2]));
                D[s][2] = 1. / ((NotlsTPlus1[s+2] - NotlsTPlus1[s])
                               *(NotlsTPlus1[s+2] - NotlsTPlus1[s+1]));
            }
        }


        /************************  1 FACTOR    ****************************/
        if (tree_data->NbFactor == 1)
        {
			/* Index into tree */
            offset = Node_Offset(1, 0, 0, t, tree_data);

			/* Find index into vector of subsequent (t+1) values. */
            for (s = 0; s < NbStates ; s++) ZeroBinaryL[s] = ZeroBinary[s] + offset;
            StepL  = Step  + offset;

			/* Loop through space dimension */
            for (i = Bottom1; i <= Top1; i ++)
            {                
                /* In this state (i), we accrete and payout a value controlled by */
				/* the variable StepL[i] which is given by a call to ladder_t     */
				/* (externally to this function)                                  */
                Accretion = 1. + ACC_FN(StepL[i], DcfSt, isSimpleSt);

                /* Iterate through states */
                for (s = 0; s < NbStates; s++)
                {
					/* Compute the notional that would we get to at the next time */
					/* step (t+1) we went from the current notional in this state */
					/* and accreted at the given step value. */
                    NewNotional = NotlsCurrent[s] * Accretion;

					/* Work backwards to the payoff */
					ZBPayoff[s] = StateVar_Interp(i, NbStates, NotlsTPlus1, ZeroBinaryL, D, NewNotional);
                } /* for state idx */

                /* Iterate through states */
                for (s = 0; s < NbStates; s++)
				{
					/* replace sticky prices with new values for each state */
					/* and accrete notional                                 */
					ZeroBinaryL[s][i] = ZBPayoff[s];
                } /* for state idx */
            } /* for i */

        } /* if NbFactor == 1 */
        /************************  2 FACTORS   ****************************/
        else if (tree_data->NbFactor == 2)
        {
			/* Loop through space dimension #1*/
            for (i = Bottom1; i <= Top1; i ++)
            {
				/* Index into tree */
                offset = Node_Offset(2, i, 0, t, tree_data);

				/* Find index into vector of subsequent (t+1) values. */
                for (s = 0; s < NbStates ; s++) ZeroBinaryL[s] = ZeroBinary[s] + offset;
                StepL  = Step  + offset;

				/* Loop through space dimension #2 */
                for (j = Bottom2[i]; j <= Top2[i]; j++)
                {
					/* In this state (i, j), we accrete and payout a value controlled */
					/* by the variable StepL[j] which is given by a call to ladder_t  */
					/* (externally to this function)                                  */
                    Accretion = 1. + ACC_FN(StepL[j], DcfSt, isSimpleSt);

	                /* Iterate through states */
                    for (s = 0; s < NbStates; s++)
                    {
						/* Compute the notional that would we get to at the next time */
						/* step (t+1) we went from the current notional in this state */
						/* and accreted at the given step value. */
                        NewNotional = NotlsCurrent[s] * Accretion;

						/* Work backwards to the payoff */
						ZBPayoff[s] = StateVar_Interp(j, NbStates, NotlsTPlus1, ZeroBinaryL, D, NewNotional);
                    } /* for state idx */

					/* Iterate through states */
					for (s = 0; s < NbStates; s++)
					{
						/* replace sticky prices with new values for each state */
						/* and accrete notional                                 */
						ZeroBinaryL[s][j] = ZBPayoff[s];
					} /* for state idx */
                }  /* for j */
            }  /* for i */
        }
        /************************  3 FACTORS   ****************************/
        else if (tree_data->NbFactor == 3)
        {
			/* Loop through space dimension #1 */
            for (i = Bottom1; i <= Top1; i ++)
				/* Loop through space dimension #2 */
                for (j = Bottom2[i]; j <= Top2[i]; j++)
                {
					/* Index into tree */
                    offset = Node_Offset(3, i, j, t, tree_data);

					/* Find index into vector of subsequent (t+1) values. */
                    for (s = 0; s < NbStates ; s++) ZeroBinaryL[s] = ZeroBinary[s] + offset;
                    StepL  = Step  + offset;

					/* Loop through space dimension #3 */
                    for (k = Bottom3[i][j]; k <= Top3[i][j]; k++)
                    {
						/* In this state (i, j, k), we accrete and payout a value controlled */
						/* by the variable StepL[k] which is given by a call to ladder_t     */
						/* (externally to this function)                                     */
                        Accretion = 1. + ACC_FN(StepL[k], DcfSt, isSimpleSt);

		                /* Iterate through states */
                        for (s = 0; s < NbStates; s++)
                        {
							/* Compute the notional that would we get to at the next time */
							/* step (t+1) we went from the current notional in this state */
							/* and accreted at the given step value. */
                            NewNotional = NotlsCurrent[s] * Accretion;

							/* Work backwards to the payoff */
							ZBPayoff[s] = StateVar_Interp(k, NbStates, NotlsTPlus1, ZeroBinaryL, D, NewNotional);
                        } /* for state idx */

						/* Iterate through states */
						for (s = 0; s < NbStates; s++)
						{
							/* replace sticky prices with new values for each state */
							/* and accrete notional                                 */
							ZeroBinaryL[s][k] = ZBPayoff[s];
						} /* for state idx */
                    }  /* for k */
                }  /* for j */
        }  /* if then else */
    }  /* if ResetFlagZeroBinary */

    status = SUCCESS;
    
RETURN:

	/* Bin the old t+1 state information */
	if (NotlsTPlus1 != NULL)
	{
		if (Free_DR_Array (NotlsTPlus1, DOUBLE, 0, NbStates-1) == FAILURE)
			goto RETURN;
	}

	/* Return the result code */
    return (status);
}  /* ZeroBinary_t */


