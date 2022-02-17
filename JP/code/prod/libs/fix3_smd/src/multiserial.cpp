/////////////////////////////////////////////////////////////////////
// Multiserial option
//
// multiserial.cpp
//////////////////////////////////////////////////////////////////////

#ifdef _MSC_VER
// For some tiresome reason you get shed loads of warnings about 
// inline functions from math.h being removed if you don't put this 
// disabler in. 
#pragma warning(disable : 4514)
#endif


#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "fix123head.h"
#include "fix123plus.h"



/*****  Multiserial_t  ****************************************************/
/*
*   Calculates the zero binary leg price, no dev, just acretes notional.
*	Note that memory is allocated here for the state variable vector
*	pStateValue.  The responsibility for freeing this at the end rests
*	with the caller.
*/
extern "C" int Fix3_Multiserial_t(
	double		**MultiSerial,	/* (O) Prices for all states			*/
	int		iResetNumber,	/* (I) Reset 0, 1 or 2                              */
	double		*ResetSlice,	/* (I) Slice containing rate at this reset          */
	double		DCF,		/* (I) Day count fraction for accrual               */
	double		dNotional,      /* (I) Notional for this reset                      */
	char		Compounding,	/* (I) Compounding is 'S'imple or 'C'omplex         */
	int		NbStates,	/* (I) Nb of notional states				*/
	double		**pStateValue,	/* (I) Values of the state variable			*/
	double		MaxState,	/* (I) Upper notional state bound, this timestep	*/
	double		MinState,	/* (I) Lower notional state bound, this timestep	*/
	int		 t,		/* (I) Current time point				*/
	FIX3_TREE_DATA	*tree_data)	/* (I) Tree data structure				*/
{
    /* Local slice pointers */
    double				*MultiSerialL[MAXNBSTATES+1];
    double				*ResetL;

	/* Local notional pointers */
	double				*StatesTPlus1 = NULL;
	double				*StatesCurrent = NULL;

    /* state-variable variables */
    double				D[MAXNBSTATES][3];  /* Precomputed quadratic coeffs    */
    int					s;                  /* State variable index            */

    /* payoff variables */
    double				OptionPayoff[MAXNBSTATES];   /* Zero binary payoff due to future ResetSlices */
    int					bSimpleCompound;	/* Compounding type */
	double				LiborArrears;		/* Value of payoff due to arrears libor */
	double				LiborPairValue;		/* Value of payoff due to libor pair */
	double				LiborPairKnownPart;	/* Value of payoff due to known part of libor pair */
	double				LiborPairFirstPart;	/* Value of payoff due to first part of libor pair */

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
        DR_Error((char*)"Multiserial_t: nb of state variables exceeds limit (200)!");
        goto RETURN;
    }

	/* What kind of compounding? */
    bSimpleCompound = (Compounding == 'S');

    /* ------------------------------------------------ */
    /*   Prepare state variables						*/
    /* ------------------------------------------------ */

	/* Store a copy of notionals from t+1 because we'll need them. */
	StatesTPlus1 = *pStateValue;

	/* Generate the notionals for this timestep */
	if ((rc = Fix3_StateVar_Generate(
					pStateValue,		/* (I/O) Notionals at current zero binary reset date	*/
					NbStates,			/* (I) Nb of notional states							*/
					MaxState,			/* (I) Upper notional state bound, this timestep		*/
					MinState,			/* (I) Lower notional state bound, this timestep		*/
					'G'					/* (I) Geometric generation of notionals				*/
					)) != SUCCESS)
	{
		DR_Error((char*)"Multiserial_t: Failed to generate the state variables!");
		goto RETURN;
	}

	/* Dereference */
	StatesCurrent = *pStateValue;

    /* ------------------------------------------------ */
    /*   Prepare prev states for interp if appropriate  */
    /* ------------------------------------------------ */
    if (!IS_EQUAL(StatesTPlus1[0], StatesTPlus1[NbStates-1]))
    {
        for (s = 0; s < NbStates-2; s++)
        {
            D[s][0] = 1. / ((StatesTPlus1[s] - StatesTPlus1[s+1])
                           *(StatesTPlus1[s] - StatesTPlus1[s+2]));
            D[s][1] = 1. / ((StatesTPlus1[s+1] - StatesTPlus1[s])
                           *(StatesTPlus1[s+1] - StatesTPlus1[s+2]));
            D[s][2] = 1. / ((StatesTPlus1[s+2] - StatesTPlus1[s])
                           *(StatesTPlus1[s+2] - StatesTPlus1[s+1]));
        }
    }

	/* What reset? */
	if (iResetNumber == 2)
	{
		/* ------------------------------------------------ */
		/* Final reset date									*/
		/* ------------------------------------------------ */

        /************************  1 FACTOR    ****************************/
        if (tree_data->NbFactor == 1)
        {
			/* Index into tree */
            offset = Fix3_Node_Offset(1, 0, 0, t, tree_data);

			/* Find index into vector of subsequent (t+1) values. */
            for (s = 0; s < NbStates ; s++) MultiSerialL[s] = MultiSerial[s] + offset;
            ResetL  = ResetSlice  + offset;

			/* Loop through space dimension */
            for (i = Bottom1; i <= Top1; i ++)
            {   
				/* In this state (i), we want the value of the arrears libor */
				LiborArrears = ACC_FN(ResetL[i], DCF, bSimpleCompound) * dNotional;

                /* Iterate through states */
                for (s = 0; s < NbStates; s++)
                {
					/* What is the value of the libor pair in this state? */
                    LiborPairValue = StatesCurrent[s];

					/* Work backwards to the payoff */
					OptionPayoff[s] = MAX(LiborArrears, LiborPairValue);
                } /* for state idx */

                /* Iterate through states */
                for (s = 0; s < NbStates; s++)
				{
					/* replace sticky prices with new values for each state */
					/* and accrete notional                                 */
					MultiSerialL[s][i] = OptionPayoff[s];
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
                offset = Fix3_Node_Offset(2, i, 0, t, tree_data);

				/* Find index into vector of subsequent (t+1) values. */
                for (s = 0; s < NbStates ; s++) MultiSerialL[s] = MultiSerial[s] + offset;
                ResetL  = ResetSlice  + offset;

				/* Loop through space dimension #2 */
                for (j = Bottom2[i]; j <= Top2[i]; j++)
                {
					/* In this state (i), we want the value of the arrears libor */
					LiborArrears = ACC_FN(ResetL[j], DCF, bSimpleCompound) * dNotional;

	                /* Iterate through states */
                    for (s = 0; s < NbStates; s++)
                    {
						/* What is the value of the libor pair in this state? */
	                    LiborPairValue = StatesCurrent[s];

						/* Work backwards to the payoff */
						OptionPayoff[s] = MAX(LiborArrears, LiborPairValue);
                    } /* for state idx */

					/* Iterate through states */
					for (s = 0; s < NbStates; s++)
					{
						/* replace sticky prices with new values for each state */
						/* and accrete notional                                 */
						MultiSerialL[s][j] = OptionPayoff[s];
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
                    offset = Fix3_Node_Offset(3, i, j, t, tree_data);

					/* Find index into vector of subsequent (t+1) values. */
                    for (s = 0; s < NbStates ; s++) MultiSerialL[s] = MultiSerial[s] + offset;
                    ResetL  = ResetSlice  + offset;

					/* Loop through space dimension #3 */
                    for (k = Bottom3[i][j]; k <= Top3[i][j]; k++)
                    {
						/* In this state (i), we want the value of the arrears libor */
						LiborArrears = ACC_FN(ResetL[k], DCF, bSimpleCompound) * dNotional;

		                /* Iterate through states */
                        for (s = 0; s < NbStates; s++)
                        {
							/* What is the value of the libor pair in this state? */
		                    LiborPairValue = StatesCurrent[s];

							/* Work backwards to the payoff */
							OptionPayoff[s] = MAX(LiborArrears, LiborPairValue);
                        } /* for state idx */

						/* Iterate through states */
						for (s = 0; s < NbStates; s++)
						{
							/* replace sticky prices with new values for each state */
							/* and accrete notional                                 */
							MultiSerialL[s][k] = OptionPayoff[s];
						} /* for state idx */
                    }  /* for k */
                }  /* for j */
        }  /* if then else */
	} /* if final reset */
	else if (iResetNumber == 1)
	{
		/* ------------------------------------------------ */
		/* Penultimate reset date							*/
	    /* ------------------------------------------------ */

        /************************  1 FACTOR    ****************************/
        if (tree_data->NbFactor == 1)
        {
			/* Index into tree */
            offset = Fix3_Node_Offset(1, 0, 0, t, tree_data);

			/* Find index into vector of subsequent (t+1) values. */
            for (s = 0; s < NbStates ; s++) MultiSerialL[s] = MultiSerial[s] + offset;
            ResetL  = ResetSlice  + offset;

			/* Loop through space dimension */
            for (i = Bottom1; i <= Top1; i ++)
            {   
				/* In this state (i), we want the value of the second libor */
				LiborPairKnownPart = ACC_FN(ResetL[i], DCF, bSimpleCompound) * dNotional;

                /* Iterate through states */
                for (s = 0; s < NbStates; s++)
                {
					/* What is the value of the libor pair in this state? */
					LiborPairValue = StatesCurrent[s] + LiborPairKnownPart;

					/* Work backwards to the payoff */
					OptionPayoff[s] = Fix3_StateVar_Interp(i, NbStates, StatesTPlus1, MultiSerialL, D, LiborPairValue);
                } /* for state idx */

                /* Iterate through states */
                for (s = 0; s < NbStates; s++)
				{
					/* replace sticky prices with new values for each state */
					/* and accrete notional                                 */
					MultiSerialL[s][i] = OptionPayoff[s];
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
                offset = Fix3_Node_Offset(2, i, 0, t, tree_data);

				/* Find index into vector of subsequent (t+1) values. */
                for (s = 0; s < NbStates ; s++) MultiSerialL[s] = MultiSerial[s] + offset;
                ResetL  = ResetSlice  + offset;

				/* Loop through space dimension #2 */
                for (j = Bottom2[i]; j <= Top2[i]; j++)
                {
					/* In this state (i, j), we want the value of the second libor */
					LiborPairKnownPart = ACC_FN(ResetL[j], DCF, bSimpleCompound) * dNotional;

					/* Iterate through states */
					for (s = 0; s < NbStates; s++)
					{
						/* What is the value of the libor pair in this state? */
						LiborPairValue = StatesCurrent[s] + LiborPairKnownPart;

						/* Work backwards to the payoff */
						OptionPayoff[s] = Fix3_StateVar_Interp(j, NbStates, StatesTPlus1, MultiSerialL, D, LiborPairValue);
					} /* for state idx */

					/* Iterate through states */
					for (s = 0; s < NbStates; s++)
					{
						/* replace sticky prices with new values for each state */
						/* and accrete notional                                 */
						MultiSerialL[s][j] = OptionPayoff[s];
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
                    offset = Fix3_Node_Offset(3, i, j, t, tree_data);

					/* Find index into vector of subsequent (t+1) values. */
                    for (s = 0; s < NbStates ; s++) MultiSerialL[s] = MultiSerial[s] + offset;
                    ResetL  = ResetSlice  + offset;

					/* Loop through space dimension #3 */
                    for (k = Bottom3[i][j]; k <= Top3[i][j]; k++)
                    {
						/* In this state (i, j, k), we want the value of the second libor */
						LiborPairKnownPart = ACC_FN(ResetL[k], DCF, bSimpleCompound) * dNotional;

						/* Iterate through states */
						for (s = 0; s < NbStates; s++)
						{
							/* What is the value of the libor pair in this state? */
							LiborPairValue = StatesCurrent[s] + LiborPairKnownPart;

							/* Work backwards to the payoff */
							OptionPayoff[s] = Fix3_StateVar_Interp(k, NbStates, StatesTPlus1, MultiSerialL, D, LiborPairValue);
						} /* for state idx */

						/* Iterate through states */
						for (s = 0; s < NbStates; s++)
						{
							/* replace sticky prices with new values for each state */
							/* and accrete notional                                 */
							MultiSerialL[s][k] = OptionPayoff[s];
						} /* for state idx */
                    }  /* for k */
                }  /* for j */
        }  /* if then else */
	}
	else if (iResetNumber == 0)
	{
		/* ------------------------------------------------ */
		/* First reset date									*/
	    /* ------------------------------------------------ */

        /************************  1 FACTOR    ****************************/
        if (tree_data->NbFactor == 1)
        {
			/* Index into tree */
            offset = Fix3_Node_Offset(1, 0, 0, t, tree_data);

			/* Find index into vector of subsequent (t+1) values. */
            for (s = 0; s < NbStates ; s++) MultiSerialL[s] = MultiSerial[s] + offset;
            ResetL  = ResetSlice  + offset;

			/* Loop through space dimension */
            for (i = Bottom1; i <= Top1; i ++)
            {   
				/* In this state (i), we want the value of the first libor */
				LiborPairFirstPart = ACC_FN(ResetL[i], DCF, bSimpleCompound) * dNotional;

				/* Given the valie of the first libor we now know the value */
				/* of the whole option.  We stick this into the slice for   */
				/* the first state variable.                                */
				MultiSerialL[0][i] = Fix3_StateVar_Interp(i, NbStates, StatesTPlus1, MultiSerialL, D, LiborPairFirstPart);
            } /* for i */

        } /* if NbFactor == 1 */
        /************************  2 FACTORS   ****************************/
        else if (tree_data->NbFactor == 2)
        {
			/* Loop through space dimension #1*/
            for (i = Bottom1; i <= Top1; i ++)
            {
				/* Index into tree */
                offset = Fix3_Node_Offset(2, i, 0, t, tree_data);

				/* Find index into vector of subsequent (t+1) values. */
                for (s = 0; s < NbStates ; s++) MultiSerialL[s] = MultiSerial[s] + offset;
                ResetL  = ResetSlice  + offset;

				/* Loop through space dimension #2 */
                for (j = Bottom2[i]; j <= Top2[i]; j++)
                {
					/* In this state (i, j), we want the value of the first libor */
					LiborPairFirstPart = ACC_FN(ResetL[j], DCF, bSimpleCompound) * dNotional;

					/* Given the valie of the first libor we now know the value */
					/* of the whole option.  We stick this into the slice for   */
					/* the first state variable.                                */
					MultiSerialL[0][j] = Fix3_StateVar_Interp(j, NbStates, StatesTPlus1, MultiSerialL, D, LiborPairFirstPart);
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
                    offset = Fix3_Node_Offset(3, i, j, t, tree_data);

					/* Find index into vector of subsequent (t+1) values. */
                    for (s = 0; s < NbStates ; s++) MultiSerialL[s] = MultiSerial[s] + offset;
                    ResetL  = ResetSlice  + offset;

					/* Loop through space dimension #3 */
                    for (k = Bottom3[i][j]; k <= Top3[i][j]; k++)
                    {
						/* In this state (i, j), we want the value of the first libor */
						LiborPairFirstPart = ACC_FN(ResetL[k], DCF, bSimpleCompound) * dNotional;

						/* Given the valie of the first libor we now know the value */
						/* of the whole option.  We stick this into the slice for   */
						/* the first state variable.                                */
						MultiSerialL[0][j] = Fix3_StateVar_Interp(k, NbStates, StatesTPlus1, MultiSerialL, D, LiborPairFirstPart);
                    }  /* for k */
                }  /* for j */
        }  /* if then else */
	}
	else
	{
		/* ------------------------------------------------ */
		/* Not a reset date                                 */
	    /* ------------------------------------------------ */
        DR_Error((char*)"Multiserial_t: Reset should be first, second or third - 0, 1 or 2!");
        goto RETURN;
	}

	/* If we get this far then it was okay! */
    status = SUCCESS;
    
RETURN:

	/* Bin the old t+1 state information */
	Free_DR_Array(StatesTPlus1, DOUBLE, 0, NbStates-1);
	StatesTPlus1=NULL;

	/* Return the result code */
    return (status);
}  /* Fix3_Multiserial_t */
