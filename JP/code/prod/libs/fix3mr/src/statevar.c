/****************************************************************************/
/*      State variable toolkit routines.                                    */
/****************************************************************************/
/*      STATEVAR.c                                                          */
/****************************************************************************/




#include <stdio.h>
#include <stdlib.h>  
#include <math.h>
#include <string.h>
#include "fix123head.h"



/*****  StateVar_Interp  *******************************************************/
/*                                                                           
 *      Interpolates (quadratically between state variables).  The matrix D 
 *		must be set up in advance to contain the necessary coefficients for
 *		the sqinterp function at each point.
 */
double StateVar_Interp(
	int					i,					/* (I) Rate index */
	int					nbStates,			/* (I) Number of states */
	double				*vStates,			/* (I) Vector of states */
	double				**Payoff,			/* (I) Sub-slice containing payoffs */
	double				D[MAXNBSTATES][3],	/* (I) Precomputed data for interpolation */
	double				dNewState)			/* (I) New state to interp for */
{
	int					sj;					/* Index */
	int					q;					
	double				dInterpPayoff;		/* Interpolated payoff */

	/* Is there more than one state at next step?  We 
	   represent that by setting all state variables equal. */
    if (IS_EQUAL(vStates[1], vStates[0])) 
    {
		/* Just one state level */
        sj = 0;
        dInterpPayoff = Payoff[0][i];
    }
    else
    {
        if (dNewState >= vStates[nbStates - 1])
        {
			/* Required value is more than highest state variable */
            dInterpPayoff = Payoff[nbStates-1][i];
        }
        else if (dNewState <= vStates[0])
        {
			/* Required value is less than least state variable */
            dInterpPayoff = Payoff[0][i];
        }
        else /* in between */
        {
			/* Search the state variable vector for the notional
			   just less than this one. */
			for (sj=1; sj<nbStates; sj++)
			{
				if (vStates[sj] > dNewState) break;
			}
			sj--;

			/* Find index of first notional state which is */
			/* necessary for the quadratic interpolation.  */
            q = MIN(sj, nbStates - 3);
    
			/* Do a quadratic interpolation */
            sqinterp(
				vStates[q],			/* (I) x-value #1 */
				vStates[q+1]		/* (I) x-value #2 */, 
				vStates[q+2]		/* (I) x-value #3 */,
				Payoff[q][i],		/* (I) y-value #1 */
				Payoff[q+1][i],		/* (I) y-value #2 */
				Payoff[q+2][i],		/* (I) y-value #3 */
				D[q][0], 			/* (I) z-value #1 */
				D[q][1], 			/* (I) z-value #2 */
				D[q][2],			/* (I) z-value #3 */
				dNewState,			/* (I) x-coordinate required */
				&dInterpPayoff);	/* (O) interpolated value */
        } /* if sj */
    }

	/* return the payoff */
	return dInterpPayoff;
}


/*****  StateVar_Generate  ****************************************************/
/*
*   Generates states (either 'G'eometrically or 'L'inearly between the given
*	bounds.  Note that the vStates vector is allocated but not freed by this 
*	function.
*/
int StateVar_Generate(
	double				**vStates,			/* (I/O) States at current reset date	*/
	int					NbStates,			/* (I) Nb of states						*/
	double				MaxState,			/* (I) Upper state bound, this timestep	*/
	double				MinState,			/* (I) Lower state bound, this timestep	*/
	char				cMethod)			/* (I) Method to use 'G'eom or 'L'inear */
{
    int					s;                  /* State variable index					*/
    double				deltaS;             /* Increment between consec states		*/
    int					status = FAILURE;   /* Error status							*/

	/* No good if we have too many states. */
    if (NbStates > MAXNBSTATES)
    {
        DR_Error("StateVar_Generate: nb of state variables exceeds limit (200)!");
        goto RETURN;
    }
	
	/* Require at least two states. */
    if (NbStates < 2)
    {
        DR_Error("StateVar_Generate: nb of state variables must be at least two.");
        goto RETURN;
    }
	
	/* The maximum state must be less than the minimum state */
    if (MaxState < MinState)
    {
        DR_Error("StateVar_Generate: max state smaller than min state!");
        goto RETURN;
    }

	/* We have to allocate space for the new state information */
	if (((*vStates) = (double *) DR_Array (DOUBLE, 0, NbStates-1)) == NULL)
    {
        DR_Error("StateVar_Generate: could not allocate memory for State!");
        goto RETURN;
    }

	/* State spacing */
    /* The line below is appropriate for linear spacing of notionals    */
	switch (toupper(cMethod))
	{
	case 'G':	/* Geometric spacing */

				/* We can't have minState = 0 because it causes a division by zero */
				if (fabs(MinState) < TINY)
				{
					DR_Error("StateVar_Generate: Minimum state cannot be zero for geometric interpolation!");
					goto RETURN;
				}

				/* Step for geometric spacing of notionals */
				deltaS = pow(MaxState / MinState, 1.0/(NbStates - 1.0));

				/* Compute state variables */
				for (s=0; s<NbStates; s++)
				{
					/* Geometric spacing */
					(*vStates)[s] = MinState * pow(deltaS, s);
				}
				break;

	case 'L':	/* Linear */

				/* Step for geometric spacing of notionals */
				deltaS = (MaxState - MinState)/(NbStates - 1);

				/* Compute state variables */
				for (s=0; s<NbStates; s++)
				{
					/* Linear spacing */
					(*vStates)[s] = MinState + deltaS * s;
				}
				break;

	default:	/* Unknown */
				DR_Error("StateVar_Generate: Generation must be either 'L'inear or 'G'eometric.");
				status=FAILURE;
				goto RETURN;
				break;
	}

	/* If we get this far then we were successful */
    status = SUCCESS;
    
RETURN:

	/* Return the result code */
    return (status);
}



