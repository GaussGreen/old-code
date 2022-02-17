//////////////////////////////////////////////////////////////////////
// Enhanced utilities for maintaining zero bank structures.
//
// zerobank_plus.cpp (based on code from zerobank.c)
//////////////////////////////////////////////////////////////////////


/* For some tiresome reason you get shed loads of warnings about 
   inline functions from math.h being removed if you don't put this 
   disabler in. */
// #pragma warning(disable : 4514)


#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "fix123head.h"


#define CASH_RATE 'C'
#define SWAP_RATE 'S'

extern "C" int Fix3_ZbkDLFromIdx_Plus(
	long				ValueDate,			// (I) Value date				
	int					NbResetDates,		// (I) size of reset date list
	long				*ResetDL,			// (I) reset date list
	long				*SwapStDL,			// (I) Rate start date list		
	char				IndexType,			// (I) Rate type (cash or swap)
	int					IndexMat,			// (I) index tenor in months
	char				IndexFreq,			// (I) index freq
	int					*NbMatDates,		// (O) Nb of zero mat dates
	long				**MatDL,			// (O) zero mats
	int					*NbUseDates,		// (O) Nb zero use dates
	long				**UseDL)			// (O) zero use dates
//////////////////////////////////////////////////////////////////////
// Given a list of reset dates and descriptions of an index generate 
// a list of zero maturity dates for index estimation and their
// corresponding usage dates.
//
// Returns SUCCESS or FAILURE
{
    int					status = FAILURE;   // Error status
    int					i, j;
    long				LastZeroMat;
    long				*tmpDL = NULL;		// tmp datelist storage
    int					NbtmpDates = 0;

	// Handle exceptions
	try {
		// Loop through the list of reset dates
		for (i=0; i<NbResetDates; i++)
		{
			// If this reset is in the past, ignore it.
			if (ResetDL[i] < ValueDate)
				continue;

			// Get the dates the right way around
			if (ResetDL[i] > SwapStDL[i])
				throw "Rate resets after it starts.";

			// Get the end date for this rate
			LastZeroMat = Nxtmth(SwapStDL[i], IndexMat, 1L);

			// Are we doing cash or swap rates?
			if (IndexType == CASH_RATE)
			{
				// Add the end date to the maturity list
				if (AddDateToList(NbMatDates, MatDL, LastZeroMat) == FAILURE)
					throw "Failed to add date to list.";

				// Add the start date to the use list
				if (AddDateToList(NbUseDates, UseDL, SwapStDL[i]) == FAILURE)
					throw "Failed to add date to list.";
			}
			else if (IndexType == SWAP_RATE)
			{
				// Generate the coupon dates for the swap
				if (DateListFromFreq(SwapStDL[i],			// (I) Start of date list
									 LastZeroMat,			// (I) End of date list
									 IndexFreq,				// (I) Frequency of list
									 'N',					// (I) Stub location Front/back
									 &NbtmpDates,			// (O) Number of dates in list
									 &tmpDL					// (O) List of dates asc order
									 ) == FAILURE)
					throw "Failed to generate coupon dates.";

				// Add them to the maturity and use lists
				j = (ResetDL[i] == SwapStDL[i]) ? 1 : 0;
				for (; j<NbtmpDates; j++)
				{
					if (AddDateToList(NbMatDates, MatDL, tmpDL[j]) == FAILURE)
						throw "Failed to add date to list.";

					if (AddDateToList(NbUseDates, UseDL, ResetDL[i]) == FAILURE)
						throw "Failed to add date to list.";
				}

				// free tmp storage for next loop
				Free_DR_Array(tmpDL,LONG,0,NbtmpDates-1);
				tmpDL = NULL;
				NbtmpDates = 0;
			}
			else
			{
				// Unknown rate type
				throw "Rate type must be cash or swap.";
			}

		} // for i

		// Successful completion
	    status = SUCCESS;
	}
	catch (char *s)
	{
		// Text exception
        DR_Error((char*) "ZbkDLFromIdx_Plus: %s\n", s);                 
		status = FAILURE;
	}
	catch (...)
	{
		// Unhandled exception
        DR_Error((char*)"ZbkDLFromIdx_Plus: An unknown exception was raised\n.");
		status = FAILURE;
	}

	// Release memory
    Free_DR_Array(tmpDL,LONG,0,NbtmpDates-1);

	// Done
    return status;
}
