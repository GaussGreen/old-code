/************************************************************************
 * Module:	PenGuin
 * File:
 * Function:    
 * Author:      D. Liu
 ************************************************************************/
/*
#include "vtlkobinary.h"
*/

#include "kvtspar.h"	// Slice parser
#include "kutilios.h"


extern "C" {
#include "drltime.h"
#include "drlmath.h"
};



const	int	nObsMax=120;


//---------------------------------------------------------------

void KVTreeKoBinary(
	KVTree          &vt,
	KRate           &obsRate,
	const KVector(TDate)  &obsDates,
	const KVector(double) &loBarriers,
	const KVector(double) &hiBarriers,
	const KVector(TDate)  &binDates,
	const String          &curveName,
	KVector(double)       &prob)
{
static	char	routine[] = "KVTreeKoBinary";

	int	nObs, idxObs;
	int	nBin, idxBin;


	KVector(KRateReset*)	mObsRateResets;
	KVector(TDate)		    mModelResetDates;
	KVector(KTSlice*)	    mValue;


	TDate			currDate;
	TDate			todayDate = vt.TPToday();

	KTSlice			*mRateTs = NULL;


	KVector(KTSlice*) mProb;
	KTSliceNode       tsNode;         // Slice node

    double            idxRate,        // Rate at each node
					  idxStep,        // Step between nodes
					  lowerStep,      // Step function at low barrier
					  upperStep,      // Step function at hight barrier
					  rangeInStep;    // double-barrier step function


#define	DEV(sliceP)	{if (sliceP != NULL) {sliceP->Dev(curveName);}}


    try {


	ASSERT_OR_THROW(obsDates.size() == loBarriers.size());
	ASSERT_OR_THROW(obsDates.size() == hiBarriers.size());
	nObs = obsDates.size();

	ASSERT_OR_THROW(nObs < nObsMax);


	// enough space to store the model reset dates and resets.
	mObsRateResets.resize(nObs);
	mModelResetDates.resize(nObs);

    nBin = binDates.size();

	for (idxBin=0; idxBin<nBin; idxBin++) 
	{
		mProb[idxBin] = NULL;
	}


	//
	// Add events in tree
	//
	for (idxObs=0; idxObs<nObs; idxObs++) 
	{
		//
		// Add events to tree.
		// Rate reset: we collect the model reset date
		// (i.e. where the rate is known in the tree)
		//
		//
		if (obsDates[idxObs] <  todayDate)
			throw KFailure("%s: obsDate %d %s < todayDate %s.\n",
				routine, idxObs+1,
				DrlTDatePrint(NULL, obsDates[idxObs]),
				DrlTDatePrint(NULL, todayDate));

		if (loBarriers[idxObs] >= hiBarriers[idxObs])
			throw KFailure("%s: loBarriers[%d] %lf < hiBarriers %lf.\n",
				routine, idxObs+1,
				loBarriers[idxObs],
				hiBarriers[idxObs]);


		mObsRateResets[idxObs] = new KRateReset(
			curveName,
			obsDates[idxObs],
			obsDates[idxObs], // !!! effective dates!
			obsRate);

		mModelResetDates[idxObs] = vt.Insert(*mObsRateResets[idxObs]);

	}

    // Binary dates
	for (idxBin=0; idxBin<nBin; idxBin++) 
	{
        vt.Insert(binDates[idxBin]);
	}




	//
	// Initialize tree timeline and calibrate tree
	//
	vt.SetUpTimeline();
	vt.Calibrate();

	//
	// Rollback
	//
	for (int tpIdx=vt.TPNum(); tpIdx >= 0; tpIdx--) 
	{
	    //
	    // Update Tree
	    //
	    vt.Update(tpIdx);


	    currDate = vt.TPDateCurrent();

	    //
	    // Dev
	    //
	    for (idxBin=0; idxBin<nBin; idxBin++) 
		{
		    DEV(mProb[idxBin]);

			// Allocate each slice and initialize to $1 on binary date
			if (currDate == binDates[idxBin])
			{
                mProb[idxBin] = new KTSlice(vt, "BinaryProb", curveName);
				*mProb[idxBin] = 1e0;
			}
				
	    }


	    //
	    // if observation date: the observation rate can be computed.
	    //
	    for (idxObs=0; idxObs<nObs; idxObs++) 
		{
            if (currDate == mModelResetDates[idxObs]) 
			{

		        //
		        // Ensure slices allocated
		        //
		        if (mRateTs == NULL)
		        	mRateTs = new KTSlice(vt, "mObsRateTs", curveName);

		        // Get rate reset
		        //
		        vt.Get(*mRateTs, *mObsRateResets[idxObs]);

                // Perform knock out calculation at each slice node
				// with single smoothing method
                for (tsNode.begin(*mRateTs); !tsNode.end(); ++tsNode)
				{
                    idxRate = (*mRateTs)[tsNode];
                    idxStep = mRateTs->GetNodeStepMax(tsNode);

			        ////////////////////////////////////////////
			        // Given two step functions:
			        // 
			        //               f1
			        // 
			        //            ------------------------
			        //            |         
			        //            |         
			        //    --------         
			        //
    			    //
    			    //               f2
    			    // 
    			    //                      --------------
    			    //                      |         
    			    //                      |         
    			    //    ------------------         
    			    //
    			    //
    			    // Knock-out range payoff function
    			    //
    			    //             (1-f2)*f1
    			    //
    			    //            -----------
    			    //            |         |
    			    //            |         |
    			    //    --------          -------------- 
    			    //
     			    //////////////////////////////////////////////

					lowerStep = DrlSmoothStepFcn(
                                    idxRate - loBarriers[idxObs],
									idxStep);
                    upperStep = DrlSmoothStepFcn(
								    idxRate - hiBarriers[idxObs],
								    idxStep);
                    rangeInStep = (1. - upperStep) * lowerStep;

					// Apply knock-out to each binary payoff
                    for (idxBin=0; idxBin<nBin; idxBin++)
					{
                        if (currDate <= binDates[idxBin])
							(*mProb[idxBin])[tsNode] *= rangeInStep;
					}

				}    // tsNode loop
                
            }    // if observation date

        }    // Observation date loop



	    //
	    // Last TP: store the result.
	    //
	    if (vt.TPIdxCurrent() == 0) 
		{
            for (idxBin=0; idxBin<nBin; idxBin++) 
			{
				prob[idxBin] = mProb[idxBin]->GetCenter();
			}
        }


	}    // t loop



	//
	// Free memory
	//
	delete mRateTs;
	for (idxObs=0; idxObs<nObs; idxObs++) 
		delete mObsRateResets[idxObs];

	for (idxBin=0; idxBin<nBin; idxBin++) 
		delete mProb[idxBin];



    }
    catch (KFailure) {
	    delete mRateTs;
	    for (idxObs=0; idxObs<nObs; idxObs++) 
	    	delete mObsRateResets[idxObs];

	    for (idxBin=0; idxBin<nBin; idxBin++) 
	    	delete mProb[idxBin];

    	throw KFailure("%s: failed.\n", routine);
    }
}

