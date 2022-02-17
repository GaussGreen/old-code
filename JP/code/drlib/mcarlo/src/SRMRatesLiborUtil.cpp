//----------------------------------------------------------------------------
//
//   Filename    : SRMRatesLiborUtil.cpp (from SRMRatesUtil)
//
//   Description : Helper for BGM model
//
//   Date        : 2 May 2006
//
//
//----------------------------------------------------------------------------

#include "edginc/config.hpp"
#include "edginc/SRMUtil.hpp"
#include "edginc/SRMRatesLiborUtil.hpp"
#include "edginc/SRMSwaption.hpp"
#include "edginc/IRCalib.hpp"
#include "edginc/SwapTool.hpp"
#include "edginc/DayCountConventionFactory.hpp"
#include "edginc/Actual365F.hpp"
#include "edginc/CompoundBasis.hpp"
#include "edginc/CriticalDateCollector.hpp"
#include "edginc/MaturityTimePeriod.hpp"
#include <cassert>

DRLIB_BEGIN_NAMESPACE


///// constructs and populates SRMRatesUtil object
SRMRatesLiborUtil::SRMRatesLiborUtil(
    const DateTime&      baseDate,
    int                  numFactors,
    const string&        modelParamsKey,
    const string&        smileParamsKey,
    CVolProcessedSP      _processedVol,
    IYieldCurveConstSP   discYC,
    IYieldCurveConstSP   diffYC,
    bool                 _skipFlag,
    double               _flatVolIr,
    const string&        cutoffChoice,
    double               constCutoffValue,
    const string&        corrSwapStart, // eg 1Y  (offset to yc spot date)
    const string&        corrSwapMat,   // eg 10Y (offset to start)
    const string&        corrSwapDCC,   // eg Act/365F
    const string&        corrSwapFreq): // eg 6M
    SRMRatesUtil(
      baseDate,
      numFactors,
      modelParamsKey,
      smileParamsKey,
      _processedVol,
      discYC,
      diffYC,
      _skipFlag,
      _flatVolIr,
      cutoffChoice,
      constCutoffValue,
      corrSwapStart, // eg 1Y  (offset to yc spot date)
      corrSwapMat,   // eg 10Y (offset to start)
      corrSwapDCC,   // eg Act/365F
      corrSwapFreq)
{
    static const string method("SRMRatesLiborUtil::SRMRatesLiborUtil");
    
	try
	{
        nbFactors = numFactors;

        // get alpha, beta, qLeft, qRight
        IRCalib::SmileRequest smileRequest(smileParamsKey);
        CVolProcessed* vol = diffYC->getProcessedVol(&smileRequest);
        smartPtr<IRCalib::VolProcessed> volData(
            &dynamic_cast<IRCalib::VolProcessed&>(*vol));

        const DoubleArray& smileParams = volData->getParams();
        
		if (smileParams.size() < 3)
		{
            // internal error
            throw ModelException(method, "Number of ir vol smile params wrong");
        }
        
		qLeft = 1.0 - smileParams[0]; 
        qRight = 1.0 - smileParams[1];

        fwdShift = smileParams[2];
        
		if (Maths::equals(fwdShift, -1.0))
		{
            throw ModelException(method, 
                                 "Pivot ratio is 0 for "+diffYC->getCcy());
        }
        
		IRCalib::ModelRequest modelRequest(modelParamsKey);
        vol = diffYC->getProcessedVol(&modelRequest);
        volData.reset(&dynamic_cast<IRCalib::VolProcessed&>(*vol));
        const DoubleArray& modelParams = volData->getParams();

       
    } 
	catch (exception& e)
	{
        throw ModelException(e, method, "For currency "+discYC->getCcy());
    }
}

/** number of factors in model */
int SRMRatesLiborUtil::numFactors() const
{
    return nbFactors;
}

void SRMRatesLiborUtil::setTimeLine(DateTimeArrayConstSP simDates) // called when we know allDates
{
    static const string method("SRMRatesLiborUtil::setTimeLine");
    
	if (initialized) 
	{
        if (0 && (simDates->size() != dates->size() || ! DateTime::isSubset(*simDates, *dates)))
            throw ModelException(method, "Re-initialized with a different timeline");
        return;
    }

    initialized = true;
    
	try 
	{
        dates = simDates;
        calcExtendedTimeLine(); // get 'extended' timeline
        
        if (!processedVol){
            // this wasn't in the original plan but it seems you can do stuff
            // even if you don't have the swaption vols
            spotVol(flatVolIr); // populate SpotVol with 1.0 at each point
            return;
        }
        else 
		{
            if (VolProcessedBSIR::TYPE->isAssignableFrom(processedVol->getClass()))
            {
                VolProcessedBSIRSP processedVol = VolProcessedBSIRSP::dynamicCast(this->processedVol);
                swapFrequency = processedVol->getSwapFrequency();
                swapDCC = processedVol->getSwapDCC();
                // get hold of swaptionExpiries, swapStartDates, swapMatDates, vols
                processedVol->getBMDetails(swaptionExpiries, swapStartDates,
                                            swapMatDates, swaptionVols);
        
                // populate LastDate - first find benchmark on/after last sim date
                int bmIdx = (*dates).back().findUpper(swaptionExpiries);
                if (bmIdx == swaptionExpiries.size())
			    {
                    bmIdx--;
                }
                LastDate = swapMatDates[bmIdx].max((*dates).back());
                spotVolBM(skipFlag); // calculate SpotVol and tau
            }
            else
            {
                // what to do here?
            }
        }
        
    } catch (exception& e){
        throw ModelException(e, method, "For currency "+discYC->getCcy());
    }
  
}

/** vol param */
double SRMRatesLiborUtil::getAlpha(int factorIdx) const
{
	return 0.;
}

void SRMRatesLiborUtil::getAlpha(vector<double>& alphaByFactor) const
{
}

void SRMRatesLiborUtil::getBeta(vector<double>& betaByFactor) const
{
}

/** Returns beta for specified factor */
double SRMRatesLiborUtil::getBeta(int factorIdx) const
{
	return 0.;
    //return model->factors[factorIdx].beta;
}

/** Returns rho[rhoIndex] */
double SRMRatesLiborUtil::getRho(int rhoIndex) const
{
	return 0.;
    //return model->rho[rhoIndex];
}

/** Returns the model rho parameters - length corresponds to off diagonal
    elements of a symmetric matrix (size of which is number of factors) */
const vector<double>& SRMRatesLiborUtil::getRho() const
{
	static vector<double> dummy(nbFactors, 0 );

	return (dummy);
    //return model->rho;
}

/* Calibration routine : populates SpotVol and tau array */
void SRMRatesLiborUtil::spotVolBM(bool skipFlag) // From swapvol::IR_SpotVol_BM
{
    try{
        swaptionSpotVol.resize(swaptionExpiries.size());
        vector<double> swaptionSpotVolNew(swaptionSpotVol.begin(), swaptionSpotVol.end());
        swaptionSpotVolAtSimDates = SRMUtil::extendVol(swaptionExpiries,
                                              swaptionSpotVolNew, 
                                              (*dates));

        /* Interpolate spot vol throughout the timeline */
        extendSpotVol(swaptionSpotVol);
    } catch (exception& e){
        throw ModelException(e, "SRMRatesLiborUtil::spotVolBM");
    }
}

/** Calls bFactor for the swap used for correlation purposes */
vector<double> SRMRatesLiborUtil::bFactor() const
{
    return bFactor(corrSwapStart,
                   corrSwapMat,
                   corrSwapDCC,
                   MaturityPeriodSP(new MaturityPeriod(corrSwapFreq)));
}

/** Calls bFactor for the swap given by the specified index */
vector<double> SRMRatesLiborUtil::bFactor(int swapIndex) const
{
    return bFactor(swapStartDates[swapIndex], 
                   swapMatDates[swapIndex],
                   swapDCC, 
                   swapFrequency);
}

/*****  From swapvol.c: BFactor    *****************************************/
/*
*       Determine the value of the B coefficient (see Vladimir's memo)
*       this function is Flat Forward ready
*/
vector<double> SRMRatesLiborUtil::bFactor(const DateTime&      swapStart,
                                     const DateTime&      swapMat,
                                     DayCountConventionSP swapDayCC,
                                     MaturityPeriodSP     swapFreq) const
{
    /* MAR: Note that in the SRM3 code the zero curve passed in thinks that
       today is in fact value date (ie spot date) so all pv's etc are wrt
       spot date */
    static const string method("SRMRatesLiborUtil::bFactor");
    try{

    } catch (exception& e){
        throw ModelException(e, method);
    }
    return vector<double>();
}  /* bFactor */

// returns the instantaneous vol by factor in the form (factor, time point)
void SRMRatesLiborUtil::instFactorVol(
			vector< vector<double> >& vol,				  
	        const vector<double>& DeltaTime,  // (I)  passed for convenience 
	        const DateTimeArray& TPDate,      // (I)
            const vector<double>& IrFwdRate,  // (I)  pre-computed for performance
	        int expiryIndex,                  // (I)  index in TPDate
			int fwdMatIndex) const            // (I)  in case, the underlying has mat > expiry
{
	static const string method("SRMRatesLiborUtil::instFactorVol");

	try
	{
	    // checks 
		if ( ((int)DeltaTime.size() < TPDate.size() - 1) )
		{
			throw ModelException(method, "DeltaTime vector is too short");
		}

		int NbFac = this->numFactors();
		int NbTP = TPDate.size();

		if ((int)vol.size() < NbFac)
		{
			throw ModelException(method, "Factor dimension of vol is too low");
		}

		return;

	}
	catch (exception& e)
	{
		throw ModelException(e, method, "For currency "+discYC->getCcy());
	}
}

// Calculates a 'factor variance' as used by equity/FX vol calibration 
void SRMRatesLiborUtil::irAssetVariance(
    vector<double>& irVariance,        // (O)  variance over [T(i-1),  T(i)] for all i <= N
    vector<double>& irAssetCovar,      // (O)  covariance in [T(i-1),  T(i)] for all i <= N
	const vector<double>& rhoEqIr,     // (I)
	const vector<double>& DeltaTime,   // (I)  passed for convenience 
	const DateTimeArray& TPDate,       // (I)  needed for simpleFwdCurve
    const vector<double>& IrFwdRate,   // (I)  pre-computed for performance
	const vector<int>& periodIndex,    // (I)  indices in TPDate
 	int fwdMatIndex) const             // (I)  indices in TPDate
{
	static const string method("SRMRatesLiborUtil::irAssetVariance");

	try
	{
	    // checks 
		if (irVariance.size() != irAssetCovar.size())
		{
			throw ModelException(method, "irVariance and irAssetCovar have different lengths");
		}

		if (irVariance.size() < periodIndex.size())
		{
			throw ModelException(method, "irVariance vector is too short");
		}

		if ((int) DeltaTime.size() < TPDate.size() - 1)
		{
			throw ModelException(method, "DeltaTime vector is too short");
		}

		int IrNbFac = this->numFactors(); 
		int NbNoEqInt = IrNbFac + (IrNbFac - 1) * IrNbFac / 2;
		int NbWithEqInt = IrNbFac;

		if ((int)rhoEqIr.size() < IrNbFac)
		{
			throw ModelException(method, "Too few correlations between EQ and IR factors");
		}

		// extend IR vols 
		int NbTP = TPDate.size();
		vector<double> irExtSpotVol(NbTP-1);
		const DateTimeArray& irDates = this->getExtendedTimeLine();
		int irNbPoint = irDates.size();
		const vector<double>& irSpotVol = this->getSpotVols();
		int i, k, l;

		for (i = 0, k = 1 ; i <= NbTP - 2; i++) 
		{
			if (k < (irNbPoint - 1)) 
			{
				if (TPDate[i] >= irDates[k]) 
				{
					k++;
				}
			}
			irExtSpotVol[i] = irSpotVol[k-1];
		}

		vector<double> IrA(IrNbFac);

		int NbIndices = periodIndex.size();
		int beginIndex = periodIndex[NbIndices-1];
		int endIndex = beginIndex;

		if (fwdMatIndex < endIndex)
		{
			throw ModelException(method, "forward maturity < expiry");
		}

		//for (i = fwdMatIndex - 1; i >= endIndex; i--) 
		//{
		//	this->aFactor(IrFwdRate[i], DeltaTime[i], IrA);
		//}

		// start calculation of variances
		for (k = NbIndices - 1; k >= 0; k--)  
		{   
			endIndex = beginIndex;
			beginIndex = (k>0) ? periodIndex[k-1] : 0; 

			// reset
			vector<double> IrVolOnly(NbNoEqInt);  // NbNoEqInt    = 1 (IR=1F), 3 (IR=2F), 6 (IR=3F) (s1,s2,s3,s12,s13,s23)
			vector<double> IrAsset(NbWithEqInt);  // NbWithEqInt  = 1 (IR=1F), 2 (IR=2F), 3 (IR=3F) (s1e,s2e,s3e)

			for (i = endIndex - 1; i >= beginIndex; i--) 
			{
				double deltaTime = DeltaTime[i]; // for ease
				double irVol = irExtSpotVol[i];

				// NB. SRM3	has	some adjustments due to	dividends, but also
				// requires	that there be no dollar	divs. This actually
				// means the adjustment	is 0. So, here we ignore the
				// adjustment code entirely.
				//this->aFactor(IrFwdRate[i], deltaTime, IrA);

				// variance of IR factors
				//vector<double>::iterator irVar = IrVolOnly.begin();
				//irVar = 
				//	this->factorVariance(irVar, IrA,
				//	irVol, deltaTime);
				//// covariance between IR factors
				//irVar = 
				//	this->factorCovariance(irVar, IrA,
				//	irVol, deltaTime);
				//// covariance between EQ/FX and IR
				//vector<double>::iterator irEqCoVar = IrAsset.begin();
				//irEqCoVar = 
				//	this->factorFXCovariance(true, // "+" and not "-" contribution
				//	irEqCoVar, IrA, rhoEqIr,
				//	irVol, deltaTime);
			}   /* for i*/

			irVariance[k] = 0.;
			irAssetCovar[k] = 0.;

			for(l = 0; l < NbNoEqInt; l++)
			{
				irVariance[k] += IrVolOnly[l];
			}

			for(l = 0; l < NbWithEqInt; l++)
			{
				irAssetCovar[k] += IrAsset[l];
			}
		}

		return;

	}
	catch (exception& e)
	{
		throw ModelException(e, method, "For currency "+discYC->getCcy());
	}
}

// TODO
/*************** From util_s.c:Triangulation ******************************
* Produces the "usual" Lower triangular matrix for orthogonalising corrolated
* factors.  
* NOTE: If Nbfac < 3 then unused matrix elements are set to zero.
*
****************************************************************************/
DoubleMatrix SRMRatesLiborUtil::triangulation() const
{
    static const string method("SRMRatesLiborUtil::triangulation");
    // for ease
    int Nbfac = numFactors();   /* Number of factors */
    /* initialise matrix (return value) */
    DoubleMatrix TriangMtx(Nbfac, Nbfac);

    return TriangMtx;
}/* Triangulation */


/** this function returns product of spot vols and forward rates
    "interest rate basis point vol". The array returned is of the same
    length as the supplied array. 
    From CMLib:SpreadCurveAlgorithms.cpp:BasisPointVol */
DoubleArray SRMRatesLiborUtil::basisPointVol(
    const DateTimeArray& dates) const /* excludes today */
{
    const DateTimeArray& volDates = swaptionExpiries; // make port easier

    // currentVol is in effect until volEndDate
    // after that volEndDate is updated with volDates[nextVolIndex]
    // we check that the dates array contains all dates from volDates
    ASSERT(DateTime::isSubset(dates, volDates));

    DateTime volEndDate = baseDate;
    int nextVolIndex = 0;

    int numDates = dates.size(); // for ease
    DoubleArray rBpVol(numDates);
    double currentVol=0;
    Actual365F act365F;
    for (int n = 0; n < numDates; n++) {
        const DateTime& startDate = n == 0? baseDate: dates[n-1];
        if (startDate == volEndDate) {
            bool behindLastVol = nextVolIndex == volDates.size();

            // for dates past the last vol date we take overnight forward vol
            volEndDate = behindLastVol? 
                (volDates.back().rollDate(1)): volDates[nextVolIndex];
            double spotVol = behindLastVol? 
                swaptionSpotVol.back(): swaptionSpotVol[nextVolIndex++];
            currentVol = spotVol *
                diffYC->fwd(startDate, volEndDate, 
                            &act365F, CompoundBasis::ANNUAL); // why ANNUAL?

            if (behindLastVol){
                // if we are past the last vol 
                // use current vol for the rest of the timeline
                volEndDate = dates.back();
            }
        }
        ASSERT(startDate < volEndDate);
        
        rBpVol[n] = currentVol; // save the vol
    }
    return rBpVol;
}

DoubleMatrix SRMRatesLiborUtil::get3A() const
{
   return DoubleMatrix();
}

DRLIB_END_NAMESPACE