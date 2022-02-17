//----------------------------------------------------------------------------
//
//   Group       : EDG Derivatives Research
//
//   Filename    : IQMCDiffusibleAsset.cpp 
//
//   Description : An interface of hooks for asset diffusion
//

//
//----------------------------------------------------------------------------

#include "edginc/config.hpp"
#include "edginc/IQMCDiffusibleAsset.hpp"
#include "edginc/SVQmcImplemented.hpp"

USING_DRLIB_NAMESPACE


/** creating a state variable to access the DiscFactor information */
SVQmcDiscFactor* IQMCDiffusibleInterestRate::createDFSV(const DateTimeArray& measurementDates, bool mm)
{ 
    if (mm) // moment-matching ON
        return new SVQmcDiscFactorMM(this, measurementDates);
    else    // moment matching OFF
        return new SVQmcDiscFactor(this, measurementDates);
}

/** creating a state variable to access the ExpDiscFactor information */
SVExpectedDiscFactor* IQMCDiffusibleInterestRate::createExpDiscFactorSV(
                            const DateTime&         measurementDate,
                            const DateTimeArray&    futureDates,
                            bool                    doLog,
                            YieldCurveConstSP       yc,
                            bool mm)
{ 
    if (mm) // moment-matching ON
        return new SVQmcExpectedDiscFactorMM( this, 
                                            measurementDate, 
                                            futureDates, 
                                            doLog,
                                            yc);
    else    // moment matching OFF
        return new SVQmcExpectedDiscFactor( this, 
                                            measurementDate, 
                                            futureDates, 
                                            doLog,
                                            yc);
}


/** creating a state variable to access the SurvivalDiscFactor information */
SVSurvivalDiscFactor* IQMCDiffusibleCreditSpreadBase::createSDFSV(
                            const DateTimeArray& measurementDates, bool mm)
{ 
    if (mm) // moment-matching ON
        return new SVQmcSurvivalDiscFactorMM(this, measurementDates);
    else
        return new SVQmcSurvivalDiscFactor(this, measurementDates);
}

/** creating a state variable to access the ExpDiscFactor information */
SVExpSurvDiscFactor* IQMCDiffusibleCreditSpreadBase::createExpSurvDiscFactorSV(
    const DateTime&         measurementDate,
    const DateTimeArray&    futureDates,
    bool                    doLog,
    bool                    mm)    
{ 
    if (mm) // moment-matching ON
        return new SVQmcExpSurvDiscFactorMM(this, measurementDate, 
                                            futureDates, doLog);
    else
        return new SVQmcExpSurvDiscFactor(this, measurementDate, 
                                            futureDates, doLog);
}


/** creating a state variable to gain compact access to a set of SDFs and ESDFs */
SVAggregatedSurvDiscFactor* IQMCDiffusibleCreditSpreadBase::createAggregatedSurvDiscFactorSV(
            DateTimeArraySP   sdfDates,            // set of dates for discount factors
            SpotIdxArraySP    sdfIdxSP,
            DateTimeArraySP   esdfRequestedDates,  //  union of {t_i} for all expected factors
            FwdIdxArraySP     esdfReqIdxSP,
            DateTimeArraySP   esdfForwardDates,    // union of all {T_j} for all expected factors
            FwdIdxArraySP     esdfForIdxSP,
            const DateTime&   maxDiffDate,
            const DateTime&   maxCurveMat,
            bool              doLog,
            bool              mm)
{
    QLIB_VERIFY(mm==false, "The moment-matching is not yet implemented for createAggregatedSurvDiscFactorSV");

    return new SVQmcAggregatedSurvDiscFactor(
                    this,
                    sdfDates,            // set of dates for discount factors
                    sdfIdxSP,
                    esdfRequestedDates,  //  union of {t_i} for all expected factors
                    esdfReqIdxSP,
                    esdfForwardDates,    // union of all {T_j} for all expected factors
                    esdfForIdxSP,
                    maxDiffDate,
                    maxCurveMat,
                    doLog);
}


/** creating a state variable to access the DateOfDefault information */
SVDateOfDefault* IQMCDiffusibleDefaultableCreditSpread::createSVDateOfDefault()
{
    return new SVQmcDateOfDefault(this);
}


/** creating a state variable to access the SurvivalDiscFactor information */
SVSpotFX* IQMCDiffusibleFX::createSpotFXSV(const DateTimeArray& measurementDates, bool mm)
{ 
    if (mm)
        return new SVSpotFXMM(this, measurementDates);
    else
        return new SVSpotFX(this, measurementDates);
}

/** creating a state variable to access the ExpDiscFactor information */
SVExpectedFX* IQMCDiffusibleFX::createExpectedFXSV(
    const DateTime&                measurementDate,
    const DateTimeArray&    futureDates,
    bool                    doLog, 
    bool                    mm)    
{ 
    if (mm)
        return new SVQmcExpectedFXMM(this, measurementDate, futureDates, doLog);
    else
        return new SVQmcExpectedFX(this, measurementDate, futureDates, doLog);
}




/** creating a state variable to access the SurvivalDiscFactor information */
SVSpotEQ* IQMCDiffusibleEQ::createSpotEQSV(const DateTimeArray& measurementDates, bool mm)
{ 
    if (mm)
        return new SVSpotEQMM(this, measurementDates);
    else
        return new SVSpotEQ(this, measurementDates);
}

/** creating a state variable to access the ExpDiscFactor information */
SVExpectedEQ* IQMCDiffusibleEQ::createExpectedEQSV(
    const DateTime&         measurementDate,
    const DateTimeArray&    futureDates,
    bool                    doLog, 
    bool                    mm)    
{ 
    if (mm)
        return new SVQmcExpectedEQMM(this, measurementDate, futureDates, doLog);
    else
        return new SVQmcExpectedEQ(this, measurementDate, futureDates, doLog);
}



/************************************************************************/
/* Energy stuff                                                         */
/************************************************************************/
// creating a state variable to access the expected energy future
SVExpEnergyFuture * IQMCDiffusibleEnergy::createExpEnergyFutureSV(
    const DateTime&         measurementDate,
    const DateTimeArray&    futureDates,
    bool                    doLog,
    bool                    mm)    
{ 
    // disregard the mm key for now -- not implemented for energies yet
    return new SVQmcExpEnergyFuture(this, measurementDate, futureDates, doLog);
}


//////////////////////////////////////////////////////////////////////////
// Basis Spread stuff
//////////////////////////////////////////////////////////////////////////

/** creating a state variable to access the ExpectedBasisFwdSpread information */
SVExpectedBasisFwdSpread* IQMCDiffusibleBasisIndex::createExpBasisFwdSpreadSV(
    const DateTime& measurementDate,
    const DateTimeArray& resetDates )
{
    return new SVQmcExpectedBasisFwdSpread(this, measurementDate, resetDates);
}



