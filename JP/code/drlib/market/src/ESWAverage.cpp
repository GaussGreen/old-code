//----------------------------------------------------------------------------
//
//   Group       : EDG Derivatives Research
//
//   Filename    : ESWAverage.cpp
//
//   Description   averaging classes for equity swap.
//
//
//----------------------------------------------------------------------------

#include "edginc/config.hpp"
#include "edginc/Maths.hpp" //isZero//
#include "edginc/ESWEquity.hpp" //NumHistoricalDates

DRLIB_BEGIN_NAMESPACE

// to do: put this in a header file, have SampleList.cpp use it as well
#define EDG_CLOSE_ENOUGH_TO_ONE 1e-9

/////////////////////////////////////////////
////////////////// ESWAvIn /////////
/////////////////////////////////////////////
// helpers
void ESWAvIn::load(CClassSP& clazz){
    clazz->setPublic(); // make visible to EAS/spreadsheet
    REGISTER(ESWAvIn, clazz);
    SUPERCLASS(CObject);
    EMPTY_SHELL_METHOD(defaultESWAvIn);
    FIELD(sampleDates,  "average in sample dates");
    FIELD(sampleLevels,  "spot sample levels");
    FIELD_MAKE_OPTIONAL(sampleLevels);
    FIELD(weights,  "percentage sample weighting");
    FIELD(liborFixingDates,  "libor fixing dates for non consolidated averaging in");
    FIELD(liborFixing,  "libor fixing for non consolidated averaging in");
    FIELD_MAKE_OPTIONAL(liborFixing);
    FIELD(rateTypes,  "rate types (e.g. 3M) for non consolidated averaging in");
    FIELD_MAKE_OPTIONAL(rateTypes);
//    FIELD(spreads,  "spreads over libor for non consolidated averaging in");
//    FIELD_MAKE_OPTIONAL(spreads);
    FIELD(fxFixing,  "fx fixing levels");
    FIELD_MAKE_OPTIONAL(fxFixing);
    FIELD(fwdLevels,  "agreed upon forward prices to reference date in consolidation mode");
    FIELD_MAKE_OPTIONAL(fwdLevels);
    FIELD(isConsolidated,  "if fixings are consolidated on one refix date");
    FIELD(consoliDate,  "refix date for consolidation");
    FIELD_MAKE_OPTIONAL(consoliDate);
    FIELD(ccyTreatment, "ccyTreatment");
    FIELD_MAKE_TRANSIENT(ccyTreatment);
}

CClassConstSP const ESWAvIn::TYPE = CClass::registerClassLoadMethod(
    "ESWAvIn", typeid(ESWAvIn), load);

// constructor
ESWAvIn::ESWAvIn(): CObject(TYPE){}
int ESWAvIn::getLength() const
{
    return sampleDates.size();
}



void ESWAvIn::validate(const DateTime& beginDate, 
                       const DateTime& endDate,
                       const DateTime& firstEquityAccrueEndDate,
                       const DateTime& valDate,
                       bool                    isFloating) const
{
    int numDates = getLength();
    if (isConsolidated)
    {
        if (ESWEquity::NumHistoricalDates(sampleDates, consoliDate,true) < numDates)
        {
            throw ModelException("Equity Swap: Average In, consolidated case: all sample dates"  
                                 " should lie on or before consolidation date.");
        }

        // consolidation date must be in first period.
        if (consoliDate > firstEquityAccrueEndDate)                                                     
        {
            throw ModelException("Equity Swap: Average In, consolidated case: consolidated date"  
                                 " canoot lie after first equity accrue end date.");
        }
    }

    int i;
    for ( i = 1 ; i < numDates; ++ i ) 
    {
        if (sampleDates[i]<sampleDates[i-1])
        {
            throw ModelException("Equity Swap: Average in leg", 
                                 "Date list must be in increasing order.\nDate " + Format::toString(i-1) + " = "
                                 + sampleDates[i-1].toString() + "\nDate " + Format::toString(i) + " = " + 
                                 sampleDates[i].toString()    );
        }
    }

    if (numDates != weights.size())
    {
        throw ModelException("Equity Swap: Average in schedule", 
                             "sampleDates and weights must all have the same size.");
    }
    

    double sum = 0;
    for ( i = 0 ; i < numDates; ++ i ) 
    {
        sum += weights[i];
        if (weights[i] < 0)
        {
            throw ModelException("Equity Swap: Average In:", "Avg in weights cannot be negative!");
        }
    }
        
    if (!Maths::areEqualWithinTol(sum, 1.0, EDG_CLOSE_ENOUGH_TO_ONE))
    {
        throw ModelException("Equity Swap: sum of averaging in weights must equal 1, \n"
                             "but it equals " + Format::toString(sum));
    }


    if ((sampleDates[0]<beginDate)||(sampleDates[numDates-1]>=endDate))
    {
        throw ModelException("Equity Swap: Average in schedule", 
                             "Averaging in schedule must be entirely contained in the first equity refixing period."   );
    }

    // to do: should be if isFloating and non Consolidated
    if (!isConsolidated)
    {
        if (isFloating)
        {
            if (rateTypes.size() != sampleDates.size())
            {
                throw ModelException("Equity Swap: Average in schedule", 
                                     "Rate types array size must equal sample dates array size."   );
            }

            if (liborFixingDates.size() != sampleDates.size())
            {
                throw ModelException("Equity Swap: Average in schedule", 
                                     "Libor Fixing Dates array size must equal sample dates array size."   );
            }
        }
                
        if (liborFixing.size() < ESWEquity::NumHistoricalDates(liborFixingDates, valDate, true))
        {
            throw ModelException("Equity Swap: Average in schedule", 
                                 "Not enought libor fixings!  (Require at least the number of historical libor fixing dates.)");
        }
    }
}

double ESWAvIn::getLevel(const DateTime& valueDate, int iAveSample, CAssetSP asset) const
{       
    // to do: code cleanup.  Can reuse this code, and call the same func for both in & out

    /** if in past or present, or positive sample later today (i.e EOD sample on SOD pricing), use recorded sample Level */
    const DateTime& fixDate = sampleDates[iAveSample];

    double fx = 1.0;
    if (ccyTreatment == CAsset::CCY_TREATMENT_STRUCK)
    {
        if (
            valueDate >= fixDate ||
            (fixDate.equals(valueDate, false) && Maths::isPositive(fxFixing[iAveSample])) 
            )
        {
            fx = fxFixing[iAveSample];
        }
        else
        {
            fx = dynamic_cast<Asset::IStruck&>(*asset).getFXSpot();
        }
    }

    if (
        valueDate >= fixDate ||
        (fixDate.equals(valueDate, false) && Maths::isPositive(sampleLevels[iAveSample])) 
        )
                
    {
        return fx*sampleLevels[iAveSample];
    }
    else /** in future, need to compute forward price */
    {
        return asset->fwdValue( sampleDates[iAveSample] );
    }
        
}


double ESWAvIn::numDollarsPerShare(const DateTime& valueDate, CAssetSP asset, const DateTime& refDate, bool soFar) const
{/** Compute number of dollars */


    int numSamples;

    if (soFar)
    {
        numSamples = ESWEquity::NumHistoricalDates(sampleDates, refDate);
    }
    else
    {
        numSamples = getLength();
    }

/** Compute inLevel */

    double sum = 0;
// to do: EOD sample treatment here
    if (isConsolidated)
    {
        for (int iAvDate=0;iAvDate<numSamples;iAvDate++)
        {
            /** if in past, use recorded fwdValue */
            if(valueDate > sampleDates[iAvDate])
            {
                sum+=weights[iAvDate]*fwdLevels[iAvDate];
            }
                
            else /** in future, need to compute forward price */
            {
                sum+=
                    weights[iAvDate]*(asset->fwdValue(
                        consoliDate));
            }
        }
    }

    else
    {
        for (int iAvDate=0;iAvDate<numSamples;iAvDate++)
        {
            sum+=weights[iAvDate]*getLevel(valueDate, iAvDate,asset);
                        
        }
                
    }
    return sum;
}



double ESWAvIn::numSharesPerDollar(const DateTime& valueDate, CAssetSP asset, const DateTime& refDate, bool soFar) const
{
    int numSamples;

    if (soFar)
    {
        numSamples = ESWEquity::NumHistoricalDates(sampleDates, refDate);
    }
    else
    {
        numSamples = getLength();
    }
                
    double sum = 0;
    // to do: EOD sample treatment here.
    if (isConsolidated)
    {
        for (int iAvDate=0;iAvDate<numSamples;iAvDate++)
        {
            /** if in past, use recorded fwdValue */
            if(valueDate > sampleDates[iAvDate])
            {
                sum+=weights[iAvDate]/fwdLevels[iAvDate];
            }
                        
            else /** in future, need to compute E[1/forward price] */
            {
                sum+=
                    weights[iAvDate]*(asset->expNumberShares(
                        valueDate, consoliDate, false));
            }
        }
    }
        
    else /** not Consolidated */
    {
                
        for (int iAvDate=0;iAvDate<numSamples;iAvDate++)
        {
            sum+=weights[iAvDate]/getLevel(valueDate, iAvDate,asset);
        }
                
    }
    return sum;
}

double ESWAvIn::fractionOfShares(const DateTime& valueDate, CAssetSP asset, const DateTime& refDate) const
{

    int numSamples = ESWEquity::NumHistoricalDates(sampleDates, refDate);

    double sum= 0;
    for (int iAvDate=0;iAvDate<numSamples;iAvDate++)
    {
        sum+=weights[iAvDate];
                
    }
                
    return sum;

}


void ESWAvIn::setFixing(const DateTime& valDate, const DateTime& rollDate, CAssetWrapper assetWrapper, const Theta* shift,
                        string ccyTreatment)
{
    // to do: (code cleanup) replace the nasty condition (which is repeated 3x below) to determine whether or not to update history with a func.
    double fx = 1.0;
    for (int i=0; i<sampleDates.size(); i++)
    {
        if (rollDate >= sampleDates[i] && valDate <= sampleDates[i])
        {
                                // currency struck case: need to update fx levels
            if(ccyTreatment == CAsset::CCY_TREATMENT_STRUCK)
            {
                fx = dynamic_cast<Asset::IStruck&>(
                    *assetWrapper.get()).getFXSpot();
            
                // Only overwrite today's zero or blank levels.
                // Overwrite any level in the future after today.
                // (only compare the dates since we don't want to
                // override EOD sample today)
                if (valDate.getDate() < sampleDates[i].getDate() ||                     // always overwrite if >= tomorrow
                    (fxFixing.size() > i && Maths::isZero(fxFixing[i])) ||          // zero level
                    (fxFixing.size()<=i))                                                                           // blank level
                {
                    if (fxFixing.size()<=i)
                    {
                        fxFixing.resize(i+1);
                    }
                    fxFixing[i] = fx;
                }
            }

                                //sample levels

                                // Only overwrite today's zero or blank levels.  Overwrite any level in the future after today.
                                // (only compare the dates since we don't want to override EOD sample today)
            if (valDate.getDate() < sampleDates[i].getDate() ||                             // always overwrite if >= tomorrow
                (sampleLevels.size()>i && Maths::isZero(sampleLevels[i]) ||             // zero level
                 sampleLevels.size()<=i))                                                                                // blank level
            {
                if (sampleLevels.size()<=i)
                {
                    sampleLevels.resize(i+1);
                }
                sampleLevels[i] = assetWrapper->getThetaSpotOnDate(shift,sampleDates[i])/fx;
            }
        
                                // consolidated case: need to update forward levels as well             
            if (isConsolidated)
            {
                if (valDate.getDate() < sampleDates[i].getDate() ||                             // always overwrite if >= tomorrow
                    (fwdLevels.size()>i && Maths::isZero(fwdLevels[i]) ||                   // zero level
                     fwdLevels.size()<=i))                                                                                   // blank level
                {
                    if (fwdLevels.size() <= i) { fwdLevels.resize(i+1); }
                    fwdLevels[i] = assetWrapper->fwdFwd(sampleDates[i], sampleLevels[i], this->consoliDate);
                }
            }
        } // if we're rolled over this date
    } // loop through dates

    // need to do: libor fixings
}


void ESWAvIn::init()
{

    int num = sampleDates.size();
    liborFixing.resize(num);
    sampleLevels.resize(num);
    fxFixing.resize(num);

}

void ESWAvIn::validateHistoricalSamples(const DateTime& valDate, bool isFloating) const
{
    int numHistoricalLevels = ESWEquity::NumHistoricalDates(sampleDates, valDate, false);   // false, so exclusive.  Today's data will be 
    if (sampleLevels.size() < numHistoricalLevels)                                                                                      // filled by setFixing.
    {
        throw ModelException("Equity Swap: Average In: there are " + 
                             Format::toString(numHistoricalLevels) +
                             " historical sample date(s), but only " +
                             Format::toString(sampleLevels.size()) +
                             " historical level(s) supplied.");
    }

    // if consolidated, must verify that we have enough fwd levels
    if (isConsolidated)
    {        
        if(fwdLevels.size() < numHistoricalLevels)
        {
            throw ModelException("Equity Swap: Average In, consolidated case: there are " + 
                                 Format::toString(numHistoricalLevels) +
                                 " historical sample date(s), but only " +
                                 Format::toString(fwdLevels.size()) +
                                 " forward level(s) supplied.");
        }


    }

    // if ccyStruck, need fx levels
    if (ccyTreatment == CAsset::CCY_TREATMENT_STRUCK &&
        fxFixing.size() <  numHistoricalLevels)
    {
        throw ModelException("Equity Swap: Average In, struck currency treatment: there are " + 
                             Format::toString(numHistoricalLevels) +
                             " historical sample date(s), but only " +
                             Format::toString(fxFixing.size()) +
                             " historical fx level(s) supplied.");
    }

    if (isFloating && !isConsolidated)
    {
        int numHistLibRefDates = ESWEquity::NumHistoricalDates(liborFixingDates, valDate, false);               // false, so exclusive.  Today's data will be
        if(liborFixing.size() < numHistLibRefDates)                                                                                                             // filled by setFixing.
        {
            throw ModelException("Equity Swap: Average In, non consolidated case: there is(are) " + 
                                 Format::toString(numHistLibRefDates) +
                                 " historical refix date(s), but only " +
                                 Format::toString(liborFixing.size()) +
                                 " libor fixing(s) supplied.");
        }
                                                                                                                                                                                                                
    }

}
    



//double ESWAvIn::getNotional(int iAveSample) const
//{
//      return weights[iAveSample]/sampleLevels[iAveSample];
//      
//      // to do: fx //
//}



/////////////////////////////////////////////
////////////////// ESWAvOut /////////
/////////////////////////////////////////////
// helpers
void ESWAvOut::load(CClassSP& clazz){
    clazz->setPublic(); // make visible to EAS/spreadsheet
    REGISTER(ESWAvOut, clazz);
    SUPERCLASS(CObject);
    EMPTY_SHELL_METHOD(defaultESWAvOut);
    FIELD(sampleDates,  "average out sample dates");
    FIELD(sampleLevels,  "spot sample levels");
    FIELD_MAKE_OPTIONAL(sampleLevels);
    FIELD(weights,  "percentage sample weighting");
    FIELD(discountFactors,  "if accrueAvOutType is fixed, required factors.");
    FIELD_MAKE_OPTIONAL(discountFactors);
    FIELD(accrueIRonFullNotional,  "determines if IR leg accrues on full or amortizing notional in avg out.");
    FIELD(accrueAvOutType,  
                 "determines how avging out accrues from settlement to pay date. \n"
                 "0=no accrual, 1=libor w/o spread"
                 "2=libor with spread, 3=use DiscountFactor input");
    FIELD(fxFixing,  "fx fixing levels");
    FIELD_MAKE_OPTIONAL(fxFixing);
    FIELD(payDates,  "pay date list, one for each average out date.  These dates are also the corresponding IR pay dates.");
    FIELD(ccyTreatment, "ccyTreatment");
    FIELD_MAKE_TRANSIENT(ccyTreatment);

}

CClassConstSP const ESWAvOut::TYPE = CClass::registerClassLoadMethod(
    "ESWAvOut", typeid(ESWAvOut), load);

// constructor
ESWAvOut::ESWAvOut(): CObject(TYPE){}

int ESWAvOut::getLength() const
{
    return sampleDates.size();
}

void ESWAvOut::validate(const DateTime& beginDate, const DateTime& endDate, const DateTime& valDate) const
{
    const string method = "ESWAvOut:validate";
    int numDates = getLength();
    int i;
    for ( i = 1 ; i < numDates; ++ i ) 
    {
        if ( sampleDates[i]<sampleDates[i-1])
        {
            throw ModelException("Equity Swap: Average out schedule", 
                                 "Date list must be in increasing order.\nDate " + Format::toString(i-1) + " = "
                                 + sampleDates[i-1].toString() + "\nDate " + Format::toString(i) + " = " + 
                                 sampleDates[i].toString()    );
        }
    }
        
    if (numDates != weights.size() || numDates != payDates.size())
    {
        throw ModelException("Equity Swap: Average out schedule", 
                             "sampleDates, weights, and pay dates arrays must all have the same size.");
    }
    
    // historical levels

        
    // check sum is 1, and that pay dates are on or after sample dates.
    // do two together for efficiency.
    double sum = 0;
    for ( i = 0 ; i < numDates; ++ i ) 
    {
        sum += weights[i];
        if (weights[i] < 0)
        {
            throw ModelException("Equity Swap: Average Out:", "Avg out weights cannot be negative!");
        }
        if (payDates[i] < sampleDates[i])
        {
            throw ModelException("Equity Swap: Average out schedule: \n" 
                                 "Pay Date " + Format::toString(i) 
                                 + " is " + payDates[i].toString() + " which is before \nsample date"
                                 " which is " + sampleDates[i].toString());
        }
    }

    if (!Maths::areEqualWithinTol(sum, 1.0, EDG_CLOSE_ENOUGH_TO_ONE))
    {
        throw ModelException("Equity Swap: sum of averaging out weights must equal 1, \n"
                             "but it equals " + Format::toString(sum));
    }

    if ( sampleDates[0] < beginDate || 
         sampleDates[numDates-1] > endDate )
    {
        throw ModelException("Equity Swap: Average out schedule", 
                             "Averaging out schedule must be entirely contained in last equity refixing period."   );
    }

    if (accrueAvOutType == ESWAvOut::FIXED)  
    {
        if (discountFactors.size() < sampleDates.size())
        {
            throw ModelException("Equity Swap: Average Out, accrual style = FIXED: there are " + 
                                 Format::toString(sampleDates.size()) +
                                 " sample date(s), but only " +
                                 Format::toString(discountFactors.size()) +
                                 " discount factor(s) supplied.");
        }

    }


    
    


    for (i=0; i<discountFactors.size(); i++)
    {   
        if (discountFactors[i] <= 0)
        {
            throw ModelException(method," Discount Factor on averOut date " 
                                 + sampleDates[i].toString() + " must be positive!");
        }
    } 
    

}


double ESWAvOut::getLevel(const DateTime& valueDate, int iAveSample, CAssetSP asset) const
{       
    
    // to do: code cleanup.  Can reuse this code, and call the same func for both in & out
    
    /** if in past or present, or positive sample later today (i.e EOD sample on SOD pricing), use recorded sample Level */
    const DateTime& fixDate = sampleDates[iAveSample];

    double fx = 1.0;
    if (ccyTreatment == CAsset::CCY_TREATMENT_STRUCK)
    {
        if (
            valueDate >= fixDate ||
            (fixDate.equals(valueDate, false) && Maths::isPositive(fxFixing[iAveSample])) 
            )
        {
            fx = fxFixing[iAveSample];
        }
        else
        {
            fx = dynamic_cast<Asset::IStruck&>(*asset).fxFwdValue(fixDate);
        }
    }

    if (
        valueDate >= fixDate ||
        (fixDate.equals(valueDate, false) && Maths::isPositive(sampleLevels[iAveSample])) 
        )
                
    {
        return fx*sampleLevels[iAveSample];
    }
    else /** in future, need to compute forward price */
    {
        return asset->fwdValue( sampleDates[iAveSample] );
    }
}

double ESWAvOut::fractionOfShares(const DateTime& valueDate, CAssetSP asset, const DateTime& refDate) const
{

    int numSamples = sampleDates.size() - ESWEquity::NumHistoricalDates(sampleDates, refDate);

    double sum= 0;
    for (int iAvDate=getLength()-1;iAvDate>getLength()-numSamples-1;iAvDate--)
    {
        sum+=weights[iAvDate];
                
    }
                
    return sum;

}

double ESWAvOut::numDollarsPerShare(const DateTime& valueDate, CAssetSP asset, const DateTime& refDate, bool soFar) const
{/** Compute number of dollars */


    int numSamples;

    if (soFar)
    {
        numSamples = sampleDates.size() - ESWEquity::NumHistoricalDates(sampleDates, refDate);
    }
    else
    {
        numSamples = getLength();
    }
        
/** Compute inLevel */

    double sum= 0;
    for (int iAvDate=getLength();iAvDate>getLength()-numSamples;iAvDate--)
    {
        sum+=weights[iAvDate]*getLevel(valueDate, iAvDate,asset);
                
    }
                

    return sum;
}


double ESWAvOut::numSharesPerDollar(const DateTime& valueDate, CAssetSP asset, const DateTime& refDate, bool soFar) const
{
    int numSamples;

    if (soFar)
    {
        numSamples = sampleDates.size() - ESWEquity::NumHistoricalDates(sampleDates, refDate);
    }
    else
    {
        numSamples = getLength();
    }
                

    double sum = 0;
    for (int iAvDate=getLength();iAvDate>getLength()-numSamples;iAvDate--)
    {
        sum+=weights[iAvDate]/getLevel(valueDate, iAvDate,asset);
    }
                
        
    return sum;
}

void ESWAvOut::setFixing(const DateTime& valDate, const DateTime& rollDate, CAssetWrapper assetWrapper, const Theta* shift,
                         string ccyTreatment)
{
    // to do: code cleanup.  Can reuse this code, and call the same func for both in & out
    double fx = 1.0;
    for (int i=0; i<sampleDates.size(); i++)
    {
        if (rollDate >= sampleDates[i] && valDate <= sampleDates[i])
        {
        
            if(ccyTreatment == CAsset::CCY_TREATMENT_STRUCK)
            {

                fx = dynamic_cast<Asset::IStruck&>(
                    *assetWrapper.get()).getFXSpot();
                // Only overwrite today's zero or blank levels.
                // Overwrite any level in the future after today.
                // (only compare the dates since we don't want to
                // override EOD sample today)
                if (valDate.getDate() < sampleDates[i].getDate() ||                     // always overwrite if >= tomorrow
                    (fxFixing.size() > i && Maths::isZero(fxFixing[i])) ||          // zero level
                    (fxFixing.size()<=i))                                                                           // blank level
                {
                    if (fxFixing.size() <= i)
                    {
                        fxFixing.resize(sampleDates.size());
                    }
                                                
                    fxFixing[i] = fx;
                }
            }
                                

            if ( valDate.getDate() < sampleDates[i].getDate() ||               // don't override EOD levels 
                 (sampleLevels.size() > i && Maths::isZero(sampleLevels[i])) || 
                 sampleLevels.size() <= i)
            {
                if (sampleLevels.size() <= i)
                {
                    sampleLevels.resize(sampleDates.size());
                }
                sampleLevels[i] = assetWrapper->getThetaSpotOnDate(shift,sampleDates[i])/fx;
            }
                                
        }
    }
}

void ESWAvOut::init()
{
    int numLevels = sampleDates.size();
    sampleLevels.resize(numLevels);
    fxFixing.resize(numLevels);
    eqGrowthFactor.resize(numLevels);
}

void ESWAvOut::validateHistoricalSamples(const DateTime& valDate)
                                                                                 
{
    int numHistoricalLevels = ESWEquity::NumHistoricalDates(sampleDates, valDate, false);   // false, so exclusive.  Today's data will be 
    // filled by setFixing.
    
    if (sampleLevels.size() < numHistoricalLevels)
    {
        throw ModelException("Equity Swap: Average Out: there are " + 
                             Format::toString(numHistoricalLevels) +
                             " historical sample date(s), but only \n" +
                             Format::toString(sampleLevels.size()) +
                             " historical level(s) supplied.");
    }
        
    // if ccyStruck, need fx levels
    if (ccyTreatment == CAsset::CCY_TREATMENT_STRUCK &&
        fxFixing.size() <  numHistoricalLevels)
    {
        throw ModelException("Equity Swap: Average Out, struck currency treatment: there are " + 
                             Format::toString(numHistoricalLevels) +
                             " historical sample date(s), but only " +
                             Format::toString(fxFixing.size()) +
                             " historical fx level(s) supplied.");
    }
        
    if (accrueAvOutType == ESWAvOut::LIBOR_WITHOUT_SPREAD ||
        accrueAvOutType == ESWAvOut::LIBOR_WITH_SPREAD)
    {
        if (discountFactors.size() < numHistoricalLevels)
        {
            throw ModelException("Equity Swap: Average Out, accrual style = "
                                 "LIBOR_WITHOUT_SPREAD, LIBOR_WITH_SPREAD: \nthere are " + 
                                 Format::toString(numHistoricalLevels) +
                                 " historical sample date(s), but only " +
                                 Format::toString(discountFactors.size()) +
                                 " discount factor(s) supplied.");
        }
                
    }

}       

bool ESWAverageLoad() {
    return ESWAvIn::TYPE != NULL && ESWAvOut::TYPE != NULL;
}

DRLIB_END_NAMESPACE
