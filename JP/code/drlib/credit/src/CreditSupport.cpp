//----------------------------------------------------------------------------
//
//   Group       : EDG Derivatives Research
//
//   Filename    : CreditSupport.cpp
//
//   Description : Credit support class
//
//   Author      : Ning Shen
//
//   Date        : 10 June 2002
//
//
//----------------------------------------------------------------------------

#include "edginc/config.hpp"
#include "edginc/Maths.hpp"
#include "edginc/CreditSupport.hpp"
#include "edginc/StruckEquity.hpp"
#include "edginc/ThetaFwdRate.hpp"

DRLIB_BEGIN_NAMESPACE

//////////////////////// CreditUnd ////////////////////////////////////////

AssetUnd::AssetUnd( IGeneralAsset * asset ) :
    m_asset( asset )
{
}
string AssetUnd::getName() const
{
    return m_asset->getName();
}
string AssetUnd::getTrueName() const
{
    return m_asset->getTrueName();
}
CVolProcessed * AssetUnd::getProcessedVol(
    const CVolRequest * volRequest ) const
{
    return m_asset->getProcessedVol( volRequest );
}
void AssetUnd::fwdValue(
    const DateTimeArray & dateList,
    CDoubleArray        & result ) const
{
    m_asset->fwdValue( dateList, result );
}
double AssetUnd::getSpot() const
{
    return m_asset->getSpot();
}
void AssetUnd::setSpot(double spot)
{
    SpotLevel newSpot(spot, true);  // set floorNegFwdPrice to true;
 
    OutputNameConstSP name(new OutputName(getTrueName()));
    if (!newSpot.findAndShift(IObjectSP::attachToRef(m_asset), name)){
        throw ModelException( "AssetUnd::setSpot", "Could not change spot of '" + getTrueName() + "', spot level may not be tweakable" );
    }
}
/** shift value date to newDate */
void CreditSupport::shiftValueDate(IObjectSP& inst, const DateTime& valDate, const DateTime& newDate)
{
    // blank holiday for shifting value date
    static HolidaySP nohols(Holiday::noHolidays());
    if (useThetaFwdRate)
    {
        ThetaFwdRate thetaShift(newDate.getDate() - valDate.getDate() , nohols);
        thetaShift.applyScenario(inst);
    }
    else
    {
        Theta thetaShift(newDate.getDate() - valDate.getDate() , nohols);
        thetaShift.applyScenario(inst);
    }
}

/** put fixng dates into groups according to sample date intervals */
void CreditSupport::groupFixingDates(const DateTimeArray& fixingDatesInput, 
                                     const DateTimeArray& pathDates,
                                     vector<DateTimeArray>& fixDates)
{
    fixDates.resize(pathDates.size());
    
    int j = 0;
    // group fixing dates for each path dates interval
    for (int i=0; i<fixingDatesInput.size(); i++)
    {
        // locate next path date >= fixing date
        while (j<pathDates.size() && pathDates[j] < fixingDatesInput[i])
            j ++;
        
        if (j==pathDates.size())
            break; // done
        
        if (j > 0)
        {
            DateTimeArray dates;
            while (pathDates[j] >= fixingDatesInput[i])
            {
                fixDates[j-1].push_back(fixingDatesInput[i]);
                i++;
                if (i==fixingDatesInput.size())
                    break; // done
            }
        }
    }
}

/** compute prices and store in priceCache */
void CreditSupport::computePriceCache
(
	IObjectSP& inst,	
    const DateTimeArray& dates,
    CreditSupport *creditSupport,
    const IModelSP &model,
    DoubleArray &priceCache					// output
    )
{
    const string method = "CreditSupport::computePriceAndDeltaGrids";
    int i;
    
    // keep value date
    CInstrument* instPtr = dynamic_cast<CInstrument*>(inst.get());
    if (!instPtr)
    {
        throw ModelException(method, "first parameter must be a smart pointer to an instrument.");
    }
    // pre-calc forwards for linear interpolation
    priceCache.resize(dates.size());
    //    deltaCache.resize(dates.size());
    
    CResults results;
    // loop through dates */
    for (i = 0; i < dates.size(); i++)
    {
        CreditSupport::shiftValueDate(inst, instPtr->getValueDate(), dates[i]);
        model->Price(instPtr,smartPtr<Control>(new Control(SensitivityArrayConstSP(),
                        OutputRequestArrayConstSP(   ),0,"")).get(), &results);
        priceCache[i] = results.retrievePrice();
    }
}

/** compute delta's and store in deltaCache 
    priceCache must already be calculated !!! */
void CreditSupport::computeDeltaCache
(
     IObjectSP& inst,	
     const DateTimeArray& dates,
     CreditSupport *creditSupport,
     const IModelSP &model,
     double tweakSize,
     const DoubleArray &priceCache,
     DoubleArray &deltaCache					// output
	)
{
    const string method = "CreditSupport::computePriceAndDeltaGrids";
    int i;
    
    // keep value date
    CInstrument* instPtr = dynamic_cast<CInstrument*>(inst.get());
    if (!instPtr)
    {
        throw ModelException(method, "first parameter must be a smart pointer to an instrument.");
    }
    
    deltaCache.resize(dates.size());

    CResults results;
    CreditUndSP und = creditSupport->getUnderlier();
    double spot = und->getSpot();
    
    for (i = 0; i < dates.size(); i++)
    {
        // restore original spot. useful for setting spot-at-forward-start
        und->setSpot( spot );
        CreditSupport::shiftValueDate(inst, instPtr->getValueDate(), dates[i]);

        // shift spot for delta
        und->setSpot( spot*(1.0 + tweakSize) );
        model->Price(instPtr, smartPtr<Control>(new Control(SensitivityArrayConstSP(   ),
                                                            OutputRequestArrayConstSP(   ),0,"")).get(), &results);
        deltaCache[i] = (results.retrievePrice() - priceCache[i])/(tweakSize*spot);
    }
}

/** compute scaling factor for forward starting. it's ratio of spotRef over spotAtStart */
double CreditSupport::computeScaleFactor
(
    DateTime			startDate,
    const DateTimeArray	&dates,
    const double		*spots,
    DateTime			valueDate,
    double				spotRef)
{
    DateTime date0;
    double spot0, frac, spotAtStart, scaleFactor = 1.0;

    for(int i=0; i<dates.size(); i++)
    {
        if( dates[i] >= startDate )
        {
            if( i==0 ) // in case, dates[0] is already after start, use value date for dates[i-1]
            {
                ASSERT( valueDate <= startDate );
                date0 = valueDate;
                spot0 = spotRef;
            }
            else
            {
                date0 = dates[i-1];
                spot0 = spots[i-1];
            }
			
            frac = (double)(startDate.getDate() - dates[i].getDate()) /
                (double)(date0.getDate() - dates[i].getDate());
            spotAtStart = spots[i] + (spot0 - spots[i])*frac;
            spotAtStart = Maths::max( spotAtStart, MIN_EQUITY_PRICE);
            scaleFactor = spotRef/spotAtStart;
            break;
        }
    }

    return scaleFactor;
}

/** set up spot grid spanning around the equity atm fwd price */
void CreditSupport::setupSpotGrid(
		const DateTimeArray	&dates, 
		const DoubleArray	&atmFwd,
		const DoubleArray	&atmVar,
		int					numFactor,
		const double		*factor,
		DoubleArrayArray	&spotGrid)	// output
{
    int i, j;
    double var = 0.0, v_dt, center;

    spotGrid.resize(dates.size());
    for (i=0; i<dates.size(); i++)
    {
        spotGrid[i].resize(numFactor);

        // grid scaling according to vol
        // accumulate var up to t
        var += atmVar[i];
        v_dt = sqrt(var);
		center = atmFwd[i] * exp( -0.5 * var );
        for (j=0; j<numFactor; j++)
        {
            spotGrid[i][j] = center * exp( factor[j]*v_dt );
		}
	}

	return;
}


/** returns linearly interpolated value on otherDate given 
function goes through (date1, value1) and (date2, value2) */
double CreditSupport::interpolateDates(const DateTime& date1, double value1, 
                                       const DateTime& date2, double value2,
                                       const DateTime& otherDate)
{	
    if (date2.getDate() == date1.getDate())
    {
        throw ModelException("CreditSupport::interpolateDates",
            "Requires two distinct dates.");
    }
    
    return value1 + 
        (double)(otherDate.getDate() - date1.getDate()) / (double)(date2.getDate() - date1.getDate()) *
        (value2 - value1);
}

/** populates sampleLevels using linear interpolation.  Similar to rollLinear, 
only applied to flat data, for example as seen in the Equity Swap. */
void CreditSupport::populateSampleArray(const DateTimeArray& sampleDates,
                                        DoubleArray&   sampleLevels,
                                        const DateTime& date1, double value1, 
                                        const DateTime& date2, double value2)
{
    if (sampleLevels.size() != sampleDates.size())
    {
        sampleLevels.resize(sampleDates.size());
    }
    for (int i = 0; i < sampleDates.size(); i++)
    {
        if( (sampleDates[i] >= date1) == (sampleDates[i] <= date2) )
            sampleLevels[i] = interpolateDates(date1, value1, date2, value2, sampleDates[i]);
    }
}

void CreditSupport::setValueDateShift(bool useFwdRate)
{
    useThetaFwdRate = useFwdRate;
}

DRLIB_END_NAMESPACE

