//----------------------------------------------------------------------------
//
//   Group       : Quantitative Research & Development
//
//   Filename    : ExpectedLossSurface.cpp
//
//   Description : Class to store an surface of base expected losses
//
//   Author      : Matthias Arnsdorf
//
//   Date        : October 2005
//
//
//----------------------------------------------------------------------------
#include "edginc/config.hpp"
#include "edginc/Format.hpp"
#include "edginc/ExpectedLossSurface.hpp"
#include "edginc/Maths.hpp"

DRLIB_BEGIN_NAMESPACE


const double ExpectedLossSurface::TINY = 1.0e-10; // A small number.


void ExpectedLossSurface::load(CClassSP& clazz){
    clazz->setPublic();
    REGISTER(ExpectedLossSurface, clazz);
    SUPERCLASS(CObject);
    EMPTY_SHELL_METHOD(defaultExpectedLossSurface);

	FIELD(dates, "Dates for expected losses");
	FIELD(strikes, "Strikes for expected losses");
	FIELD(losses, "Expected losses at grid points");
	FIELD(RR, "Recovery rate");
        
}

CClassConstSP const ExpectedLossSurface::TYPE = CClass::registerClassLoadMethod(
    "ExpectedLossSurface", typeid(ExpectedLossSurface), load);

bool  ExpectedLossSurfaceLoad() {
    return (ExpectedLossSurface::TYPE != 0);
}
 
ExpectedLossSurface::ExpectedLossSurface(): 
    CObject(TYPE)
{
}

IObject* ExpectedLossSurface::defaultExpectedLossSurface(){
    return new ExpectedLossSurface();
}

/************************************************************************************
 validation
 ***********************************************************************************/
void ExpectedLossSurface::validatePop2Object()
{
	static const string method ="ExpectedLossSurface::validatePop2Object";
    try
    {
		 int  i;
        // validation ------------------------------------------------------------------
        if(losses.get()==NULL) throw ModelException("losses array is empty");
		if(strikes.get() == NULL) throw ModelException("strikes array is empty");
        if(dates.get() == NULL) throw ModelException("dates array is empty");
        
		if(strikes->size() < 1) throw ModelException("Need at least one strike");
		if(dates->size() <1) throw ModelException("Need at least one date");

        if (losses->size() != strikes->size())
        {
            throw ModelException("Num Strikes("+ Format::toString(strikes->size()) +
                                 ") different to num cols ("+
                                 Format::toString((*losses).size()) +
                                 ") in losses matrix");
        }

		for(i = 0; i<losses->size();i++)
		{
			if ((*losses)[i].size() != dates->size())
			{
				throw ModelException("Num dates ("+
									 Format::toString(dates->size()) +
									 ") different to num rows ("+
									 Format::toString((*losses)[i].size()) +
									 ") in losses matrix at column "+Format::toString(i));
			}
		}

                
        // check that strike grid is increasing
       
        for(i = 1;i<strikes->size();i++)
        {
            if((*strikes)[i] < (*strikes)[i-1]) 
            {
                throw ModelException("Strike grid not increasing");
            }
        }

		// check that dates are increasing
		DateTime::ensureStrictlyIncreasing(*dates,
                                           string("dates are not strictly increasing"),
                                           true);

        if(RR >1 || RR < 0)
		{
			throw ModelException("Recovery rate ("+Format::toString(RR) +") is out of bounds [0,1]");
		}

    } catch(exception & e)
    {
        throw ModelException(e,method);
    }
}

/***************************************************************************************
public constructor 
***************************************************************************************/
ExpectedLossSurface::ExpectedLossSurface(DateTimeArraySP dates,
                                         DoubleArraySP strikes,
                                         DoubleArrayArraySP losses,
                                         double RR) :

    CObject(TYPE),
    dates(dates),
    strikes(strikes),
    losses(losses),
    RR(RR)
{
    static const string method ="ExpectedLossSurface::ExpectedLossSurface";
    try
    {  
		validatePop2Object();
    } catch(exception & e)
    {
        throw ModelException(e,method);
    }

}

/*******************************************************************************************
 get expected losses for base tranches 
******************************************************************************************/
double ExpectedLossSurface::getBaseEL(double          strike, 
                                      const DateTime& date) const
{
    static const string method = "ExpectedLossSurface::getBaseEL";
    try{
                

        // TODO should check that dates are increasing
        int tLow = date.findLower((*dates));

        if(tLow<0) throw ModelException("date < (*dates)[0]");                  
        int tHigh = tLow; // first time index after date (if date is beyond then use constant extrapolation)
        if(date != (*dates)[tLow] && (tLow+1) < dates->size()) tHigh = tLow+1; 


        double tLowEL = strikeInterp(strike, tLow);
        double tHighEL = strikeInterp(strike, tHigh);
                
        // linear interpolation in time dimension ---------------------------------------
        double loss = tLowEL;
        double dt = (*dates)[tLow].yearFrac((*dates)[tHigh]);
        if(dt > TINY)
        {
            loss += (tHighEL-tLowEL)*((*dates)[tLow].yearFrac(date))/dt;
        }

        return loss;
    }
    catch (exception& e){
        throw ModelException(e,method);
    }
}


/********************************************************************************************************

  get contingent leg effective curve for given strikes from expected loss surface

********************************************************************************************************/
DoubleArraySP ExpectedLossSurface::getEffectiveCurve(double lowStrike, double highStrike) const
{
        
    static const string method = "ExpectedLossSurface::getEffectiveCurve";
    try{
                
        /** validation */
        if(lowStrike >1 || lowStrike < 0) throw ModelException("Low strike is out of bounds [0,1]");
        if(highStrike >1 || highStrike < 0) throw ModelException("High strike is out of bounds [0,1]");
        if(lowStrike>=highStrike) throw ModelException("Low Strike >= high Strike");

        int numPoints = dates->size();
                
        DoubleArraySP effCurve(new DoubleArray(numPoints));
        for(int i = 0;i < numPoints; i++)
        {
            double loss = strikeInterp(highStrike,i) - strikeInterp(lowStrike,i);
            (*effCurve)[i] = 1 - loss/(highStrike-lowStrike);
        }
        return effCurve;
    }
    catch (exception& e){
        throw ModelException(e,method);
    }
}

/********************************************************************************************************

  get fee leg effective curve for given strikes from expected loss surface

********************************************************************************************************/
DoubleArraySP ExpectedLossSurface::getFeeLegEffectiveCurve(double lowStrike, double highStrike) const
{
        
    static const string method = "ExpectedLossSurface::getFeeLegEffectiveCurve";
    try{
                
        if(RR > 0)
        {
            /** validation */
            if(lowStrike >1 || lowStrike < 0) throw ModelException("Low strike is out of bounds [0,1]");
            if(highStrike >1 || highStrike < 0) throw ModelException("High strike is out of bounds [0,1]");
            if(lowStrike>=highStrike) throw ModelException("Low Strike >= high Strike");

            int numPoints = dates->size();
                        
            DoubleArraySP effCurve(new DoubleArray(numPoints));
            for(int i = 0;i < numPoints; i++)
            {
                /** unscaled expected loss of tranche [lowStrike, highStrike] */
                double loss = strikeInterp(highStrike,i) - strikeInterp(lowStrike,i);
                                
                /** expected loss for tranche [0,A] where A = (1-RR)/RR*(1-low) */
                double baseELLow = strikeInterp(Maths::min(1.0,(1-RR)/RR*(1-lowStrike)),i);

                /** expected loss for tranche [0,B] where B = (1-RR)/RR*(1-high) */
                double baseELHigh = strikeInterp(Maths::min(1.0,(1-RR)/RR*(1-highStrike)),i);
                                
                                
                (*effCurve)[i] = feeLegEffCurve(
                    loss,
                    highStrike - lowStrike,  // size of tranche
                    baseELLow,
                    baseELHigh,
                    RR
                    );
            }
            return effCurve;
        }
        else // for 0 RR cont leg and fee leg eff curve are the same
        {
            return getEffectiveCurve(lowStrike,highStrike);
        }
    }
    catch (exception& e){
        throw ModelException(e,method);
    }
}

/*****************************************************************************************************
fee leg effCurve

 no validation
***************************************************************************************************/
double ExpectedLossSurface::feeLegEffCurve(
    double expectedLoss,  // unscaled
    double trancheSize,
    double baseELLow,
    double baseELHigh,
    double R
    )
{
    double out = 1;
    out -= expectedLoss/trancheSize;
    out -= (baseELLow-baseELHigh)*R/(1-R)/trancheSize;
    return out;
}
/*******************************************************************************************

 linear interpolation of exp loss by strike 

*******************************************************************************************/
double ExpectedLossSurface::strikeInterp(double strike, int timeIndex) const
{
    static const string method = "ExpectedLossSurface::strikeInterp";
    try
    {
		/** output loss */
		double loss;

        if(strike < 0.0 || strike > 1.0)
        {
            throw ModelException("strike ("+Format::toString(strike)+") is out of bounds [0,1]");
        }
		if(timeIndex < 0 || timeIndex >= dates->size())
		{
			throw ModelException("timeIndex ("+Format::toString(timeIndex)+
				") is out of bounds [0, " + Format::toString(dates->size()-1) + "]");
		}

		if(strike >= strikes->back()) // lie between last strike and 1. Assume constant EL extrapolation
		{
			// return last expected loss
			loss = (*losses)[strikes->size()-1][timeIndex];	
		}
		else if(strike < strikes->front()) // lie between first strike an 0. Loss at 0 is 0.
		{
			// interpolate
            loss = (*losses)[0][timeIndex]/(*strikes)[0]*strike;
		}
		else // lie in between lowest and highest strike
		{
			/** lowest index such that strikes[highI] >= strike */
			int highI = 0;
			while(highI < strikes->size() && (*strikes)[highI] < strike)
			{
				highI++;
			}
                
			if(highI>=strikes->size())
			{
				throw ModelException("strike ("+Format::toString(strike)+") too high"); // this can't happen
			}

			loss =  (*losses)[highI][timeIndex];
                
			// check if strike lies on strikes grid in which case we do nothing
			// otherwise interploate
			if((*strikes)[highI] -strike > TINY)
			{
				if(highI-1< 0 ) throw ModelException("highI too low"); // this can't happen

				/** distance between strikes */
				double ds = (*strikes)[highI]-(*strikes)[highI-1];

				loss -= ((*losses)[highI][timeIndex] - (*losses)[highI-1][timeIndex])*((*strikes)[highI]-strike)/ds;
			}
		}
        return loss; 
    }
    catch (exception & e)
    {
        throw ModelException(e, method);
    }

}


/*****************************************************************************************
Get
***************************************************************************************/
DoubleArrayArraySP ExpectedLossSurface::getLosses() const
{
	// don't return const array because this simplifies things elsewhere
	// hence we copy the internal field
    return DoubleArrayArraySP(new DoubleArrayArray(*losses));
}

DateTimeArrayConstSP ExpectedLossSurface::getDates() const
{
    return dates;
}

DoubleArrayConstSP ExpectedLossSurface::getStrikes() const
{
	return strikes;
}
double ExpectedLossSurface::getRecovery() const
{
	return RR;
}

DRLIB_END_NAMESPACE
