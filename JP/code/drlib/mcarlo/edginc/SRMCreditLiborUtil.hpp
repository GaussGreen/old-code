//----------------------------------------------------------------------------
//
//   Group       : QR Cross Asset
//
//   Filename    : SRMCreditLiborUtil.hpp
//
//   Description : Derived SRMCreditUtil class for CreditLibor model
//
//
//----------------------------------------------------------------------------

#ifndef SRMCREDITLIBORUTIL_HPP
#define SRMCREDITLIBORUTIL_HPP

#include "edginc/SRMCreditUtil.hpp"
#include "edginc/CDSParSpreads.hpp"
#include "edginc/VirtualDestructorBase.hpp"
#include "edginc/DECLARE.hpp"

#include <cassert>


DRLIB_BEGIN_NAMESPACE


class SRMCreditLiborUtil : public SRMCreditUtil
{
public:
    // constructor
    SRMCreditLiborUtil(const DateTime&       baseDate,
                 const string&         smileParams,
                 ICDSParSpreadsConstSP stochCDSCurve); 
           
    // destructor
    virtual ~SRMCreditLiborUtil(){};

    // accessors
	vector<double> getResetDates() {return mResetDates;}
	vector<double> getAccruedFrac() {return mAccruedFrac;}
	vector<double> getSpotVols() {return mSpotVols;}
	vector<double> getInitialIntensities() {return mInitialIntensities;}
	vector<double> getSpotSurvProb() {return mSpotSurvProb;}
	double         getYearFracToFirstReset(){return mYearFracToFirstReset;}
	double		   getNumInt(){return mNumInt;}	

    double getQRight() const{ return qRight;}
    double getQLeft() const{ return qLeft;}

     // finishes initialization
    virtual void setTimeLine(DateTimeArrayConstSP simDates);

private:
	int						mNumInt,					// number of intensities, i.e. state variables  
							mIntensityInterval; 	    // in months


	double                  qLeft,						// smile parameter
                            qRight,                     // smile parameter
	                        mLastSimTime,               // in years
							mYearFracToFirstReset;      // time elapsed from today to first reset time.

    vector<double>			mResetDates,
							mAccruedFrac,
							mSpotVols,
							mInitialIntensities,		
							mSpotSurvProb;				// Q(0,T_i), T_i = discreteIntensityDate
};

//declare smartPtr versions
DECLARE(SRMCreditLiborUtil);

DRLIB_END_NAMESPACE

#endif // SRMCREDITLIBORUTIL_HPP
