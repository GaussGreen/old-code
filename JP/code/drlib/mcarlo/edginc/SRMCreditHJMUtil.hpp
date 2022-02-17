//----------------------------------------------------------------------------
//
//   Group       : QR Cross Asset
//
//   Filename    : SRMCreditHJMUtil.hpp
//
//   Description : Derived SRMCreditUtil class for Markovian HJM / RS model
//
//
//----------------------------------------------------------------------------

#ifndef SRMCREDITHJMUTIL_HPP
#define SRMCREDITHJMUTIL_HPP

#include "edginc/SRMCreditUtil.hpp"
#include "edginc/CDSParSpreads.hpp"
#include "edginc/VirtualDestructorBase.hpp"
#include "edginc/DECLARE.hpp"

#include <cassert>

DRLIB_BEGIN_NAMESPACE

class SRMCreditHJMUtil : public SRMCreditUtil
{
public:

    /** constructor */
    SRMCreditHJMUtil(
                 const DateTime&       baseDate,
                 const string&         smileParams,
                 ICDSParSpreadsConstSP stochCDSCurve); 

    /** destructor */
    virtual ~SRMCreditHJMUtil();

    /** computes 'KFactors' on srmUtil 'dates' 
        kFactor param is initialised to correct length. */
    void computeKFactor(vector<double>& kFactor) const;

    /** computes 'GFactors' on srmUtil 'dates' 
        gFactor param  is initialised to correct length. */
    void computeGFactor(vector<double>& gFactor) const;
        
    /** returns GfactorCR */
    const vector<double>& getGfactorCR() const;

    /** returns the spot vols */
    const vector<double>& getSpotVols() const;

    /** returns beta */
    double getBeta() const;

    /** returns qRight */
    double getQRight() const;

    /** returns qLeft */
    double getQLeft() const;

    /** returns fwdShift */
    double getFwdShift() const;

    /** for calculating 'lbar' when moving FROM a date in dates whose index
        is given by futureDatesIndex */
    double lBar(int futureDatesIndex) const;

    /** for calculating 'Lbar' on supplied date */
    double lBar(const DateTime& date) const;

    /** computes KFactor between the 2 dates supplied */
    double kFactor(const DateTime& dateFrom,
                   const DateTime& dateTo) const; 

    /** computes GFactor between the 2 dates supplied */
    double gFactor(const DateTime& dateFrom,
                   const DateTime& dateTo) const; 

    /** Utility routine to calculate the cutoff rate and populate
        the cutoff rate arrays */
    void calcEffRateLimit(double            NbSigmasMax,    // (I) Number or sigmas to cut at
                          double            NbSigmasMin,    // (I) Number or sigmas to cut at
                          vector<double>&   MaxRate,        // (O) Cutoff forward rates
                          vector<double>&   MinRate) const; // (O) Cutoff forward rates

    /** this function is called on the edfForwardDates dates */
    void populatePartialIntCR(const DateTime& today,
                              const DateTimeArray& myDatesIn, 
                              vector<double>& partialIntegralOut);

    /** this function is called on the requestedEDFdates */
    void populatePartialZeta(const DateTime& today,
                             const DateTimeArray& myDates,
                             vector<double>& zeta);

    /** finishes initialization */
    virtual void setTimeLine(DateTimeArrayConstSP simDates);

private:

    /** populates FwdIntensity field on extendedTimeLine. 
    requires 'dates' field */
    void calcFwdIntensity();

    /** Calculates and populates the GfactorCR[].
    This array is used in Gfactor() to streamline the calculation. From
    irdiffuse::PopulateGfactorCR */
    void populateGfactorCR();

    /** Utility routine to calculate the cutoff rate and populate the cutoff rate arrays 
        see crdiffuse::CalcEffRateLimitCR for details */
    void calcEffRateLimitCR(double          NbSigmasMax,    // (I) Number or sigmas to cut at
                            double          NbSigmasMin,    // (I) Number or sigmas to cut at
                            vector<double>& MaxRate,        // (O) Cutoff forward rates
                            vector<double>& MinRate) const; // (O) Cutoff forward rates

    double                  beta; 
    vector<double>          GfactorCR;          /** across dates (size = extendedTimeLine.size)
                                                    populated in SRMCreditHJMUtil::populateGfactorCR
                                                    used in computeGFactor as well as gFactor */
    vector<double>          FwdIntensity;       // on extendedTimeLine
    vector<double>          SpotVol;            // LOGNORMAL vols, same size as FwdIntensity
    double                  qLeft;
    double                  qRight;
    double                  fwdShift;
};

/** declare smartPtr versions */
DECLARE(SRMCreditHJMUtil);

DRLIB_END_NAMESPACE

#endif // SRMCREDITHJMUTIL_HPP
