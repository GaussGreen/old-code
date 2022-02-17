//----------------------------------------------------------------------------
//
//   Group       : Cross Asset QR
//
//   Filename    : SRMBasisSpreadHJMUtil.hpp
//
//   Description : Rates SP model diffusion
//
//   Date        :
//
//----------------------------------------------------------------------------

#ifndef QLIB_SRM_SP_UTIL_HPP
#define QLIB_SRM_SP_UTIL_HPP

#include "edginc/DECLARE.hpp"
#include "edginc/SRMRatesUtil.hpp"
#include "edginc/BasisIndexCurve.hpp"
#include <cassert>


DRLIB_BEGIN_NAMESPACE


/** Helper class for models to get data out of market cache */
class SRMBasisSpreadUtil:   public virtual VirtualDestructorBase
{
public:

    ///// constructs and populates SRMCRUtil object
    SRMBasisSpreadUtil(const DateTime&       _baseDate,
                       IBasisIndexCurveConstSP _basisIndexCurve);

    ~SRMBasisSpreadUtil();

    /** Utility routine to calculate the cutoff rate and populate
        the cutoff rate arrays */
    void calcEffRateLimit(double            NbSigmasMax,    // (I) Number or sigmas to cut at
                          double            NbSigmasMin,    // (I) Number or sigmas to cut at
                          vector<double>&   MaxRate,        // (O) Cutoff forward spreads
                          vector<double>&   MinRate) const; // (O) Cutoff forward spreads

    void setTimeLine(DateTimeArrayConstSP allDates); // finishes initialization

    vector<double> getSpotSpreads() const { return fwdSpreads; }; // for the TimeLine points
    vector<double> getSpotVols() const { return spotVols; }// for the TimeLine points
    vector<double> computeSqrtYearFrac() const; // for the TimeLine points

    double hFactor(const DateTime& d1, const DateTime& d2) const;

    vector<double> computePartialH(const DateTime& today, const DateTimeArray& forwardDates) const;
    vector<double> computeFwdParSpreads(const DateTimeArray& forwardDates) const;

    inline double getQ() const{ return q; }

    /** returns the basis curve */
    IBasisIndexCurveConstSP getBasisCurve() const { return basisCurve; }

private:

    double calcFwdSpread(const DateTime& resetDate) const;



    double                  q;
    double                  beta;
    DateTime                baseDate;           // today
    DateTimeArrayConstSP    dates;              // sim start + dates after sim start

    IBasisIndexCurveConstSP basisCurve;

    vector<double>          spotVols;
    vector<double>          fwdSpreads;

    bool                    initialized;



//-----------------------
// uncertain if we'd need anything below this line
    /* some vol parameters ... */
    double getBeta() const { assert(initialized); return beta; }

    /** returns all the simulation dates excluding today */
    const DateTimeArray& getSimDates() const { assert(initialized); return *dates; }
};
DECLARE(SRMBasisSpreadUtil);

DRLIB_END_NAMESPACE
#endif


