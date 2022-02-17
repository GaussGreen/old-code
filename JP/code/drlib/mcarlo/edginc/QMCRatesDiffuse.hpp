//----------------------------------------------------------------------------
//
//   Group       : QR
//
//   Filename    : QMCRatesDiffuse.hpp (code previously in MCPathConfigSRM.cpp)
//
//   Description : base class for SRM IR diffusions
//
//
//----------------------------------------------------------------------------

#ifndef QMCRATESDIFF_HPP
#define QMCRATESDIFF_HPP

#include <cstdio>
#include <cassert>
#include <set>

#include "edginc/FORWARD_DECLARE.hpp"
#include "edginc/IQMCDiffusibleAsset.hpp"
#include "edginc/Maths.hpp"

#include "edginc/QMCHelperBoundedDiffusion.hpp"
#include "edginc/QMCFXBaseDiffuse.hpp"


DRLIB_BEGIN_NAMESPACE

class IQMCHelperTimeLogic;
FORWARD_DECLARE(QMCRatesUtil);

//// base class for low level IR path generator
class QMCRatesDiffuse: public IQMCDiffusibleInterestRate
{

public:

    virtual ~QMCRatesDiffuse();

	//// constructor
    QMCRatesDiffuse();

    /** called after origSVol externally recalibrated via ICE */
    virtual void recalibrate(QMCRatesUtilSP thisSRMRatesUtilSP) = 0;

	// ir->detachFX() will be used to detach shared fx pointer
    //       for resolving the circular dependency of pointers
    // (!!!) not using detachFX() will result in a memory leak
    void detachFX() { fx = QMCFXBaseDiffuseSP(NULL); }

    void        setFXBase(QMCFXBaseDiffuseSP _fx);/* { fx=_fx; }*/
    QMCFXBaseDiffuseSP getFXBase() const { return fx; };
    const vector<double>* getSigmaFX();
    const vector<double>* getSpotFXFullPath();

    /** Allow access to domLnMONEY variable */
    virtual vector<double>& getDomLnMONEY(){ return domLnMONEY; }

    /** Allow access to sigmaR variable */
    vector<double>& getSigmaR(){ return sigmaR; }

    /** allows modification to origSVol vector */
    vector<double>& getOrigSVol(){ return origSVol; }

    /** getting the simulation start date */
    virtual DateTime getBaseDate() { return today; }

    /** this is for efficient access to DiscFactor values, but optional */
    virtual double* getInternalDFArray(); // this asset stores DFs

    double getDiscFactor(SpotIdx i) {
        return df[i];
    }

	size_t registerYCFlavor(IYieldCurveConstSP yc); // returns index in the internal YCForwardsDB

    virtual size_t getDiscYCIdx(void) = 0; // returns Idx of Discount YC
    virtual size_t getDiffYCIdx(void) = 0; // returns Idx of Diffusion YC

    virtual QMCRatesUtilSP getRatesUtil(void) = 0;

    double  getSqrtYearFrac(size_t i) { return sqrtYearFrac[i];}
    virtual QMCHelperBoundedDiffusion* getDiffusionBound() { return & diffBound;}
 
    virtual double getOriginalLnDiscFactor(SpotIdx spotDateIdx)
    {
        return originalDFs.at(spotDateIdx);
    }

    virtual double getOriginalLnExpectedDiscFactor(size_t ycIdx, FwdIdx measurementDateIdx, FwdIdx futureDateIdx)
    {
        return getYCForward(ycIdx, measurementDateIdx, futureDateIdx);
    }

protected:

    int             randomIndex; // index into random numbers
    QMCFXBaseDiffuseSP     fx;          // link to FX, may be null
    vector<double>  domLnMONEY;  // per time pt: populated when explicitly requested to [Nall] or [0]
    vector<double>  sigmaR;  // per time pt: populated when needed for credit [Nall]
    //// cached variables
    vector<int>     dfIndexes; // indexes for when we save discount factors [Ndf]
    vector<int>     expDFIndexes; // indexes for when we compute expected  [Nedf]
                                  // discount factors
    vector<double>  logDiscFactor; // determinstic from discount zero curve[Nall]
    vector<double>  effRSVol; // for zero q, effR * svol[Nall]

    vector<double>  svol; /* contains SpotVol initially. If zeroQ is false
                             then processed in calcSigmaRParams */ // [Nall]
    vector<double>  MaxEffRateIR; // the cutoff limit at each ZDate time point [Nall]
    vector<double>  MinEffRateIR; // the cutoff limit at each ZDate time point[Nall]
    vector<double>  sqrtYearFrac; // sqrt of yearFrac between sim dates[Nall]
    vector<double>  rbar;         // [Nall]

    DateTime        today;     /* number of historic sim dates */
    // we used to have
    // int todayIdx = today.findLower(allDates); where allDates = pastDates + simDates;
    int             todayIndex; // index of today in dfRequestedDates array

    // computed discount factors at dfIndexes
    vector<double>  df; // [Ndf]
    vector<double>  origSVol; // unmodified svol [Nall]
    bool            calibrateAgainstSwaptionVols;
    int             todayIdx; // index in allDates

    bool            initialized;
    bool            assetDatesFixed;

	// Calculates fwd(i, j) given SV's YC index.
	double getYCForward(size_t idx, FwdIdx i, FwdIdx j) {
       return ycForwardsDB[idx].second[j] - ycForwardsDB[idx].second[i];
    }

    SpotIdx         getFwdIdx2EdfIdx(FwdIdx idx); // translates FwdIdx into RequestedIdx

    bool           saveDomLnMoney;
    bool           saveSigmaR;
    double         NbSigmasMax;
    double         NbSigmasMin;

	YCForwardsDB  ycForwardsDB;
    vector<double> originalDFs; // [Ndf], used only by moment matching as reference values

    DateTimeArrayConstSP timelineSP; // remember our timeline
    QMCIRBoundedDiffusion diffBound;
};

//declare smartPtr versions
DECLARE(QMCRatesDiffuse);

DRLIB_END_NAMESPACE
#endif // QMCRATESDIFF_HPP

