//----------------------------------------------------------------------------
//
//   Group       : QR Credit Hybrids
//
//   Filename    : QMCDetermRates.hpp
//
//   Description : deterministic IRs asset
//
//----------------------------------------------------------------------------

#ifndef QMCDETERMRATES_HPP
#define QMCDETERMRATES_HPP

#include "edginc/DECLARE.hpp"
#include "edginc/IQMCDiffusibleAsset.hpp"
#include "edginc/QMCHelperBoundedDiffusion.hpp"

DRLIB_BEGIN_NAMESPACE

/******************************************************************************
*******************************************************************************
** Deterministic IRs asset                                                   **
*******************************************************************************
******************************************************************************/
class MCARLO_DLL QMCDetermRates : public IQMCDiffusibleInterestRate
{
public:
///////////////////////////////////////////////////////////////////////////////
// CONSTRUCTION AND DESTRUCTION ///////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////

    /** 
      * Constructor. It is empty **********************************************
    **/
    QMCDetermRates();
    /** 
      * Quasi-constructor. Does all work **************************************
    **/
    void setQMCDetermRates(DateTime today, IYieldCurveConstSP gYC);
    /** 
      * Destructor ************************************************************
    **/
    virtual	~QMCDetermRates();

///////////////////////////////////////////////////////////////////////////////
// MONTE-CARLO METHODS ////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////

    /** 
      * Finalize the timelines, allocate necessary memory *********************
    **/
    virtual void finalize(DateTimeArrayConstSP allDates);
    /** 
      * Generate path across all dates ****************************************
    **/
    virtual void generatePath(IQMCRNGManagerSP rngMgr, QMCStrataConstSP /*strata*/) ;

///////////////////////////////////////////////////////////////////////////////
// COMPUTATION ////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////

    /** 
      * Returns the DF 4 date getTimeLogic()->getDFDates()[measurementDateIdx]*
    **/
    virtual double  getDiscFactor(SpotIdx measurementDateIdx);
    /** 
      * Returns forward DF for dates 
      *                        getTimeLogic()->getDFDates()[measurementDateIdx]
      *                    and getTimeLogic()->getDFDates()[futureDateIdx] ****
    **/
    virtual double  getExpectedDiscFactor(size_t ycIdx,
                                          FwdIdx measurementDateIdx,
                                          FwdIdx futureDateIdx);
    /** 
      * Returns the natural log of forward DF for dates 
      *                        getTimeLogic()->getDFDates()[measurementDateIdx]
      *                    and getTimeLogic()->getDFDates()[futureDateIdx] ****
    **/
    virtual double getLnExpectedDiscFactor(size_t ycIdx,
                                           FwdIdx measurementDateIdx,
                                           FwdIdx futureDateIdx);
    /** 
      * Moment-matching related stuff.
      * Here implementation is returns ln(getDiscFactor(...)) *****************
    **/
    virtual double getOriginalLnDiscFactor(SpotIdx futureDateIdx);
    /** 
      * Moment-matching related stuff.
      * Here implementation is identical to getLnExpectedDiscFactor(...) ******
    **/
    virtual double getOriginalLnExpectedDiscFactor(size_t ycIdx,
                                                   FwdIdx measurementDateIdx,
                                                   FwdIdx futureDateIdx) ;

///////////////////////////////////////////////////////////////////////////////
// ADVANCED ///////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////

    /** 
      * Getting indicator of the YC type. 
      * Returns 0 because here we do not need complex manipulations 
      * with the different YC flavors
      * since it is quick to compute the exact result *************************
    **/
    virtual size_t registerYCFlavor(IYieldCurveConstSP );
    /** 
      * Getting indicator of the Discount YC. 
      * Returns 0 because here we do not need complex manipulations 
      * with the different YC flavors
      * since it is quick to compute the exact result *************************
    **/
    virtual size_t getDiscYCIdx(void);
    /** 
      * Getting indicator of the Diffusion YC. 
      * Returns 0 because here we do not need complex manipulations 
      * with the different YC flavors
      * since it is quick to compute the exact result *************************
    **/
    virtual size_t getDiffYCIdx(void);

///////////////////////////////////////////////////////////////////////////////
// UNDER CONSTRUCTION /////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////

    /** 
      * Allow access to domLnMONEY variable. Throws an error. 
      * YV: I do not understand what this function means **********************
    **/
    virtual vector<double>& getDomLnMONEY();
    /** 
      * Returns FX rates. Check that this method should really be here ********
    **/
    virtual const vector<double>* getSigmaFX();

    /** Returns full FX path. */
    virtual const vector<double>* getSpotFXFullPath();

    /**
      * Getting the simulation start date = 'today'. 
      * Currently throws an error since - as far as I have understood -
      * this method is likely to become as the one of TimeLogic.
      * Hence, it might be retired from the asset.
      * Same for the bound on the latest needed date **************************
    **/
    virtual DateTime getBaseDate() ;

    /** getting information about asset specific diffusion bounds */
    virtual QMCHelperBoundedDiffusion * getDiffusionBound() { return &pBoundedDiffusion; }

private:

///////////////////////////////////////////////////////////////////////////////
// AUXILIARY //////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////

    void operator=(const QMCDetermRates & g2Clone); // not defined
    QMCDetermRates(const QMCDetermRates & g2Clone); // not defined

    /** 
      * Fills the array of discount factors (gResult) given dates *************
    **/
    void finalizeForGivenDates(const DateTimeArray & gDates,
                                     DoubleArray   & gResult);
////////////////////////////////////////////////////////////////////////////////
// FIELDS //////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////

    /** 
      * This wrapped market object contains all numbers ***********************
    **/
    IYieldCurveConstSP yc;
    /** 
      * Discount factors (== present values, ZCs) 
      * pulled out of the YieldCurveWrapper in 'finalize()'.
      * It is the info source for state variables. ****************************
    **/
    DoubleArray       cachedSpotIdxDFs; // for getTimeLogic()->getDFDates()
    DoubleArray       cachedFwdIdxDFs;  // for getTimeLogic()->getFwdEDFDates()

    DateTime          baseDate;
    QMCHelperBoundedDiffusion pBoundedDiffusion;
};
// declare smartPtr versions //////////////////////////////////////////////////
DECLARE(QMCDetermRates);
///////////////////////////////////////////////////////////////////////////////
DRLIB_END_NAMESPACE
#endif // QMCDETERMRATES_HPP
