//----------------------------------------------------------------------------
//
//   Group       : EDG Derivatives Research
//
//   Filename    : VolProcessedDVFParam.hpp
//
//   Description : Processed DVF parameterised vols. 
//                 Previously, was part of VolParam.hpp
//
//   Author      : Regis Guichard
//
//   Date        : 02 Mai 2001
//
//
//----------------------------------------------------------------------------

#ifndef EDR_VOL_PROCESSED_DVF_HPP
#define EDR_VOL_PROCESSED_DVF_HPP

#include "edginc/VolParam.hpp"
#include "edginc/VolProcessedDVF.hpp"
#include "edginc/VolRequestDVF.hpp"

DRLIB_BEGIN_NAMESPACE

/** Default implementation CVolProcessedDVF for parameterised vols */
class MARKET_DLL CVolProcessedDVFParam: public CVolProcessedDVF{
    friend class IVolCalculatorDVF;
public:
    static CClassConstSP const TYPE;
    static void load(CClassSP& clazz);

    ~CVolProcessedDVFParam();

    /** create a volatility calculator using a cache of forward values */
    virtual CVolProcessedDVF::IVolCalculator* CreateVolCalculator(
        const DateTimeArray&   maturity, // the time axis of the lattice
        bool                   isIntraDayInterp = true) const;
    
    virtual void CalcLocVol(
        CLatticeDouble*     strikes,
        DateTimeArray&      maturities, // the time axis of the lattice above
        CLatticeDouble*     locVol,
        bool                isIntraDayInterp=true) const;
  
    virtual void CalcLocVar(
        CLatticeDouble*     strikes,
        DateTimeArray&      maturities, // the time axis of the lattice above
        CLatticeDouble*     locVar,
        bool                isIntraDayInterp=true) const;

    void CalcLocVar(const DateTimeArray&  maturity,
                    const double*         strike,
                    const int             NStr,
                    vector<double>&       locVar,
                    bool		  isIntraDayInterp=true) const;
    
    /** necessary, or at least should be private */
    virtual void CalcLocVol(const DateTimeArray&  maturity,
                            const double*         strike,
                            const int             NStr,
                            vector<double>&       locVol,
                            bool                  isIntraDayInterp=true) const;
    
    virtual string getName() const;
    /** calculates the trading time between two dates */
    virtual double calcTradingTime(const DateTime &date1, 
                                   const DateTime &date2) const;
    /** retieve time measure for the vol */
    virtual TimeMetricConstSP GetTimeMetric()const;

    /** Given a maturity date and a strike, returns the implied vol */
    virtual double computeImpVol(const DateTime& maturity,
                                 double          strike) const;

    /** Given a maturity date and a strike, returns the local vol */ 
    virtual double CalcLocVol(const DateTime&  maturity,
                              double    strike,
							  bool		isIntraDayInterp=true) const;

    //// construct a vol processed DVF for parameterised vol
    CVolProcessedDVFParam(const CVolBase*         vol,
                          const CVolParamConstSP&  volParam,
                          const CVolRequestDVF*    volRequest, 
                          const CAsset*            asset,
                          const VolSurface*        VolSurf);


    //// construct a vol processed DVF for parameterised vol when struck
   CVolProcessedDVFParam(const CVolBase*         vol,
                          const CVolParamConstSP&  volParam,
                          const CVolRequestDVF*    volRequest, 
                          const CAsset*            asset,
                          const VolSurface*        VolSurf,
                          const FXAsset*           fxAsset,
                          const Correlation*       eqFXCorr); 

protected:
   CVolProcessedDVFParam(const CClassConstSP&    clazz, 
                          const CVolBase*          vol,
                          const CVolParamConstSP&  volParam,
                          const CVolRequestDVF*    volRequest, 
                          const CAsset*            asset,
                          const VolSurface*        VolSurf,
                          const FXAsset*           fxAsset = 0,
                          const Correlation*       eqFXCorr = 0);  
	// for child class
   CVolProcessedDVFParam(const CClassConstSP&    clazz) : CVolProcessedDVF(clazz){};


    /** Given a  maturity dates and a lettice of strikes,
        returns the implied vol, its time derivative 
        and its first and second strike derivatives. 
        Time derivs are computed by tweaking 
        using timeTweakUnscaled. 
        Strike derivs are calculated by tweaking using strikeTweakUnscaled. */
    virtual void computeLocVol(const DateTime&  maturity,
                               const double*         strike,
                               const int             NStr,
                               vector<SImpV>&        impV,
                               vector<double>&           locV,
							   bool				    isIntraDayInterp) const;  

    /** Given an array of maturity dates and a lettice of strikes,
        returns the implied vol, its time deriv its first and 
        second strike derivs as well as the local vol. 
        Growth rate is computed by tweaking using timeTweakUnscaled .
        Time and strike derivs are computed  by tweaking 
        using timeTweakUnscaled and strikeTweakUnscaled, respectively. */
    void computeLocVNextStep(const DateTimeArray&  maturity,
                             const double*         strike,
                             const int             NStr,
                             vector<SImpV>&        impV,
                             vector<double>&       locV,
                             bool                  reqVar,
                             bool                  isIntraDayInterp) const;

    /** Given an array of maturity dates and a lettice of strikes,
        returns the implied vol, its time deriv its first and 
        second strike derivs as well as the local vol at (t[i] + t[i+1])/2.
        Time derivs computed by tweaking using timeTweakUnscaledand. 
        Grwoth rate is computed by finite diffs
        using next time step. Strike derivs are calculated by tweaking 
        using strikeTweakUnscaled.*/
    void computeLocVarNextStep(const DateTimeArray&  maturity,
                               const CLatticeDouble&   strike,
                               ImpVSurf&             impV,
                               CLatticeDouble&       locV) const;

    virtual void computeImpVolAndDerivsNextStep(const double*         strike,
                                                const int             NStr,
                                                const DateTime&       maturity,
                                                const DateTime&       maturityup,
                                                const double          timeTweakyrs,
                                                vector<SImpV>&        impV) const;

    static const double SMALL_YEARS_TO_MATURITY;
    static const double MIN_STRIKE_TWEAK;
    static const double MIN_TIME_TWEAK;
    static const double MIN_NEXT_DATE;

private:
    void Init(  const CVolBase*         vol,
                const CVolParamConstSP& volParam,
                const CVolRequestDVF* volRequest, 
                const CAsset* asset,
                const VolSurface* VolSurf,
                const FXAsset*   fxAsset = 0, 
                const  Correlation* eqFXCorr = 0);

    inline DateTime tweakDate(const DateTime& date, double& timeTweak) const;
 
    void computeImpVol(const CLatticeDouble&    strikes,
                       const DateTimeArray&     maturities,
                       CLatticeDouble&          impV) const;

    // result vol^2 is in v2. returns true if v2 update, false if v2 is not computed
    virtual bool computeLocV2(double  yrsToMat,
                              double  strike,
                              double  forward,
                              double  growthRate,
                              SImpV&  impV,
                              double  &v2) const;

    CVolBaseConstSP         myVol; // $unregistered
    CVolParamConstSP        myParamVol; // $unregistered
    CVolRequestDVFConstSP   myVolRequest; // $unregistered
    AssetConstSP            asset; // $unregistered
    VolSurfaceConstSP       myVolSurf; // $unregistered
    // the two below are potentially NULL
    FXAssetConstSP          myFXAsset; // $unregistered
    CorrelationConstSP      myEqFXCorr; // $unregistered
    CVolParam::FwdStart     fwdStart; // $unregistered

    DateTime startDate; // $unregistered
    
    // used to cache of forwad values and growth rates
    DoubleArray forwardmid; // $unregistered
    DoubleArray growthRate; // $unregistered
    DateTimeArray maturity; // $unregistered
    
};

typedef smartConstPtr<CVolProcessedDVFParam> CVolProcessedDVFParamSP;

DRLIB_END_NAMESPACE
#endif 
