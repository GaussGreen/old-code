//----------------------------------------------------------------------------
//
//   Group       : QR Cross Asset
//
//   Filename    : SRMFXUtil.hpp
//
//   Description : Utility class for FX SRM
//
//
//----------------------------------------------------------------------------

#ifndef EDR_SRMFXUTIL_HPP
#define EDR_SRMFXUTIL_HPP

#include "edginc/VirtualDestructorBase.hpp"
#include "edginc/DECLARE.hpp"

#include "edginc/SRMRatesUtil.hpp"
#include "edginc/SRMUtil.hpp"
#include "edginc/SRMFXVol.hpp"

DRLIB_BEGIN_NAMESPACE    

class MCARLO_DLL SRMFXUtil : public virtual VirtualDestructorBase
{
public:
    ~SRMFXUtil();

    SRMFXUtil(
        SRMRatesUtilSP       forIR,
        SRMRatesUtilSP       domIR,
        FXAssetConstSP       fxAsset,    
        double               InpCorrFxFor, // correlation between fx and foreign IR
        double               InpCorrFxDom, // correlation between fx and domestic IR 
        double               InpCorrDomFor,/* correlation between domestic
                                           and foreign IR */
        const string&        VolBootstrapMode, // choice for bootstrapping vol
        double               FxCutOffLevel,
        const string&        volType,
        bool                 forceTimeOfDay = false);  /* true - set of times of
                                                       day to today's */
    
    /** Returns the dates that were used in the calibration - essentially
        option expiry dates */
    const DateTimeArray& calibratedSpotVolDates();

    /** Returns the spot vols along the getCalibratedSpotVolDates() */
    const vector<double>& calibratedSpotVols();

    //// Returns name of fx asset
    const string& getName() const;

    // accessors
    SRMRatesUtilSP getDomIR() const { return domIR; }
    /** returns bootstrapped spot vol */
    // FIXME remove calls to initialize() after we fix util/timeline logic
    const vector<double>& getSpotVol()     { initialize(); return spotVol; }
    const DoubleArray& getLnFwdFx()        { initialize(); return LnFwdFx; }
    const vector<double>& getFXSmile_a1()  { initialize(); return FXSmile_a1; }
    const vector<double>& getFXSmile_a2()  { initialize(); return FXSmile_a2; }
    const vector<double>& getFXSmile_a3()  { initialize(); return FXSmile_a3; }
    double getLogSpotFX() const { return logSpotFX; }

    void  setTimeLine(DateTimeArrayConstSP simDates); // called when we know allDates
    void setMomentMatchingFlag(bool mm) {momentMatching = mm;}
    bool getMomentMatchingFlag() const {return momentMatching;}
    void computeOriginalSpotPrices(const DateTimeArray& dates, DoubleArray & spots) const;

private:
    SRMFXUtil(const SRMFXUtil& rhs);
    SRMFXUtil& operator=(const SRMFXUtil& rhs);

    /** Calculate the spot vol curve from the input vol curves.
        From fxdiffuse:FX_SpotVol. Returns timeline and the spot
        vols along it which can then be fed into extendSpotVol below */
    void calcSpotVol(DateTimeArray&  SpotVolDate,
                     vector<double>& SpotFxVol);

    /** given a timeline for the spot vols extends over the actual simulation
        dates. Inputs ignored for VolBootstrapMode == CONSTANT_SPOT_VOL */
    vector<double> extendSpotVol(const DateTimeArray&  SpotVolDate,
                                 const vector<double>& SpotFxVol);

    /**********  From fxvol.c:CalcPerturb  ********************************
     *
     *  Returns integral of
     *         vol^2(s) * volA(s,T_n) * volA_xx(s,T_n) / vol^2(T_{n-1})
     *  over (t,T_{n-1})
     *  output[i] contains the value of integral over [T_i,T_{n-1}]
     *
     *  Notes:
     *
     *  1) uses two approximations,
     *
     *                [ \tau \Sigma_{A,x}^2 ]_{\tau} \approx \Sigma_{A,x}^2
     *                     sigma \approx sigma(t)  instead of sigma(t,T)
     *
     *  2) t[n] used
     *
     *****************************************************************/
    void calcPerturb (
        vector<double>&        output,      /* (O) indexed [0,n-1] */
        int                    n,           /* (I) index of last time
                                               point accessed, t[n] etc */
        const vector<double>&  vol,         /* (I) spot fx vol, [0,n-1] */
        const vector<double>&  loc_vol_x,   /* (I) slope of local vol wrt
                                               log-money, [0,n]          */
        const vector<double>&  loc_vol_xx,  /* (I) curvature of local vol wrt
                                               log-money, [0,n]      */
        const vector<double>&  t,           /* (I) time grid,
                                               Time[0] = 0, [0,n] */
        const vector<double>&  dt);

    /********** From fxvol.c:MultiFac_Spot_FXVol2 ***************************/
    /* Outputs as many spot fx vols as there are Vol Expiry dates (i.e NbVol).
       The correlation factors are with respect to the exponential factors

       NOTES:

       1) Assumes that TPDate[0]=fx_data.Today.
       2) Assumes and checks that each VolDate and VolIntegrationDate is
       on the time
       line.
       3) resetDate and expiryDate must be entered in a
       strictly ascending order.

    ***************************************************************/
    void multiFacSpotFxVol2(
        vector<double>&       SpotFxVol, /* (M) output  */
        const vector<double>& RhoFxDomFac, /* (I) correl fx/dom. curve */
        const vector<double>& RhoFxForFac, /* (I) correl fx/for. curve */
        const vector<double>& RhoDomFacForFac, /* (I) correl dom./for. curves */
        const DateTimeArray&  expiryDate,  //(I) Modified vol->compVolDate
        const vector<double>& A_Fx, /* (I) FX smile param   */
        const vector<double>& B_Fx, /* (I) FX smile param   */
        const vector<double>& C_Fx, /* (I) FX smile param   */
        const DateTimeArray&  TPDate); /* (I) date of each pt.              */

    /*****  From fxdiffuse.c:FX_SmileParamExtend  *******************
     *
     *  Extend the input smile parameters (smilea1, smileA2 and smileA3)
     *      using the flat interp method to generate the extended FXSmile_a1,
     *      FXSmile_a2, FXSmile_a3, and modifies the date indexing of
     *      the extended curve.
     *
     *      Cases where a period in the extended curve covers multiple periods
     *      in the source curve are NOT ALLOWED. It means that the set of source
     *      dates up to the last time point must be a subset of the output dates
     *
     *      On input:
     *      ---------
     *      fxInp->FXSmile_xx[t] covers between fxInp->FXSmileDate[t-1] to
     *      fxInp->FXSmileDate[t] (when t=0, it will be between Today and
     *      fxInp->FXSmileDate[0])
     *
     *      Note: SrcDates must be in strictly ascending order
     *
     *      On output:
     *      ----------
     *      fx->FXSmile_xx[t] is the extended param values applicable between
     *      fx->SpotVolDate[t] and fx->SpotVolDate[t+1].
     *
     *      Note: 1) We assume TimePt[0] to be Today
     *            2) Size of fx->FXSmile_xx[] must be the same as
     *               that of TimePt[]
     *            3) We assume fx->FXSmile_xx[NbTimePts-1] is never used because
     *               there is no need for the smile parameters after
     *               fx->SpotVolDate[NbTimePts-1]
     *            4) TimePts must be in strictly ascending order
     *            5) all memories must be alloc'ed on entry
     *            6) {fxInp->FXSmileDate[j]} must be a subset of
     *               {fx->SpotVolDate[j]}
     *
     */
    void smileParamExtend();

    /*****  From fxdiffuse.c:FX_SmileParamTrans  ************************
     *
     *      Transform the input smile parameters
     *      It is assumed that the inputs are:
     *      Skew                            = FXSmile_a1
     *      Curvature                       = FXSmile_a2
     *      M (variation around ATM vol)    = FXSmile_a3
     *
     *      Call this function after calling FX_SmileParamExtend
     */
    void smileParamTrans();

    static const double SPOT_CUTOFF_R; // max of (CompFX Vol/SpotFX Vol) allowed



    //// fields ////
    SRMRatesUtilSP         domIR;
    SRMRatesUtilSP         forIR;
    SRMFXVolBaseSP         vol;    // fx vol
    double                 InpCorrFxDom; // correlation between fx and domestic IR
    double                 InpCorrFxFor; // correlation between fx and foreign IR
    double                 InpCorrDomFor;// correlation between domestic and foreign IR
    string                 VolBootstrapMode; // choice for bootstrapping vol
    double                 FxCutOffLevel;
    vector<double>         FXSmile_a1; // processed smile parameters
    vector<double>         FXSmile_a2;
    vector<double>         FXSmile_a3;
    vector<double>         spotVol; /* processed spot Vol (do we need this as a field?)
                                    along timeline */
    DoubleArray            LnFwdFx; // log of determinstic forward fx rate
    double                 logSpotFX; // log of spot fx today
    string                 fxName; // name of fx asset
    DateTimeArray          calibSpotVolDates; /* essentially option expiries (but with
                                              hard coded spot vol dates as well) */
    vector<double>         calibSpotVols;     /* spot vols on calibSpotVols */

	FXAssetConstSP         fxAsset;

    bool                   momentMatching;
	
	bool initialized; // true when allDates are known
    void initialize(); // called when we know domIR->simDates()
	bool isInitialized() const {return initialized;}
	void setInitialized(bool val=true) {initialized=val;}
	
};

DECLARE(SRMFXUtil);

DRLIB_END_NAMESPACE    

#endif // EDR_SRMFXUTIL_HPP
