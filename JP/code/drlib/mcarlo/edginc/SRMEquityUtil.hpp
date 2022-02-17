//----------------------------------------------------------------------------
//
//   Group       : QR Cross Asset
//
//   Filename    : SRMQuityUtil.hpp
//
//   Description : Utility class for Equity SRM
//
//
//----------------------------------------------------------------------------

#ifndef EDR_SRMEQUITYUTIL_HPP
#define EDR_SRMEQUITYUTIL_HPP

#include "edginc/SRMRatesUtil.hpp"
#include "edginc/SRMEQVol.hpp"
#include "edginc/Asset.hpp"
#include "edginc/BorrowCurve.hpp"
#include "edginc/VirtualDestructorBase.hpp"
#include "edginc/DECLARE.hpp"

DRLIB_BEGIN_NAMESPACE    

typedef enum
{
    ccyVanilla=0,       //payoff in same ccy as underlying
    ccyProtected,       //payoff taken in domestic currency without multiplication by FX
    ccyStruck,          //payoff multiplied by FX
} CurrencyTreatment;

class SRMEquityUtil : public virtual VirtualDestructorBase
{
public:
    ~SRMEquityUtil() {};

    SRMEquityUtil(
        SRMRatesUtilSP  baseIR,
        CAssetConstSP   eqAsset,    
        double          InpCorrEqIr, // correlation between eq and IR
        CurrencyTreatment ccyTreatment,
        double          InpCorrEqFX,
        const string&   VolBootstrapMode, // choice for bootstrapping vol
        double          EqCutOffLevel);
                                           
    //// Returns name of asset
    const string& getName() const;

    
    // accessors
    CurrencyTreatment getCurrencyTreatment() const { return ccyTreatment; }
    SRMRatesUtilSP getBaseIR() const { return baseIR; }
    /** returns bootstrapped spot vol */
    const vector<double>& getSpotVol()  { initialize(); return spotVol; }
    const DoubleArray& getLnFwdEQ() { initialize(); return LnFwdEq; }
    const vector<double>& getEQSmile_a1() { initialize(); return EqSmile_a1; }
    const vector<double>& getEQSmile_a2() { initialize(); return EqSmile_a2; }
    const vector<double>& getEQSmile_a3() { initialize(); return EqSmile_a3; }
    const DateTimeArray& getSmileDate() { initialize(); return vol->getSmileDate(); }

    double getLogSpotEQ() const { return logSpotEQ; }
    const DateTimeArray& getDivExDates() const { return divExDates; }
    DividendListSP getDivs() const { return divs; }
    BorrowCurveConstSP getBorrowCurve() const { return borrowCurve; }
    double getInpCorrEqFX() const { return InpCorrEqFX; }
    
    void  setTimeLine(DateTimeArrayConstSP simDates); // called when we know allDates

    void setMomentMatchingFlag(bool mm) {momentMatching = mm;}
    bool getMomentMatchingFlag() const {return momentMatching;}

    void computeOriginalSpotPrices(const DateTimeArray& dates, DoubleArray & spots) const;

private:
   
    /** Calculate the spot vol curve from the input vol curves.
        From EQdiffuse:EQ_SpotVol */
    vector<double> calcSpotVol();

    /**********  From eqvol.c:CalcPerturb  ********************************
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
    void calcEqPerturb( 
        vector<double>&        output,      /* (O) indexed [0,n-1] */
        int                    n,           /* (I) index of last time 
                                               point accessed, t[n] etc */
        const vector<double>&  vol,         /* (I) spot EQ vol, [0,n-1] */
        const vector<double>&  loc_vol_x,   /* (I) slope of local vol wrt 
                                               log-money, [0,n]          */
        const vector<double>&  loc_vol_xx,  /* (I) curvature of local vol wrt 
                                               log-money, [0,n]      */
        const vector<double>&  t,           /* (I) time grid,
                                               Time[0] = 0, [0,n] */
        const vector<double>&  dt);         /* (I) intervals in grid, [0,n-1] */


    /********** From EQvol.c:MultiFac_Spot_EQVol2 ***************************/
    /* Outputs as many spot EQ vols as there are Vol Expiry dates (i.e NbVol).
       The correlation factors are with respect to the exponential factors
       
       NOTES:
       
       1) Assumes that TPDate[0]=EQ_data.Today.
       2) Assumes and checks that each VolDate and VolIntegrationDate is 
       on the time 
       line.
       3) VolDate and VolIntegrationDate must be entered in a 
       strictly ascending order.
        
    ***************************************************************/
    void multiFacSpotEqVol(
        vector<double>&       SpotEqVol,          /* (M) output  */
        const vector<double>& RhoEqIr,            /* (I) correl EQ/IR curve */
        const DateTimeArray&  ImpVolDate,         //(I) Modified vol->compVolDate
        const vector<double>& A_Eq,               /* (I) EQ smile param   */
        const vector<double>& B_Eq,               /* (I) EQ smile param   */
        const vector<double>& C_Eq,               /* (I) EQ smile param   */
        const DoubleArray*    fwdsToDivPmtDates,  // (I) not currently used
        DividendListConstSP   divList,            // (I) 
        const DateTimeArray&  TPDate);            /* (I) date of each pt.              */

    /*****  From EQdiffuse.c:EQ_SmileParamExtend  *******************
     *
     *  Extend the input smile parameters (smilea1, smileA2 and smileA3)
     *      using the flat interp method to generate the extended EQSmile_a1, 
     *      EQSmile_a2, EQSmile_a3, and modifies the date indexing of 
     *      the extended curve. 
     *
     *      Cases where a period in the extended curve covers multiple periods 
     *      in the source curve are NOT ALLOWED. It means that the set of source
     *      dates up to the last time point must be a subset of the output dates
     *
     *      On input:
     *      ---------
     *      EQInp->EQSmile_xx[t] covers between EQInp->EQSmileDate[t-1] to 
     *      EQInp->EQSmileDate[t] (when t=0, it will be between Today and 
     *      EQInp->EQSmileDate[0])
     *
     *      Note: SrcDates must be in strictly ascending order
     *
     *      On output:
     *      ----------
     *      EQ->EQSmile_xx[t] is the extended param values applicable between 
     *      EQ->SpotVolDate[t] and EQ->SpotVolDate[t+1].
     *
     *      Note: 1) We assume TimePt[0] to be Today
     *            2) Size of EQ->EQSmile_xx[] must be the same as
     *               that of TimePt[]
     *            3) We assume EQ->EQSmile_xx[NbTimePts-1] is never used because
     *               there is no need for the smile parameters after
     *               EQ->SpotVolDate[NbTimePts-1]
     *            4) TimePts must be in strictly ascending order
     *            5) all memories must be alloc'ed on entry
     *            6) {EQInp->EQSmileDate[j]} must be a subset of
     *               {EQ->SpotVolDate[j]} 
     *
     */
    // XXX exactly same as fx
    void smileParamExtend();

    /*****  From EQdiffuse.c:EQ_SmileParamTrans  ************************
     *
     *      Transform the input smile parameters 
     *      It is assumed that the inputs are:
     *      Skew                            = EQSmile_a1 
     *      Curvature                       = EQSmile_a2 
     *      M (variation around ATM vol)    = EQSmile_a3 
     *       
     *      Call this function after calling EQ_SmileParamExtend 
     */
    // XXX exactly same as fx
    void smileParamTrans();

    static const double SPOT_CUTOFF_R; // max of (CompEQ Vol/SpotEQ Vol) allowed
    //// fields ////
    SRMRatesUtilSP      baseIR;  // base IR 
    CAssetConstSP       eqAsset;
    SRMEQVolConstSP     vol;    // EQ vol
    double         InpCorrEqIr; // correlation between EQ and base IR 
    CurrencyTreatment ccyTreatment;
    double         InpCorrEqFX;
    string         VolBootstrapMode; // choice for bootstrapping vol
    double         EqCutOffLevel;
    vector<double> EqSmile_a1; // processed smile parameters
    vector<double> EqSmile_a2;
    vector<double> EqSmile_a3;
    vector<double> spotVol; // processed spot Vol (do we need this as a field?)
    DoubleArray    LnFwdEq; // log of determinstic forward EQ rate
    double         logSpotEQ; // log of spot EQ today
    DividendListSP divs;       
    DateTimeArray  divExDates; // with adjusted time so as to fit in timeline XXX
    BorrowCurveConstSP  borrowCurve; 
    string         eqName; // name of asset
    bool           momentMatching;
	bool initialized; // true when allDates are known
    void initialize(void); // called when we know domIR->simDates()
	bool isInitialized(void) const {return initialized;}
	void setInitialized(bool val=true) {initialized=val;}
    DateTimeArrayConstSP  simDates; // diffusion dates   
};

DECLARE(SRMEquityUtil);

DRLIB_END_NAMESPACE    

#endif // EDR_SRMEQUITYUTIL_HPP
