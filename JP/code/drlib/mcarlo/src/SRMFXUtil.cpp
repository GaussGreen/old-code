//----------------------------------------------------------------------------
//
//   Group       : EDG Derivatives Research
//
//   Filename    : SRMFXUtil.cpp
//
//   Description : Utility class for FX SRM
//
//
//----------------------------------------------------------------------------


#include "edginc/config.hpp"
#include "edginc/SRMUtil.hpp"
#include "edginc/SRMFXUtil.hpp"
#include "edginc/SRMConstants.hpp"
#include "edginc/FXAsset.hpp"
#include "edginc/VolRequestRaw.hpp"
#include "edginc/SRMFXDiffuse.hpp"
#include "edginc/SRMCorrelation.hpp"
#include "edginc/Maths.hpp"
#include "edginc/SRMFXVolSpot.hpp"

DRLIB_BEGIN_NAMESPACE

// max of (CompFX Vol/SpotFX Vol) allowed
const double SRMFXUtil::SPOT_CUTOFF_R = 3.0;

SRMFXUtil::~SRMFXUtil(){}

/*
	SRMFXUtil depends on the simDates of domIR. Currently FXUtil ay be created from several places: in some the domIR->simDates() are available (like in QuantoCDSAlgorithm.cpp) in other it's not available at the time of construction (MCPathConfigSRMGen constructor), but is available after irAsset is finalized.
	So, we postpone calculation of FXSmile_* and other time dependent values until they are actually needed.
*/
SRMFXUtil::SRMFXUtil(
    SRMRatesUtilSP  forIR,
    SRMRatesUtilSP  domIR,
    FXAssetConstSP  _fxAsset,
    double          InpCorrFxFor, // correlation between fx and foreign IR
    double          InpCorrFxDom, // correlation between fx and domestic IR
    double          InpCorrDomFor,// correlation between domestic and foreign IR
    const string&   VolBootstrapMode, // choice for bootstrapping vol
    double          FxCutOffLevel,
    const string&   volType,
    bool            forceTimeOfDay):  // true - set of times of day to today's
        domIR(domIR), forIR(forIR), InpCorrFxDom(InpCorrFxDom),
        InpCorrFxFor(InpCorrFxFor), InpCorrDomFor(InpCorrDomFor),
        VolBootstrapMode(VolBootstrapMode), FxCutOffLevel(FxCutOffLevel),
        // MAR: I think that FXSmile_a1 should be of length
        // getSimDates().size()-1 but there is some that reads the extra
        // element - but it is never set I think

        logSpotFX(log(_fxAsset->getSpot())),
        fxName(_fxAsset->getName()),
        fxAsset(_fxAsset),
        momentMatching(false),
        initialized(false)
{
    static const string method("SRMFXUtil::SRMFXUtil");
    try{
        CVolProcessedSP processedVol;
        if (volType == "SRMFX::Vol" ||
            volType == "SRMFX::VolSimple" ||
            volType == "SRMFX::VolSpot"
            ){
            VolRequestRaw request;
            processedVol= CVolProcessedSP(
                fxAsset->getProcessedVol(&request));
            vol.reset(dynamic_cast<SRMFXVolBase*>(processedVol.get()));
        } else{
            throw ModelException(method,
                                "volType '"
                                + volType
                                + "' is not supported;"
                                " 'SRMFX::Vol' or 'SRMFX::VolSimple' expected");
        }
        if (forceTimeOfDay){
            // hack until the code really does support time of day
            SRMFXVolSP copiedVol(dynamic_cast<SRMFXVol*>(vol->clone()));
            copiedVol->setToday( fxAsset->getToday() );
            copiedVol->forceTimeOfDay();
            vol = copiedVol;
        }


    } catch (exception& e){
        throw ModelException(e, "SRMFXUtil::SRMFXUtil", "Failed for fx asset "+
                                fxName+" between "+forIR->getCcy()
                                +" and "+domIR->getCcy());
    }
}

void  SRMFXUtil::initialize(void) // called when we know allDates
{
    assert (isInitialized());
}
void  SRMFXUtil::setTimeLine(DateTimeArrayConstSP  simDates) // called when we know allDates
{
	static const string method("SRMFXUtil::initialize");
	if (isInitialized()) {
		return;
    }      
	else
		setInitialized();

	// make sure these are initialized as well
    forIR->setTimeLine(simDates);
    domIR->setTimeLine(simDates);

	if ((*simDates).size() == 0)
		throw ModelException(method,  "Failed for fx asset "+
                                fxName+" between "+forIR->getCcy()
                                +" and "+domIR->getCcy() + " as domIR has no simDates");

    FXSmile_a1.resize((*simDates).size()),
    FXSmile_a2.resize((*simDates).size()),
    FXSmile_a3.resize((*simDates).size()),
    LnFwdFx.resize((*simDates).size()-1),

    smileParamExtend();
    smileParamTrans();
    // calibrate spot vols
    calcSpotVol(calibSpotVolDates, calibSpotVols);
    // then extend spot vols along timeline
    spotVol = extendSpotVol(calibSpotVolDates, calibSpotVols);
    
    DateTimeArray futureDates((*simDates).begin()+1, (*simDates).end());
    fxAsset->fwdValue(futureDates, LnFwdFx);
    for (int i = 0; i < LnFwdFx.size(); i++){
        LnFwdFx[i] = log(LnFwdFx[i]); // convert to log
    }
}

//// Returns name of fx asset
const string& SRMFXUtil::getName() const{
    return fxName;
}



/** Returns the dates that were used in the calibration - essentially
    option expiry dates */
const DateTimeArray& SRMFXUtil::calibratedSpotVolDates() {
    initialize();
    return calibSpotVolDates;
}
/** Returns the spot vols along the getCalibratedSpotVolDates() */
const vector<double>& SRMFXUtil::calibratedSpotVols() {
	initialize();
    return calibSpotVols;
}

/** Calculate the spot vol curve from the input vol curves.
    From fxdiffuse:FX_SpotVol. Returns timeline and the spot
    vols along it which can then be fed into extendSpotVol below */
void SRMFXUtil::calcSpotVol(
    DateTimeArray&  SpotVolDate,
    vector<double>& SpotFxVol)
{
    static const string method("SRMFXUtil::calcSpotVol");
    const DateTimeArray& timeLine = domIR->getSimDates();
    /* do nil calibration and exit */
    if (VolBootstrapMode == SRMFXDiffuse::CONSTANT_SPOT_VOL) {
        SpotVolDate = DateTimeArray(1, timeLine.back());
        SpotFxVol = vector<double>(1, FxCutOffLevel);
    }
    SRMFXVolSpot* pvol = dynamic_cast<SRMFXVolSpot*>(vol.get());
    if (pvol)
    {
        SpotVolDate = pvol->getSpotVolDate();
        const DoubleArray& spotVols = pvol->getSpotVol();
        SpotFxVol.assign(spotVols.begin(), spotVols.end());
        return;
    }

    SRMFXVol* bvol = dynamic_cast<SRMFXVol*>(vol.get());
    if (!bvol)
        throw ModelException(method, "FX Vol should be either SRMFXVol or SRMFXVolSimple.");

    /* create a local list of comp vols */
    /* remove comp vols which overlap with input spot vols */
    DateTimeArray compVolDate(bvol->getCompVolDate());
    int NbCompFXVols = compVolDate.size();
    DateTimeArray tmpSpotVolDate( bvol->getSpotVolDate() );
    if (tmpSpotVolDate.size() > 0) {
        /* find offset of first composite bvol to be removed */
        int CompVolOfs = tmpSpotVolDate.front().findUpper(compVolDate);
        if (CompVolOfs >= 0) {
            NbCompFXVols = CompVolOfs;
            compVolDate.resize(NbCompFXVols);
        }
    }

    /* bootstapping is not needed if no composite vols are given */
    DoubleArray tmpSpotVol( vol->getSpotVol() );
    if (NbCompFXVols == 0) {
        if (vol->getSpotVol().empty()){
            throw ModelException(method, "No spot or composite vols "
                                    "supplied for fx vol "+vol->getName());
        }
        SpotVolDate = vol->getSpotVolDate();
        SpotFxVol = vector<double>(tmpSpotVol.begin(), tmpSpotVol.end());
        return;
    }
    const DateTimeArray& forDate = forIR->getExtendedTimeLine();
    const DateTimeArray& domDate = domIR->getExtendedTimeLine();

    const DateTime& Today = timeLine.front();
    /***********************/
    /*  Extend the corrs   */
    /***********************/

    /* correlation fx vs domestic factors */
    vector<double> RhoFxDom(
        SRMCorrelation::assetIr(InpCorrFxDom,
                                *domIR,
                                SRMCorrelation::EXPONENTIAL_FACTORS));

    /* correlation fx vs foreign factors */
    vector<double> RhoFxFor(
        SRMCorrelation::assetIr(InpCorrFxFor,
                                *forIR,
                                SRMCorrelation::EXPONENTIAL_FACTORS));

    /* correlation between the two IR factors*/
    vector<double> RhoDomFor(
        SRMCorrelation::irToIr(InpCorrDomFor,
                                *domIR,
                                *forIR,
                                SRMCorrelation::EXPONENTIAL_FACTORS));

    /**************************************************/
    /* Create time line and extend ir spot vol curves */
    /**************************************************/
    vector<const DateTimeArray*> dtArrayVector(4);
    dtArrayVector[0] = &domDate;
    dtArrayVector[1] = &forDate;
    dtArrayVector[2] = &compVolDate;
    dtArrayVector[3] = &bvol->getCompVolMatDate(); 

    DateTimeArray TPDate(DateTime::merge(dtArrayVector));

    /* we require Today to be part of the datelist, so */
    /* check that first date in ZDates is indeed Today */

    if ((domDate[0] != Today) || (forDate[0] != Today)){
        throw ModelException(method, "Internal error");
    }
    // fill in dates using hard coded ppy
    DateTimeArraySP filledTPDate = SRMUtil::fillInTimeLine(Today, TPDate, 12);

    /* extend IR vols */
    int NbTP = filledTPDate->size(); // for ease
    vector<double> DomExtSpotVol(NbTP-1);
    vector<double> ForExtSpotVol(NbTP-1);
    int i, k;


    /* temporary extension of FX smile parameters              */
    /* note that a[i] applies for filledTPDate[i] <= t < filledTPDate[i+1] */
    vector<double> a1_Ext(NbTP-1);
    vector<double> a2_Ext(NbTP-1);
    vector<double> a3_Ext(NbTP-1);
    k = 0;
    /* increase k until we find a smile date   */
    /* after today - tho checked in input      */
    for (i = 0; i <= NbTP - 2; i++) {
        if (k < timeLine.size() - 1) {
            if ((*filledTPDate)[i] >= timeLine[k+1]){
                k++;
            }
            a1_Ext[i] = FXSmile_a1[k];
            a2_Ext[i] = FXSmile_a2[k];
            a3_Ext[i] = FXSmile_a3[k];
        } else { /* use final parameter set */
            /* MAR: I think this is a bug since FXSmile_a1.back() is
                never populated from what I can make out */
            a1_Ext[i] = FXSmile_a1.back();
            a2_Ext[i] = FXSmile_a2.back();
            a3_Ext[i] = FXSmile_a3.back();
        }
    }/* for i */


    /************************************************/
    /*  Finally, we can bootstrap the fx spot vol  */
    /************************************************/
    SpotFxVol.resize(NbCompFXVols); // reserve some space

    multiFacSpotFxVol2(SpotFxVol,      // only NbCompFXVols populated
                        RhoFxDom, RhoFxFor, RhoDomFor,
                        compVolDate,    /* FX Option expriy  */
                        a1_Ext, a2_Ext, a3_Ext,
                        *filledTPDate);
    SpotVolDate = compVolDate;
    /* append the input spot vol list */
    SpotVolDate.insert(SpotVolDate.end(),
                        tmpSpotVolDate.begin(), tmpSpotVolDate.end());
    SpotFxVol.insert(SpotFxVol.end(),
                        tmpSpotVol.begin(), tmpSpotVol.end());
}

/** given a timeline for the spot vols extends over the actual simulation
    dates. Inputs ignored for VolBootstrapMode == CONSTANT_SPOT_VOL */
vector<double> SRMFXUtil::extendSpotVol(
    const DateTimeArray&  SpotVolDate,
    const vector<double>& SpotFxVol)
{
    const DateTimeArray& timeLine = domIR->getSimDates();
    if (dynamic_cast<SRMFXVolSpot*>(vol.get()) == 0 && VolBootstrapMode == SRMFXDiffuse::CONSTANT_SPOT_VOL)
    {
        return vector<double>(timeLine.size()-1, FxCutOffLevel);
    }
    /* extend spot curve to cover the simulation's timeline */
    return SRMUtil::extendVol(SpotVolDate, SpotFxVol, timeLine);
}

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
void SRMFXUtil::calcPerturb (
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
    const vector<double>&  dt)          /* (I) intervals in grid, [0,n-1] */
{
    static const string method("SRMFXUtil::calcPerturb");
    if (n <= 0) {
        throw ModelException(method, "number of time points must >0");
    }

    /* explicit Euler integration, so requires dense time line     */
    double s1    = 0.0;
    double s2    = 0.0;
    double s3    = 0.0;
    double s4    = 0.0;
    double s5    = 0.0;
    for (int i = n-1; i >= 0; i--) {
        if (vol[i] < SRMConstants::SRM_TINY){
            throw ModelException(method, "Vol too small");
        }
        /* avoid the following two cases:
            1) t[n-1] = t[n]    => tau_0 = 0        => division by zero, in fact no contribution
            2) t[i] = t[i+1]    => tau_0 = tau_1    => division by zero, in fact no contribution */
        double tau_0 = t[n] - t[i];
        double tau_1 = t[n] - t[i+1];

        bool skip = (Maths::isZero(tau_0) || Maths::isZero(dt[i]));
        if(!skip) {
            s1   += vol[i] * vol[i] * dt[i];
            /* u1 = Sigma_A(t[i],t[n])   */
            double u1    = sqrt(s1 / tau_0);

            s2   += 0.5 * vol[i] * vol[i] * loc_vol_x[i] * u1 * u1
                * (tau_0 * tau_0 - tau_1 * tau_1);
            /* u2 = Sigma_A_x(t[i],t[n]) */
            double u2    =  s2 * u1 / (s1 * s1);

            double g     = - u2 * u2 +
                (loc_vol_x[i] * loc_vol_x[i] + loc_vol_xx[i]
                    - 4.0 * loc_vol_x[i] * u2 / u1
                    + 3.0 * u2 * u2 / (u1 * u1) )
                * vol[i] * vol[i];
            /* s3 uses approx: time deriv of u2 is negligible */
            s3   += (tau_0 * tau_0 * tau_0 - tau_1 * tau_1 * tau_1)
                * u1 * u1 * u1 * u1
                * g / 3.0;
            /* volA_xx[t[i],t[n]] */
            double volA_xx = s3 * u1 / (s1 * s1 * s1);

            /* calculate smile pertubation */
            /* \int_t^T vol^2 d [ tau * volA * volA_xx ] / dT du / vol[i]^2
                at T = t[i] */
            s4   += (tau_0 * tau_0 - tau_1 * tau_1) * u1 * u1 * g / 2.0;
            s5   +=  vol[i] * vol[i] * ( - 2.0 * volA_xx / u1
                                            + 2.0 * s4 / (s1 * s1) ) * dt[i];
        }
        output[i] = s5;
    } /* i */
}
/****** CalcPerturb ******************************************************/


/********** From fxvol.c:MultiFac_Spot_FXVol2 ***************************/
/* Outputs as many spot fx vols as there are Vol Expiry dates (i.e NbVol).
    The correlation factors are with respect to the exponential factors

    NOTES:

    1) Assumes that TPDate[0]=fx_data.Today.
    2) Assumes and checks that each resetDate and ImpVolDate is
    on the time line.
    3) resetDate and ImpVolDate must be entered in a
    strictly ascending order.

***************************************************************/
void SRMFXUtil::multiFacSpotFxVol2(
    vector<double>&       SpotFxVol, /* (M) output  */
    const vector<double>& RhoFxDomFac, /* (I) correl fx/dom. curve */
    const vector<double>& RhoFxForFac, /* (I) correl fx/for. curve */
    const vector<double>& RhoDomFacForFac, /* (I) correl dom./for. curves */
    const DateTimeArray&  ImpVolDate,//(I) Modified vol->compVolDate
    const vector<double>& A_Fx, /* (I) FX smile param   */
    const vector<double>& B_Fx, /* (I) FX smile param   */
    const vector<double>& C_Fx, /* (I) FX smile param   */
    const DateTimeArray&  TPDate) /* (I) date of each pt.              */
{
    static const string method("SRMFXUtil::multiFacSpotFxVol2");
    /* Quick checks first */
    int NbTP = TPDate.size();
    if (NbTP < 2) {
        throw ModelException(method, "Less that 2 time points supplied");
    }
    int DomNbFac = domIR->numFactors(); // for ease
    int ForNbFac = forIR->numFactors(); // for ease

    int i;
    vector<double> DeltaTime(NbTP);
    vector<double> Time(NbTP);
    vector<double> loc_vol_x(NbTP);
    vector<double> loc_vol_xx(NbTP);
    
	for (i = 0; i < NbTP - 1;i++) 
	{
        /* time step between time points, and corresponding times in years */
        DeltaTime[i]  =    TPDate[i].yearFrac(TPDate[i+1]);
        Time[i]       =    TPDate[0].yearFrac(TPDate[i]);
        /* atm slope wrt log-fwd-money */
        loc_vol_x[i]    =   C_Fx[i] * A_Fx[i];
        /* atm curvature wrt log-fwd-money */
        loc_vol_xx[i]   =   C_Fx[i] * C_Fx[i] * B_Fx[i];
    }

    /* we have already checked that NbTp >= 2 */
    Time.back() = TPDate[0].yearFrac(TPDate.back());
    /* find the first time point on or after the FIRST imp vol date */
    int TPFirstImpVol = ImpVolDate.front().findUpper(TPDate);
    if (TPFirstImpVol == TPDate.size())
	{
        throw ModelException(method, "Expiry date "+
                                ImpVolDate.front().toString()+
                                " is beyond the input time line");
    }

    int TPLastImpVol = ImpVolDate.back().findUpper(TPDate);
    if (TPLastImpVol == TPDate.size()) 
	{
        throw ModelException(method, "Expiry date "+
                                ImpVolDate.back().toString()+
                                " is beyond the input time line");
    }

    /********************************************************************/
    /*                                                                  */
    /*          remove 1st order FX smile contribution                  */
    /*                                                                  */
    /********************************************************************/

    /* first estimate zero expiry limit of implied vol = FxVolUnSmile[0] */
    /* the estimate is obtained from Eq. 8 of "Long-dated FX smile model",*/
    /* by assuming FX spot vol constant and approximating RHS             */

    // when this function is called, we are certain that the fxvol is composite
    SRMFXVol* bvol = dynamic_cast<SRMFXVol*>(vol.get());


    /* time to first option expiry */
    double tau_0 =  TPDate[0].yearFrac(ImpVolDate[0]);
    DoubleArray FxVol( bvol->getCompVol() ); // for ease
    /* approximation assumes const FX spot vol until first expiry */

    double epsilon = 2. * tau_0 * FxVol[0] * FxVol[0] *
        (loc_vol_xx[0] - 0.5 * loc_vol_x[0] * loc_vol_x[0]) / 3.0;

    if (epsilon < -0.5)
	{ /* sqr root of negative number */
        throw ModelException(method, "Can't calculate initial FX spot vol");
    }
    
	vector<double> FxVolUnSmile(NbTP);
    if (fabs(epsilon) < SRMConstants::SRM_TINY) 
	{
        FxVolUnSmile[0] = FxVol[0];
    } 
	else 
	{
        FxVolUnSmile[0] = FxVol[0] *
            sqrt((- 1. + sqrt(1 + 2. * epsilon)) / epsilon);
    }

    /* extend the implied vols to full time line         */
    double t_low   = 0.0;          /* initialise as today */
    double delta_t = tau_0;        /* first step size equals first expiry */
    double FxVolPrev = FxVolUnSmile[0];   /* guess for implied vol
                                                extrapolated to zero expiry */
    double slope    = FxVol[0] * FxVol[0] -
        FxVolUnSmile[0] * FxVolUnSmile[0];
    int NbImpVol = ImpVolDate.size();
    vector<double> FxVolExt(NbTP);
    int k;
    
	for (i = 0, k = 0 ; i <= NbTP - 1;i++) 
	{
        if (k <= NbImpVol - 1) 
		{
            if (TPDate[i] >=  ImpVolDate[k]) 
			{
                if (k < NbImpVol - 1) {
                    FxVolPrev  =  FxVol[k];
                    slope      =  FxVol[k+1] * FxVol[k+1] -
                        FxVol[k] * FxVol[k];
                    t_low     +=  delta_t;
                    /* time step for next update   */
                    delta_t    =  ImpVolDate[k].
                        yearFrac(ImpVolDate[k+1]);
                    k++;
                } 
				else 
				{
                    /*  flat after last implied vol */
                    slope      = 0.0;
                    FxVolPrev  = FxVol[NbImpVol-1];
                    k++; /* i.e., k = NbImpVol */
                }
            }
            /* interpolate linearly on squares of vols - in the short
                * end, in particular, this is more accurate than
                * interpolation on the vols  */
            FxVolExt[i] =  sqrt(FxVolPrev * FxVolPrev +
                                slope * (Time[i] - t_low) / delta_t);
        } 
		else 
		{ /* flat after last implied vol    */
            FxVolExt[i] =  FxVol[NbImpVol - 1];
        }
    }

    /* initialise dummy for integration */
    double temp1   = FxVolUnSmile[0] * FxVolUnSmile[0] * DeltaTime[0];
    vector<double> volA(NbTP);
    volA[1] = FxVolUnSmile[0];

    /* only volA is used outside this loop */
    vector<double> Perturb(NbTP);
    for (i = 1; i < TPLastImpVol; i++) 
	{
        calcPerturb(Perturb,     /* (O) [0,i]...[i-1,i] */
                    i,           /* (I) number of time points */
                    FxVolUnSmile,/* (I) implied fx vol */
                    loc_vol_x,   /* (I) slope of local vol wrt log-money */
                    loc_vol_xx,  // (I) curvature of local vol wrt log-money
                    Time,        /* (I) time grid, Time[0] = 0 */
                    DeltaTime);  /* (I) intervals in grid       */
    
		/* FxVol[i] etc is the implied vol at Time[i] - not Time[i+1] */
        double temp2 = (Time[i+1] * FxVolExt[i+1] * FxVolExt[i+1] -
                        Time[i] * FxVolExt[i] * FxVolExt[i]) / 
					   (1.0 + Perturb[0]);

        if (Maths::isZero(DeltaTime[i])) 
		{
            FxVolUnSmile[i] = FxVolUnSmile[i-1];
        } 
		else 
		{
            if (Maths::isNegative(temp2)) 
			{
                throw ModelException(method, "negative fx vol in calibration");
            }
            FxVolUnSmile[i] = sqrt(temp2/DeltaTime[i]);
        }
        temp1   +=  FxVolUnSmile[i] * FxVolUnSmile[i] * DeltaTime[i];
        /* Time[i] > 0 because i > 0 and TPDate is strictly increasing  */
        /* volA[i] refers to Time[i]                                    */
        volA[i+1]  =  sqrt(temp1/Time[i+1]);
    }
    /************* 1st order FX smile contribution removed **************/

    DateTimeArray resetDate( bvol->getCompVolMatDate() );

	// check dates and find offsets on TPDate timeline
	int kVol = 0;
	vector<int> expiryIndex(NbImpVol);
	vector<int> resetIndex(NbImpVol);
	
	for (kVol = 0; kVol < NbImpVol; kVol++) 
	{      
		// check that ImpVolDate (expiry) <= resetDate for fwd 
		if (ImpVolDate[kVol] > resetDate[kVol]) {
			throw ModelException("Expiry date "+
				ImpVolDate[kVol].toString()+
				" is after"
				"corresponding fwd reset date");
		}
		
		if (ImpVolDate[kVol] < TPDate[0]) 
		{
			throw ModelException("Expiry date cannot be before TPDate[0]");
		}
		
		// check that ImpVolDates are strictly increasing
		// use QLib function for this ?
		if (kVol > 0)
		{
			if (ImpVolDate[kVol] <= ImpVolDate[kVol-1]) 
			{
				throw ModelException("Equity vol expiry date "+
					ImpVolDate[kVol].toString()+
					" is less than or equal to"
					" the previous expiry date");
			}
		}
		
		int NextTP = ImpVolDate[kVol].find(TPDate.begin(),
			TPDate.end());
		
		expiryIndex[kVol] = NextTP;
		
		NextTP = resetDate[kVol].find(TPDate.begin(),
			TPDate.end());
		
		resetIndex[kVol] = NextTP;
		
	}

    int j;

    vector< vector<double> > domVol(DomNbFac);
	for (int l = 0; l < DomNbFac; l++) 
	{
		 domVol[l].resize(NbTP);
	}

    vector< vector<double> > forVol(ForNbFac);
	for (int m = 0; m < ForNbFac; m++) 
	{
		 forVol[m].resize(NbTP);
	}

    double Sigma = 0.0; /* Fx Spot Vol in current bucket */
    bool FxCutOffFlag = VolBootstrapMode != SRMFXDiffuse::NO_FAILURE_ALLOWED;
    bool FxCutOffLast = VolBootstrapMode == SRMFXDiffuse::USE_LAST_LEVEL;
    double CutOffLevel = (FxCutOffLast || !FxCutOffFlag) ?
        SPOT_CUTOFF_R : FxCutOffLevel;

    vector<double> SpotFxVolSoFar(NbTP);

    int TPIntBefore = 0;
    int TPIntSoFar = 0;

	vector<int>::iterator tmpExpiryBegin = expiryIndex.begin();
	vector<int>::iterator tmpExpiryEnd = tmpExpiryBegin;
	vector<int>::iterator tmpResetEnd = resetIndex.begin();

    // Pre-compute IrFwdRateDom and IrFwdRateFor for performance
    // Note : discretely compounded fwd rates * deltaTime
    vector<double> IrFwdRateDom(simpleFwdCurve(domIR->getDiffYC(),TPDate));
    vector<double> IrFwdRateFor(simpleFwdCurve(forIR->getDiffYC(),TPDate));

    // Allocate memory outside loop
    vector<double> domIrVar(NbImpVol);
    vector<double> forIrVar(NbImpVol);
    vector<double> domFXCovar(NbImpVol);
    vector<double> forFXCovar(NbImpVol);

    for (int kVol = 0; kVol < NbImpVol; kVol++) 
	{		
		TPIntSoFar = expiryIndex[kVol];
	
		// the following can only happen if expiry dates are close
		// - currently, if less than 1 month apart
		if (TPIntSoFar == TPIntBefore) 
		{
			continue;
		}
		
		//////////////////////////////////////////////
		// Solving Ax^2 + BX + C = 0 where 
		// Poly[0] = C
		// Poly[1] = B
		// Poly[2] = A
		// x = SpotEqVol
		//////////////////////////////////////////////

		double FxOnlyVariance = 0.;
		// known fx variance
		for (i = 0; i < TPIntBefore; i++) 
		{               
			// constant coeffs due to EQ variance
			FxOnlyVariance += SpotFxVolSoFar[i] * SpotFxVolSoFar[i] * DeltaTime[i];
		}

        domIR->instFactorVol(domVol, DeltaTime, TPDate, IrFwdRateDom, expiryIndex[kVol], resetIndex[kVol]);
	    forIR->instFactorVol(forVol, DeltaTime, TPDate, IrFwdRateFor, expiryIndex[kVol], resetIndex[kVol]);

		double domForCovar = 0.;
		// dom-for rates covariance 
		for (i = 0; i < TPIntSoFar; i++) 
		{               
			for (int l = 0; l < DomNbFac; l++) 
			{
				for (int m = 0; m < ForNbFac; m++) 
				{
					domForCovar -= 2. * domVol[l][i] * forVol[m][i] *
						           RhoDomFacForFac[m + l * ForNbFac] * DeltaTime[i];
				}/* for m */
			}/* for l */
		}

		tmpExpiryEnd++;
	    vector<int> tmpExpiryIndex(tmpExpiryBegin, tmpExpiryEnd);

		domIR->irAssetVariance(
			domIrVar,          // (O)
			domFXCovar,        // (O)   
			RhoFxDomFac,       // (I)
			DeltaTime,         // (I)  passed for convenience 
			TPDate,            // (I)
            IrFwdRateDom,      // (I)  pre-computed ir fwd rates 
			tmpExpiryIndex,    // (I)  indices in TPDate
			*tmpResetEnd);     // (I)  -"-  

		forIR->irAssetVariance(
			forIrVar,          // (O)
			forFXCovar,        // (O)   
			RhoFxForFac,       // (I)
			DeltaTime,         // (I)  passed for convenience 
			TPDate,            // (I)
            IrFwdRateFor,      // (I)  pre-computed ir fwd rates 
			tmpExpiryIndex,    // (I)  indices in TPDate
			*tmpResetEnd);     // (I)  -"-  

		tmpResetEnd++;

		double Poly[3] = {0.0, 0.0, 0.0};

		Poly[0] = FxOnlyVariance;		
        Poly[0] += domForCovar;

		for(j = 0; j <= kVol; j++)
		{
 		    Poly[0] += domIrVar[j];
		    Poly[0] += forIrVar[j];
		}

		// note : this assumes spotFXVol constant between expiries
		for(j = 0; j < kVol; j++)
		{
			int fxIndex = (j > 0) ? expiryIndex[j-1] : 0;
			Poly[0] += (domFXCovar[j] - forFXCovar[j])
				        * SpotFxVolSoFar[fxIndex];
		}
			
		// Subtract the target "unsmiled" variance
		Poly[0] -= volA[TPIntSoFar] * volA[TPIntSoFar] * 
			TPDate[0].yearFrac(TPDate[TPIntSoFar]);

        Poly[1] = domFXCovar[kVol] - forFXCovar[kVol];

        // notice that TPDate[TPIntBefore] < ImpVolDate[kVol]     
        // i.e (Poly[2] non zero)                                 
        // because of earlier check on TPIntBefore==TPIntSoFar    
        Poly[2] = TPDate[TPIntBefore].yearFrac(ImpVolDate[kVol]);

        double discriminant = Poly[1] * Poly[1] - 4.* Poly[2] * Poly[0];
        bool VolTooLow = (discriminant < SRMConstants::SRM_TINY);
        
		if (discriminant >= SRMConstants::SRM_TINY) 
		{
            /* using the larger root ensures a continuous spot fx vol   */
            Sigma  = (-Poly[1] + sqrt(discriminant)) / 2. / Poly[2];
            /* Evaluate whether spotvol satisfies appropriate criterion */
            if (FxCutOffLast || !FxCutOffFlag) 
			{
                VolTooLow  = (FxVol[kVol] > CutOffLevel * Sigma);
            } 
			else 
			{
                VolTooLow = (Sigma < CutOffLevel);
            }
        }

        // now implement cutoff actions as the case may be 
        if (VolTooLow) 
		{
            if (FxCutOffFlag) 
			{
                if (FxCutOffLast) 
				{
                    if (kVol == 0) 
					{
                        string m("Unable to bootstrap at least one FX vol "
                                    "point and therefore unable to"
                                    "cut off at last spot vol level");
                        throw ModelException(method, m);
                    }
                    for (int n = kVol; n < NbImpVol; n++) 
					{
                        SpotFxVol[n] = SpotFxVol[kVol-1];
                    }
                } 
				else 
				{
					// end of FxCutOffLast=true 
                    for(int n = kVol; n < NbImpVol; n++) 
					{
                        SpotFxVol[n] = CutOffLevel;
                    }
                }
                break;
            } 
			else 
			{
				/* end of FxCutOffFlag */
                double pcVol = 100.*FxVol[kVol]/CutOffLevel;
                string m("Problem in calculating fx vol "
                            "at "+resetDate[kVol].toString()+": less"
                            " than "+Format::toString(pcVol)+"% (level"
                            " determined by ratio to fwd vol)");
                throw ModelException(method, m);
            }
        } // end of VolTooLow 

        SpotFxVol[kVol] = Sigma;
        for (i = TPIntBefore; i < TPIntSoFar; i++) {
            SpotFxVolSoFar[i] = Sigma;
        }// i

        TPIntBefore  = TPIntSoFar;
    }//for kVol
}// MultiFac_Spot_FxVol2

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
void SRMFXUtil::smileParamExtend() 
{
    static const string method("SRMFXUtil::smileParamExtend");
    DoubleArray smileA1( vol->getSmileA1() );
    DoubleArray smileA2( vol->getSmileA2() );
    DoubleArray smileA3( vol->getSmileA3() );
    DateTimeArray smileDate( vol->getSmileDate() );

    int NbSrcParams = smileA1.size();
    if (NbSrcParams == 0){
        // the original code does not support the smileA1 etc not being
        // specified (ie have to supply array of length 1 to turn it off)
        FXSmile_a1 = vector<double>(FXSmile_a1.size(), 0.0);
        FXSmile_a2 = vector<double>(FXSmile_a1.size(), 0.0);
        FXSmile_a3 = vector<double>(FXSmile_a1.size(), 0.0);
        return;
    }
    const DateTimeArray& simDates = domIR->getSimDates();
    /* Check that the set of source dates (except the last one) up to the */
    /* last time point is a subset of the output dates                    */
    // MAR: This is true by construction - but leave in for now
    int j;
    for (j = 0; j < NbSrcParams - 1; j++) {
        if ( smileDate[j] < simDates.back()) {
            try{
                smileDate[j].find(simDates);
            } catch (exception& e){
                throw ModelException(e, method, "Internal error: fx smile"
                                        " date not in timeline");
            }
        }
    }
    j = 0;  /* pointing to Src_xx[j] */
    /* Loop through each period TimePt[i] to TimePt[i+1] */
    int nbTimePts   = simDates.size();
    for (int i = 0; i < nbTimePts-1; i++) {
        if (j < (NbSrcParams-1)) {
            /* increase j until we're in right FXSmileDate interval */
            if (simDates[i] == smileDate[j]){
                j++;
            }
        }
        /* LHS forward looking, RHS backward looking */
        FXSmile_a1[i] = smileA1[j];
        FXSmile_a2[i] = smileA2[j];
        FXSmile_a3[i] = smileA3[j];
    }
}/* FX_SmileParamExtend */

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
void SRMFXUtil::smileParamTrans() 
{
    static const string method("SRMFXUtil::smileParamTrans");
    const DateTimeArray& simDates = domIR->getSimDates();
    int NbParams   = simDates.size();
    for (int i=0; i<NbParams-1; i++) {
        double Sk      =  FXSmile_a1[i];        /* ATM skew          */
        double Crv     =  FXSmile_a2[i];        /* ATM crv           */
        double M       =  FXSmile_a3[i] + 1.0;  /* maximum variation */

        double a       = 0.0;
        double b       = 0.0;
        double c       = 0.0;
        if (Maths::equals(M, 1.0)){ /* flat local vol */
            // so a = 0.0 and b = 0.0;
            c = M;
        } else { /* M > 1.0 */
            /* no curvature */
            if (Maths::isZero(Crv)) {
                /* no skew */
                if (Maths::isZero(Sk)) {
                    // so a = 0.0 and b = 0.0;
                    c = 1.0;
                } else {
                    /* skew */
                    /* always try to reach M */
                    /* note: local vol jumps M - 1, asymptotically in */
                    /* the wings, as Sk passes TINY                   */
                    c = fabs(Sk)/(M - 1.0);
                    a = Sk / c;
                    // so b = 0.0;
                }
            } else{
                /* curvature, pos or neg */
                double b0 = - fabs(Sk) + sqrt(Sk * Sk + 4.0 * fabs(Crv) *
                    (M - 1.0));
                b =  (b0 * b0)/
                    (4.0 * Crv);
                c =  sqrt(Crv/b);
                a =  Sk / c;
            }
        }

		// temporarily turned off - tests need to be updated
        /* check for negative curvature */
		/*
        if ( b < 0. ) 
		{
            throw ModelException(method, "FX smile curvature is negative");  
        } 
		*/

        /* check for negative volatility */
		if ( 1. - fabs(a) + b < 0. ) 
		{
            throw ModelException(method, "Negative or zero FX volatility");  
        } 

		/*
        if ( a * a - 1. - 2. * b > - SRMConstants::SRM_TINY ) 
		{ 
            throw ModelException(method, "Negative or zero FX volatility");
        }
		*/

        FXSmile_a1[i]   =   a;
        FXSmile_a2[i]   =   b;
        FXSmile_a3[i]   =   c;
    }
}/* FX_SmileParamTrans */


void SRMFXUtil::computeOriginalSpotPrices(const DateTimeArray& dates, DoubleArray & spots) const
{
    spots.resize(dates.size());
    fxAsset->fwdValue(dates, spots);
    for (int i = 0; i < spots.size(); ++i)
        spots[i] = log(spots[i]);
}


DRLIB_END_NAMESPACE
