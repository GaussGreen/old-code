//----------------------------------------------------------------------------
//
//   Group       : EDG Derivatives Research
//
//   Filename    : SRMFXDiffuse.cpp
//
//   Description : A generator of paths using stochastic rates
//                 for FX Assets
//
//   Date        : 13 Aug 2004
//
//
//----------------------------------------------------------------------------


#include "edginc/config.hpp"

#include "edginc/SRMConstants.hpp"
#include "edginc/SRMUtil.hpp"
#include "edginc/VolRequestRaw.hpp"
#include "edginc/VolProcessed.hpp"
#include "edginc/EquityBase.hpp"
#include "edginc/SRMCorrelation.hpp"
#include "edginc/DividendCollector.hpp"
#include "edginc/Maths.hpp"
#include "edginc/SRMEquityDiffuse.hpp"
#include "edginc/SRMEquityUtil.hpp"

DRLIB_BEGIN_NAMESPACE

// max of (CompEQ Vol/SpotEQ Vol) allowed
const double SRMEquityUtil::SPOT_CUTOFF_R = 3.0;

SRMEquityUtil::SRMEquityUtil(
    SRMRatesUtilSP  baseIR,
    CAssetConstSP   eqAsset,    
    double          InpCorrEqIr, // correlation between EQ and  IR
    CurrencyTreatment _ccyTreatment,
    double          InpCorrEqFX,
    const string&   VolBootstrapMode, // choice for bootstrapping vol
    double          EqCutOffLevel): 
        baseIR(baseIR), 
        eqAsset(eqAsset), 
        InpCorrEqIr(InpCorrEqIr),
        ccyTreatment(_ccyTreatment),
        InpCorrEqFX(InpCorrEqFX),
        VolBootstrapMode(VolBootstrapMode), 
        EqCutOffLevel(EqCutOffLevel),
        // MAR: I think that FXSmile_a1 should be of length 
        // getSimDates().size()-1 but there is some that reads the extra
        // element - but it is never set I think
        logSpotEQ(log(eqAsset->getSpot())),
        eqName(eqAsset->getName()),
        momentMatching(false),
        initialized(false)
{
}

/* Initialize the needed smiles when we have baseIR->simDates */
void SRMEquityUtil::initialize(void)
{
    assert(isInitialized());
}

void  SRMEquityUtil::setTimeLine(DateTimeArrayConstSP  _simDates) // called when we know allDates
{
	static const string method("SRMEquityUtil::initialize");
	if (isInitialized()) {
        assert(DateTime::equals(*_simDates, *simDates));
		return;
    }      
	else
		setInitialized();

    simDates = _simDates;
	
	if ((*simDates).size() == 0)
		throw ModelException(method,  "Failed for equity "+
                                getName() + " as baseIR has no simDates");
    
    baseIR->setTimeLine(simDates);
    
    EqSmile_a1.resize((*simDates).size());
    EqSmile_a2.resize((*simDates).size());
    EqSmile_a3.resize((*simDates).size());
    LnFwdEq.resize((*simDates).size()-1);

    
    try{
        // similar as in SRMFXDiffuse
        VolRequestRaw request;
        CVolProcessedSP processedVol(eqAsset->getProcessedVol(&request));
        vol.reset(dynamic_cast<SRMEQVol*>(processedVol.get()));
        if (!vol) {
            throw ModelException("SRMEquityUtil", "MCPathConfigSRM requires SRMEQVol as vol for EQ");
        }
        smileParamExtend();
        smileParamTrans();
        spotVol = calcSpotVol();
        // calculate determinstic fwds
        // const DateTimeArray& simDates = baseIR->getSimDates();
        DateTimeArray futureDates((*simDates).begin()+1, (*simDates).end());
        eqAsset->fwdValue(futureDates, LnFwdEq);
        for (int i = 0; i < LnFwdEq.size(); i++){
            LnFwdEq[i] = log(LnFwdEq[i]); // convert to log
        }

        // Extract discrete divs - for now we look only at ex-dates.
        // Pay dates are accounted for by a determinsitic PV using "baseIR->getDiscYC()"
        // Also need levels -> hence the divList for now...
        SRMEquityDiffuse::getDivData(eqAsset,
                          baseIR->getSimDates().front(), // start,
                          baseIR->getSimDates().back(), // end,
                          divExDates,
                          true, // withDivList -> populate divs
                          divs,
                          baseIR->getDiscYC());
        
        // Need cts payments - only borrow for now
        if (const EquityBase* eqb = dynamic_cast<const EquityBase*>(eqAsset.get())) {
            borrowCurve = eqb->getEquity()->getBorrow();
        }

    } catch (exception& e){
        throw ModelException(e, "SRMEquityUtil", "Failed for EQ asset "+
                             eqAsset->getName());
    }
}

//// Returns name of asset
const string& SRMEquityUtil::getName() const{
    return eqName;
}

/** Calculate the spot vol curve from the input vol curves.
    From EQdiffuse:EQ_SpotVol */
vector<double> SRMEquityUtil::calcSpotVol()
{
    static const string method("SRMEquityUtil::calcSpotVol");
    try {
        const DateTimeArray& timeLine = baseIR->getSimDates();
        /* do nil calibration and exit */
        if (VolBootstrapMode == SRMEquityDiffuse::CONSTANT_SPOT_VOL) {
            return vector<double>(timeLine.size()-1, EqCutOffLevel);
        }
        /* create a local list of comp vols */
        /* remove comp vols which overlap with input spot vols */
        DateTimeArray compVolDate( vol->getCompVolDate() );
        int NbCompEqVols = compVolDate.size();
        DateTimeArray tmpSpotVolDate( vol->getSpotVolDate() );
        if (tmpSpotVolDate.size() > 0) {
            /* find offset of first composite vol to be removed */ 
            int CompVolOfs = 
                tmpSpotVolDate.front().findUpper(compVolDate);
            if (CompVolOfs >= 0) {
                NbCompEqVols = CompVolOfs;
                compVolDate.resize(NbCompEqVols);
            }
        }

        /* bootstapping is not needed if no composite vols are given */
        DoubleArray tmpSpotVol( vol->getSpotVol() );
        if (NbCompEqVols == 0) {
            if (tmpSpotVol.empty()){
                throw ModelException(method, "No spot or composite vols "
                                     "supplied for EQ vol "+vol->getName());
            }
            /* extend spot curve cover the simulation's timeline */
            return SRMUtil::extendVol(
                tmpSpotVolDate,
                vector<double>(tmpSpotVol.begin(), tmpSpotVol.end()),
                timeLine);
        }
        const DateTimeArray& baseIRDate = baseIR->getExtendedTimeLine();

        const DateTime& Today = timeLine.front();
        /***********************/
        /*  Extend the corrs   */
        /***********************/

        /* correlation EQ vs base IR factors */
        vector<double> RhoEqIr(
            SRMCorrelation::assetIr(InpCorrEqIr,
                                    *baseIR,
                                    SRMCorrelation::EXPONENTIAL_FACTORS));
    
        /******************************************/
        /*  Calc eq fowards to div payment dates  */
        /******************************************/
        // For each future discrete dividend need the fwd equity
        // price at that date (BEFORE the div is paid)
        // Exploiting the rather convenient facility to convert divs from dollar to yield
        // avoids the restriction (and error trap below) implicit in the SRM model.
        DividendListSP divList = 
            DividendCollector::divsBetweenDates(
                eqAsset.get(),
                Today,
                Today, // start
                baseIRDate.back(), // end. Seems sensible
                                // Note "timeLine" is not long enough to match SRM3!
                DividendCollector::DOLLAR_TO_YIELD);
#if 0
        // NOT CURRENTLY USED - not this causes problems if the payment dates are not
        // in increasing order.

        DateTimeArrayConstSP discreteDivPmtDates = divList->getPayDates();
        DateTimeArray eqFwdDates(discreteDivPmtDates->size());
        for(int j=0; j<eqFwdDates.size(); j++) {
            // Need the fwd price BEFORE DISCRETE DIVS ARE PAID
            // which is achieved by setting the time for the discrete div dates
            eqFwdDates[j] = DateTime((*discreteDivPmtDates)[j].getDate(),
                                     DateTime::BEFORE_EX_DIV_TIME);
        }
        DoubleArray fwdsToDivPmtDates(eqFwdDates.size(), 0.0);
        eqAsset->fwdValue(eqFwdDates, fwdsToDivPmtDates);
#endif

        /**************************************************/
        /* Create time line and extend ir spot vol curves */
        /**************************************************/
        vector<const DateTimeArray*> dtArrayVector(4);
        dtArrayVector[0] = &baseIRDate;
        dtArrayVector[1] = &divExDates; 
        dtArrayVector[2] = &compVolDate;
        dtArrayVector[3] = &vol->getCompVolResetDate();

        DateTimeArray TPDate(DateTime::merge(dtArrayVector));

        /* we require Today to be part of the datelist, so */
        /* check that first date in ZDates is indeed Today */

        if (baseIRDate[0] != Today) {
			throw ModelException(method, "Internal error : base IR date != today");
        }

        // fill in dates using hard coded ppy - really ?
        DateTimeArraySP filledTPDate = SRMUtil::fillInTimeLine(Today, TPDate, 12);
	    int NbTP = filledTPDate->size();

        /* temporary extension of EQ smile parameters              */
        /* note that a[i] applies for filledTPDate[i] <= t < filledTPDate[i+1] */
        vector<double> a1_Ext(NbTP-1);
        vector<double> a2_Ext(NbTP-1);
        vector<double> a3_Ext(NbTP-1);
        int i, k = 0;
        
		/* increase k until we find a smile date   */
        /* after today - tho checked in input      */
        for (i = 0; i <= NbTP - 2; i++) 
		{
            if (k < timeLine.size() - 1) 
			{
                if ((*filledTPDate)[i] >= timeLine[k+1])
				{
                    k++;
                }
                a1_Ext[i] = EqSmile_a1[k];
                a2_Ext[i] = EqSmile_a2[k];
                a3_Ext[i] = EqSmile_a3[k];
            } 
			else 
			{ /* use final parameter set */
                a1_Ext[i] = EqSmile_a1.back();
                a2_Ext[i] = EqSmile_a2.back();
                a3_Ext[i] = EqSmile_a3.back();
            }
        }// for i 

        /************************************************/
        /*  Finally, we can bootstrap the EQ spot vol  */
        /************************************************/
        vector<double> SpotEqVol(NbCompEqVols); // reserve some space

        multiFacSpotEqVol(SpotEqVol,      // only NbCompEqVols populated   
                          RhoEqIr, 
                          compVolDate,    /* EQ Option expriy  */
                          a1_Ext, a2_Ext, a3_Ext,       
                          0, // was fwdsToDivPmtDates but not currently used
                          divList,
                          *filledTPDate);
    
		DateTimeArray SpotVolDate(compVolDate);
        /* append the input spot vol list */
        SpotVolDate.insert(SpotVolDate.end(),
                           tmpSpotVolDate.begin(), tmpSpotVolDate.end());
        SpotEqVol.insert(SpotEqVol.end(),
                         tmpSpotVol.begin(), tmpSpotVol.end());

        /* extend spot curve to cover the simulation's timeline */
        return SRMUtil::extendVol(SpotVolDate, SpotEqVol, timeLine);
    } 
	catch (exception& e)
	{
        throw ModelException(e, method);
    }

}

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
void SRMEquityUtil::calcEqPerturb( 
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
    const vector<double>&  dt)          /* (I) intervals in grid, [0,n-1] */
{
    static const string method("SRMEquityUtil::calcEqPerturb");
    try {
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
    } catch (exception& e){
        throw ModelException(e, method);
    }
}
/****** CalcPerturb ******************************************************/


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
void SRMEquityUtil::multiFacSpotEqVol(
    vector<double>&       SpotEqVol,          /* (M) output  */
    const vector<double>& RhoEqIr,            /* (I) correl EQ/IR curve */
    const DateTimeArray&  ImpVolDate,         //(I) Modified vol->compVolDate
    const vector<double>& A_Eq,               /* (I) EQ smile param   */
    const vector<double>& B_Eq,               /* (I) EQ smile param   */
    const vector<double>& C_Eq,               /* (I) EQ smile param   */
    const DoubleArray*    fwdsToDivPmtDates,  // (I) not currently used
    DividendListConstSP   divList,            // (I) 
    const DateTimeArray&  TPDate)             /* (I) date of each pt.              */
{
    static const string method("SRMEquityUtil::multiFacSpotEqVol");
    try 
	{
        /* Quick checks first */
        int NbTP = TPDate.size();
        if (NbTP < 2) {
            throw ModelException(method, "Fewer than 2 time points supplied");
        }
        int IrNbFac = baseIR->numFactors(); // for ease
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
            loc_vol_x[i]    =   C_Eq[i] * A_Eq[i];  
            /* atm curvature wrt log-fwd-money */  
            loc_vol_xx[i]   =   C_Eq[i] * C_Eq[i] * B_Eq[i];
        }

        /* we have already checked that NbTp >= 2 */
        Time.back() = TPDate[0].yearFrac(TPDate.back());
        /* find the first time point on or after the FIRST imp vol date */
        int TPFirstImpVol = ImpVolDate.front().findUpper(TPDate);
        
		if (TPFirstImpVol == TPDate.size())
		{
            throw ModelException(method, "Vol Integration date "+
                                 ImpVolDate.front().toString()+
                                 " is beyond the input time line");
        }

        int TPLastImpVol = ImpVolDate.back().findUpper(TPDate);
        if (TPLastImpVol == TPDate.size()) 
		{
            throw ModelException(method, "Vol Integration date "+
                                 ImpVolDate.back().toString()+
                                 " is beyond the input time line");
        }

        /************************************************************/
        /* postion dividend payment dates (if any) on the time line */
        /************************************************************/
        // to do XXX
    
        /******************************************************/
    	/* presence of discrete dividends or smile parameters */
    	/******************************************************/
        /* detect if there are discrete dividends or non-zero smile parameters fail       */
        /* in the case where discrete dividends and non-zero smile parameters are inputed */
        bool isSmile = false; // used further below also
        for (i = 0; !isSmile && i < NbTP - 1;i++) 
		{
            isSmile = (!Maths::isZero(A_Eq[i]) ||
                       !Maths::isZero(B_Eq[i]) ||
                       (Maths::isPositive(C_Eq[i] - 1.0)));
        }
        // hasDollarDividend() will return true even if div = $0, but that is surely rare
        if (isSmile && divList->hasDollarDividend()) 
		{
            // XXX This should never trip since the divList is built with
            // XXX a conversion of all dollar divs to yield.
            throw ModelException(method, "Can't have non-zero smile parameters "
                                 " and discrete dividends\n");
        }

        /********************************************************************/
        /*                                                                  */
        /*          remove 1st order EQ smile contribution                  */
        /*                                                                  */
        /********************************************************************/

        /* first estimate zero expiry limit of implied vol = EQVolUnSmile[0] */ 
        /* the estimate is obtained from Eq. 8 of "Long-dated FX smile model",*/
        /* by assuming EQ spot vol constant and approximating RHS             */

        /* time to first option expiry */
        double tau_0 =  TPDate[0].yearFrac(ImpVolDate[0]);
        const DoubleArray& EqVol = vol->getCompVol(); // for ease
        /* approximation assumes const EQ spot vol until first expiry */
        double epsilon = 2. * tau_0 * EqVol[0] * EqVol[0] * 
            (loc_vol_xx[0] - 0.5 * loc_vol_x[0] * loc_vol_x[0]) / 3.0;

        if (epsilon < -0.5)
		{ /* sqr root of negative number */
            throw ModelException(method, "Can't calculate initial EQ spot vol");
        }  
        vector<double> EqVolUnSmile(NbTP);
        if (fabs(epsilon) < SRMConstants::SRM_TINY) 
		{
            EqVolUnSmile[0] = EqVol[0];
        } 
		else 
		{
            EqVolUnSmile[0] = EqVol[0] * 
                sqrt((- 1. + sqrt(1 + 2. * epsilon)) / epsilon);
        }

        /* extend the implied vols to full time line         */    
        double t_low   = 0.0;          /* initialise as today */
        double delta_t = tau_0;        /* first step size equals first expiry */
        double EqVolPrev = EqVolUnSmile[0];   /* guess for implied vol
                                                 extrapolated to zero expiry */
        double slope = EqVol[0] * EqVol[0] - 
                       EqVolUnSmile[0] * EqVolUnSmile[0];
        int NbImpVol = ImpVolDate.size();
        vector<double> EqVolExt(NbTP);
        int k;

        for (i = 0, k = 0 ; i <= NbTP - 1;i++) 
		{   
            if (k <= NbImpVol - 1) 
			{
                if (TPDate[i] >=  ImpVolDate[k]) 
				{
                    if (k < NbImpVol - 1) 
					{
                        EqVolPrev  =  EqVol[k];
                        slope      =  EqVol[k+1] * EqVol[k+1] - 
                                      EqVol[k] * EqVol[k];
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
                        EqVolPrev  = EqVol[NbImpVol-1];
                        k++; /* i.e., k = NbImpVol */
                    }       
                }
                /* interpolate linearly on squares of vols - in the short
                 * end, in particular, this is more accurate than
                 * interpolation on the vols  */
                EqVolExt[i] =  sqrt(EqVolPrev * EqVolPrev + 
                                    slope * (Time[i] - t_low) / delta_t);
            } 
			else 
			{ /* flat after last implied vol    */
                EqVolExt[i] =  EqVol[NbImpVol - 1];
            }
        }

        /* initialise dummy for integration */
        double temp1 = EqVolUnSmile[0] * EqVolUnSmile[0] * DeltaTime[0];
        vector<double> volA(NbTP);
        volA[1] = EqVolUnSmile[0];

        /* only volA is used outside this loop */
        vector<double> Perturb(NbTP);
        for (i = 1; i < TPLastImpVol; i++) 
		{   
            /* EqVol[i] etc is the implied vol at Time[i] - not Time[i+1] */
            double temp2 = (Time[i+1] * EqVolExt[i+1] * EqVolExt[i+1] - 
                            Time[i] * EqVolExt[i] * EqVolExt[i]);
        
            /* if smile then use the bootstrapping routine described in "Long-dated FX" */
            if (isSmile) 
			{
                calcEqPerturb(Perturb,     /* (O) [0,i]...[i-1,i] */
                              i,           /* (I) number of time points */
                              EqVolUnSmile,/* (I) implied EQ vol */
                              loc_vol_x,   /* (I) slope of local vol wrt log-money */
                              loc_vol_xx,  // (I) curvature of local vol wrt log-money
                              Time,        /* (I) time grid, Time[0] = 0 */
                              DeltaTime);  /* (I) intervals in grid       */
                temp2 /= ( 1.0 + Perturb[0]);
            } 
            
			if (Maths::isZero(DeltaTime[i])) 
			{
                EqVolUnSmile[i] = EqVolUnSmile[i-1];
            } 
			else 
			{
                if (Maths::isNegative(temp2)) 
				{
                    throw ModelException(method, "negative Eq vol in calibration");
                }
                EqVolUnSmile[i] = sqrt(temp2/DeltaTime[i]);
            }
            temp1 += EqVolUnSmile[i] * EqVolUnSmile[i] * DeltaTime[i];
            /* Time[i] > 0 because i > 0 and TPDate is strictly increasing  */
            /* volA[i] refers to Time[i]                                    */
            volA[i+1]  =  sqrt(temp1/Time[i+1]); 
        }
        /************* 1st order EQ smile contribution removed **************/

	    /******************************************/
	    /* bootstrap equity spot volatilities     */
	    /******************************************/
                
		// check dates and find offsets on TPDate timeline
		int kVol = 0;
		vector<int> expiryIndex(NbImpVol);
		vector<int> resetIndex(NbImpVol);
        
		// reset dates are kept only for backwards compatibility 
		const DateTimeArray& resetDate = vol->getCompVolResetDate(); 

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

        double Sigma = 0.0; // Eq Spot Vol in current bucket 
        bool EqCutOffFlag = VolBootstrapMode != SRMEquityDiffuse::NO_FAILURE_ALLOWED;
        bool EqCutOffLast = VolBootstrapMode == SRMEquityDiffuse::USE_LAST_LEVEL;
        double CutOffLevel = (EqCutOffLast || !EqCutOffFlag) ?
            SPOT_CUTOFF_R : EqCutOffLevel;

		int j;
        vector<double> SpotEqVolSoFar(NbTP);

        int TPIntBefore = 0;
        int TPIntSoFar = 0;

		vector<int>::iterator tmpExpiryBegin = expiryIndex.begin();
		vector<int>::iterator tmpExpiryEnd = tmpExpiryBegin;
		vector<int>::iterator tmpResetEnd = resetIndex.begin();

        // Pre-compute IrFwdRate for performance
        // Note : discretely compounded fwd rates * deltaTime
        vector<double> IrFwdRate(simpleFwdCurve(baseIR->getDiffYC(),TPDate));

        // Allocate memory outside loop
        vector<double> irVariance(NbImpVol);
        vector<double> irEqCovariance(NbImpVol);

        for (kVol = 0; kVol < NbImpVol; kVol++) 
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
            
			tmpExpiryEnd++; // because it need's to point one past last element
			vector<int> tmpExpiryIndex(tmpExpiryBegin, tmpExpiryEnd);

			baseIR->irAssetVariance(
				irVariance,        // (O)
				irEqCovariance,    // (O)  
				RhoEqIr,           // (I)
				DeltaTime,         // (I)  passed for convenience 
				TPDate,            // (I)
				IrFwdRate,         // (I)  pre-computed ir fwd rates
                tmpExpiryIndex,    // (I)  indices in TPDate
				*tmpResetEnd);     // (I)    

			tmpResetEnd++;

			double Poly[3] = {0.0, 0.0, 0.0};
        
 		    double EqOnlyVariance = 0.;
			// update EqOnlyVariance for next round
			for (i = 0; i < TPIntBefore; i++) 
			{               
                // constant coeffs due to EQ variance
                EqOnlyVariance += SpotEqVolSoFar[i] * SpotEqVolSoFar[i] * DeltaTime[i];
            }
			
			Poly[0] = EqOnlyVariance;

			for(j = 0; j <= kVol; j++)
			{
		    	Poly[0] += irVariance[j];
			}

			// note : this assumes spotEqVol constant between expiries
			for(j = 0; j < kVol; j++)
			{
				int eqIndex = (j > 0) ? expiryIndex[j-1] : 0;
                Poly[0] += irEqCovariance[j] * SpotEqVolSoFar[eqIndex];
			}
            
			// Subtract the target "unsmiled" variance
            Poly[0] -= volA[TPIntSoFar] * volA[TPIntSoFar] * 
                TPDate[0].yearFrac(TPDate[TPIntSoFar]);

            Poly[1] = irEqCovariance[kVol];
    
            // notice that TPDate[TPIntBefore] < ImpVolDate[kVol]     
            // i.e (Poly[2] non zero)                                         
            // because of earlier check on TPIntBefore==TPIntSoFar            
            Poly[2] = TPDate[TPIntBefore].yearFrac(ImpVolDate[kVol]);

            double discriminant = Poly[1] * Poly[1] - 4.* Poly[2] * Poly[0];
            bool VolTooLow = (discriminant < SRMConstants::SRM_TINY);
            if (discriminant >= SRMConstants::SRM_TINY) 
			{       
                /* using the larger root ensures a continuous spot EQ vol   */ 
                Sigma  = (-Poly[1] + sqrt(discriminant)) / 2. / Poly[2];
                /* Evaluate whether spotvol satisfies appropriate criterion */
                if (EqCutOffLast || !EqCutOffFlag) 
				{
                    VolTooLow  = (EqVol[kVol] > CutOffLevel * Sigma);
                } 
				else 
				{
                    VolTooLow = (Sigma < CutOffLevel);
                }
            }
            /* now implement cutoff actions as the case may be */
            if (VolTooLow) 
			{
                if (EqCutOffFlag) 
				{
                    if (EqCutOffLast) 
					{
                        if (kVol == 0) 
						{
                            string m("Unable to bootstrap at least one EQ vol "
                                     "point and therefore unable to"
                                     "cut off at last spot vol level");
                            throw ModelException(method, m);
                        }
                        for (int n = kVol; n < NbImpVol; n++) 
						{
                            SpotEqVol[n] = SpotEqVol[kVol-1];
                        }
                    } 
					else 
					{/* end of EqCutOffLast=true */
                        for(int n = kVol; n < NbImpVol; n++) 
						{
                            SpotEqVol[n] = CutOffLevel;
                        }
                    }
                    break;
                } 
				else 
				{/* end of EqCutOffFlag */
                    double pcVol = 100.*EqVol[kVol]/CutOffLevel;
                    string m("Problem in calculating EQ vol "
                             "at "+ImpVolDate[kVol].toString()+": less"
                             " than "+Format::toString(pcVol)+"% (level"
                             " determined by ratio to fwd vol)");
                    throw ModelException(method, m);
                }
            } /* end of VolTooLow */

            SpotEqVol[kVol] = Sigma;
            for (i = TPIntBefore; i < TPIntSoFar; i++) {
                SpotEqVolSoFar[i] = Sigma;
            }/* for i*/

            TPIntBefore = TPIntSoFar;
        }/*for kVol*/
    } 
	catch (exception& e)
	{
        throw ModelException(e, method);
    }
}/* multiFacSpotEqVol*/

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
void SRMEquityUtil::smileParamExtend() 
{
    static const string method("SRMEquityDiffuse::Utils::smileParamExtend");
    DoubleArray smileA1( vol->getSmileA1() );
    DoubleArray smileA2( vol->getSmileA2() );
    DoubleArray smileA3( vol->getSmileA3() );
    DateTimeArray smileDate( vol->getSmileDate() );
    try {
        int NbSrcParams = smileA1.size();
        const DateTimeArray& simDates = baseIR->getSimDates();
        /* Check that the set of source dates (except the last one) up to the */
        /* last time point is a subset of the output dates                    */
        // MAR: This is true by construction - but leave in for now
        int j;
        for (j = 0; j < NbSrcParams - 1; j++) {
            if (smileDate[j] < simDates.back()) {
                try{
                    smileDate[j].find(simDates);
                } catch (exception& e){
                    throw ModelException(e, method, "Internal error: EQ smile"
                                         " date not in timeline for equity " +vol->getName());
                }
            }
        }
        j = 0;  /* pointing to Src_xx[j] */ 
        /* Loop through each period TimePt[i] to TimePt[i+1] */
        int nbTimePts   = simDates.size();
        for (int i = 0; i < nbTimePts-1; i++) {
            if (j < (NbSrcParams-1)) {
                /* increase j until we're in right EQSmileDate interval */
                if (simDates[i] == smileDate[j]){
                    j++;
                }
            }
            /* LHS forward looking, RHS backward looking */
            EqSmile_a1[i] = smileA1[j];
            EqSmile_a2[i] = smileA2[j];
            EqSmile_a3[i] = smileA3[j];
        }        
    } catch (exception& e){
        throw ModelException(e, method);
    }
}/* EQ_SmileParamExtend */

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
void SRMEquityUtil::smileParamTrans() 
{
    static const string method("SRMEquityUtil::smileParamTrans");
    try {
        const DateTimeArray& simDates = baseIR->getSimDates();
        int NbParams   = simDates.size();
        for (int i=0; i<NbParams-1; i++) {
            double Sk      =  EqSmile_a1[i];        /* ATM skew          */
            double Crv     =  EqSmile_a2[i];        /* ATM crv           */
            double M       =  EqSmile_a3[i] + 1.0;  /* maximum variation */

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

                    double b0 = - fabs(Sk) + sqrt(Sk * Sk + 4.0 * fabs(Crv) *
                        (M - 1.0));
                    /* curvature, pos or neg */
                    b =  (b0 * b0)/
                        (4.0 * Crv); 
                    c =  sqrt(Crv/b);
                    a =  Sk / c;
                }
            }

			// Comment : checks on negative vols have been turned off because vols are now floored
			// temnporarily turned off - tests need to be updated
            /* check for negative curvature */
			/*
            if ( b < 0. ) 
			{
                throw ModelException(method, "Equity smile curvature is negative");  
            } 
			*/

            /* check for negative volatility */
			
			/*
            if ( 1. - fabs(a) + b < 0. ) 
			{
                throw ModelException(method, "Negative or zero equity volatility");  
            } 
			*/

            /*
            if ( a * a - 1. - 2. * b > - SRMConstants::SRM_TINY ) 
			{ 
                throw ModelException(method, "Negative or zero equity volatility");
            }
			*/

            EqSmile_a1[i]   =   a;
            EqSmile_a2[i]   =   b;
            EqSmile_a3[i]   =   c;
        }        
    } catch (exception& e){
        throw ModelException(e, method);
    }
}/* EQ_SmileParamTrans */

void SRMEquityUtil::computeOriginalSpotPrices(const DateTimeArray& dates, DoubleArray & spots) const
{
    spots.resize(dates.size());
    eqAsset->fwdValue(dates, spots);
    for (int i = 0; i < spots.size(); ++i)
        spots[i] = log(spots[i]);
}

DRLIB_END_NAMESPACE
