//----------------------------------------------------------------------------
//
//   Group       : EDG Derivatives Research
//
//   Filename    : SRMInitialGuess.cpp
//
//   Description : Initial Guess for SRM
//
//----------------------------------------------------------------------------

#include "edginc/config.hpp"
#include "edginc/Function.hpp"
#include "edginc/Asset.hpp"
#include "edginc/MarketDataFetcher.hpp"
#include "edginc/MarketDataFetcherLN.hpp"
#include "edginc/NonPricingModel.hpp"
#include "edginc/VolProcessedBS.hpp"
#include "edginc/LinearStrikeTSVolRequest.hpp"
#include "edginc/LocVolRequest.hpp"
#include "edginc/VolProcessedDVF.hpp"
#include "edginc/Optimizer.hpp"
#include "edginc/Addin.hpp"
#include "edginc/Maths.hpp"
#include "edginc/VanillaGrid.hpp"
#include "edginc/Calibrator.hpp"
#include "edginc/Spline.hpp"
#include "edginc/IRVol.hpp"
#include "edginc/MarketDataFetcherSRM.hpp"
#include "edginc/SwapMaturityVolRequest.hpp"
#include "edginc/SRMInitialGuess.hpp"
#include "edginc/SRMUtil.hpp"
#include "edginc/SRMEQVol.hpp"
#include "edginc/SRMRatesUtil.hpp"
#include "edginc/SRMRatesHJMUtil.hpp"
#include "edginc/SRMEquityUtil.hpp"

DRLIB_BEGIN_NAMESPACE

// Class for objective function

// This function calculates the square difference between input function and 
// hyperbolic local volatility model
class SRMIV_ObjFuncDupire : public MFunctionND{
    
    vector<double> logMoneyFwd,
		           normTargetVol; // normalised target vol

	double t, stdDevCutoff; // don't calibrate in the wings 

public:

    virtual void operator()(const CDoubleArray&  inParams,  // in
                            CDoubleArray&        f) const   // out
	{

        static const string method = "SRMIV_ObjFuncDupire::operator()";

        if (inParams.size() != 3) 
		{
			throw ModelException(method, "There must be exactly 3 parameter inputs");
		}

        if (f.size() != 1) 
		{
			throw ModelException(method, "Function must have size = 1");
		}

		int nbStrikes = logMoneyFwd.size(),
			nbUsedStrikes = 0,
			i = 0;

		double aSmile = 0.,
			   bSmile = 0.,
			   cSmile = 0.,
			   x = 0., y = 0., z = 0., w = 0., 
			   sigma = 0., 
			   sumError = 0., weight = 0.;

        aSmile = inParams[0];
	    bSmile = inParams[1];
	    cSmile = inParams[2];

		if (cSmile < 0.000001)
		{
			throw ModelException(method, "Internal error - scale parameter is too small");
		}

        const double SMILE_CUTOFF = 100.0;

        if (nbStrikes < 3) 
		{
			throw ModelException(method, "Need at least 3 strikes");
		}
		
		const vector<double> *logMoney = &logMoneyFwd;
		
		for (i = 0; i < nbStrikes; i++)
		{
			// note : different sign convention than for diffusion model 
		    // - hence the "-"
			x = - logMoneyFwd[i];
			w = cSmile * x;
			
			if (w < SMILE_CUTOFF && w > - SMILE_CUTOFF) 
			{
				y = exp(w);
				z = 1.0/y;
				/* loc vol = spotvol * (1 + A * Tanh(C * x) + 
				*                   B * ( 1 - 1/Cosh(C * x)) ) */
				sigma = 1.0 + aSmile * (y - z)/(y + z)  
					+ bSmile * (1.0 - 2.0/(y + z));
			} 
			else if (w <= - SMILE_CUTOFF) 
			{
				/* y = 0 */
				sigma = 1.0 - aSmile + bSmile;
			} 
			else 
			{
				/* z = 0 */
				sigma = 1.0 + aSmile + bSmile;
			}

            sumError += Maths::square(sigma - normTargetVol[i]);
			nbUsedStrikes++;
		} // else do nothing

		if (nbUsedStrikes > 0)
		{
		    f[0] = sumError / nbUsedStrikes;
		}
		else
		{
			f[0] = 0.;
		}
    }

    SRMIV_ObjFuncDupire(const RangeArray& defRanges,
                        int nbVars,
						int nbFuncs,
						vector<double> logMoneyFwd,
						vector<double> normTargetVol,
						double stdDevCutoff,
						double t):  
	MFunctionND(nbVars, nbFuncs, defRanges),
	logMoneyFwd(logMoneyFwd),
	normTargetVol(normTargetVol),
	stdDevCutoff(fabs(stdDevCutoff)),
	t(t){}

};

// interpolation using a second-order polynimial 
DoubleArray localLeastSquareFit( DoubleArray x,  /* logMoneyFwd */
                                 DoubleArray y) /* impVar */
{
    static const string method = "localLeastSquareFit";

    int nbStrike = x.size();


    int j = 0;


    if( nbStrike < 3)
    {
        throw ModelException(method, " we need at least 3 points " );
    }

    double A = 0.,
           B = 0.,
           C = 0.,
           D = 0.,
           E = 0.,
           F = 0.,
           G = 0.,
           H = 0.;

   
    for( j = 0 ; j < nbStrike ; j++)
    {
        A += ::pow(x[j], 4) ;//slow
        B += ::pow(x[j], 3) ;
        C += Maths::square(x[j]);
        D += x[j] ;
        E += 1.;
        F += y[j];
        G += x[j] * y[j];
        H += Maths::square(x[j]) * y[j];
    }
   

    double a = 0.,
           b = 0.,
           c = 0.;

   
    double denom = (B*E - C*D)*(A*C - B*B) -
           (A*D - B*C)*(B*D - C*C);

    double num = (B*F - C*G)*(A*C - B*B) -
           (A*G - B*H)*(B*D - C*C);

    if(denom == 0.)
    {
        if( num != 0.)
        {
            throw ModelException(method , " can't fit the data with a second-order polynomial ");
        }
        else
        {
            throw ModelException(method , " there is not a unique solution to the least square fit problem ");
        }
    }
    else
        {
            c = num/denom;
            b = ((A*G - B*H) - (A*D - B*C)*c)/(A*C - B*B);
            a = (H - b*B - c*C)/A;
        }
    

    DoubleArray result( 3);

   
    result[0] = a;
    result[1] = b;
    result[2] = c;


    return result;
}


CDoubleMatrixSP SRMInitialGuess::smileDupireRates() 
{
	static const string method = "SRMInitialGuess::smileDupireRates";
		
	try 
	{
		
		DateTime 
			refDate = market->GetReferenceDate(),
			lastDate = expiries->back()->toDate(refDate);
		
		int nbMaturities = expiries->size(),
			nbVar = 3;   // nb parameters per expiry 
		// - if you change this, then x[] and guess[] below need to change too
		
        const int nbStrikes = 7; 
        double relMoney[nbStrikes] = {-0.25,-0.15, -0.05, 0., 0.05, 0.15, 0.25};

		// outParams contains spotvol + (a1,a2,a3) smile parameters
		CDoubleMatrixSP outParams(new CDoubleMatrix(nbMaturities,nbVar+1));

		DoubleMatrix  
			impVar(nbMaturities, nbStrikes),    // implied variance
			impVar_T(nbMaturities, nbStrikes),  // derivative w.r.t. expiry
			logMoneyFwd(nbMaturities, nbStrikes),
			effLocalVol(nbMaturities, nbStrikes),
			smileAdjust(nbMaturities, nbStrikes);
		
		// pick our vol from the market
		// MarketDataFetcherSP mdf(new MarketDataFetcherLN(volType)); 
		// bool dummy = mdf->setStochasticYieldCurveMode(true);

		// copied from SRMFXVolSimple
        StringArray fxVolType(1, "SRMFX::VolSimple");
        MarketDataFetcherSP mdfSRM(
			new MarketDataFetcherSRM(
			                "IRCalib::Smile2Q",  // irCalibSmileType
							"IRCalib::Model1FL", // irCalibModelType
							true,                // getSwaptionVols?
							false,               // useIRVolPair?
							volType,             // already have that one
							fxVolType,           // fxVolType - not needed  
							"CRCalib::Smile",    // credit smile - not needed
                            "FlatCDSSpotVol"));  // cds vol type - not needed


		NonPricingModel model(mdfSRM);
		
		// equity data
		asset.getData(&model, market);

		// now look for the yield curve
		const string& ycName = asset->getYCName();

		// iso code just in case
        if (!ycName.empty()){
            const string& YCIsoCode = 
                market->getYieldCurveISOCode(ycName);
        }

		// look in market
		// ? MarketObjectSP mo = mdfSRM->fetch(&(*market), YC, ycType, &model);
        MarketObjectSP moYC (market->GetData(&model,ycName,IYieldCurve::TYPE)); 
		moYC->getMarket(&model,&(*market));

		// now get the yield curve
		IYieldCurveConstSP domYC = IYieldCurveConstSP::dynamicCast(moYC);

        // get EQ-IR correlation
	    const CClassConstSP assetType = CClass::forName(volType);
		const CClassConstSP ycType = CClass::forName("IYieldCurve");

        bool found = market->hasCorrelationData(asset->getName(), ycName);
		
		if ( !found ) 
		{
			throw ModelException(method, "No equity-IR correlation available");
		}

        const string& corrName = market->getCorrelationName(asset->getName(), ycName);
        CorrelationBaseSP corrBase(model.getCorrelation(
			                                     corrName,
                                                 assetType,
                                                 ycType,
                                                 Correlation::TYPE, 
                                                 &(*market)));
	    
		CorrelationSP corr = CorrelationSP::dynamicCast(corrBase);
		// must be vector for comptability with factorFXCovariance()
	    vector<double> rho(1,corr->getCorrelation());

        // rates calibration 
		const string calibrationStyle = "FIX";
        const string calibrationMaturity = "10Y";

        ExpirySP calibExpiry(new MaturityPeriod(calibrationMaturity));
        CVolRequestSP request(new SwapMaturityVolRequest(calibExpiry.get(), 
                                                         calibrationStyle));

        // get hold of processed vol, domestic
        CVolProcessedSP volProcessed(domYC->getProcessedVol(request.get()));
        VolProcessedBSIRSP domProcessedVol(VolProcessedBSIRSP::dynamicCast(
                                               volProcessed));
		FXAssetConstSP fxDummy;

		// VNFM type constructor
        SRMRatesHJMUtilSP domIR (new SRMRatesHJMUtil(
			                             refDate, 
			                             fxDummy,
                                         "1F-STANDARD", 
									     domProcessedVol,
                                         domYC,
                                         domYC, 
										 true,       // skip bad vols 
										 "1Y",       // corrSwapStart 
                                         "10Y",      // corrSwapMat
										 "Act/365F", // corrSwapDCC 
										 "6M"));     // corrSwapFreq

        const DateTimeArray irDates = domIR->getSwaptionExpiryDates();

		if (irDates.size() < 1)
		{
			throw ModelException(method, "irDates is empty");
		}
		
		DateTimeArray matDates(nbMaturities); // period begin
		DateTimeArray ratesDates(nbMaturities+1);

		// calculate f and its derivatives
		int i = 0, 
			t = 0;  

		vector<double> deltaTime(nbMaturities);
		vector<double> tNow(nbMaturities);

		vector<double> fwdPrice(nbMaturities);

        double spotPrice = asset->fwdValue(refDate);

		ratesDates[0] = refDate;
		for (t = 0; t < nbMaturities; t++)
		{
			ratesDates[t+1] = (*expiries)[t]->toDate(refDate); 
			matDates[t] = (*expiries)[t]->toDate(refDate);
			tNow[t] = refDate.yearFrac(matDates[t]);
			deltaTime[t] = (t > 0) ? tNow[t] - tNow[t-1]: tNow[0];
			
			if (deltaTime[t] < SRMConstants::SRM_TINY)
			{
				throw ModelException(method, "Expiries are equal or too close");
			}
			
			fwdPrice[t] = asset->fwdValue(matDates[t]);
			
			if (fwdPrice[t] < SRMConstants::SRM_TINY)
			{
				throw ModelException(method, "forward price too small");
			}
			
		}

		// generate forward rates
        // attention: IrFwdRate is a fwdRate multiplied by deltaT 
        vector<double> fwdRate(simpleFwdCurve(domIR->getDiffYC(), ratesDates));

        const DateTimeArray& domIRDate = domIR->getExtendedTimeLine();
		vector<double> irSpotVol = domIR->getSpotVols();

		// not needed here but works
		/*
        // now build our combined timeline
        vector<const DateTimeArray*> dtArrayVector(2);
        dtArrayVector[0] = &irDates;
        dtArrayVector[1] = &matDates;
        DateTimeArray timeLine(DateTime::merge(dtArrayVector));

		// from PVDiv
        CAssetConstSP assetPtr(asset.getSP());
        SRMEquityUtil srmEquity(domIR,
			            assetPtr,
						rho,      // InpCorrEqIr
						true,     // isDomesticEquity
						0.,       // InpCorrEqFX
						0.,       // fxVol
						SRMEquityDiffuse::USE_LAST_LEVEL, // VolBootstrapMode
						0.3);     // EqCutOffLevel
        */

		// controls how many points from the vol surface are used 
		const double NUM_STD_DEV = 5.;  
		
		// expiry-dependent cutoff
		vector<double> stdDevCutoff(nbMaturities);

		// get implied volatilities 
		for (i = 0; i < nbStrikes; i++)
		{
			for (t = 0; t < nbMaturities; t++)
			{ 
                double scale = 0.3333 * (1. + 5. * t/nbMaturities + 3.0 * Maths::square(t/nbMaturities));
                double fwdStrike = (1. + scale * relMoney[i]) * fwdPrice[t];

				// Move to Rel. strike - is StrikeTS needed ?
				LinearStrikeTSVolRequest req(fwdStrike, refDate, lastDate, false);
				CVolProcessedBSSP volproc(asset->getProcessedVol(&req));
				logMoneyFwd[t][i] = log(fwdPrice[t]/fwdStrike);// Move to Rel. strike
				double interpVol = volproc->CalcVol(refDate, matDates[t]);
                impVar[t][i] = tNow[t] * ::pow(interpVol,2);
			}
		}

		// builds NRSpline with left and right 2nd derivatives = 0
		NRSpline spline(2, 0, 2, 0);
        //CubicShapePresSplineInterpECnd splineShapePre(false);

		double x;	
		
		vector<double> estEqVol(nbMaturities);
        double integral = 0.;

		vector<double> omega(nbMaturities);

		double 
			temp0 = 0., temp1 = 0., temp2 = 0., temp3 = 0., temp4 = 0.,
			temp5 = 0., temp6 = 0., ratio = 0., fac1 = 0., fac2 = 0., 
			spotVol = 0., a = 0., b = 0., c = 0., skewDrift = 0.;

		for (t = 0; t < nbMaturities; t++)
		{	
			ratio = (t>0) ? tNow[t-1]/tNow[t] : 0.;
			fac1  = 0.5 * (1. + ratio); // factors appearing after averaging over period
			fac2  = (1. + ratio + ratio * ratio)/3.;

			// SPLINE interpolation between [xData[n]-xData[n-1]], 
			// using {(xData[i], yData[i])} points (i between 1 and n-1)			
			Interpolator::InterpolantConstSP impVarInterp = 
					spline.computeInterp(logMoneyFwd[t], impVar[t], nbStrikes);   

			Interpolator::InterpolantConstSP impVarInterpPrevious = (t == 0) ? impVarInterp :
                    spline.computeInterp(logMoneyFwd[t-1], impVar[t-1], nbStrikes);


			for (i = 0; i < nbStrikes; i++)
			{
                DoubleArray slidingLogMoneyFwd(3);
                DoubleArray slidingImpVar(3);

                DoubleArray slidingLogMoneyFwdPrevious(3);
                DoubleArray slidingImpVarPrevious(3);

                if(i == 0)
                {
                    slidingLogMoneyFwd[0] = logMoneyFwd[t][i];
                    slidingLogMoneyFwd[1] = logMoneyFwd[t][i+1];
                    slidingLogMoneyFwd[2] = logMoneyFwd[t][i+2];

                    slidingLogMoneyFwdPrevious[0] = (t == 0) ? logMoneyFwd[t][i] : logMoneyFwd[t-1][i];
                    slidingLogMoneyFwdPrevious[1] = (t == 0) ? logMoneyFwd[t][i+1] : logMoneyFwd[t-1][i+1];
                    slidingLogMoneyFwdPrevious[2] = (t == 0) ? logMoneyFwd[t][i+2] : logMoneyFwd[t-1][i+2];

                    slidingImpVar[0] = impVar[t][i];
                    slidingImpVar[1] = impVar[t][i+1];
                    slidingImpVar[2] = impVar[t][i+2];

                    slidingImpVarPrevious[0] = (t == 0) ? impVar[t][i] : impVar[t-1][i];
                    slidingImpVarPrevious[1] = (t == 0) ? impVar[t][i+1] : impVar[t-1][i+1];
                    slidingImpVarPrevious[2] = (t == 0) ? impVar[t][i+2] : impVar[t-1][i+2];
                }

                if(i == nbStrikes-1)
                {
                    slidingLogMoneyFwd[0] = logMoneyFwd[t][i-2];
                    slidingLogMoneyFwd[1] = logMoneyFwd[t][i-1];
                    slidingLogMoneyFwd[2] = logMoneyFwd[t][i];

                    slidingLogMoneyFwdPrevious[0] = (t == 0) ? logMoneyFwd[t][i-2] : logMoneyFwd[t-1][i-2];
                    slidingLogMoneyFwdPrevious[1] = (t == 0) ? logMoneyFwd[t][i-1] : logMoneyFwd[t-1][i-1];
                    slidingLogMoneyFwdPrevious[2] = (t == 0) ? logMoneyFwd[t][i] : logMoneyFwd[t-1][i];

                    slidingImpVar[0] = impVar[t][i-2];
                    slidingImpVar[1] = impVar[t][i-1];
                    slidingImpVar[2] = impVar[t][i];

                    slidingImpVarPrevious[0] = (t == 0) ? impVar[t][i-2] : impVar[t-1][i-2];
                    slidingImpVarPrevious[1] = (t == 0) ? impVar[t][i-1] : impVar[t-1][i-1];
                    slidingImpVarPrevious[2] = (t == 0) ? impVar[t][i] : impVar[t-1][i];
                }

                if(i > 0 && i < nbStrikes-1)
                {
                    slidingLogMoneyFwd[0] = logMoneyFwd[t][i-1];
                    slidingLogMoneyFwd[1] = logMoneyFwd[t][i];
                    slidingLogMoneyFwd[2] = logMoneyFwd[t][i+1];

                    slidingLogMoneyFwdPrevious[0] = (t == 0) ? logMoneyFwd[t][i-1] : logMoneyFwd[t-1][i-1];
                    slidingLogMoneyFwdPrevious[1] = (t == 0) ? logMoneyFwd[t][i] : logMoneyFwd[t-1][i];
                    slidingLogMoneyFwdPrevious[2] = (t == 0) ? logMoneyFwd[t][i+1] : logMoneyFwd[t-1][i+1];

                    slidingImpVar[0] = impVar[t][i-1];
                    slidingImpVar[1] = impVar[t][i];
                    slidingImpVar[2] = impVar[t][i+1];

                    slidingImpVarPrevious[0] = (t == 0) ? impVar[t][i-1] : impVar[t-1][i-1];
                    slidingImpVarPrevious[1] = (t == 0) ? impVar[t][i] : impVar[t-1][i];
                    slidingImpVarPrevious[2] = (t == 0) ? impVar[t][i+1] : impVar[t-1][i+1];

                    
                }


				x = logMoneyFwd[t][i];
				if (t==0) 
				{   
					impVar_T[t][i] = impVar[t][i] / deltaTime[t];
				}
				else 
				{
					// slope along lines of constant log-moneyness 
					impVar_T[t][i] = (impVar[t][i] - impVarInterpPrevious->value(x,0)) / deltaTime[t];
				}

                temp0 = impVar[t][i];
                try{
                    DoubleArray  temp =  localLeastSquareFit(slidingLogMoneyFwd, slidingImpVar);
                    temp1 = 2. * temp[0] * x + temp[1];
                    temp2 = 2. * temp[0];
                }catch (exception&) {
                    temp1 = impVarInterp->value(x,1); // f_x
                    temp2 = impVarInterp->value(x,2); // f_xx
                }
               
				
				// time t
				if (t>0)
				{	
					temp3 = ( -1./temp0 + x * x / (temp0 * temp0) - 0.25) * temp1 * temp1;
					temp4 = 1. - x * temp1 / temp0 + 0.25 * temp3 + 0.5 * temp2;
					
					// previous
					temp0 = impVarInterpPrevious->value(x,0);
                    try{
                         DoubleArray  tempPrevious =  localLeastSquareFit(slidingLogMoneyFwdPrevious, slidingImpVarPrevious);
                        temp1 = 2. * tempPrevious[0] * x + tempPrevious[1];
                        temp2 = 2. * tempPrevious[0];
                    }catch (exception&) {
                        temp1 = impVarInterpPrevious->value(x,1); // f_x
                        temp2 = impVarInterpPrevious->value(x,2); // f_xx
                    }
					temp3 = ( -1./temp0 + x * x / (temp0 * temp0) - 0.25) * temp1 * temp1;
					temp4 += 1. - x * temp1 / temp0 + 0.25 * temp3 + 0.5 * temp2;
					
					//average : amounts to linear interpolation for integrand
					temp4 /= 2.;
				}
				else // t = 0
				{
					
					temp3 = ( -fac1/temp0 + x * x / (temp0 * temp0) - 0.25 * fac2) 
						    * temp1 * temp1;
					temp4 = 1. - x * temp1 / temp0 + 0.25 * temp3 + 0.5 * fac1 * temp2;
					
				}

				// this may be needed if surface allows for arbitrage
				temp4 = Maths::max(temp4,0.1); 
			    temp4 = 1./temp4;

				smileAdjust[t][i] = temp4;
				temp4 *= impVar_T[t][i];
				temp4 = sqrt(Maths::max(temp4, 0.));	 
				// we override this below
				// applies in [ T[t-1],T[t] []
				effLocalVol[t][i] = temp4; 
			}

			// function value at x = 0	
			stdDevCutoff[t] = NUM_STD_DEV * impVarInterp->value(0.,0); 
		
     	    vector<double> NoEqInt(1);
		    vector<double> WithEqInt(1);
		    vector<double> WithEqIntCoeff(1);
		    vector<double> irA(1);
		
			// crude integration but will do for now
			int u;
		
			//DateTime maxDate = refDate.rollDateInMonths(3650);
			DateTime dateToUse;
			
			// in the first step, equity spot vol is unknown, so we calculate only the 
			// coefficient WithEqIntCoeff

			// from SRMEquityUtil
			u = t; 
			dateToUse = matDates[u]; 
			int irVolIndex = Maths::max(0L, dateToUse.findLower(domIRDate));
			double spotVol = irSpotVol[irVolIndex]; 
			
			// for accuracy, add fine-grained loop
			int loops = 50;
			double loopDt = deltaTime[u] / loops;
			double loopFwd = fwdRate[u] / loops;
			
			int k;
			for (k = loops - 1; k>=0; k--)
			{	
				// from SRMEquityUtil
				// note that this gives us value of AFactor at previous expiry date
				// CAREFUL : fwdRate is scaled by dt, effectively
				domIR->aFactor(loopFwd, loopDt, irA);
				
				// variance of IR factors
				vector<double>::iterator currentNbInt = NoEqInt.begin();
				currentNbInt = domIR->factorVariance(currentNbInt, irA,
					spotVol, loopDt);
				
				vector<double>::iterator withEqIntPos = WithEqIntCoeff.begin();
				// covariance between EQ and IR - I know the routine says "FX" but it does the same thing! 
				// Note : Rasmussen (2006) uses opposite sign convention of code
				withEqIntPos = domIR->factorFXCovariance(true, // "+" and not "-" contribution
					withEqIntPos, irA, rho, spotVol, loopDt);
			}
			

			// note that the AFactor calculation forces us to step backwards 
			for(u = t - 1; u >= 0; u--)
			{
				// from SRMEquityUtil
				dateToUse = matDates[u]; 
				int irVolIndex = Maths::max(0L, dateToUse.findLower(domIRDate));
				double spotVol = irSpotVol[irVolIndex]; 
				
				// for accuracy, add fine-grained loop
				int loops = 50;
				double loopDt = deltaTime[u] / loops;
				double loopFwd = fwdRate[u] / loops;

				int k;
				for (k = loops - 1; k>=0; k--)
				{	
					// from SRMEquityUtil
					// note that this gives us value of AFactor at previous expiry date
					// CAREFUL : fwdRate is scaled by dt, effectively
					domIR->aFactor(loopFwd, loopDt, irA);
					
					// variance of IR factors
					vector<double>::iterator currentNbInt = NoEqInt.begin();
					currentNbInt = domIR->factorVariance(currentNbInt, irA,
						spotVol, loopDt);
					
					vector<double>::iterator withEqIntPos = WithEqInt.begin();
					// covariance between EQ and IR - I know the routine says "FX" but it does the same thing! 
					// Note : Rasmussen (2006) uses opposite sign convention of code
					withEqIntPos = domIR->factorFXCovariance(true, // "+" and not "-" contribution
						withEqIntPos, irA, rho,
						spotVol * estEqVol[u], loopDt);
				}
			}

			// actually 2 Omega 
			omega[t] = - integral;
			
			integral = NoEqInt[0] + 2. * WithEqInt[0];
			
			omega[t] += integral;
			omega[t] /= deltaTime[t];
			
			// then make better estimate of vol
			for (i = 0; i < nbStrikes; i++)
			{
				temp0 = impVar_T[t][i];
				temp1 = 0.;
                // set to zero because we move along x=constant
				temp3 = omega[t]; // written like this for testing
				temp4 = temp0 - temp1 - temp3; 
                temp5 = smileAdjust[t][i]; 
				double C = - Maths::max(temp4 * temp5, 0.);
				double B = - 2. * WithEqIntCoeff[0] * temp5 / deltaTime[t];
				double A = 1.;
				double discr = B * B - 4. * A * C; // positive by construction
                double vol = (-B+sqrt(discr))/ (2. * A);	 
				effLocalVol[t][i] = vol; // applies in [ T[t-1],T[t] []
			}

			// calculate spot vol for next round
			Interpolator::InterpolantConstSP volInterpolator = 
				spline.computeInterp(logMoneyFwd[t], effLocalVol[t], nbStrikes);   
			
			// estimate spot vol as effLocalVol for log moneyness x=0 
			estEqVol[t] = volInterpolator->value(0.,0); 

		}
		
		DoubleArray params(nbVar); // parameters
		DoubleArray guess(nbVar);
		
		vector<double> tempLogFwd(nbStrikes);
		vector<double> tempLogSpot(nbStrikes);
		vector<double> normTargetVol(nbStrikes);
		
		for (t = 0; t < nbMaturities; t++)
		{
			// use only available vols
			Interpolator::InterpolantConstSP volInterpolator = 
				spline.computeInterp(logMoneyFwd[t], effLocalVol[t], nbStrikes);   
			
			// estimate spot vol as effLocalVol for log moneyness x=0 
			spotVol = volInterpolator->value(0.,0); 

			spotVol = Maths::collar(spotVol,.75,0.05);

			for (i = 0; i < nbStrikes; i++)
			{
				tempLogFwd[i] = logMoneyFwd[t][i];
				normTargetVol[i] = effLocalVol[t][i]/spotVol;
			}

			// make sure that less curvature is allowed in long end
			double tExp = tNow[t];
			RangeArray ranges(nbVar);
			ranges[0] = RangeSP(new Range (OpenBoundary(-10.), OpenBoundary(10.))); // A
			ranges[1] = RangeSP(new Range (OpenBoundary(0.0),  OpenBoundary(10.))); // B
			ranges[2] = RangeSP(new Range (OpenBoundary(0.1),  OpenBoundary(50.))); // C

			SRMIV_ObjFuncDupire errorFunc(
				ranges,
				nbVar,
				1,
				tempLogFwd,
				normTargetVol,
				stdDevCutoff[t],
				tNow[t]);
			
			QuasiNewton qn;
            		
			// we use previous parameter values as starting point
			// note sign convention for skew
			double midVal = volInterpolator->value(0.0,0);
			double upVal = volInterpolator->value(0.05,0); 
			double downVal = volInterpolator->value(-0.05,0); 

			c = Maths::max(1., 20./(1.+tExp));  
			guess[2] = c; 
			guess[0] = - (upVal - downVal) / (0.1 * c); // note sign convention 	 
			guess[1] = Maths::max(0.5, (upVal + downVal - 2. * midVal) / (0.01 * c * c));  

			qn.minimize(errorFunc,guess,params);
			
			c = params[2];
			b = params[1]; // curvature
			a = params[0]; // skew
			
			// floor b
			if (b < 0.)
			{
				b = 0.;
			}
			
			(*outParams)[t][0] = spotVol; // pseudo comp vol
			(*outParams)[t][1] = a; 
			(*outParams)[t][2] = b; 
			(*outParams)[t][3] = c; 
		}	
		
		return outParams;
		
	}
	catch(exception& e)
	{
		throw ModelException(e, method);
	}
}

///////////////////////////////////////////////////////////////////////////////////////////////
// author : henrik rasmussen
//          copied from initial guess addin
//          unlike latter, this doesn't take "strike" as input, these are generated internally 
///////////////////////////////////////////////////////////////////////////////////////////////

class SRMEquityCalib : public CObject
{

public:

    static CClassConstSP const TYPE;

    MarketDataSP market;
    CAssetWrapper asset;
    ExpiryArraySP expiries;

private:

    CDoubleMatrixSP calibrate()
	{
        static const string method = "SRMEquityCalib::calibrate";
		
		try
		{		
			DoubleArray dummyStrike(1); // not used
			dummyStrike[0] = 1.;

			SRMInitialGuess guess(market,
								  asset,
								  expiries,
								  dummyStrike, 
								  SRMInitialGuess::VOLPREF,  // volType
								  SRMInitialGuess::DUPIREMET,
								  SRMInitialGuess::SPOTVOL); // volMap			

			int nbMat = guess.getNbExpiries();
			
			CDoubleMatrixSP result(new CDoubleMatrix(nbMat,4));
			
			DoubleArraySP atmVol  = guess.getAtmVol();
			DoubleArraySP smileA1 = guess.getSmileA1();
			DoubleArraySP smileA2 = guess.getSmileA2();
			DoubleArraySP smileA3 = guess.getSmileA3();

			for (int iMat=0; iMat<nbMat;iMat++)
			{
				(*result)[iMat][0] = (*atmVol)[iMat];
				(*result)[iMat][1] = (*smileA1)[iMat];
				(*result)[iMat][2] = (*smileA2)[iMat];
				(*result)[iMat][3] = (*smileA3)[iMat];
			}

			return result;
		} 
		catch(exception& e)
		{
			throw ModelException(e, method);
		}
    }

    static IObjectSP callCalibrate(SRMEquityCalib* params) 
	{
        return params->calibrate();
    }

    // for reflection
    SRMEquityCalib():
        CObject(TYPE){}

    // Invoked when Class is 'loaded'
    static void load(CClassSP& clazz)
    {
        REGISTER(SRMEquityCalib, clazz);
        SUPERCLASS(CObject);
        EMPTY_SHELL_METHOD(defaultSRMEquityCalib);
        FIELD(asset, "");
        FIELD(market, "");
        FIELD(expiries, "");

        Addin::registerClassObjectMethod(
            "SRMEquityCalib",
            Addin::MARKET,
            "compute initial values for SRM calibration",
            TYPE,
            false,
            Addin::expandMulti,
            (Addin::ObjMethod*)callCalibrate);
    }

    static IObject* defaultSRMEquityCalib(){
        return new SRMEquityCalib();
    }
};

CClassConstSP const SRMEquityCalib::TYPE =
CClass::registerClassLoadMethod("SRMEquityCalib", typeid(SRMEquityCalib), load);


bool SRMEquityCalibLoad() {
    return (SRMEquityCalib::TYPE != 0);
}

DRLIB_END_NAMESPACE
