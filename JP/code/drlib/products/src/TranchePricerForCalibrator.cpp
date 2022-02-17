//----------------------------------------------------------------------------
//
//   Group       : Credit QR
//
//   Filename    : TranchePricerForCalibrator.cpp
//
//   Description : special tranche pricer used in the calibration methodology for SCID model 
//
//   Author      : Adrian Bozdog
//
//----------------------------------------------------------------------------


#include "edginc/config.hpp"
#include "edginc/TranchePricerForCalibrator.hpp"
#include "edginc/SCIDparameters.hpp"
#include "edginc/DoubleMatrix.hpp"
#include "edginc/SwapTool.hpp"
#include "edginc/QuadraticProg.hpp"
#include "edginc/Maths.hpp"
#include <numeric> 
#include <algorithm> 

DRLIB_BEGIN_NAMESPACE

TranchePricerForCalibrator::~TranchePricerForCalibrator(){}

/* Default constructor */
TranchePricerForCalibrator::TranchePricerForCalibrator() : CInstrument(TYPE){
};

string TranchePricerForCalibrator::discountYieldCurveName() const {
    return discount.getName();
}

/** Do some asset specific validation */
void TranchePricerForCalibrator::Validate() {
    static const string method("TranchePricerForCalibrator::Validate");
    /*=========================================================================
     * TODO: something to be added later ????? 
     *=======================================================================*/   
}

/** Get the asset , par curve and discount market data */
void TranchePricerForCalibrator::GetMarket(const IModel*        model, 
                          const CMarketDataSP  market)
{
    market->GetReferenceDate(valueDate);

    /*=========================================================================
     * GET THE DISCOUNT CURVE 
     *=======================================================================*/
    discount.getData(model, market);
    rfCrv = YieldCurveSP(discount.get());
    
	sCIDparam.getData(model, market); // get the sCID parameters used for the evaluation
	sCIDCalibParam.getData(model, market); // get the sCID calibration parameters 
}

DateTime TranchePricerForCalibrator::getValueDate() const {
    return valueDate;
}


void TranchePricerForCalibrator::priceClosedForm(CResults* results, Control* control, SCID* model) {
    static const string method = "TranchePricer::priceClosedForm";

    try {
        if (! SCID::TYPE->isInstance(model))
            throw ModelException(method, "sCID instrument did not receive an SCID model");
        const string& ccy = discount->getCcy();


		// obtain numerical datas from the model
		lossCalculationInterval = model->getFreqFastMC();
		seed = model->getSeed();
		timeSteps = model->getCFtimeSteps();
		nbPathsNoJump = model->getNbPathsFastNoJump();
		nbPathsAtLeastOneJump = model->getNbPathsFastAtLeastOneJump();
		convolutionNoJump = model->getConvolutionNoJump();
		convolutionAtLeastOneJump = model->getConvolutionAtLeastOneJump();
		int i, j, m;


		// Do first calibration, i.e. compute the weights without changing 
		DoubleArraySP imslWeights;
		vector<DoubleArrayArray> riskyAnnuity, defaultLeg;
		int calibratedNrWorlds = ComputeWeights(imslWeights,riskyAnnuity, defaultLeg);

		// improve single name calibration!
		DoubleArray pelWeights = *sCIDCalibParam->getPELCalibrationWeights();
		DateTimeArray pelDates = sCIDCalibParam->getPELCalibrationMaturities();
		DateTimeArray perfectELcalibrationDates;
		for (i=0; i<pelDates.size(); i++) if (pelWeights[i]<0) perfectELcalibrationDates.push_back(pelDates[i]);

		DateTimeArray singleNameDates;
		sCIDparam->getSingleNameCalibrationDates(singleNameDates);
		DoubleMatrix spMarket1, spMarket2;
		// get market survival proba at singleNameDates
		sCIDparam->MarketSurvProb(singleNameDates, spMarket1);
		// get market survival proba at perfectELcalibrationDates
		sCIDparam->MarketSurvProb(perfectELcalibrationDates, spMarket2);
		DoubleMatrix spNewMarket = spMarket1;
		if (perfectELcalibrationDates.size()>0)
		{
			int nbNames = sCIDparam->getNbNames();
			DoubleMatrix spCalib2;

			for (int loops=0; loops<nbLoops; loops++)
			{

				// get model survival proba at perfectELcalibrationDates
				sCIDparam->changeInitialWeights(*imslWeights);
				sCIDparam->NamesSurvProb(perfectELcalibrationDates, spCalib2);

				double tweak;
				int idxPerfect, idxCalib;
				for (m=0; m<nbNames; m++)
				{
					idxCalib = 0;
					for (idxPerfect = 0; idxPerfect<perfectELcalibrationDates.size(); ++idxPerfect)
					{
						tweak = log(spCalib2[idxPerfect][m])/log(spMarket2[idxPerfect][m]);
						while ( (idxCalib<singleNameDates.size()) && (singleNameDates[idxCalib]<=perfectELcalibrationDates[idxPerfect]) )
						{
							// we change the survival probabilities used to calibrate the base world
							spNewMarket[idxCalib][m] = pow(spMarket1[idxCalib][m], 1.0/tweak);
							++idxCalib;
						}
					}
					while (idxCalib<singleNameDates.size())
					{
						spNewMarket[idxCalib][m] = pow(spMarket1[idxCalib][m], 1.0/tweak);
						++idxCalib;
					}
				}

				sCIDparam->reInitializeSurvProba(spNewMarket);
				// recompute weights
				calibratedNrWorlds = ComputeWeights(imslWeights,riskyAnnuity, defaultLeg);
			}
		}


        DoubleArray multiShiftParameters(calibratedNrWorlds);
        DoubleArray parallelShiftParameters(calibratedNrWorlds);
        DoubleArray calibWeights(calibratedNrWorlds);

        DoubleArray initialMultiShiftPars, initialParallelShiftPars;
        sCIDparam->getMultShifts(initialMultiShiftPars);
        sCIDparam->getParallelShifts(initialParallelShiftPars);
        for (i=0, j=0;i<sCIDparam->getNbWorlds();i++)
            if ((*imslWeights)[i] == 0) continue;
            else {
                calibWeights[j] = (*imslWeights)[i];
                parallelShiftParameters[j] = initialParallelShiftPars[i];
                multiShiftParameters[j] = initialMultiShiftPars[i];
                j++;
            }
        // store the nber of vars
		results->storePrice(0, ccy);
        results->storeScalarGreek(
            calibratedNrWorlds, 
            "SCIDCalibrator", 
            OutputNameSP(new OutputName("nrWorlds")));
        
        results->storeGreek(
            DoubleArraySP(new DoubleArray(multiShiftParameters)),    
            "SCIDCalibrator", 
            OutputNameSP(new OutputName("multShiftParameters")));
                        
        results->storeGreek(
            DoubleArraySP(new DoubleArray(parallelShiftParameters)),    
            "SCIDCalibrator", 
            OutputNameSP(new OutputName("parallelShiftParameters")));
                        
        results->storeGreek(
            DoubleArraySP(new DoubleArray(calibWeights)),    
            "SCIDCalibrator", 
            OutputNameSP(new OutputName("weights")));

        results->storeGreek(
            CDoubleMatrixSP(new DoubleMatrix(spNewMarket)),    
            "SCIDCalibrator", 
            OutputNameSP(new OutputName("baseWorldCalibSurvProba")));


        //introduce calibrated values
        results->storeGreek(
            getCalibratedPEL(imslWeights),    
            "SCIDCalibrator", 
            OutputNameSP(new OutputName("calibratedPEL")));

        results->storeGreek(
            sCIDparam->getMarketPortfolioEL(sCIDCalibParam->getPELCalibrationMaturities()),    
            "SCIDCalibrator", 
            OutputNameSP(new OutputName("marketPEL")));
        
        //get aggregated risky annuity and default Legs
        DoubleArrayArraySP resultRA(new DoubleArrayArray()), resultDL(new DoubleArrayArray());
	    getCalibratedLegs(imslWeights, resultRA, resultDL, riskyAnnuity, defaultLeg);   
            
        results->storeGreek(
            resultRA, 
            "SCIDCalibrator", 
            OutputNameSP(new OutputName("calibratedRA")));
        
        results->storeGreek(
            resultDL, 
            "SCIDCalibrator", 
            OutputNameSP(new OutputName("calibratedDL")));


// OUTPUT SURVIVAL PROBABILITIES

		DoubleMatrix survProba;
		sCIDparam->changeInitialWeights(*imslWeights);
		sCIDparam->NamesSurvProb(singleNameDates, survProba);
		survProba.scale(-1.0);
		survProba.scalarAdd(1.0);

        results->storeGreek(
            CDoubleMatrixSP(new DoubleMatrix(survProba)), 
            "SCIDCalibrator", 
            OutputNameSP(new OutputName("survProba")));

        
    } catch (exception& e) {
        throw ModelException(&e, method);
    }
}

int TranchePricerForCalibrator::ComputeMarketConstraints(DoubleMatrix & EqIneq, 
                                                        DoubleArray &EqIneqVal, 
                                                        DoubleMatrix &Min,
                                                        vector<DoubleArrayArray> & riskyAnnuity,
                                                        vector<DoubleArrayArray> & defaultLeg)
{
    SCIDparametersSP sCIDparamSP = sCIDparam.getSP();
    DoubleMatrix ETL;
	DoubleArray kmin, kmax;
    //get attachement points
	sCIDCalibParam->getTrancheAttachementPoints(kmin,kmax);
	int i,j,k,l;
    //get all maturities from tranche quotes and build the aggregate datetime array
    DateTimeArraySP trancheMaturities;

	DateTime lastMaturity = valueDate;
    int size = kmin.size();
    //TODO: change it if CDO Tranche Quotes will have a start date too
    for ( i = 0; i < size; ++i) {
       trancheMaturities =  sCIDCalibParam->getTrancheMaturities(kmin[i], kmax[i]);
       if (lastMaturity < trancheMaturities->back())
           lastMaturity = trancheMaturities->back();
    }

    //get other tranche specification
	DayCountConventionSP dcc = sCIDCalibParam->getTrancheDCC();
	int TrancheCoupon = sCIDCalibParam->getTrancheCoupon();

    int nrWorlds = sCIDparamSP->getNbWorlds(), nrTELMaturities , nrTranches = kmin.size();
    // Compute Dates at which we are going to compute the TEL
    DateTimeArray tDate;
	sCIDparamSP->push_backTimeLine(tDate, valueDate, lastMaturity, lossCalculationInterval, true);
    array<double> t = sCIDparamSP->DateAsDouble(tDate);

    ETL.resize(t.size(),nrTranches);
    ETL.fill(0.0);
    vector< DoubleMatrix > TELallWorlds(nrWorlds, ETL);

    // Compute Tranche Expected Losses at all times t
    sCIDparamSP->setFastMC(seed, timeSteps, int(t.back()/timeSteps), nbPathsNoJump, nbPathsAtLeastOneJump);
    sCIDparamSP->setConvolution(kmin, kmax);
    sCIDparamSP->ComputeTELinAllWorlds(t, TELallWorlds, convolutionNoJump, convolutionAtLeastOneJump);

    //Prepare the input variables
	int nbEquality = 1;
	int nbMinimization = 0;
    DoubleArrayConstSP indexWeights = sCIDCalibParam->getPELCalibrationWeights();
    size = indexWeights->size();
	for (i=0; i<size; i++)
		if ((*indexWeights)[i]>1e-12) nbMinimization++;
			else if ((*indexWeights)[i]<-1e-12) nbEquality++;
    
    //set storage structures
    riskyAnnuity.resize(nrWorlds);
    defaultLeg.resize(nrWorlds);
    for (i=0; i < nrWorlds; i++)
	{
        riskyAnnuity[i].resize(nrTranches);
        defaultLeg[i].resize(nrTranches);
	}
   
    DoubleArraySP tranchePrices, trancheWeights, trancheUpfrontPremiums;
    DoubleArray::const_iterator weightIterator;
    for ( i = 0; i < nrTranches; ++i) {
        trancheWeights = sCIDCalibParam->getTrancheWeights(kmin[i], kmax[i]);       
        for(weightIterator = trancheWeights->begin(); weightIterator != trancheWeights->end(); ++weightIterator)
			if ((*weightIterator)>1e-12) nbMinimization++;
				else if ((*weightIterator)<-1e-12) nbEquality++;
        //set-up storage structures for tranche legs
        for (j=0; j < nrWorlds; j++) {
            riskyAnnuity[j][i].resize(trancheWeights->size());
            defaultLeg[j][i].resize(trancheWeights->size());
        }
    }
    
	EqIneq.resize(nrWorlds, nrWorlds+nbEquality);
	EqIneq.fill(0);
	EqIneqVal.resize(nbEquality+nrWorlds);
	Min.resize(nrWorlds, nbMinimization);
	Min.fill(0);
    

    DateTimeArray pelMaturities = sCIDCalibParam->getPELCalibrationMaturities(), trancheStartDates;
    pel.resize(pelMaturities.size(), nrWorlds);
    
	double DL, RA;
	double DLT, DLt, RAT, RAt;
	int minCount = 0, eqCount = 1;
    
    for(l=0; l < nrWorlds; l++) 
	{
        array<double> riskyDiscount(tDate.size());
        minCount = 0; eqCount = 1;        
        for (i=0; i<nrTranches; i++)
        {
            double div = 1.0 / ( kmax[i] - kmin[i] );
            riskyDiscount[0] = 1.0;
            for (j=0; j<tDate.size(); j++) riskyDiscount[j] = Maths::max(1e-12,1 - TELallWorlds[l][j][i]*div);
            
            EffectiveCurve trancheCurve(valueDate, discount.getSP(), tDate, riskyDiscount, EffectiveCurve::FLAT_FORWARD);
            
            trancheMaturities =  sCIDCalibParam->getTrancheMaturities(kmin[i], kmax[i]);
            nrTELMaturities = trancheMaturities->size();
            trancheWeights = sCIDCalibParam->getTrancheWeights(kmin[i], kmax[i]);       
            tranchePrices = sCIDCalibParam->getTranchePrices(kmin[i], kmax[i]);       
            trancheUpfrontPremiums= sCIDCalibParam->getTrancheUpFrontPremiums(kmin[i], kmax[i]);       
            trancheStartDates =  sCIDCalibParam->getTrancheStartDates(kmin[i], kmax[i]);
            for (j=0; j<nrTELMaturities; j++)
            {               
                if(valueDate.equals(trancheStartDates[j])) {
                    CashFlowArray cashflows = SwapTool::cashflows(valueDate, (*trancheMaturities)[j], false, 1.0, TrancheCoupon, "M", &(*dcc));
                    cashflows[cashflows.size()-1].amount -= 1.0; // ugly
                    DL = trancheCurve.protectionPV(valueDate,valueDate,(*trancheMaturities)[j],IDiscountCurveRisky::RECOVER_1);
                    RA = trancheCurve.annuityPV(cashflows,valueDate,IDiscountCurveRisky::RECOVER_1);
                }
                else {
                    DLT = trancheCurve.protectionPV(valueDate,valueDate,(*trancheMaturities)[j],IDiscountCurveRisky::RECOVER_1);
                    DLt = trancheCurve.protectionPV(valueDate,valueDate,trancheStartDates[j],IDiscountCurveRisky::RECOVER_1);
                    DL = DLT - DLt; 
                    CashFlowArray cashflows = SwapTool::cashflows(valueDate, (*trancheMaturities)[j], false, 1.0, TrancheCoupon, "M", &(*dcc));
                    cashflows[cashflows.size()-1].amount -= 1.0; // ugly
                    RAT = trancheCurve.annuityPV(cashflows,valueDate,IDiscountCurveRisky::RECOVER_1);
                    cashflows = SwapTool::cashflows(valueDate, trancheStartDates[j], false, 1.0, TrancheCoupon, "M", &(*dcc));
                    cashflows[cashflows.size()-1].amount -= 1.0; // ugly
                    RAt = trancheCurve.annuityPV(cashflows,valueDate,IDiscountCurveRisky::RECOVER_1);
                    RA = RAT - RAt; 
                }

                riskyAnnuity[l][i][j] = RA;
                defaultLeg[l][i][j] = DL;
                if ((*trancheWeights)[j]>1e-12)
                {
					if ((*tranchePrices)[j]>0)
					{
	                    if ((*trancheUpfrontPremiums)[j]==0.0)
		                    Min[l][minCount] = (*trancheWeights)[j]*(-RA + DL / (*tranchePrices)[j]);
			            else
				            Min[l][minCount] = (*trancheWeights)[j]*(-RA*(*tranchePrices)[j])+ DL - (*trancheUpfrontPremiums)[j] ;
					}
				    minCount++;
                }
                else if ((*trancheWeights)[j]<-1e-12)  
                {
					if ((*tranchePrices)[j]>0)
					{
	                    if ((*trancheUpfrontPremiums)[j]==0.0)
		                    EqIneq[l][eqCount] = (*trancheWeights)[j]*(-RA + DL / (*tranchePrices)[j]);
			            else
				            EqIneq[l][eqCount] = (*trancheWeights)[j]*(-RA*(*tranchePrices)[j])+ DL - (*trancheUpfrontPremiums)[j] ;
					    EqIneqVal[eqCount]=0;
					}
                    eqCount++;
                }
            }
        }
        DoubleArray pelWorld(pelMaturities.size());
        sCIDparamSP->PortfolioELinGivenWorld(l, pelMaturities, &pelWorld[0]);
        size = indexWeights->size();
        for (k=0; k<size; k++)
        {
            if ((*indexWeights)[k]>1e-12)
            {
                Min[l][minCount] = (*indexWeights)[k]*(pelWorld[k]/sCIDparam->getMarketPortfolioEL(pelMaturities[k]) -1);
                minCount++;
            }
            else if ((*indexWeights)[k]<-1e-12)
            {
                EqIneq[l][eqCount] = sCIDparam->getMarketPortfolioEL(pelMaturities[k]) - pelWorld[k];
                EqIneqVal[eqCount]=0;
                eqCount++;
            }
            //store pel
            pel[k][l] = pelWorld[k];
        }
    }


	// finally constraints to make sure the weights are positive and add up to one	
	for (i=0; i<nrWorlds; i++) EqIneq[i][0]=1; 
	EqIneqVal[0]=1;											// the weights are probability...
	for (i=0; i<nrWorlds; i++) 
	{
		EqIneq[i][nbEquality+i] = 1;  
		EqIneqVal[nbEquality+i] = 1e-25;
	}
	return nbEquality;
}

int TranchePricerForCalibrator::ComputeWeights(DoubleArraySP &imslWeights, vector<DoubleArrayArray> &riskyAnnuity, vector<DoubleArrayArray> &defaultLeg)
{
    int i,j,k, meq, m, n = sCIDparam->getNbWorlds();

    DoubleMatrix sqrtHessian, a;
    DoubleArray b;
    DoubleArray g(n, 0.0);   // same notations used as in IMSL solver...
    DoubleArray h(n*n,0.0);

    meq = ComputeMarketConstraints(a, b, sqrtHessian, riskyAnnuity, defaultLeg);
    m = meq + n;

//	h = sqrtHessian * sqrtHessian^T +1e-10 Id;
    for (i=0; i<n; i++)
        for (j=i; j<n; j++)
            for (k=0; k<sqrtHessian.numRows(); k++)
                h[i*n+j] += sqrtHessian[i][k]*sqrtHessian[j][k];
    for (i=0; i<n; i++)
        for (j=0; j<i; j++)
            h[i*n+j] = h[j*n+i];
    for (i=0; i<n; i++)
        h[i*(n+1)] += 1e-10;

    DoubleArray vectorEqIneq(a.numRows()*a.numCols());
    for(i=0; i < a.numCols();i++)   
        for(j=0; j < a.numRows();j++)   
            vectorEqIneq[j*a.numCols()+i] = a[i][j];        

    //call the imsl quadratic programming method
    imslWeights = QuadraticProg::minimize(m, n, meq, vectorEqIneq, b ,g ,h);

    int calibratedNrWorlds = 0;
    for (i=0; i<n; i++){
        if ((*imslWeights)[i]<1e-12) (*imslWeights)[i] = 0;
        else calibratedNrWorlds ++;
    }

    double sum = accumulate(imslWeights->begin(), imslWeights->end(), 0.0);
    for (i=0; i<n; i++) (*imslWeights)[i] = (*imslWeights)[i]/sum;
	return calibratedNrWorlds;
}

void TranchePricerForCalibrator::getCalibratedLegs(DoubleArraySP& imslWeights, 
                                                    DoubleArrayArraySP& resultRA, 
                                                    DoubleArrayArraySP& resultDL, 
                                                    vector<DoubleArrayArray> &riskyAnnuity, 
                                                    vector<DoubleArrayArray> &defaultLeg)
{
    int nrWorlds = riskyAnnuity.size(), nrTranches = riskyAnnuity[0].size(), nrMat;
    double RA, DL;
    resultRA->resize(nrTranches);
    resultDL->resize(nrTranches);
    for (int i = 0; i< nrTranches; i++) {
        nrMat = riskyAnnuity[0][i].size();
        for (int j=0; j < nrMat; j++){
            RA = DL = 0.0;
            for (int k=0; k < nrWorlds; k++){
                RA+=riskyAnnuity[k][i][j]*(*imslWeights)[k];
                DL+=defaultLeg[k][i][j]*(*imslWeights)[k];
            }
            (*resultRA)[i].push_back(RA);
            (*resultDL)[i].push_back(DL);
        }
    }
}

DoubleArraySP TranchePricerForCalibrator::getCalibratedPEL(DoubleArraySP &initialWeights)
{
    DateTimeArray pelMaturities = sCIDCalibParam->getPELCalibrationMaturities();
    int nrMat = pelMaturities.size(), nrWorlds = pel.numRows();
    DoubleArraySP result(new DoubleArray(nrMat,0.0));
    for (int i=0; i< nrMat; i++){
        for (int j=0; j < nrWorlds; j++)
            (*result)[i]+=pel[i][j]*(*initialWeights)[j];
    }
    return result; 
}



/*=============================================================================
 * Reflection, loading, etc.
 *===========================================================================*/
class TranchePricerForCalibratorHelper {
public:
    static IObject* defaultTranchePricerForCalibrator();
    /** Invoked when Class is 'loaded' */
    static void load(CClassSP& clazz){
        clazz->setPublic(); // make visible to EAS/spreadsheet
        REGISTER(TranchePricerForCalibrator, clazz);
        SUPERCLASS(CInstrument);
        IMPLEMENTS(SCID::IIntoProduct);
        EMPTY_SHELL_METHOD(defaultTranchePricerForCalibrator);

        FIELD(valueDate,                  "Valuation Date - default from market.");
        FIELD(discount,                   "Discount curve");        
 		FIELD(sCIDparam,                  "sCID parameters");
 		FIELD(sCIDCalibParam,             "sCID calibration parameters");
 		FIELD(nbLoops,                    "nbLoops done in the calibration to improve single name Calibration");

        FIELD(pel,                               "the vector with portfolio ELs per world for each index maturity"); 
        FIELD_MAKE_TRANSIENT(pel);
		FIELD(lossCalculationInterval,    "how often do we compute the loss distribution");
        FIELD_MAKE_TRANSIENT(lossCalculationInterval);
		FIELD(seed,						 "seeds: one for Common Factor, and two for Poisson processes");
        FIELD_MAKE_TRANSIENT(seed);
		FIELD(timeSteps,                  "time Steps used in the discretization of the CIR process");
        FIELD_MAKE_TRANSIENT(timeSteps);
 		FIELD(nbPathsNoJump,              "nbPaths conditional on no jumps");
        FIELD_MAKE_TRANSIENT(nbPathsNoJump);
 		FIELD(nbPathsAtLeastOneJump,      "nbPaths conditional on at least one jump");
        FIELD_MAKE_TRANSIENT(nbPathsAtLeastOneJump);
 		FIELD(convolutionNoJump,          "Convolution Method conditional on no jumps ");
        FIELD_MAKE_TRANSIENT(convolutionNoJump);
 		FIELD(convolutionAtLeastOneJump,  "Convolution Method conditional on at least one jump");
        FIELD_MAKE_TRANSIENT(convolutionAtLeastOneJump);
        
	}
};

IObject* TranchePricerForCalibratorHelper::defaultTranchePricerForCalibrator() {
    return new TranchePricerForCalibrator();
}

CClassConstSP const TranchePricerForCalibrator::TYPE = 
    CClass::registerClassLoadMethod("TranchePricerForCalibrator", typeid(TranchePricerForCalibrator),TranchePricerForCalibratorHelper::load);

/*=============================================================================
 * Pricing Models
 *===========================================================================*/
class TranchePricerForCalibratorSCID : public SCID::IProduct {
private:
    const TranchePricerForCalibrator* instr; // a reference
    SCID* model;

public:
    TranchePricerForCalibratorSCID(const TranchePricerForCalibrator* instr, SCID* model): instr(instr), model(model){}
    void price(SCID* model,
               Control*         control, 
               CResults*        results) const {
        const_cast<TranchePricerForCalibrator*>(instr)->priceClosedForm(results, control, model);
    }
};

/** Implementation of ClosedFormBSImpliedSmile::IntoProduct interface */
SCID::IProduct* TranchePricerForCalibrator::createProduct(SCID* model) const {
    return new TranchePricerForCalibratorSCID(this, model);
}

/** Included in ProductsLib to force the linker to include this file */
bool TranchePricerForCalibratorLoad() {
    return (TranchePricerForCalibrator::TYPE != 0);
}

DRLIB_END_NAMESPACE

