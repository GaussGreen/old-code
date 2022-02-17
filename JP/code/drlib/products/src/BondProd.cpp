//----------------------------------------------------------------------------
//
//   Group       : Convertibles DR
//
//   Filename    : BondProd.cpp
//
//   Description : Bond Product Code
//
//   Author      : André Segger
//
//   Date        : 08 December 2003
//
//
//----------------------------------------------------------------------------


#include "edginc/config.hpp"
#include "edginc/CorporateBond.hpp"
#include "edginc/LinearStrikeVolRequest.hpp"
#include "edginc/FD1FLNGeneric.hpp"
#include "edginc/FD1FLNResettable.hpp"

DRLIB_BEGIN_NAMESPACE

/////////////////////////////////////////////////////////
//           fd product
/////////////////////////////////////////////////////////


/** private class */
class BondFDProd: virtual public FD1F::IProduct,
                  virtual public FD1FGeneric::IProduct {
    
public:

    BondFDProd(const CorporateBond* instr): bond(instr) {
        coupon              = 0;
        addCouponBool       = 0;
        callLevel           = 0;
        callBool            = 0;
        callIsAdjForAccrued = 0;
        putLevel            = 0;
        putBool             = 0;
        defaultPayment      = 0;

        control             = 0;
        isE2C               = false;

        // make the risky curve
        // make sure that we have a BootstrappedYieldCurve otherwise we can't add the spreadCurve
        if (!BootstrappedYieldCurve::TYPE->isInstance(bond->discount.get()))
        {
    		throw ModelException("BondFDProd::BondFDProd", 
                                 "Yield curve must be a cash swap curve");
	    }

        DateTime tmp(bond->bond->getMaturityDate());
        DateTimeSP matDate(new DateTime(tmp));

        asset = CAssetWrapper(copy(bond->asset.get()));
    }

    virtual ~BondFDProd();

    virtual CAssetConstSP GetAssetRef() {
        return CAssetConstSP::dynamicCast((IObjectConstSP)asset.getSP());
    }

    virtual bool GetFwdStartLV()
    {
           return false;
    }

    virtual DateTime GetFwdStartDateLV()
    {
           return bond->getValueDate();
    }
    
    virtual YieldCurveConstSP GetDiscCurveRef();

    virtual CVolRequestConstSP GetLNRequest();

    void InitGridSize();

    void InitFD(Control* control);
    
    void InitProd();

    virtual bool Positive()
    {
        // option payoff
        return false;
    }

    virtual string getCcyTreatment(){ return bond->ccyTreatment;}

    /** make price refinement - control variate */
    double RefinePrice(double basePrice, double discFactor, bool)
    {
       // scale premium for forward starting.
       double price = discFactor*scalePremium(basePrice);
       return price;
    }

    virtual double scalePremium(const double& fairValue);

    /** extra output requests */
    virtual void recordOutputRequests(Control* control, Results* results, double fairValue);


    virtual void preCalcFD(int step, int idx, int pStart, int pEnd);
    
    void PayoffAtMat(const double * s, int step, int bot, int top, int pStart, int pEnd,
                     double * const * price);

    void PayoffBeforeMat(const double * s, int step, int bot, int top, int pStart, int pEnd,
                         double * const * price);

    virtual double getCoupon(int step, const double* s, int start, int end);

    virtual bool hasEquityLayer();

    virtual void postProcessFD(const double *price, const int numLayers);

    bool                isE2C;

protected:

private:

    const CorporateBond* bond;

    double *coupon;
    bool   *addCouponBool;
    double *callLevel;
    bool   *callBool;
    bool   *callIsAdjForAccrued;
    double *putLevel;
    bool   *putBool;
    double *defaultPayment;

    CAssetWrapper       asset;

    // members required for put probability calculation
    const Control* control;

    // processed vol is required for Black on call
    CVolProcessedBSSP processedVol;

    // some members used to speed up performance;
    double          faceValue;
    double          recovery;
    double          defaultBarrier;
};

BondFDProd::~BondFDProd() {
    if (coupon != 0) {
        delete [] coupon;
        coupon = 0;}
    if (addCouponBool != 0) {
        delete [] addCouponBool;
        addCouponBool = 0;}
    if (callLevel != 0) {
        delete [] callLevel;
        callLevel = 0;}
    if (callBool != 0) {
        delete [] callBool;
        callBool = 0;}
    if (callIsAdjForAccrued != 0) {
        delete [] callIsAdjForAccrued;
        callIsAdjForAccrued = 0;}
    if (putLevel != 0) {
        delete [] putLevel;
        putLevel = 0;}
    if (putBool != 0) {
        delete [] putBool;
        putBool = 0;}
    if (defaultPayment != 0) {
        delete [] defaultPayment;
        defaultPayment = 0;}
}


YieldCurveConstSP BondFDProd::GetDiscCurveRef() {
    return YieldCurveConstSP::dynamicCast((IObjectConstSP)bond->discount.getSP());
}

void BondFDProd::InitGridSize()
{
    if ( genericFDModel ) {
        if (genericFDModel->stockStepsToUse == -1) { // grid hasn't been set yet
            FD1FLNGeneric *lnModel;

            lnModel = dynamic_cast<FD1FLNGeneric*>(genericFDModel);

            if (genericFDModel->stepsPerYear == -1 ) {
                if (lnModel != 0) {
                    int BaseStepAmount = 300;
                    double var;
                    int totalSteps;
                
                    processedVol = CVolProcessedBSSP(dynamic_cast<CVolProcessedBS*>(asset->getProcessedVol(GetLNRequest().get())));
                    var = processedVol->CalcVar(bond->getValueDate(), bond->bond->getMaturityDate());

                    totalSteps = int(BaseStepAmount*Maths::max(1., var));

                    genericFDModel->stepsPerYearToUse = int(ceil((double)totalSteps/bond->getValueDate().yearFrac(bond->bond->getMaturityDate())));
            
                    // always use at least 50 a year
                    genericFDModel->stepsPerYearToUse = Maths::max(genericFDModel->stepsPerYearToUse, 50);
                } else {
                    genericFDModel->stepsPerYearToUse = 100;
                }
            } else {
                genericFDModel->stepsPerYearToUse = genericFDModel->stepsPerYear;
            }
        
            if (genericFDModel->stockSteps == -1) {
                if (lnModel != 0) {
                    int totalSteps = int(genericFDModel->stepsPerYearToUse * bond->getValueDate().yearFrac(bond->bond->getMaturityDate()));
                    genericFDModel->stockStepsToUse = int(2*model1F->TruncationStd *sqrt(2.*totalSteps));
                    // make it even
                    genericFDModel->stockStepsToUse = int((2*genericFDModel->stockStepsToUse + 1)/2);
                } else {
                    genericFDModel->stockStepsToUse = 200;
                }
            } else {
                genericFDModel->stockStepsToUse = genericFDModel->stockSteps;
            }
        }
    } else {
        if (fdModel->stockStepsToUse == -1) { // grid hasn't been set yet
            FD1FLN* lnModel;

            lnModel = dynamic_cast<FD1FLN*>(fdModel);

            if (fdModel->stepsPerYear == -1 ) {
                if (lnModel != 0) {
                    int BaseStepAmount = 300;
                    double var;
                    int totalSteps;
                
                    processedVol = CVolProcessedBSSP(dynamic_cast<CVolProcessedBS*>(asset->getProcessedVol(GetLNRequest().get())));
                    var = processedVol->CalcVar(bond->getValueDate(), bond->bond->getMaturityDate());

                    totalSteps = int(BaseStepAmount*Maths::max(1., var));

                    fdModel->stepsPerYearToUse = int(ceil((double)totalSteps/bond->getValueDate().yearFrac(bond->bond->getMaturityDate())));
            
                    // always use at least 50 a year
                    fdModel->stepsPerYearToUse = Maths::max(fdModel->stepsPerYearToUse, 50);
                } else {
                    fdModel->stepsPerYearToUse = 100;
                }
            } else {
                fdModel->stepsPerYearToUse = fdModel->stepsPerYear;
            }
        
            if (fdModel->stockSteps == -1) {
                if (lnModel != 0) {
                    int totalSteps = int(fdModel->stepsPerYearToUse * bond->getValueDate().yearFrac(bond->bond->getMaturityDate()));
                    fdModel->stockStepsToUse = int(2*model1F->TruncationStd *sqrt(2.*totalSteps));
                    // make it even
                    fdModel->stockStepsToUse = int((2*fdModel->stockStepsToUse + 1)/2);
                } else {
                    fdModel->stockStepsToUse = 200;
                }
            } else {
                fdModel->stockStepsToUse = fdModel->stockSteps;
            }
        }
    }
}

void BondFDProd::InitFD(Control* control)
{
    static const string method = "BondFDProd::InitFD";
    try {

        // required to calculate put probability as part of the PayoffBeforeMat callback
        this->control = control;

        // critical dates
        DateTimeArray critDates;

        int i;
        CashFlowArraySP myCFs = bond->bond->getExAdjCashFlows(bond->getValueDate(), bond->discount.getSP()); 
        for (i=0; i<myCFs->size(); i++) {
            critDates.push_back((*myCFs)[i].date);
            critDates.push_back((*myCFs)[i].date.rollDate(-1));
        }

        // add all dividend pass through dates to the critical dates
        if ( bond->putSchedule.get() ) {
            const DateTimeArray& putDates = bond->putSchedule->getDateArray();
            for (i=0; i<putDates.size(); i++)
                critDates.push_back(putDates[i]);  
        }

        double minGap = -1; // default value
        
        DateTimeArray hardCallDates;
        if (bond->callSchedule.get()) {
            hardCallDates = bond->callSchedule->getDates();
        }

        // All call dates are critical dates
        for (i=0; i<hardCallDates.size(); i++)
            critDates.push_back(hardCallDates[i]);
   
        // Hard calls are segment dates but soft calls are standard critical dates

        TimeMetricConstSP timeMetric;
        if (fdModel) {
            timeMetric = fdModel->GetTimeMetric();
        } else {
            timeMetric = genericFDModel->GetTimeMetric();
        }

        // only the next call date is a segment date
        DateTimeArray segDates;
        segDates.resize(1);
        segDates[0] = bond->getValueDate();
        for (i=0; i<hardCallDates.size(); i++) {
            if (hardCallDates[i] > bond->getValueDate() && hardCallDates[i] < bond->bond->getMaturityDate()) {
                if ( i == 0 || Maths::isPositive(timeMetric->yearFrac(segDates[segDates.size()-1], hardCallDates[i]))) {
                    segDates.push_back(hardCallDates[i]);
                    break;
                }
            }
        }

        segDates.push_back(bond->bond->getMaturityDate());

        // sort critical dates
        DateTime minDate = bond->getValueDate().rollDate(-1);
        DateTime maxDate = bond->endDate(0).rollDate(+1);
        timeMetric->SortDate(minDate, maxDate, true, segDates);

        vector<int> density;
        density.resize(segDates.size()-1);
        for (i=0; i<segDates.size()-1; i++) {
            density[i] = 1;
        }
        bool useEqualTime = true;

        int numOfPriceArray;

        isE2C           = false;
        numOfPriceArray =  1;

        if (fdModel) {
            fdModel->TimePts.SetTreeStepBase(500);
            // the try/catch block is here to avoid a mysterious crash in gcc on Solaris optimised
            try { 
                fdModel->Setup(bond->getValueDate(), segDates, density, &critDates,
                               minGap, useEqualTime, numOfPriceArray);
            } catch (exception&) {
                throw;
            }
        } else {
            genericFDModel->TimePts.SetTreeStepBase(500);
            // the try/catch block is here to avoid a mysterious crash in gcc on Solaris optimised
            try { 
                genericFDModel->Setup(bond->getValueDate(), segDates, density, &critDates,
                               minGap, useEqualTime, numOfPriceArray);
            } catch (exception&) {
                throw;
            }

        }
        

    } catch (exception& e) {
        throw ModelException(&e, method);
    }
}

void BondFDProd::InitProd()
{
  	static const string method = "BondFDProd::InitProd";
    int i;
    CashFlowArraySP myCFs = bond->bond->getExAdjCashFlows(bond->getValueDate(), bond->discount.getSP()); 
    int numStep           = model1F->TimePts.NumOfStep;
    coupon                = new double[numStep+1];
    addCouponBool         = new bool[numStep+1];
    int cfIdx             = 0;

    // for performance

    FirmAsset* firmAsset;
    if ( FirmAsset::TYPE->isInstance(asset.get())) {
        CAsset*    ncAsset   = const_cast<CAsset*>(asset.get());
        firmAsset = dynamic_cast<FirmAsset*>(ncAsset);
    }


    faceValue      = bond->bond->getFaceValue();
    recovery       = bond->getBondRecovery();
    defaultBarrier = firmAsset->getDefaultBarrier();

    // assume the cashFlows are in order
    for (i=0; i<=numStep; i++) {
        if ( cfIdx <  myCFs->size() && 
            model1F->TimePts.StepDates[i].getDate() == (*myCFs)[cfIdx].date.getDate() ) {
            if ( i == 0 || 
                 model1F->TimePts.StepDates[i].getDate() != 
                 model1F->TimePts.StepDates[i-1].getDate() ) {
                addCouponBool[i] = true;
                // I might have a problem here with coupons at maturity
                if (cfIdx == myCFs->size()-1) {
                    addCouponBool[i] = true;
                    coupon[i] = (*myCFs)[cfIdx].amount - bond->bond->getRedemption();
                } else {
                    coupon[i] = (*myCFs)[cfIdx].amount;
                }
                cfIdx++;
            }
        } else {
            addCouponBool[i] = false;
            coupon[i] = 0.;
        }
    }
    if (cfIdx != myCFs->size()) {
        throw ModelException("InitTree","Can't find all coupon dates in timeLine");
    }

    callLevel           = new double[numStep+1];
    callBool            = new bool[numStep+1];
    callIsAdjForAccrued = new bool[numStep+1];
    
    double dummy = asset->getSpot();
    vector<double> fwdVol;
    for (i=0; i<=numStep; i++) {
        model1F->GetStepVol(i, fwdVol, &dummy,0,0);
        callBool[i] = bond->getCallLevel(model1F->TimePts.StepDates[i], callLevel[i]);
    } 

    putLevel   = new double[numStep+1];
    putBool    = new bool[numStep+1];
    
    for (i=0; i<=numStep; i++) {
        putBool[i] = bond->getPutLevel(model1F->TimePts.StepDates[i], putLevel[i]);
    }

    defaultPayment      = new double[numStep+1];
    for (i=0; i<=numStep; i++) {
        double claim = bond->bond->getNotional(model1F->TimePts.StepDates[i]);
        if (bond->payAccruedUponDefault) {
            claim += bond->bond->getAccruedAtDate(model1F->TimePts.StepDates[i]);
        }
        defaultPayment[i] = claim * recovery;
    } 

    if (model1F->isTree() == false) {
        // set the max stock level for each sement
        int        j;
        double     maxStock;
        bool       newStockMax;
        FirmAsset* firmAsset;

        if ( FirmAsset::TYPE->isInstance(asset.get())) {
            CAsset*    ncAsset   = const_cast<CAsset*>(asset.get());
            firmAsset = dynamic_cast<FirmAsset*>(ncAsset);
        }

        for (j=0; j<int(model1F->TimePts.SegmentEnd.size()); j++) {
            maxStock = 0.;
            newStockMax = true;
            for (i=model1F->TimePts.SegmentStart[j]; i<=model1F->TimePts.SegmentEnd[j]; i++) {
                newStockMax = false;
                break; // keep the default truncation level
            }

            genericFDModel->stockMinSeg[j] = defaultBarrier;
        }
    }
}

/** calculate barriers if needed */
void BondFDProd::preCalcFD(int step, int idx, int pStart, int pEnd)
{
    int i;
    int layer;

    double conversionRatio;
    double conversionRatioNext;

    if (fdModel) {
        for (layer=0 ; layer<= pEnd ; ++layer) {
            i = 0;
            conversionRatio     = 0.0;
            conversionRatioNext = 0.0;

            if (step < fdModel->TimePts.NumOfStep) {
                fdModel->fdEngine.upBarrierNew[i]       = -1;
                fdModel->fdEngine.upPayoutNew[i]        = 0.;
                fdModel->fdEngine.upPayoutDeltaNew[i]   = 0.;
    
                fdModel->fdEngine.upBarrierOld[i]       = -1;
                fdModel->fdEngine.upPayoutOld[i]        = 0.;
                fdModel->fdEngine.upPayoutDeltaOld[i]   = 0.;
    
                fdModel->fdEngine.upBarrierOld[i]     = -1;
                fdModel->fdEngine.upPayoutOld[i]      = 0.;
                fdModel->fdEngine.upPayoutDeltaOld[i] = 0.;
            }

            ++i;
        }

    } else {
        double accruedToday = 0.0;
        if (bond->payAccruedUponDefault) {
            accruedToday = bond->bond->getAccruedAtDate(genericFDModel->TimePts.StepDates[step]);
        }

        for (layer=0 ; layer<= pEnd ; ++layer) {
            i = 0;

            conversionRatio     = 0.0;
            conversionRatioNext = 0.0;

            if (step < genericFDModel->TimePts.NumOfStep) {
                genericFDModel->fdEngine.upBarrierNew[i]       = -1;
                genericFDModel->fdEngine.upPayoutNew[i]        = 0.;
                genericFDModel->fdEngine.upPayoutDeltaNew[i]   = 0.;
    
                genericFDModel->fdEngine.upBarrierOld[i]       = -1;
                genericFDModel->fdEngine.upPayoutOld[i]        = 0.;
                genericFDModel->fdEngine.upPayoutDeltaOld[i]   = 0.;

                genericFDModel->fdEngine.downBarrierNew[i]       = defaultBarrier;
                genericFDModel->fdEngine.downPayoutNew[i]        = (faceValue + accruedToday) * recovery;

                genericFDModel->fdEngine.downPayoutDeltaNew[i]   = 0.0;
            }
            ++i;
        }
    }
}


/** product payoff method at steps earlier than maturity */
void BondFDProd::PayoffAtMat(const double * s, int step, int bot, int top, int pStart, int pEnd,
                                    double * const * price)
{
    int     i,j;

    // the price is simply the redemption value
    for (i=pStart; i<=pEnd; i++) {
        for (j=-bot; j<=top; ++j) {
            if (putBool[step] == true) {
                price[i][j] = Maths::max(bond->bond->getRedemption(), putLevel[step]);
            } else {
                price[i][j] = bond->bond->getRedemption();
            }
            // check coupon paid today ...
            if (addCouponBool[step] == true) {
                price[i][j] += coupon[step];
            }
        }
    }
}


/** product payoff method at steps earlier than maturity */
void BondFDProd::PayoffBeforeMat(const double * s, int step, int bot, int top, int pStart, int pEnd,
                                        double * const * price)
{
    static const string method = "BondFDProd::PayoffBeforeMat";
    int         j;
    int         layer;

    int endIdx;
    endIdx = pStart;

    for (layer=pStart; layer<=endIdx; ++layer) {
        // calculate relevant stock interval for Black on call
        for (j=-bot; j<=top; ++j) {
            // check coupon paid today ...
            if (addCouponBool[step] == true) {
                price[layer][j] += coupon[step];
            }

            // check callable ...
            if (callBool[step]) {
                price[layer][j] = Maths::min(price[layer][j], callLevel[step]);
            }
            
            // check puttable ....
            if (putBool[step]) {
                price[layer][j] = Maths::max(price[layer][j], putLevel[step]);
            }
        }  
    }
}

double BondFDProd::getCoupon(int step, const double* s, int start, int end)
{
    return defaultPayment[step];
}

bool BondFDProd::hasEquityLayer()
{
    return false;
}

void BondFDProd::postProcessFD(const double *price, const int numLayers) {

    if (numLayers != 1) {
        throw ModelException("BondFDProd::postProcessFD",
            "More than one gf layer is not currently supported for the finite difference E2C bond model");
    }
    if ( !bond->e2cBasePriceCalculated ) {
        bond->e2cBasePriceCalculated = true;
        bond->e2cBasePrice           = price[0];
    }
};


/** premium scaling */
double BondFDProd::scalePremium(const double& fairValue)
{
    return model1F->PriceEnd[0];
}

/** extra output requests */
void BondFDProd::recordOutputRequests(Control* control, Results* results, double fairValue)
{
    static const string method = "BondFDProd::recordOutputRequests";
    OutputRequest* request = NULL;

    if ( control && control->isPricing() ) {
        // need to put that back in 
        // bond->recordOutputRequests(control, results, fairValue, fairValue, false);

        DateTime maturity = bond->bond->getMaturityDate();
        HolidaySP hols(Holiday::weekendsOnly());
        
        DefaultRatesSP psDefRates = bond->cdsParSpreads->defaultRates();
        CashFlowArraySP cleanSpreadCurve = psDefRates->getCleanSpreadCurve();
        
        
        IObjectSP currentSpread = bond->cdsParSpreads->getCurrentSpreadOrUntweakable(
                                      bond->valueDate, maturity);
        bond->addRequests(control, results, cleanSpreadCurve, currentSpread);

        // Indicative vol
        if (control->requestsOutput(OutputRequest::IND_VOL, request)) {
            double             indVol     = 0.0;
            CVolRequestConstSP volRequest   = GetLNRequest();
            CVolRequestSP      ncVolRequest =CVolRequestSP::constCast(volRequest);

            LinearStrikeVolRequest* lsVolRequest = dynamic_cast<LinearStrikeVolRequest*>(ncVolRequest.get());
            // interpolate the vol using our LN request
            CVolProcessedBSSP volBS(asset->getProcessedVol(lsVolRequest));
            // calculate the indicative vol
            try {
                indVol = volBS->CalcVol(bond->valueDate, bond->bond->getMaturityDate());
            }
            catch (exception& ) {
                indVol = 0.0;
            }
            results->storeRequestResult(request, indVol);
        }

		// Duration
		if ( control->requestsOutput(OutputRequest::BOND_DURATION, request) ) {
			results->storeRequestResult(request, bond->bond->duration(bond->getValueDate(),bond->discount.getSP()));
		}
		// Convexity
		if ( control->requestsOutput(OutputRequest::BOND_CONVEXITY, request) ) {
			results->storeRequestResult(request, bond->bond->convexity(bond->getValueDate(),bond->discount.getSP()));
		}
    }
}

/** returns a vol request for log-normal vol */
CVolRequestConstSP BondFDProd::GetLNRequest() {
    CVolRequestConstSP volRequest(
        new LinearStrikeVolRequest(bond->asset->getSpot(),
                                   bond->getValueDate(),
                                   bond->bond->getMaturityDate(), 
                                   false));
  
    return volRequest;
}

/** Implementation of CFDGridPass::IntoProduct interface */
FD1FGeneric::IProduct* CorporateBond::createProduct(
    FD1FGeneric* model) const {
    BondFDProd *tmpBondFDProd = new BondFDProd(this);
    tmpBondFDProd->genericFDModel = model;
    tmpBondFDProd->model1F        = model;
    return tmpBondFDProd;
}

DRLIB_END_NAMESPACE
