 //----------------------------------------------------------------------------
//
//   Group       : Convertibles DR
//
//   Filename    : ConvBondProd.cpp
//
//   Description : Convertible Bond Product Code
//
//   Author      : Tycho von Rosenvinge
//
//   Date        : January 28, 2002
//
//----------------------------------------------------------------------------

#include "edginc/config.hpp"
#include "edginc/LinearStrikeVolRequest.hpp"
#include "edginc/BootstrappedYieldCurve.hpp"
#include "edginc/Black.hpp"
#include "edginc/Untweakable.hpp"
#include "edginc/mathlib.hpp"
#include "edginc/StruckEquity.hpp"
#include "edginc/Untweakable.hpp"
#include "edginc/OptOnConvBond.hpp"
#include "edginc/FD1FLNGeneric.hpp"
#include "edginc/DividendCollector.hpp"
#include "edginc/FD1FLNResettable.hpp"
#include "edginc/FD1FDDE.hpp"
#include "edginc/FD1FE2C.hpp"
#include "edginc/FirmAsset.hpp"
#include "edginc/Tree1fLN.hpp"
#include "edginc/VegaParallel.hpp"

#include "edginc/IndexSpecEQ.hpp"
#include "edginc/LatticeProdEDR.hpp"

DRLIB_BEGIN_NAMESPACE

/////////////////////////////////////////////////////////
//           tree1f and fd product
/////////////////////////////////////////////////////////

 double conversionOption(bool isCall, const DateTime& valueDate, const DateTime& matDate, double fwd, double strike, double variance, double discFact)
{
    // call Black to value option
    double callValue = Black::price(isCall,
                                    fwd,
                                    strike,
                                    discFact,
                                    variance);

    return callValue;
}

/** private class */
class ConvBond1fProd: virtual public FD1F::IProduct,
                      virtual public FD1FGeneric::IProduct {
    
public:

    ConvBond1fProd(const ConvBond* instr): cvb(instr) {
        coupon              = 0;
        addCouponBool       = 0;
        callLevel           = 0;
        callBool            = 0;
        softPutBool         = 0;
        softPutLevel        = 0;
        softPutTrigger      = 0;
        callIsHard          = 0;
        callIsAdjForAccrued = 0;
        putLevel            = 0;
        putBool             = 0;
        convRatio           = 0;
        convCash            = 0;
        parityFloor         = 0;
        convBool            = 0;
        convBeforeCall      = 0;
        bondFloor           = 0;
        bondFloorGenericMem = 0;
        bondRiskFreeArray   = 0;
        nonEqCreditPvs        = 0;
        control             = 0;
        putStockLevel       = 0.0;
        putProbability      = 0.0;
        firstPutIndex       = -1;
        calcedFirstPutProb  = false;
        isOCB               = false;
        hasE2CLayer         = false;
        resetType           = 0;
        numCVBPriceArrays   = 0;
        pvCFbeforeStart     = 0.0; 

        ocbExerBool = 0;
        ocbStrike = 0;
        ocbRiskyAdj = 0;

        // make the risky curve
        // make sure that we have a BootstrappedYieldCurve otherwise we can't add the spreadCurve
        if (!BootstrappedYieldCurve::TYPE->isInstance(cvb->discount.get()))
        {
            throw ModelException("ConvBond1fProd::ConvBond1fProd", 
                                 "Yield curve must be a cash swap curve");
    }

        DateTime tmp(cvb->getEffMaturityDate());
        DateTimeSP matDate(new DateTime(tmp));

        risky = cvb->getRiskyCurve(matDate.get());
        asset = CAssetWrapper(copy(cvb->asset.get()));
        if (cvb->riskyGrowth && ICanBeRisky::TYPE->isInstance(getEquityAsset().get())) {
            // cast ICanBeRisky pointer
            IObject* obj            = dynamic_cast<IObject*>(getEquityAsset().get());
            ICanBeRisky* riskyAssetI = dynamic_cast<ICanBeRisky*>(obj);
            CreditCurveWrapper* cs = 0;
            ICreditCurveSP stockCurve;
            if ( cvb->useJointDefaultCorrelation ) {
                ExpiryArraySP stockExp(new ExpiryArray(1));
                DoubleArraySP stockCredit(new DoubleArray(1));

                (*stockExp)[0]    = ExpirySP(new MaturityPeriod("30Y"));
                (*stockCredit)[0] = cvb->stockCreditSpread;

                stockCurve = CreditSpreadCurveSP(new CreditSpreadCurve("STOCK_CREDIT", stockExp.get(), *(stockCredit.get())));

                riskyAssetI->makeRisky(stockCurve,matDate.get());
            } else {
                cs = (CreditCurveWrapper*)(&(cvb->creditSpreads));
                riskyAssetI->makeRisky(cs->getSP(),matDate.get());
            }
        }
    }
    
    virtual ~ConvBond1fProd();
    
    virtual CAssetConstSP GetAssetRef() {
        if (hasE2CLayer) {
            return CAssetConstSP::dynamicCast((IObjectConstSP)asset.getSP());
        } else {
            return CAssetConstSP::dynamicCast((IObjectConstSP)getEquityAsset());
        }
    }
    
    virtual bool GetFwdStartLV()
        {
            if ( !!cvb->resetSchedule && cvb->resetSchedule->length() > 0 ) {
                return cvb->resetSchedule->isFwdStart();
            }
            else
                return false;                   
        }
    
    virtual DateTime GetFwdStartDateLV()
        {
            if ( !!cvb->resetSchedule && cvb->resetSchedule->length() > 0 && cvb->resetSchedule->isFwdStart()) {
                return cvb->resetSchedule->getStartDate(cvb->valueDate);
            }
            else
                return cvb->valueDate;
        }
    
    virtual YieldCurveConstSP GetDiscCurveRef();
    
    virtual CVolRequestConstSP GetLNRequest();
    
    virtual YieldCurveConstSP GetRiskFreeCurve();
    
    void InitGridSize();
    
    void InitFD(Control* control);
    
private:
    double InitProd1();
    void InitProd2(double);
public:
    void InitProd();
    
    virtual bool Positive() {
        // option payoff may < 0 if soft call level is less than bond flr.
        return false;
    }
    
    virtual string getCcyTreatment(){ return cvb->ccyTreatment;}
    
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

    virtual void preCalcFDGeneric(int step, int idx, int pStart, int pEnd,
                                  double const * const * optionArray,
                                  double const * const * optionArrayOld);
    
    void PayoffAtMat(const double * s, int step, int bot, int top, int pStart, int pEnd,
                     double * const * price);

    void PayoffBeforeMat(const double * s, int step, int bot, int top, int pStart, int pEnd,
                         double * const * price);

    void preCalcBlackOnCall(int step, double midGrid, double &bocLowerStockPrice, double &bocUpperStockPrice, 
                            DoubleArray &offset, DoubleArray &coefficient);   

    void PenultStraight(const double* s, int step, int bot, int top, double *price);
    
    void PenultMandatory(const double* s, int step, int bot, int top, double *price);

    virtual double getCoupon(int step, const double* s, int start, int end);

    virtual bool hasEquityLayer();

    CAssetSP     getEquityAsset() 
    {
        if ( FirmAsset::TYPE->isInstance(asset.get()) ) {
            FirmAsset* firmAsset = 0;
            if ( FirmAsset::TYPE->isInstance(asset.get())) {
                CAsset*    ncAsset   = const_cast<CAsset*>(asset.get());
                firmAsset = dynamic_cast<FirmAsset*>(ncAsset);
            } else {
                throw ModelException("ConvBondProd::getEquityAsset",
                                     "Underlying asset must be a FirmAsset");
            }

            return firmAsset->getEquityAsset().getSP();
        } else {
            return asset.getSP();
        }
    }


    // following functions for spot dependent spread
    virtual bool isStaticSpread() const { return true; }

    void mergePriceForReset(const double* s, int step, int bot, int top, int pStart, int pEnd,
                            double * const * price);

    inline double adjConversionValueDECDelay(int step, double s, double convertValue, double conversionValue);
    
    inline double calcAccelerateDECCallValue(int step, double s, double bondFloor);

    void calcPutLevelNProb(int step, int bot, int top, const double *s, const double *priceLayer, double putValue, double &putStockLevel, double &putProbability);

    double calcPutToStockValue(int step);

    int calcRiskyBondNAdjPriceLayer(int step, int bot, int top, 
                    int pEnd, double * const * price, double *&bondFloorGeneric);

    void calcBondRiskyNFloorBeforeMat(int step, int bot, int top, double *bondRiskyArray, double *bondFloorArray);

    bool                isOCB;
    int                 numCVBPriceArrays;

    const OptOnConvBond *ocb;

    bool                hasE2CLayer;


protected:

    const ConvBond* cvb;

protected:

    double *coupon;
    bool   *addCouponBool;
    double *callLevel;
    bool   *callBool;
    bool   *softPutBool;
    double *softPutLevel;
    double *softPutTrigger;
    bool   *callIsAdjForAccrued;
    bool   *callIsHard;
    double *putLevel;
    bool   *putBool;
    double *convRatio;
    double *convCash;
    double *parityFloor;
    bool   *convBool;
    bool   *convBeforeCall;
    double *bondFloor;

    double *bondFloorGenericMem;
    double *bondRiskFreeArray;
    double *nonEqCreditPvs;    // extra pv between steps for the portion of clean spread curve that eq does not participate 

    bool   *ocbExerBool;
    double *ocbStrike;
    double *ocbRiskyAdj;
    
    YieldCurveSP        risky;
    CAssetWrapper       asset;

    // members for resettable convertible bond - need to get these into a proper class
    // when the dust has settled
    CDoubleArraySP   minConvPrice;
    CDoubleArraySP   maxConvPrice;
    CDoubleArraySP   resetParity;

    CBoolArraySP     resetArray;
    int              resetType;
    int              firstGridStep;
    CDoubleMatrixSP  conversionRatioGrid;
    CDoubleMatrixSP  conversionRatioGridnoReset; // for triggerReset
    double           actualCR;

    // members required for put probability calculation
    const Control* control;
    int            firstPutIndex;
    double         putStockLevel;
    double         putProbability;
    bool           calcedFirstPutProb;

    // processed vol is required for Black on call
    CVolProcessedBSSP processedVol;

    // some members used to speed up performance;
    double          faceValue;
    bool            hasReset;
    CDoubleMatrix   priceCopy;

    // Currency struck reset schedule
    ResetScheduleSP adjustedResetSchedule;

    // for fwdStarting Resettable
    bool            fwdStarting;
    double          pvCFbeforeStart;

};

ConvBond1fProd::~ConvBond1fProd() {
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
    if (softPutBool != 0) {
        delete [] softPutBool;
        softPutBool = 0;}
    if (softPutLevel != 0) {
        delete [] softPutLevel;
        softPutLevel = 0;}
    if (softPutTrigger != 0) {
        delete [] softPutTrigger;
        softPutTrigger = 0;}
    if (callIsHard != 0) {
        delete [] callIsHard;
        callIsHard = 0;}
    if (callIsAdjForAccrued != 0) {
        delete [] callIsAdjForAccrued;
        callIsAdjForAccrued = 0;}
    if (putLevel != 0) {
        delete [] putLevel;
        putLevel = 0;}
    if (putBool != 0) {
        delete [] putBool;
        putBool = 0;}
    if (convRatio != 0) {
        delete [] convRatio;
        convRatio = 0;}
    if (convCash != 0) {
        delete [] convCash;
        convCash = 0;}
    if (parityFloor != 0) {
        delete [] parityFloor;
        parityFloor = 0;}
    if (convBool != 0) {
        delete [] convBool;
        convBool = 0;}
    if (convBeforeCall != 0) {
        delete [] convBeforeCall;
        convBeforeCall = 0;}
    if (bondFloor != 0) {
        delete [] bondFloor;
        bondFloor = 0;}
    if (bondFloorGenericMem != 0) {
        delete [] bondFloorGenericMem;
        bondFloorGenericMem = 0;}
    if (bondRiskFreeArray != 0) {
        delete [] bondRiskFreeArray;
        bondRiskFreeArray = 0;}
    if (nonEqCreditPvs != 0) {
        delete [] nonEqCreditPvs;
        nonEqCreditPvs = 0;}
    if (ocbExerBool != 0) {
        delete [] ocbExerBool;
        ocbExerBool = 0;}
    if (ocbStrike != 0) {
        delete [] ocbStrike;
        ocbStrike = 0;}
    if (ocbRiskyAdj != 0) {
        delete [] ocbRiskyAdj;
        ocbRiskyAdj = 0;}
}



YieldCurveConstSP ConvBond1fProd::GetDiscCurveRef() {
    if (cvb->convertIntoIssuerShares == true && cvb->riskyGrowth == false) {
        return YieldCurveConstSP::dynamicCast((IObjectConstSP)cvb->discount.getSP());
    } else {
        // exchangeables will be discounted at the risky rate
        if ( cvb->useJointDefaultCorrelation ) {
            DateTime tmp(cvb->getEffMaturityDate());
            DateTimeSP matDate(new DateTime(tmp));
            YieldCurveConstSP riskyDiscountCurve = cvb->getRiskyCurve(matDate.get(), cvb->useJointDefaultCorrelation);
            return riskyDiscountCurve;
        } else {
            return YieldCurveConstSP::dynamicCast((IObjectConstSP)risky);
        }
    }
}


YieldCurveConstSP ConvBond1fProd::GetRiskFreeCurve() {
    return YieldCurveConstSP::dynamicCast((IObjectConstSP)cvb->discount.getSP());
}


void ConvBond1fProd::InitGridSize()
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
                    var = processedVol->CalcVar(cvb->getValueDate(), cvb->bond->getMaturityDate());

                    totalSteps = int(BaseStepAmount*Maths::max(1., var));

                    genericFDModel->stepsPerYearToUse = int(ceil((double)totalSteps/cvb->getValueDate().yearFrac(cvb->bond->getMaturityDate())));
            
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
                    int totalSteps = int(genericFDModel->stepsPerYearToUse * cvb->getValueDate().yearFrac(cvb->bond->getMaturityDate()));
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
                    var = processedVol->CalcVar(cvb->getValueDate(), cvb->bond->getMaturityDate());

                    totalSteps = int(BaseStepAmount*Maths::max(1., var));

                    fdModel->stepsPerYearToUse = int(ceil((double)totalSteps/cvb->getValueDate().yearFrac(cvb->bond->getMaturityDate())));
            
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
                    int totalSteps = int(fdModel->stepsPerYearToUse * cvb->getValueDate().yearFrac(cvb->bond->getMaturityDate()));
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

void ConvBond1fProd::InitFD(Control* control)
{
    static const string method = "ConvBond1fProd::InitFD";
    try {

        DateTime valDate = cvb->valueDate;
        hasReset  = (!!cvb->resetSchedule && cvb->resetSchedule->length() > 0 );
        fwdStarting = hasReset?cvb->resetSchedule->isFwdStart():false;

        if (fwdStarting)
            valDate = cvb->resetSchedule->getStartDate(valDate);

        // required to calculate put probability as part of the PayoffBeforeMat callback
        this->control = control;

        // critical dates
        DateTimeArray critDates;

        // some validation
        if (  !cvb->convertIntoIssuerShares && FD1FE2C::TYPE->isInstance(genericFDModel)) {
           throw ModelException("ConvBond1fProd::InitProd",
                     "The E2C model is not supported for Echangeables");
        }

        int i;
        CashFlowArraySP myCFs = cvb->bond->getExAdjCashFlows(valDate, risky); 
        for (i=0; i<myCFs->size(); i++) {
            critDates.push_back((*myCFs)[i].date);
            critDates.push_back((*myCFs)[i].date.rollDate(-1));
        }

        // add all dividend pass through dates to the critical dates
        DividendListSP divSchedule;
        if ( cvb->dividendPassThrough ) {
            divSchedule = DividendCollector::divsBetweenDates(getEquityAsset().get(),
                                                              valDate,
                                                              valDate,
                                                              cvb->bond->getMaturityDate());

            DateTimeArrayConstSP payDates = divSchedule->getPayDates();
            for (i=0; i<payDates->size(); i++) {
                critDates.push_back((*payDates)[i]);
            }
        }

        if ( cvb->putSchedule.get() ) {
            const DateTimeArray& putDates = cvb->putSchedule->getDateArray();
            for (i=0; i<putDates.size(); i++)
                critDates.push_back(putDates[i]);  
        }
        if ( cvb->conversionRatios.get() ) {
            const DateTimeArray& conversionDates = cvb->conversionRatios->getDateArray();
            for (i=0; i<conversionDates.size(); i++)
                critDates.push_back(conversionDates[i]);
        }
        // push in additional addon conversion dates for embedded warrant
        if( cvb->hasAddOnConvRatios() )
            cvb->addOnConvRatios->insertCritDates(critDates);

        // push in additional dates for contingenet conversion
        if( cvb->hasContConversion() )
            cvb->contConversion->insertCritDates(critDates);

        if (isOCB == true) {
            const DateTimeArray& exerDates = ocb->exerSched->getDateArray();
            for (i=0; i<exerDates.size(); i++)
                critDates.push_back(exerDates[i]);
        }
        
        double minGap = -1; // default value
        
        DateTimeArray softCallDates, hardCallDates, softPutDates;
        if (cvb->softCallTriggerSchedule.get()) {
            softCallDates = cvb->softCallTriggerSchedule->getDates();
        }
        if (cvb->callSchedule.get()) {
            hardCallDates = cvb->callSchedule->getDates();

            // move critical dates forward for Black-on-call
            if ( cvb->callTreatment == "B") {
               for (i=0; i<hardCallDates.size(); i++) {
                   hardCallDates[i] = hardCallDates[i].rollDate(-cvb->callNotification);
               }
            }
        }
        if (cvb->softPutTriggerSchedule.get()) {
            softPutDates = cvb->softPutTriggerSchedule->getDates();
        }

        // All call dates are critical dates
        for (i=0; i<softCallDates.size(); i++)
            critDates.push_back(softCallDates[i]);
        for (i=0; i<hardCallDates.size(); i++)
            critDates.push_back(hardCallDates[i]);
        for (i=0; i<softPutDates.size(); i++)
            critDates.push_back(softPutDates[i]);
   
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
        segDates[0] = valDate;
        if (!hasE2CLayer) {
            for (i=0; i<hardCallDates.size(); i++) {
                if (hardCallDates[i] > valDate && hardCallDates[i] < cvb->bond->getMaturityDate()) {
                    if ( i == 0 || Maths::isPositive(timeMetric->yearFrac(segDates[segDates.size()-1], hardCallDates[i]))) {
                        segDates.push_back(hardCallDates[i]);
                        break;
                    }
                }
            }
        }

        // all reset dates are critDates dates
        if ( !!cvb->resetSchedule && cvb->resetSchedule->length() > 0 ) {
            const DateTimeArray& resetDates = cvb->resetSchedule->getDateArray();
            for(i=0;i<resetDates.size();++i) {
                critDates.push_back(resetDates[i]);
            }
        }

        segDates.push_back(cvb->bond->getMaturityDate());

        // sort critical dates
        DateTime minDate = valDate.rollDate(-1);
        DateTime maxDate = cvb->getBondMaturityDate().rollDate(+1);
        timeMetric->SortDate(minDate, maxDate, true, segDates);

        vector<int> density;
        density.resize(segDates.size()-1);
        for (i=0; i<segDates.size()-1; i++) {
            density[i] = 1;
        }
        bool useEqualTime = true;

        int numOfPriceArray;

        if (isOCB) {
            if ( cvb->triggerReset ) {
                throw ModelException(method, "Cannot price asset swaps on trigger resettables");
            }
            if ( cvb->hasContConversion() ) {
                throw ModelException(method, "Cannot price asset swaps on cvb with contingenet conversion");
            }

            if ( !!cvb->resetSchedule && cvb->resetSchedule->length() > 0 ) {
                if ( FD1FLNResettable::TYPE->isInstance(fdModel)) {
                    FD1FLNResettable* resettableModel = dynamic_cast<FD1FLNResettable*>(fdModel);
                    numOfPriceArray   =  2 * resettableModel->getNum3DSteps();
                    numCVBPriceArrays =      resettableModel->getNum3DSteps();
                } else {
                    // default number of conversion ratio samples
                    numOfPriceArray   = 30;
                    numCVBPriceArrays = 15;
                }
            } else {
                // four price arrays for trigger reset, 1 otherwise
                // left trigger reset in here in case I ever get around to implement it
                if (cvb->triggerReset) {
                    numOfPriceArray   =  4;
                    numCVBPriceArrays =  2;
                } else {
                    numOfPriceArray   =  2;
                    numCVBPriceArrays =  1;
                }
            }
        } else {
            if ( !!cvb->resetSchedule && cvb->resetSchedule->length() > 0 ) {
                if ( FD1FLNResettable::TYPE->isInstance(fdModel)) {
                    FD1FLNResettable* resettableModel = dynamic_cast<FD1FLNResettable*>(fdModel);
                    numOfPriceArray =  resettableModel->getNum3DSteps();
                } else {
                    // default number of conversion ratio samples
                    numOfPriceArray = 15;
                }
            } else {
                // two price arrays for trigger reset, 1 otherwise
                if (cvb->triggerReset || 
                    // if contingent conversion and convertible for rest of life, 1 array
                    cvb->hasContConversion() )
                {
                    numOfPriceArray =  2;
                } else {
                    numOfPriceArray =  1;
                }
            }
            
        }

        // add a last price array for bondrisky
        if( !isStaticSpread() ) numOfPriceArray++;

        if (fdModel) {
            fdModel->TimePts.SetTreeStepBase(500);
            // the try/catch block is here to avoid a mysterious crash in gcc on Solaris optimised
            try { 
                fdModel->Setup(valDate, segDates, density, &critDates,
                               minGap, useEqualTime, numOfPriceArray);
            } catch (exception&) {
                throw;
            }
        } else {
            genericFDModel->TimePts.SetTreeStepBase(500);

            if (hasE2CLayer) 
                numOfPriceArray = numOfPriceArray+2;

            // the try/catch block is here to avoid a mysterious crash in gcc on Solaris optimised
            try { 
                genericFDModel->Setup(valDate, segDates, density, &critDates,
                               minGap, useEqualTime, numOfPriceArray);
            } catch (exception&) {
                throw;
            }

        }
    } catch (exception& e) {
        throw ModelException(&e, method);
    }
}

void ConvBond1fProd::InitProd() {
  // See :-( below
  InitProd2(InitProd1());
}

double ConvBond1fProd::InitProd1()
{
    static const string method = "ConvBond1fProd::InitProd";

    int numStep           = model1F->TimePts.NumOfStep;
    coupon                = new double[numStep+1];
    addCouponBool         = new bool[numStep+1];

    // for performance
    faceValue = cvb->bond->getFaceValue();
    hasReset  = (!!cvb->resetSchedule && cvb->resetSchedule->length() > 0 );

    // prepare for fwdStarting case
    fwdStarting = hasReset?cvb->resetSchedule->isFwdStart():false;
    DateTime startDate = fwdStarting ? cvb->resetSchedule->getStartDate(cvb->valueDate):cvb->valueDate;
    double fwdAtStart = cvb->asset->fwdValue(startDate);
    if (fwdStarting) {
        // :-( Under Solaris/gcc3.2, we dump core here if getStartLevel throws
        // an exception (test resettableinp/fwdRst24.xml), because the
        // compiler tries to be too clever about deconstructing the objects
        // that thereby go out of scope.

        fwdAtStart = cvb->resetSchedule->getStartLevel(cvb->valueDate, fwdAtStart);   // use Known Value
    }
    return fwdAtStart;
}

void ConvBond1fProd::InitProd2(double fwdAtStart)
{
    static const string method = "ConvBond1fProd::InitProd";
    
    int i;
    CashFlowArraySP myCFs = cvb->bond->getExAdjCashFlows(cvb->valueDate, risky); 
    int numStep           = model1F->TimePts.NumOfStep;
    int cfIdx             = 0;
    double putLevelAtMat  = 0.0;
    DateTime startDate = fwdStarting ? cvb->resetSchedule->getStartDate(cvb->valueDate):cvb->valueDate;
    
    int cfIdxFirst = 0;
    if (fwdStarting){
        while (startDate.getDate()>(*myCFs)[cfIdxFirst].date.getDate()){
            cfIdxFirst++;   // finding for the first cfIdx after TimePts.StepDates[0]
        }
    }
    cfIdx = cfIdxFirst;
                
    // generate adjusted reset schedule 
    if ( !!cvb->resetSchedule) {
        adjustedResetSchedule = ResetScheduleSP(cvb->resetSchedule.clone());

        // convert the schedule to bond currency for currency struck bonds
        if ( StruckEquity::TYPE->isInstance(getEquityAsset().get())) {
            StruckEquity* struckAsset = dynamic_cast<StruckEquity*>(getEquityAsset().get());
            double       fxSpot      = struckAsset->getFXSpot();

            // get all reset dates
            DateTimeArray resetDates = adjustedResetSchedule->getDates();
            DoubleArraySP fxFwds(new DoubleArray(resetDates.size()));
            for(int i=0 ; i<fxFwds->size() ; ++i) {
                if (resetDates[i] >= cvb->valueDate ) {
                    (*fxFwds)[i] = struckAsset->fxFwdValue(resetDates[i]);
                } else {
                    (*fxFwds)[i] = 0.0;
                }
            }

            // calculate initial conversion price
            double initialCR        = cvb->conversionRatios->firstValue();
            if (fwdStarting){
                initialCR /= fwdAtStart;
                // fwdAtStart is already in payoff currency Unit (not Underlying Ccy Unit)
                // in a few lines later, adjustedResetSchedule->preProcessSchedule(....)
                // fx rate is multiplying ResetLevels, by assuming that reset level are given by Underlying Ccy.
                // Thus, the Reset Level at heare should be Ccy Unit of underling.  Not payoff.
                adjustedResetSchedule->scaleLevels(fwdAtStart/struckAsset->fxFwdValue(startDate));
            }
            if (Maths::isZero(initialCR)) {
                initialCR = DBL_EPSILON;
            }
            double initialConvPrice = cvb->bond->getFaceValue() / initialCR;

            // convert the reset schedule into a currency struck schedule, if necessary
            adjustedResetSchedule->preProcessSchedule(true, fxSpot, fxFwds, initialConvPrice);
        } else {
            double       fxSpot      = 1.0;

            // get all reset dates
            DateTimeArray resetDates = adjustedResetSchedule->getDates();
            DoubleArraySP fxFwds(new DoubleArray(resetDates.size()));
            for(int i=0 ; i<fxFwds->size() ; ++i) {
                (*fxFwds)[i] = 1.0;
            }

            // calculate initial conversion price
            double initialCR        = cvb->conversionRatios->firstValue();
            if (fwdStarting){
                initialCR /= fwdAtStart;
                adjustedResetSchedule->scaleLevels(fwdAtStart);
            }
            if (Maths::isZero(initialCR)) {
                initialCR = DBL_EPSILON;
            }
            double initialConvPrice = cvb->bond->getFaceValue() / initialCR;

            // convert the reset schedule into a currency struck schedule, if necessary
            adjustedResetSchedule->preProcessSchedule(false, fxSpot, fxFwds, initialConvPrice);
        }
    }


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
                    if (cvb->getCouponAtMat == true) {
                        addCouponBool[i] = true;
                        coupon[i] = (*myCFs)[cfIdx].amount - cvb->bond->getRedemption();
                    } else {
                        addCouponBool[i] = false;
                        coupon[i] = 0.;
                    }
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

    // add the dividend schedule 
    if ( cvb->dividendPassThrough ) {
        DividendListSP divSchedule = DividendCollector::divsBetweenDates(getEquityAsset().get(),
                                                          cvb->valueDate,
                                                          cvb->valueDate,
                                                          cvb->bond->getMaturityDate());

        if (!!divSchedule) {
            const DividendArray& divArray = divSchedule->getArray();

            double convRatio = 0.0;
            if (cvb->DECS == true ) {
               convRatio = cvb->initialConvRatio;
            } else {
               if ( cvb->conversionRatios.get() && cvb->conversionRatios->length() > 0) {
                  convRatio = cvb->conversionRatios->firstValue();
               }
            }

            // assume the dividends are in order
            cfIdx = cfIdxFirst;
            for (i=0; i<=numStep; i++) {
                if ( cfIdx <  divArray.size() && 
                    model1F->TimePts.StepDates[i].getDate() == divArray[cfIdx].getPayDate().getDate() ) {
                    if ( i == 0 || 
                         model1F->TimePts.StepDates[i].getDate() != 
                         model1F->TimePts.StepDates[i-1].getDate() ) {
                        addCouponBool[i] = true;
                        // I might have a problem here with coupons at maturity
                        if (cfIdx == divArray.size()-1) {
                            if (cvb->getCouponAtMat == true) {
                                addCouponBool[i] = true;
                                coupon[i] += convRatio * cvb->dividendPassThroughPct * divArray[cfIdx].getDivAmount();
                            }
                        } else {
                            coupon[i] += convRatio * cvb->dividendPassThroughPct * divArray[cfIdx].getDivAmount();
                        }
                        cfIdx++;
                    }
                }
            }
        }
    }

    callLevel           = new double[numStep+1];
    callBool            = new bool[numStep+1];
    callIsAdjForAccrued = new bool[numStep+1];
    callIsHard          = new bool[numStep+1];
    
    double dummy = getEquityAsset()->getSpot();
    vector<double> fwdVol;
    for (i=0; i<=numStep; i++) {
        model1F->GetStepVol(i, fwdVol, &dummy,0,0);
        cvb->getCallLevel(model1F->TimePts.StepDates[i], fwdVol[0], &(callBool[i]), &(callIsAdjForAccrued[i]), 
            &(callIsHard[i]), &(callLevel[i]));
    } 

    putLevel   = new double[numStep+1];
    putBool    = new bool[numStep+1];
    
    for (i=0; i<=numStep; i++) {
        cvb->getPutLevel(model1F->TimePts.StepDates[i], &(putBool[i]), &(putLevel[i]));
        // first put date is required for put probability calculation
        if ( putBool[i] && firstPutIndex < 0 ) {
            firstPutIndex = i;
        }
    }

    softPutBool         = new bool[numStep+1];
    softPutLevel        = new double[numStep+1];
    softPutTrigger      = new double[numStep+1];

    for (i=0; i<=numStep; i++) {
        cvb->getSoftPutLevel(model1F->TimePts.StepDates[i], &(softPutBool[i]), &(softPutTrigger[i]), &(softPutLevel[i]));
    }
    
    convRatio   = new double[numStep+1];
    convCash    = new double[numStep+1];
    parityFloor = new double[numStep+1];
    convBool    = new bool[numStep+1];
    
    for (i=0; i<=numStep; i++) {
        cvb->getConversionInfo(model1F->TimePts.StepDates[i], 0/*dummyspot*/, &(convBool[i]), &(convRatio[i]), &(convCash[i]), true, &(parityFloor[i]));
        if ( Maths::isZero(convRatio[i])) {
            convBool[i] = false;
        }
        if (fwdStarting)
            convRatio[i] /= fwdAtStart;
    }

    convBeforeCall = new bool[numStep+1];
    
    convBeforeCall[0] = true; // handle the first point specially. no harm in setting it true. OCB will be wrong if set to false.
    // There can be some trouble with setting convBeforeCall[0] = true. The cvb value can be above the call level 
    // when parity is below. Should fix this in way that doesn't screw up OCB.
    for (i=1; i<=numStep; i++) {
        if (callBool[i] == false){
            // can't call now
            convBeforeCall[i] = true; // must be true or OCB will screw up
        } else {
            if (callBool[i-1] == false) {
                // can call now but not at previous step
                convBeforeCall[i] = true;
            } else {
                // can call now and at previous step
                if (callLevel[i] < callLevel[i-1]*.9975) { // give 25 bp leeway
                    convBeforeCall[i] = true;
                } else {
                    convBeforeCall[i] = false;
                } 
            }
        }
    }

    if (firstPutIndex > 0 && model1F->TimePts.StepDates[firstPutIndex] <= cvb->riskFreeCouponEndDate)
        throw ModelException("ConvBond1fProd::InitProd", 
                             "the model does not handle puts that are at or before the end date for risk free coupons");

    bondFloor = new double[numStep+1];
    if (cvb->getCouponAtMat == true){
        if (cvb->DECS == true || cvb->PERCS == true) {
            bondFloor[numStep] = 0.; // note that getCouponAtMat is validated to be true.
        } else {
            if (putBool[numStep] == true) {
                bondFloor[numStep] = Maths::max(cvb->bond->getRedemption(), putLevel[numStep]);
            } else {
                bondFloor[numStep] = cvb->bond->getRedemption();
            }
        }
    } else {
        if (putBool[numStep] == true) {
            // we always get the coupon if we put at mat - only in the first case is the accrued interest
            // not yet included in the put level
            if ( !cvb->putAdjustForAccrued ) {
                putLevelAtMat = putLevel[numStep] + cvb->getAccruedAtDate(model1F->TimePts.StepDates[numStep]);
            } else {
                putLevelAtMat = putLevel[numStep];
            }
            bondFloor[numStep] = Maths::max((*myCFs)[myCFs->size()-1].amount, putLevelAtMat);
        } else {
            bondFloor[numStep] = (*myCFs)[myCFs->size()-1].amount;
        }
    }

    double  bondRiskFree = 0;
    double  bondRisky = bondFloor[numStep];

    // check for stock being in Escrow account
    if (cvb->stockInEscrow) {
        double parity     = getEquityAsset()->getSpot() * convRatio[numStep];
        double difference = parity - bondRiskFree;
        bondRiskFree += difference;
        bondRisky    -= difference;
    }

    if( !isStaticSpread() )
    {
        bondRiskFreeArray = new double[numStep+1];
        bondRiskFreeArray[numStep] = bondRiskFree;
        
        // make sure to get non equity portion of spread
        if( dynamic_cast<IHaveNonEqSpread *>(model1F) )
        {
            const CleanSpreadCurve *nonEqCleanSprds = dynamic_cast<IHaveNonEqSpread *>(model1F)->getNonEqSpreads();
            if( nonEqCleanSprds == 0 )
                throw ModelException(method, "Model has non eq spread but no eq spread is returned");

            nonEqCreditPvs = new double[numStep];
            for(i=0; i<numStep; i++)
                nonEqCreditPvs[i] = nonEqCleanSprds->getDefaultPV(model1F->TimePts.StepDates[i], model1F->TimePts.StepDates[i+1]);
        }
    }

    for (i=numStep-1; i>=0; i--) {
        double accruedIR = cvb->getAccruedAtDate(model1F->TimePts.StepDates[i]);
        // We assume here that puts cannot be before riskFreeCouponEndDate
        if ( putBool[i+1] == true && 
             !( model1F->TimePts.StepDates[i].getDate() == 
                model1F->TimePts.StepDates[model1F->TimePts.StepDates.size()-1].getDate()) &&
                !((cvb->DECS == true || cvb->PERCS == true) && i == numStep-1)) {
            bondRisky = Maths::max(bondRisky, putLevel[i+1]);
        } 

        if (addCouponBool[i+1] == true) {
            if (model1F->TimePts.StepDates[i+1] <= cvb->riskFreeCouponEndDate)
                bondRiskFree += coupon[i+1];
            else
                bondRisky += coupon[i+1];
        }

        if ( (model1F->TimePts.StepDates[i+1] <= cvb->riskFreeCouponEndDate && bondRiskFree > 0) ||
             cvb->stockInEscrow ) {
            bondRiskFree *= cvb->discount.getSP()->pv(model1F->TimePts.StepDates[i], 
                                 model1F->TimePts.StepDates[i+1]);
        }

        // use PV recovery rule for preferreds
        if (cvb->DECS || cvb->PERCS) {
            bondRisky = risky->riskyPV(model1F->TimePts.StepDates[i],
                                       model1F->TimePts.StepDates[i+1],
                                       bondRisky,
                                       bondRisky,
                                       cvb->useAssetRecovery,
                                       cvb->recoveryPct);
        } else {
            bondRisky = risky->riskyPV(model1F->TimePts.StepDates[i],
                                       model1F->TimePts.StepDates[i+1],
                                       bondRisky,
                                       cvb->bond->getNotional(model1F->TimePts.StepDates[i]) + accruedIR,
                                       cvb->useAssetRecovery,
                                       cvb->recoveryPct);
        }

        if (callIsHard[i] && callLevel[i] < bondRisky + bondRiskFree)
        {
            bondRisky = Maths::max((double)callLevel[i]-bondRiskFree, 0.);
            if (bondRisky == 0)
                bondRiskFree = callLevel[i];
        }
        bondFloor[i] = bondRisky + bondRiskFree;

        if( !isStaticSpread() )
        {
            bondRiskFreeArray[i] = bondRiskFree;
        }
    }

    // add all coupons before fwdStartdate, after value date.
    if (fwdStarting){
        pvCFbeforeStart = cvb->bond->couponsPV(cvb->valueDate,startDate,risky); // pv of cpns before starting.
        double dfRiskFree = cvb->discount->pv(cvb->valueDate,startDate);
        double dfRisky    = risky->pv(cvb->valueDate,startDate);
        double riskyBondPVAfterStart = cvb->bond->presentValue(cvb->valueDate, risky) - pvCFbeforeStart;
        riskyBondPVAfterStart *= (1.0-dfRiskFree/dfRisky);
        pvCFbeforeStart += riskyBondPVAfterStart;
        pvCFbeforeStart /= dfRiskFree;
    }

    if (isOCB) {
        try {
            ((OptOnConvBond *)ocb)->strikeAtOptionMat = cvb->bond->getNotional(ocb->exerSched->lastDate());
        } catch (exception&) {
            ((OptOnConvBond *)ocb)->strikeAtOptionMat = faceValue;
        }
    }

    
    if (isOCB == true) {
        ocbExerBool  = new bool[numStep+1];
        ocbStrike    = new double[numStep+1];
        for (i=0; i<=numStep; i++) {
            ocb->getStrike(model1F->TimePts.StepDates[i], bondFloor[i], true, &(ocbExerBool[i]), &(ocbStrike[i]));  
        }

        ocbRiskyAdj = new double[numStep+1];
        ocbRiskyAdj[numStep] = 0.;
        long lastExerIdx = -1;
        for (i=numStep-1; i>=0; i--) {
            if (ocbExerBool[i+1] == true) {
                ocbRiskyAdj[i] = Maths::max(bondFloor[i+1] - ocbStrike[i+1], 0.);
                lastExerIdx = i+1;
            } else if (lastExerIdx < 0) {
                ocbRiskyAdj[i] = 0;
            } else {
                ocbRiskyAdj[i] = Maths::max(bondFloor[lastExerIdx] - ocbStrike[lastExerIdx], 0.);
            }
            if (lastExerIdx > i+1) {
                 ocbRiskyAdj[i] *= 
                    risky->pv(model1F->TimePts.StepDates[i+1], model1F->TimePts.StepDates[lastExerIdx]);
            }
            if (lastExerIdx >= 0) {
                ocbRiskyAdj[i] *= 
                    (risky->pv(model1F->TimePts.StepDates[i], model1F->TimePts.StepDates[i+1]) -
                    cvb->discount->pv(model1F->TimePts.StepDates[i], model1F->TimePts.StepDates[i+1]));
            }
        }
        // a little validation that's hard to do earlier
        bool isExerAfter = false;
        for (i=numStep-1; i>=0; i--) {
            if (isExerAfter == false && ocbExerBool[i] == true) {
                isExerAfter = true;
            }
            if (isExerAfter == true && callBool[i] == true && ocbExerBool[i] == false) {
                throw ModelException("ConvBond1fProd::InitProd", 
                                     "option on convertible must be exercisable whenever it is callable");
            }
        }
    }

    if (model1F->isTree() == false && !hasReset && isStaticSpread() ) { 
        // no segment is not static spread due to possible FD extrapolation error 
        // across segments and FD boundary conditions at call barrier 
        // Set the max stock level for each segment

        bool fwdFlag;
        if (fdModel) {
            fwdFlag = fdModel->useFwdGrid;
        } else {
            fwdFlag = genericFDModel->useFwdGrid;
        }

        // AM: Need to get forward prices here to accurately compute upper grid boundary when using forward grid   
        CDoubleArray forwards;
        if (true == fwdFlag) {
            forwards.resize(numStep+1);
            asset->fwdValue(model1F->TimePts.StepDates, forwards);
        }

        int     j;
        double  gridMax;
        double  maxFwd;    // Needed for min forward grid boundary
        double  stockBound;
        bool    newStockMax;
        double  maxStock = 0;

        for (j=0; j<int(model1F->TimePts.SegmentEnd.size()); j++) {
            gridMax = 0.;
            maxFwd = 0.;
            newStockMax = true;
            for (i=model1F->TimePts.SegmentStart[j]; i<=model1F->TimePts.SegmentEnd[j]; i++) {
                if (callBool[i] == false || convBool[i] == false || hasE2CLayer) {
                    newStockMax = false;
                    break; // keep the default truncation level
                }

                if (addCouponBool[i] == true && callIsAdjForAccrued[i] == true) {
                    stockBound = (callLevel[i]+coupon[i]-convCash[i])/convRatio[i];
                } else {
                    stockBound = (callLevel[i]-convCash[i])/convRatio[i];
                }

                if (false == fwdFlag) {
                    gridMax = Maths::max(gridMax, stockBound);
                } else {
                    maxFwd = Maths::max(forwards[i], maxFwd);
                    gridMax = Maths::max(gridMax, stockBound/forwards[i]);
                }
                
            }

            if (newStockMax == true && !hasE2CLayer) {
                if (fdModel) {
                    fdModel->stockMaxSeg[j]        = gridMax;
                } else {
                    genericFDModel->stockMaxSeg[j] = gridMax;
                }
            }

            if (hasE2CLayer) {

                FirmAsset* firmAsset = 0;
                if ( FirmAsset::TYPE->isInstance(asset.get())) {
                    CAsset*    ncAsset   = const_cast<CAsset*>(asset.get());
                    firmAsset = dynamic_cast<FirmAsset*>(ncAsset);
                } else {
                    throw ModelException("FD1FE2C::PostSetup",
                            "Underlying credit asset must be a FirmAsset");
                }

                FD1FE2C* e2cModel = 0;
                if ( FD1FE2C::TYPE->isInstance(genericFDModel)) {
                    e2cModel = dynamic_cast<FD1FE2C*>(genericFDModel);
                } else {
                    throw ModelException("ConvBond1fProd::InitProd",
                            "The model must be an E2C model");
                }


                genericFDModel->stockMinSeg[j] = firmAsset->getDefaultBarrier();

                // have to set the max equity layers 
                newStockMax = true;
                for (i=model1F->TimePts.SegmentStart[j]; i<=model1F->TimePts.SegmentEnd[j]; i++) {
                    if (callBool[i] == false || convBool[i] == false ) {
                        newStockMax = false;
                        break; // keep the default truncation level
                    }
                    if (addCouponBool[i] == true && callIsAdjForAccrued[i] == true) {
                        maxStock = Maths::max(maxStock, (callLevel[i]+coupon[i]-convCash[i])/convRatio[i]);
                    } else {
                        maxStock = Maths::max(maxStock, (callLevel[i]-convCash[i])/convRatio[i]);
                    }
                }

                if (newStockMax == true) {
                    e2cModel->equityMaxSeg[j] = maxStock;
                }

                // fwd grid not currently supported properly
                // if (true == fwdFlag) {
                //    genericFDModel->stockMinSeg[j] /= maxFwd;
                // }
            }
        }
    }


    if ( genericFDModel ) {
        if (isOCB) {
            numCVBPriceArrays = (genericFDModel->numPriceArrays - (isStaticSpread()?0:1)) / 2;
        } else {
            numCVBPriceArrays = genericFDModel->numPriceArrays - (isStaticSpread()?0:1);
        }
    } else if ( fdModel ) {
        if (isOCB) {
            numCVBPriceArrays = fdModel->numPriceArrays / 2;
        } else {
            numCVBPriceArrays = fdModel->numPriceArrays;
        }
    }

    if (hasReset) {
        
        int numStep             = model1F->TimePts.NumOfStep;
        minConvPrice            = CDoubleArraySP(new CDoubleArray(numStep + 1));
        maxConvPrice            = CDoubleArraySP(new CDoubleArray(numStep + 1));
        resetParity             = CDoubleArraySP(new CDoubleArray(numStep + 1));
        
        resetArray              = CBoolArraySP(new BoolArray(numStep + 1));
        firstGridStep           = numStep+1;
        int numOfPriceArray     = 0;
        
        
        if ( FD1FLNResettable::TYPE->isInstance(fdModel)) {
            const FD1FLNResettable* resettableModel = dynamic_cast<const FD1FLNResettable*>(fdModel);
            numOfPriceArray =  resettableModel->getNum3DSteps();
        } else {
            // default number of conversion ratio samples
            numOfPriceArray =  15;
        }

        conversionRatioGrid              = CDoubleMatrixSP(new CDoubleMatrix(numStep + 1, numOfPriceArray));
        conversionRatioGridnoReset       = CDoubleMatrixSP(new CDoubleMatrix(numStep + 1, numOfPriceArray));
        
        double initCVR = cvb->conversionRatios->firstValue();
        if (fwdStarting) {
            initCVR /=  fwdAtStart;
        }
        
        double crAtValueDate    = adjustedResetSchedule->getCurrentConversionRatio(
            initCVR,
            model1F->TimePts.StepDates[0],
            faceValue);
        
        for (i=0; i<=numStep; i++) {
            // check if there is a reset today
            if ( adjustedResetSchedule->hasReset(model1F->TimePts.StepDates[i])) {
                (*resetArray)[i]     = true;

                // maximum conversion price
                (*maxConvPrice)[i]   = adjustedResetSchedule->getMaxResetStrike(model1F->TimePts.StepDates[i]); 

                // minimum conversion price
                (*minConvPrice)[i]   = adjustedResetSchedule->getMinResetStrike(model1F->TimePts.StepDates[i]); 
               
                // minimum conversion price
                (*minConvPrice)[i]   = Maths::min((*minConvPrice)[i],(*maxConvPrice)[i]);
                (*resetParity)[i]    = adjustedResetSchedule->getParity(model1F->TimePts.StepDates[i]); 
                
                if (i < firstGridStep) {
                    firstGridStep =  i;
                }
            } 
            else {
                if ( i > 0 ) {
                    (*resetArray)[i]     = false;
                    (*maxConvPrice)[i]   = (*maxConvPrice)[i-1];
                    (*minConvPrice)[i]   = (*minConvPrice)[i-1];
                    (*resetParity)[i]    = (*resetParity)[i-1];
                } 
                else {
                    (*resetArray)[i]     = false;
                    double initialCR;
                    if ( model1F->TimePts.StepDates[i] < cvb->conversionRatios->firstDate()) {
                        initialCR = cvb->conversionRatios->firstValue();
                        if (fwdStarting)
                            initialCR /= fwdAtStart;
                    } 
                    else {
                        initialCR = convRatio[i];
                    }

                    double currentCR     = adjustedResetSchedule->getCurrentConversionRatio(
                        initialCR,
                        model1F->TimePts.StepDates[i],
                        faceValue);
                    (*maxConvPrice)[i]   = faceValue / currentCR;
                    (*minConvPrice)[i]   = faceValue / currentCR;
                    (*resetParity)[i]    = 1.0;
                    actualCR=adjustedResetSchedule->getCurrentConversionRatio(
                        initialCR, model1F->TimePts.StepDates[i], faceValue);
                }
            }

            // set up the conversion ratio grid
            int j;
            double minCR = faceValue / (*maxConvPrice)[i];
            double maxCR = faceValue / (*minConvPrice)[i];

            string resetTypeStr = adjustedResetSchedule->getResetType(model1F->TimePts.StepDates[i]);

            if ( resetTypeStr == "UP" && adjustedResetSchedule->isFlat(model1F->TimePts.StepDates[i])) {
                // up reset only
                minCR = Maths::max(minCR, crAtValueDate);
                maxCR = Maths::max(minCR, maxCR);
            }

            for (j=0;j<numOfPriceArray;++j) {
                (*conversionRatioGrid)[i][j] = exp(log(minCR) + (log(maxCR) - log(minCR)) /( numOfPriceArray-1) * j);
                //(*conversionRatioGridnoReset)[i][j] = exp(log(minCR) + (log(maxCR) - log(minCR)) /( numOfPriceArray-1) * j);
            }
            
            if (hasE2CLayer) {
                FirmAsset* firmAsset = 0;
                if ( FirmAsset::TYPE->isInstance(asset.get())) {
                    CAsset*    ncAsset   = const_cast<CAsset*>(asset.get());
                    firmAsset = dynamic_cast<FirmAsset*>(ncAsset);
                } else {
                    throw ModelException("FD1FE2C::PostSetup",
                            "Underlying credit asset must be a FirmAsset");
                }

                FD1FE2C* e2cModel = 0;
                if ( FD1FE2C::TYPE->isInstance(genericFDModel)) {
                    e2cModel = dynamic_cast<FD1FE2C*>(genericFDModel);
                } else {
                    throw ModelException("ConvBond1fProd::InitProd",
                            "The model must be an E2C model");
                }

                for (int j=0; j<int(model1F->TimePts.SegmentEnd.size()); j++) {
                    genericFDModel->stockMinSeg[j] = firmAsset->getDefaultBarrier();
                }
            }
        }
    } 
    else {
        firstGridStep = 0;
    }
    
    // set vol for trigger adjustment in contingenet conversion
    if( cvb->hasContConversion() ) {
        CVolRequestConstSP volRequest = GetLNRequest();
        CVolProcessedSP  procVol(cvb->asset->getProcessedVol(volRequest.get()));
        CVolProcessedBSSP volBS = CVolProcessedBSSP::dynamicCast(procVol);
        cvb->contConversion->setVol(volBS);
    }
}

/** calculate barriers if needed */
void ConvBond1fProd::preCalcFD(int step, int idx, int pStart, int pEnd)
{
    double couponToAdd;
    int i;
    int layer;

    double conversionRatio;
    double conversionRatioNext;
    double conversionCash = convCash[step];
    double conversionCashNext;
    bool isConvertible; // dummy, value not used
    
    if (fdModel) {
        conversionCashNext = (step<fdModel->TimePts.NumOfStep?convCash[step+1]:0);
        for (layer=0 ; layer<= pEnd ; ++layer) {
            if( cvb->hasContConversion()) {
                i = layer;
            }
            else {
                i = 0;
            }
            
            conversionRatio     = convRatio[step];
            conversionRatioNext = (step < fdModel->TimePts.NumOfStep)?convRatio[step+1]:0.0;
          
            if (step < fdModel->TimePts.NumOfStep) {
                // no barrier if contConversion
                if ( callBool[step] && convBool[step]  && !cvb->hasContConversion() ) {
                    if( cvb->hasAddOnConvRatios() )
                        cvb->getConversionInfo(fdModel->TimePts.StepDates[step], callLevel[step], 
                                               &isConvertible, &conversionRatio, &conversionCash, false/*is not spot*/);
                        
                    fdModel->fdEngine.upBarrierNew[i]       = (callLevel[step] - conversionCash)/conversionRatio;
                    fdModel->fdEngine.upPayoutNew[i]        = conversionCash - bondFloor[step];
                    fdModel->fdEngine.upPayoutDeltaNew[i]   = conversionRatio;
                } 
                else {
                    fdModel->fdEngine.upBarrierNew[i]       = -1;
                    fdModel->fdEngine.upPayoutNew[i]        = 0.;
                    fdModel->fdEngine.upPayoutDeltaNew[i]   = 0.;
                    
                    fdModel->fdEngine.upBarrierOld[i]       = -1;
                    fdModel->fdEngine.upPayoutOld[i]        = 0.;
                    fdModel->fdEngine.upPayoutDeltaOld[i]   = 0.;
                }    
                
                // no barrier if contConversion
                if (callBool[step+1] && convBool[step+1] && !cvb->hasContConversion() ){ 
                    if (callIsAdjForAccrued[step+1] && addCouponBool[step+1]) {
                        couponToAdd = coupon[step+1];
                    }
                    else {
                        couponToAdd = 0.;
                    }
                    
                    if( cvb->hasAddOnConvRatios() ) {
                        cvb->getConversionInfo(fdModel->TimePts.StepDates[step+1], callLevel[step+1], 
                                               &isConvertible, &conversionRatioNext, &conversionCashNext, false/*is not spot*/);
                    }
                        
                    fdModel->fdEngine.upBarrierOld[i]     = (callLevel[step+1] + couponToAdd- conversionCashNext)/conversionRatioNext;
                    fdModel->fdEngine.upPayoutOld[i]      = conversionCashNext - Maths::max(bondFloor[step+1], putLevel[step+1]);
                    fdModel->fdEngine.upPayoutDeltaOld[i] = conversionRatioNext;
                } 
                else {
                    fdModel->fdEngine.upBarrierOld[i]     = -1;
                    fdModel->fdEngine.upPayoutOld[i]      = 0.;
                    fdModel->fdEngine.upPayoutDeltaOld[i] = 0.;
                }    
            }
                
            if ( cvb->triggerReset == true ) {
                i = 1;
                if (callBool[step] == true && convBool[step] == true) {
                    fdModel->fdEngine.upBarrierNew[i]       = (callLevel[step] - convCash[step])/cvb->resetConversionRatio;
                    fdModel->fdEngine.upPayoutNew[i]        = convCash[step] - bondFloor[step];
                    fdModel->fdEngine.upPayoutDeltaNew[i]   = convRatio[step];
                }
                else {
                    fdModel->fdEngine.upBarrierNew[i]       = -1;
                    fdModel->fdEngine.upPayoutNew[i]        = 0.;
                    fdModel->fdEngine.upPayoutDeltaNew[i]   = 0.;
                    
                    fdModel->fdEngine.upBarrierOld[i]       = -1;
                    fdModel->fdEngine.upPayoutOld[i]        = 0.;
                    fdModel->fdEngine.upPayoutDeltaOld[i]   = 0.;
                }    
                
                if (callBool[step+1] == true && convBool[step+1] == true) {
                    if (callIsAdjForAccrued[step+1] == true && addCouponBool[step+1]) {
                        couponToAdd = coupon[step+1];
                    } 
                    else {
                        couponToAdd = 0.;
                    }
                    
                    fdModel->fdEngine.upBarrierOld[i]     = (callLevel[step+1] + couponToAdd- convCash[step+1])/cvb->resetConversionRatio;
                    fdModel->fdEngine.upPayoutOld[i]      = convCash[step+1] - Maths::max(bondFloor[step+1], putLevel[step+1]);
                    fdModel->fdEngine.upPayoutDeltaOld[i] = convRatio[step+1];
                    } 
                else {
                    fdModel->fdEngine.upBarrierOld[i]     = -1;
                    fdModel->fdEngine.upPayoutOld[i]      = 0.;
                    fdModel->fdEngine.upPayoutDeltaOld[i] = 0.;
                }    
            }
            ++i;
        }
        
        if (isOCB == true) {
            for (i=numCVBPriceArrays ; i<= pEnd ; ++i) {
                if (step < fdModel->TimePts.NumOfStep) {
                    if (callBool[step] == true && convBool[step] == true && ocbExerBool[step] == true) {
                        fdModel->fdEngine.upBarrierNew[i]     = (callLevel[step] - convCash[step])/conversionRatio;
                        fdModel->fdEngine.upPayoutNew[i]      = convCash[step] - ocbStrike[step];
                        fdModel->fdEngine.upPayoutDeltaNew[i] = convRatio[step];
                    } 
                    else {
                        fdModel->fdEngine.upBarrierNew[i]     = -1;
                        fdModel->fdEngine.upPayoutNew[i]      = 0.;
                        fdModel->fdEngine.upPayoutDeltaNew[i] = 0.;
                    }    
                        
                    if (callBool[step+1] == true && convBool[step+1] == true && ocbExerBool[step+1] == true) {
                        if (callIsAdjForAccrued[step+1] == true && addCouponBool[step+1]) {
                            couponToAdd = coupon[step+1];
                            }
                        else {
                            couponToAdd = 0.;
                        }
                        
                        fdModel->fdEngine.upBarrierOld[i]     = (callLevel[step+1] + couponToAdd- convCash[step+1])/conversionRatioNext;
                        fdModel->fdEngine.upPayoutOld[i]      = convCash[step+1] - ocbStrike[step+1];
                        fdModel->fdEngine.upPayoutDeltaOld[i] = convRatio[step+1];
                        } 
                    else {
                        fdModel->fdEngine.upBarrierOld[i]     = -1;
                        fdModel->fdEngine.upPayoutOld[i]      = 0.;
                        fdModel->fdEngine.upPayoutDeltaOld[i] = 0.;
                    }    
                }
            }
        }
    } else {
        throw ModelException("ConvBond1fProd::preCalcFD",
                             "Internal error - fdModel must exist in function preCalcFD");
    }
}

/** calculate barriers if needed */
void ConvBond1fProd::preCalcFDGeneric(int step, int idx, int pStart, int pEnd, 
                                      double const * const * optionArray,
                                      double const * const * optionArrayOld)
{
    double couponToAdd;
    int i;
    int layer;

    double conversionRatio;
    double conversionRatioNext;

    if ( !isStaticSpread() ) {
        
        for (i=pStart; i<= pEnd ; ++i) {
            if (step < genericFDModel->TimePts.NumOfStep) {
                genericFDModel->fdEngine.upBarrierNew[i]       = -1;
                genericFDModel->fdEngine.upPayoutNew[i]        = 0.;
                genericFDModel->fdEngine.upPayoutDeltaNew[i]   = 0.;
                
                genericFDModel->fdEngine.upBarrierOld[i]       = -1;
                genericFDModel->fdEngine.upPayoutOld[i]        = 0.;
                genericFDModel->fdEngine.upPayoutDeltaOld[i]   = 0.;
            }
        }
    }
    else {
        conversionRatio     = convRatio[step];
        conversionRatioNext = (step < genericFDModel->TimePts.NumOfStep)?convRatio[step+1]:0.0;
        
        for (layer=0 ; layer<= pEnd ; ++layer) {
            if( cvb->hasContConversion() ) {
                i = layer;
            }
            else {
                i = 0;
            }

            if (step < genericFDModel->TimePts.NumOfStep) {
                if (callBool[step] == true && convBool[step] == true && !cvb->hasContConversion() ) // no barrier if contConversion
                {
                    if (hasE2CLayer) {
                        genericFDModel->fdEngine.upBarrierNew[i]       = (callLevel[step] - convCash[step])/conversionRatio;
                        if (step == genericFDModel->TimePts.NumOfStep-1) {
                            genericFDModel->fdEngine.upPayoutNew[i]        = convCash[step] - bondFloor[step];
                        }
                        else {
                            genericFDModel->fdEngine.upPayoutNew[i]        = convCash[step] - 
                                optionArray[pEnd+2][genericFDModel->stockSteps];
                        }
                        genericFDModel->fdEngine.upPayoutDeltaNew[i]   = conversionRatio;
                    } 
                    else {
                        genericFDModel->fdEngine.upBarrierNew[i]       = (callLevel[step] - convCash[step])/conversionRatio;
                        genericFDModel->fdEngine.upPayoutNew[i]        = convCash[step] - bondFloor[step];
                        genericFDModel->fdEngine.upPayoutDeltaNew[i]   = conversionRatio;
                    }
                } 
                else {
                    genericFDModel->fdEngine.upBarrierNew[i]       = -1;
                    genericFDModel->fdEngine.upPayoutNew[i]        = 0.;
                    genericFDModel->fdEngine.upPayoutDeltaNew[i]   = 0.;
                    
                    genericFDModel->fdEngine.upBarrierOld[i]       = -1;
                    genericFDModel->fdEngine.upPayoutOld[i]        = 0.;
                    genericFDModel->fdEngine.upPayoutDeltaOld[i]   = 0.;
                }    
                
                if (callBool[step+1] == true && convBool[step+1] == true && !cvb->hasContConversion() ) // no barrier if contConversion
                {
                    if (callIsAdjForAccrued[step+1] == true && addCouponBool[step+1]) {
                        couponToAdd = coupon[step+1];
                    } 
                    else {
                        couponToAdd = 0.;
                    }
                    
                    if (hasE2CLayer) {
                        genericFDModel->fdEngine.upBarrierOld[i]     = (callLevel[step+1] + couponToAdd -
                                                                        convCash[step+1])/conversionRatioNext;
                        genericFDModel->fdEngine.upPayoutOld[i]      = convCash[step+1] - 
                            Maths::max(optionArrayOld[pEnd+2][genericFDModel->stockSteps], putLevel[step+1]);
                        genericFDModel->fdEngine.upPayoutDeltaOld[i] = conversionRatioNext;
                    } 
                    else {
                        genericFDModel->fdEngine.upBarrierOld[i]     = (callLevel[step+1] + couponToAdd -
                                                                        convCash[step+1])/conversionRatioNext;
                        genericFDModel->fdEngine.upPayoutOld[i]      = convCash[step+1] - Maths::max(bondFloor[step+1], putLevel[step+1]);
                        genericFDModel->fdEngine.upPayoutDeltaOld[i] = conversionRatioNext;
                    }
                } 
                else {
                    genericFDModel->fdEngine.upBarrierOld[i]     = -1;
                    genericFDModel->fdEngine.upPayoutOld[i]      = 0.;
                    genericFDModel->fdEngine.upPayoutDeltaOld[i] = 0.;
                }    
            }
            
            if ( cvb->triggerReset == true ) {
                i = 1;
                if (callBool[step] == true && convBool[step] == true) {
                    genericFDModel->fdEngine.upBarrierNew[i]       = (callLevel[step] - convCash[step])/cvb->resetConversionRatio;
                    genericFDModel->fdEngine.upPayoutNew[i]        = convCash[step] - bondFloor[step];
                    genericFDModel->fdEngine.upPayoutDeltaNew[i]   = convRatio[step];
                    
                } 
                else {
                    genericFDModel->fdEngine.upBarrierNew[i]       = -1;
                    genericFDModel->fdEngine.upPayoutNew[i]        = 0.;
                    genericFDModel->fdEngine.upPayoutDeltaNew[i]   = 0.;
                    
                    genericFDModel->fdEngine.upBarrierOld[i]       = -1;
                    genericFDModel->fdEngine.upPayoutOld[i]        = 0.;
                    genericFDModel->fdEngine.upPayoutDeltaOld[i]   = 0.;
                }    
                
                if (callBool[step+1] == true && convBool[step+1] == true) {
                    if (callIsAdjForAccrued[step+1] == true && addCouponBool[step+1]) {
                        couponToAdd = coupon[step+1];
                    } 
                    else {
                        couponToAdd = 0.;
                    }
                    
                    genericFDModel->fdEngine.upBarrierOld[i]     = (callLevel[step+1] + couponToAdd
                                                                    - convCash[step+1])/cvb->resetConversionRatio;
                    genericFDModel->fdEngine.upPayoutOld[i]      = convCash[step+1] - Maths::max(bondFloor[step+1], putLevel[step+1]);
                    genericFDModel->fdEngine.upPayoutDeltaOld[i] = convRatio[step+1];
                } 
                else {
                    genericFDModel->fdEngine.upBarrierOld[i]     = -1;
                    genericFDModel->fdEngine.upPayoutOld[i]      = 0.;
                    genericFDModel->fdEngine.upPayoutDeltaOld[i] = 0.;
                }    
            }
            ++i;
        }
        
        
        if (hasE2CLayer) {
            
            FirmAsset* firmAsset = 0;
            if ( FirmAsset::TYPE->isInstance(asset.get())) {
                CAsset*    ncAsset   = const_cast<CAsset*>(asset.get());
                firmAsset = dynamic_cast<FirmAsset*>(ncAsset);
            } 
            else {
                throw ModelException("FD1FE2C::PostSetup",
                                     "Underlying credit asset must be a FirmAsset");
            }

            // don't like explicit cast here - cleanup with an interface when finished 
            FDEngine1FGeneric* genericEngine = dynamic_cast<FDEngine1FGeneric*>(&genericFDModel->fdEngine);
            if (genericEngine) {
                double accruedToday = cvb->getAccruedAtDate(model1F->TimePts.StepDates[step]);
                double recovery;
                
                if (cvb->payAccruedUponDefault) {
                    recovery = cvb->useAssetRecovery?((faceValue + accruedToday) * cvb->recoveryPct):0.0;
                } 
                else {
                    recovery = cvb->useAssetRecovery?(faceValue * cvb->recoveryPct):0.0;
                }
                
                genericEngine->downBarrierNew[pEnd+1]     = firmAsset->getDefaultBarrier();
                genericEngine->downPayoutNew[pEnd+1]      = recovery;
                genericEngine->downPayoutDeltaNew[pEnd+1] = 0.0;
                genericEngine->downBarrierOld[pEnd+1]     = firmAsset->getDefaultBarrier();
                genericEngine->downPayoutOld[pEnd+1]      = recovery;
                genericEngine->downPayoutDeltaOld[pEnd+1] = 0.0;
            }
        }
        
        
        if (isOCB == true) {
            for (i=numCVBPriceArrays ; i<= pEnd ; ++i) {
                if (step < genericFDModel->TimePts.NumOfStep) {
                    if (callBool[step] == true && convBool[step] == true && ocbExerBool[step] == true) {
                        genericFDModel->fdEngine.upBarrierNew[i]     = (callLevel[step] - convCash[step])/conversionRatio;
                        genericFDModel->fdEngine.upPayoutNew[i]      = convCash[step] - ocbStrike[step];
                        genericFDModel->fdEngine.upPayoutDeltaNew[i] = convRatio[step];
                        
                    } 
                    else {
                        genericFDModel->fdEngine.upBarrierNew[i]     = -1;
                        genericFDModel->fdEngine.upPayoutNew[i]      = 0.;
                        genericFDModel->fdEngine.upPayoutDeltaNew[i] = 0.;
                    }    
    
                    if (callBool[step+1] == true && convBool[step+1] == true && ocbExerBool[step+1] == true) {
                        if (callIsAdjForAccrued[step+1] == true && addCouponBool[step+1]) {
                            couponToAdd = coupon[step+1];
                        }
                        else {
                            couponToAdd = 0.;
                        }
                        
                        genericFDModel->fdEngine.upBarrierOld[i]     = (callLevel[step+1] + couponToAdd
                                                                        - convCash[step+1])/conversionRatioNext;
                        genericFDModel->fdEngine.upPayoutOld[i]      = convCash[step+1] - ocbStrike[step+1];
                        genericFDModel->fdEngine.upPayoutDeltaOld[i] = convRatio[step+1];
                    } 
                    else {
                        genericFDModel->fdEngine.upBarrierOld[i]     = -1;
                        genericFDModel->fdEngine.upPayoutOld[i]      = 0.;
                        genericFDModel->fdEngine.upPayoutDeltaOld[i] = 0.;
                    }    
                }
            }
        }
    }
}



/** product payoff method at steps earlier than maturity */
void ConvBond1fProd::PayoffAtMat(const double * s, int step, int bot, int top, int pStart, int pEnd,
                                    double * const * price)
{
    int     i,j;
    double  strike = 0.;
    
    // allocate local memory for bondFloorGeneric. assume that tree/fd need max bond floor memory at maturity
    if( bondFloorGenericMem == 0 && (!isStaticSpread() || hasE2CLayer ) ) {
        bondFloorGenericMem = new double[bot + top + 1];
    }
    
    // if not static spread, need to initiate the risky bond prices and get the non-eq related spread curve
    if (!isStaticSpread()) {
        int i = step;
        double  bondRiskFree = 0;
        double  bondRisky = bondFloor[step];
        
        // ********** escrow code is suspicious since bondRiskFree is always 0 prior to this calculation!
        // check for stock being in Escrow account
        if (cvb->stockInEscrow) {
            double parity     = getEquityAsset()->getSpot() * convRatio[i];
            double difference = parity - bondRiskFree;
            bondRisky    -= difference;
        }

        // We assume here that puts cannot be before riskFreeCouponEndDate
        if ( putBool[i] == true && 
             !( model1F->TimePts.StepDates[i-1].getDate() == 
                model1F->TimePts.StepDates[model1F->TimePts.StepDates.size()-1].getDate()) &&
             !(cvb->DECS == true || cvb->PERCS == true) ) {
            bondRisky = Maths::max(bondRisky, putLevel[i]);
        } 
        
        if (addCouponBool[i] == true) {
            if ( !( model1F->TimePts.StepDates[i] <= cvb->riskFreeCouponEndDate ) ) {
                bondRisky += coupon[i];
            }
        }
        for (j=-bot; j<=top; j++) {
            price[pEnd][j] = bondRisky;
        }
        
        pEnd--;
    }
    
    if (hasE2CLayer) {
        int     j;
        // the price is simply the redemption value
        // price[pEnd+2][0] = price[pEnd+1][0] = faceValue * cvb->recoveryPct;
        price[pEnd+2][0] = price[pEnd+1][0] = Maths::max(cvb->bond->getRedemption(), putLevel[step]);
        for (j=-bot+1; j<=top; ++j) {
            if (putBool[step] == true) {
                price[pEnd+2][j] = price[pEnd+1][j] = Maths::max(cvb->bond->getRedemption(), putLevel[step]);
            } 
            else {
                if (callBool[step] && callIsHard[step]) {
                    price[pEnd+2][j] = price[pEnd+1][j] = Maths::min(cvb->bond->getRedemption(), callLevel[step]);
                } 
                else {
                    price[pEnd+2][j] = price[pEnd+1][j] = cvb->bond->getRedemption();
                }
            }
        }
    }

    
    // set numCVBPriceArrays if we're pricin on a tree - this should possibly be in preCalcTree, which
    // would require touching the interface in many payouts.
    if ( !genericFDModel && !fdModel) {
        if (isOCB) {
            numCVBPriceArrays = pEnd / 2;
        } 
        else {
            numCVBPriceArrays = pEnd;
        }
    }

    int endIdx = (isOCB && hasReset)?numCVBPriceArrays-1:pEnd;
    
    // when we put on the maturity date we always get the coupon
    double put = putLevel[step];
    if ( putBool[step] == true && cvb->putAdjustForAccrued == false ) {
        put += cvb->getAccruedAtDate(cvb->bond->getMaturityDate());
    }

    for (i=pStart; i<=endIdx; i++) {
        if (cvb->DECS == true) {
            for (j=-bot; j<=top; j++) {
                if (s[j] < cvb->initialPrice) {
                    price[i][j] = s[j]*cvb->initialConvRatio;
                } 
                else if (s[j] < cvb->convPrice){
                    price[i][j] = cvb->initialPrice * cvb->initialConvRatio;
                } 
                else if ( cvb->decsHasCutoff && s[j] > cvb->decsCutoffLevel ) {
                    price[i][j] = cvb->initialPrice * cvb->initialConvRatio + (cvb->decsCutoffLevel - cvb->convPrice)*cvb->minConvRatio;
                } 
                else {
                    price[i][j] = cvb->initialPrice * cvb->initialConvRatio + (s[j] - cvb->convPrice)*cvb->minConvRatio;
                }
                
            }        
        } 
        else if (cvb->PERCS == true) {
            for (j=-bot; j<=top; j++) {
                if (s[j] < cvb->initialPrice) {
                    price[i][j] = s[j]*cvb->initialConvRatio;
                } 
                else {
                    price[i][j] = cvb->initialPrice*cvb->initialConvRatio;
                } 
            }
        } 
        else {
            strike = bondFloor[step];
            double conversionRatio, conversionCash = convCash[step];
            for (j=-bot; j<=top; j++) {
                if (hasE2CLayer) {
                    strike = price[pEnd+2][j];
                }
                
                if (hasReset && (endIdx - pStart > 0)) {
                    // trigger no reset code to has here
                    conversionRatio = (*conversionRatioGrid)[step][i];
                }
                else {
                    if( convBool[step] && cvb->hasAddOnConvRatios() ) {
                        bool isConvertible;
                        cvb->getConversionInfo(model1F->TimePts.StepDates[step], s[j], &isConvertible, &conversionRatio, &conversionCash);
            } 
                    else {
                        conversionRatio = convRatio[step];
                    }
                }
                
                double convertValue = 0.0;

                if (cvb->triggerReset                                                                      &&
                    cvb->resetType == "UP_RESET"                                                           &&
                    model1F->TimePts.StepDates[step].getDate() >= cvb->resetObservationStartDate.getDate() &&
                    model1F->TimePts.StepDates[step].getDate() <= cvb->resetObservationEndDate.getDate()   &&
                    cvb->resetTrigger <=  s[j] ) {
                    conversionRatio = cvb->resetConversionRatio;
                } 
                else if ( 
                    cvb->triggerReset                                                                       &&
                    cvb->resetType == "DOWN_RESET"                                                          &&
                    model1F->TimePts.StepDates[step].getDate() >= cvb->resetObservationStartDate.getDate()  &&
                    model1F->TimePts.StepDates[step].getDate() <= cvb->resetObservationEndDate.getDate()    &&
                    cvb->resetTrigger >=  s[j]  ) {
                    conversionRatio = cvb->resetConversionRatio;
                }
                
                if( cvb->hasContConversion() && cvb->contConversion->isTrigActiveAtMat() &&
                    i==(pStart+cvb->contConversion->isHistEnabled(cvb->getValueDate())) ) {
                    convertValue = 0.0; // non-convertible price array
                } 
                else if ( cvb->contingentConversion && cvb->triggerActiveAtMat && 
                          conversionRatio*s[j] < cvb->contingentConversionTrigger ) {
                    convertValue = 0.0;
                } 
                else {
                    convertValue = conversionRatio*s[j]+conversionCash;
                }
                double holderValue = Maths::max(convertValue, put);
                
                double issuerValue = strike;
                
                price[i][j] = Maths::max((Maths::max(holderValue, issuerValue) - strike), 0.0);
                
                if ( cvb->triggerReset ) {
                    convertValue = cvb->resetConversionRatio*s[j]+convCash[step];
                    double holderValue = Maths::max(convertValue, put);
                    double issuerValue = bondFloor[step];
                    
                    price[i+1][j] = Maths::max((Maths::max(holderValue, issuerValue) - strike), 0.0);
                }
            }    
        }
    }

    // does not handle options on resettables yet
    if (isOCB && pEnd > pStart) { // we're an OCB
        for (i=numCVBPriceArrays ; i<=pEnd ; ++i) {
            for (j=-bot; j<=top; j++) {
                if (ocbExerBool[step] == true) {
                    price[i][j] = Maths::max(price[i-numCVBPriceArrays][j] + strike - ocbStrike[step], 0.);
                } 
                else {
                    price[i][j] = 0.;
                }
            }
        }
    }
}


/** product payoff method at steps earlier than maturity */
void ConvBond1fProd::PayoffBeforeMat(const double * s, int step, int bot, int top, int pStart, int pEnd,
                                     double * const * price)
{
    static const string method = "ConvBond1fProd::PayoffBeforeMat";
    int         j;
    int         layer;
    double      convertValue;
    bool        aboveCall;
    double      slope;
    double      conversionValue;
    double      conversionValueForCall; // distinguish between conversion value from convertibility and from right to convert when called
    DoubleArray offset;
    DoubleArray coefficient;
    double      conversionRatioForCall, conversionCash = convCash[step];
    double      conversionRatio = convRatio[step];

    int         layerCoCoCV=-1, layerCoCoNCV=-1, layerTrigResetRatio=-1;
    
    //////////////////////////////////////////////////////////////
    //                                                          //
    //  For debt equity model, need to calc spot dependent      //
    //  bond floor. also adjust pEnd if needed                  //
    //                                                          //
    //////////////////////////////////////////////////////////////
    
    bool        isBondFlrGeneric = ( !isStaticSpread() || hasE2CLayer );
    double      bondFloorTmp = bondFloor[step], *bondFloorGeneric=0;
    if( isBondFlrGeneric ) {
        pEnd = calcRiskyBondNAdjPriceLayer(step, bot, top, pEnd, price, bondFloorGeneric);
    }
    
    int endIdx = pStart;
    if ( hasReset || cvb->triggerReset || cvb->hasContConversion() ) {
        if ( step >= firstGridStep ) {
            endIdx = (isOCB)?numCVBPriceArrays-1:pEnd;
        } 
        else {
            // only need single grid if there is no reset before the current step
            endIdx = pStart;
        }
        
        if( cvb->triggerReset ) {
            layerTrigResetRatio = pStart + 1;
            if( endIdx != pStart + 1 ) {
                throw ModelException(method, "Internal error. 2 cvb pricing grids needed for trigger reset");
            }
        }
        else if ( cvb->hasContConversion() ) {   
            if( endIdx != pStart + 1 ) {
                throw ModelException(method, "Internal error. 2 cvb pricing grids needed for contingent conversion");
            }
            
            bool idxFlag = cvb->contConversion->isHistEnabled(cvb->getValueDate());
            layerCoCoCV = pStart + !idxFlag;
            layerCoCoNCV = pStart + idxFlag;
        }
    }
    
    // first do Penult adjustment if required
    if (step == model1F->TimePts.NumOfStep-1) {
        for (layer=pStart; layer<=endIdx; ++layer) {
            // disable smoothing for spot dependent spread
            if( isStaticSpread() && 
                !( cvb->hasContConversion() && cvb->contConversion->isTrigActiveAtMat() && layer==layerCoCoNCV ) && 
                // disable smoothing for no-convert price slice if contingent conversion and trigger active at maturity
                ( !adjustedResetSchedule || adjustedResetSchedule->length() == 0) ) {
                // no penultimate step smoothing for resettables yet, I'm afraid
                if (cvb->DECS == false && cvb->PERCS == false)
                    PenultStraight(s, step, bot, top, price[layer]);
                else
                    PenultMandatory(s, step, bot, top, price[layer]);
            }
        }
    }

    // calculate relevant stock interval for Black on call and some coefficients
    // initial stock interval to be not-do-Black-on-call-adjustments
    double bocLowerStockPrice = -2, bocUpperStockPrice = -1;
    if ( cvb->callTreatment == "B" && callBool[step] )
        preCalcBlackOnCall(step, s[(top-bot)/2], bocLowerStockPrice, bocUpperStockPrice, offset, coefficient);

    
    //////////////////////////////////////////////////////////////
    //                                                          //
    //  Merge price layers of diff conv ratios to get cvb       //
    //  option values for conversion/call/put decision making   //
    //  Needed for: reset, triggerReset                         //
    //  Ideally should do CoCo (done at end of function)        //
    //                                                          //
    //////////////////////////////////////////////////////////////
    if ( hasReset && (pEnd - pStart) > 0 && (*resetArray)[step] ) {
        mergePriceForReset(s, step, bot, top, pStart, endIdx, price);
    }
    
    if (cvb->triggerReset &&
        model1F->TimePts.StepDates[step].getDate() >= cvb->resetObservationStartDate.getDate() &&
        model1F->TimePts.StepDates[step].getDate() <= cvb->resetObservationEndDate.getDate() ) {
        bool isUp = (cvb->resetType == "UP_RESET"); 
        for (j=-bot; j<=top; ++j) {
            if( isUp == (cvb->resetTrigger <=  s[j]) ) price[pStart][j] = price[layerTrigResetRatio][j];
        }
    }
     
    
    //////////////////////////////////////////////////////////////
    //                                                          //
    //  Pricing for each layer                                  //
    //                                                          //
    //////////////////////////////////////////////////////////////

    // normally, FD calculation is up to callLevel and need to linear extend beyond
    // for certain cases, no upBarrier is used and should not use slope due to potential inaccuracy/instability
    // therefore, ideally should be in-synch with boundary treatment in preCalcFD()
    bool useSlopeAboveCall = true;
    if ( hasE2CLayer || cvb->hasContConversion() ) useSlopeAboveCall = false;

    for (layer=pStart; layer<=endIdx; ++layer) {
        
        aboveCall = false;
        slope   = 0.0;

        // initialize for each loop over spots. used be addOnConvRatios
        double nextAddOnConvTrigger = cvb->hasAddOnConvRatios()?(-1):(s[top]*1.1); 
        
        for (j=-bot; j<=top; ++j) {
            
            // if not static spread, we replace the bond floor by the appropriate bond floor per node
            if ( isBondFlrGeneric ) {
                bondFloorTmp = bondFloorGeneric[j];
            }

            convertValue = price[layer][j] + bondFloorTmp;
            
            //////////////////////////////////////////////////////////////
            //                                                          //
            //  I: conversionValue and conversionValue/Ratio            //
            //  in case of call                                         //
            //                                                          //
            //////////////////////////////////////////////////////////////
            
            if ( !convBool[step] ) {
                conversionValue = 0.0;                
                conversionRatioForCall = 0.0;
                conversionValueForCall = 0.0;
            } 
            // if DEC but not continuousReset, use convRatio[step] (see following logic)
            else if ( cvb->DECS && cvb->continuousReset ) {
                if (s[j] < cvb->initialPrice)
                    conversionValue = s[j] * cvb->initialConvRatio;
                else if (s[j] < cvb->convPrice)
                    conversionValue = cvb->initialPrice * cvb->initialConvRatio;
                else if( !cvb->decsHasCutoff || s[j] < cvb->decsCutoffLevel)
                    conversionValue = cvb->initialPrice * cvb->initialConvRatio + (s[j] - cvb->convPrice)*cvb->minConvRatio;
                else
                    conversionValue = cvb->initialPrice * cvb->initialConvRatio + (cvb->decsCutoffLevel - cvb->convPrice)*cvb->minConvRatio;
                
                // calculate conversion option
                if ( cvb->delayReset > 0 )
                    conversionValue = adjConversionValueDECDelay(step, s[j], convertValue, conversionValue);
                
                conversionValue += conversionCash;
                
                conversionRatioForCall = (conversionValue - conversionCash)/s[j];
                conversionValueForCall = conversionValue;
                
            } 
            else if ( hasReset && (pEnd - pStart > 0) ) {
                // currently using linearly interpolated conversion ratios
                conversionRatio = (*conversionRatioGrid)[step][layer];
                conversionValue = conversionRatio*s[j]+conversionCash;
                
                if ( cvb->maxWithParity == true ) {
                    double parityCR, parityValue;
                    double  minStrike   = (*minConvPrice)[step];
                    if ( s[j] * cvb->manualResetConvPricePct < minStrike) {
                        parityCR = faceValue / minStrike;
                    } 
                    else {
                        parityCR = faceValue / (s[j] * cvb->manualResetConvPricePct);
                    }
                    parityValue = parityCR * s[j];
                    conversionValue = Maths::max(parityValue,conversionValue);
                }
                
                conversionRatioForCall = conversionRatio;
                conversionValueForCall = conversionValue;
            } 
            else if( cvb->triggerReset && layer==layerTrigResetRatio )  {
                conversionValue = cvb->resetConversionRatio*s[j]+conversionCash;
                conversionRatioForCall = cvb->resetConversionRatio;
                conversionValueForCall = conversionValue;
            }
            else {
                // if has addOnConvRatio, update conversionRatio/Cash across addOnConvRatio
                // trigger level. we know that s[j] >0 and isConvertible
                if( cvb->hasAddOnConvRatios() && s[j] >= nextAddOnConvTrigger ) {
                    bool isConvertible;
                    cvb->getConversionInfo(model1F->TimePts.StepDates[step], s[j], &isConvertible, &conversionRatio, &conversionCash);
                    nextAddOnConvTrigger = cvb->addOnConvRatios->getNextTrigger(s[j], s[top] * 1.1);
                }
                conversionValue = conversionRatio*s[j]+conversionCash;
                conversionRatioForCall = conversionRatio;
                conversionValueForCall = conversionValue;
                
                // apply contingent conversion trigger
                // does not affect conversionValue/Ratio for call purpose
                if( (cvb->hasContConversion() && layer==layerCoCoNCV ) ||
                    (cvb->contingentConversion && conversionValue < cvb->contingentConversionTrigger) )
                {
                    conversionValue = 0.0; 
                }
            }
            
            //////////////////////////////////////////////////////////////
            //                                                          //
            //  II: call level and DEC's call-like features             //
            //                                                          //
            //////////////////////////////////////////////////////////////
            
            if (callBool[step]) {
                if (convBeforeCall[step]) 
                {
                    convertValue = Maths::min(convertValue, callLevel[step]);
                    // first make sure convertValue is below parityFloor
                    if( convertValue < parityFloor[step] * conversionValue )
                        convertValue = Maths::max(conversionValue, convertValue);
                } 
                else 
                {
                    // convert value now given conversionValue and convertValue-if-not-conversion-now
                    // first make sure convertValue is below parityFloor
                    if( convertValue < parityFloor[step] * conversionValue )
                        convertValue = Maths::max(conversionValue, convertValue);
                    
                    // notice that convertValue may be less than conversionValueForCall, eg for certain ContConversion
                    // notice that (p[layer])[j-1] is already sum of (bondfloor + option)
                    if (useSlopeAboveCall && (!aboveCall) && !Maths::isNegative( convertValue - conversionValueForCall ) &&
                        ((cvb->notResetCallTrigger? cvb->conversionRatios->firstValue()*
                          s[j]:conversionValueForCall) >= callLevel[step]) && s[j] >= bocUpperStockPrice ) {
                        aboveCall = true;
                        if (j > -bot)
                            slope = (callLevel[step] - price[layer][j-1])/((callLevel[step]-conversionCash)/conversionRatioForCall-s[j-1]);
                        else
                            slope = conversionRatioForCall;
                    }
                    
                    if (!aboveCall) 
                    { 
                        // still need to min with call
                        convertValue = Maths::min(convertValue, callLevel[step]); 
                    }
                    else 
                    {
                        convertValue = callLevel[step] + slope*(s[j] - (callLevel[step]-conversionCash)/conversionRatioForCall);
                    }
                }
                
            }
            // no calls. but may have call-like feature for DECS
            else if ( cvb->DECS && cvb->accelerateDECS && model1F->TimePts.StepDates[step] >= cvb->accelerateDECSStartDate && 
                    model1F->TimePts.StepDates[step] <= cvb->accelerateDECSEndDate && s[j] >= cvb->accelerateDECSTriggerLevel )
            {
                double callValue = calcAccelerateDECCallValue(step, s[j], bondFloorTmp);
                convertValue = Maths::max(conversionValue, Maths::min(convertValue, callValue));
            }   
            else if ( cvb->DECS && cvb->decsHasCutoff && !cvb->triggerStartDate.empty() && model1F->TimePts.StepDates[step] >= cvb->triggerStartDate &&
                s[j] > cvb->decsCutoffLevel * cvb->cappedDECSTrigger )
            {
                convertValue = cvb->initialPrice * cvb->initialConvRatio + (cvb->decsCutoffLevel - cvb->convPrice)*cvb->minConvRatio;
            } 
            else 
            {
                // first make sure convertValue is below parityFloor
                if( convertValue < parityFloor[step] * conversionValue )
                    convertValue = Maths::max(conversionValue, convertValue);
            }
            
            // add soft put condition
            if ( softPutBool[step] && conversionValue <= softPutTrigger[step] ) {
                convertValue = Maths::max(convertValue, softPutLevel[step]);
            }
            
            //////////////////////////////////////////////////////////////
            // Notice that price array holds (convert = floor + option) //
            //////////////////////////////////////////////////////////////
            price[layer][j] = convertValue;
        }  
    }
    
    
    //////////////////////////////////////////////////////////////
    //                                                          //
    //  III: put strike/prob calc and adj price for put level   //
    //                                                          //
    //////////////////////////////////////////////////////////////
    
    // if put to stock, need to calc put value, which does not affect bond floor calculation.
    double putValue = putLevel[step];
    if( putBool[step] == true && !!cvb->putToStock ) putValue = calcPutToStockValue(step);

    // check if we want to calculate the put stock levels and put probs here
    // we should really check the control pointer here to avoid unnecessarily duplicating work, but this requires
    // some changes to callbacks - still need to be done.
    // Don't calc if puts a  re American style i.e. if you can put at next step.
    // if contingent conversion, the put level/prob is not unique. so we just use the cv/ncv price depending on whether today is cv or ncv      
    if ((pStart == pEnd || cvb->hasContConversion() ) && 
        firstPutIndex == step && !Maths::isPositive(putStockLevel) && putBool[step+1] == false)
        calcPutLevelNProb(step, bot, top, s, price[pStart], putValue, putStockLevel, putProbability);
    
    // max with put - this has to be done after the put probability has been calculated
    // also adjust the bond floor by put level
    int i;
    if (putBool[step] == true ) {
        bondFloorTmp = Maths::max(bondFloorTmp, putLevel[step]);
        if ( isBondFlrGeneric )
        {
            for (j=-bot; j<=top; ++j)
                bondFloorGeneric[j] = Maths::max(bondFloorGeneric[j], putLevel[step]);
        }
        
        for (i=endIdx; i>=pStart; --i) {
            for (j=-bot; j<=top; ++j)
            {
                price[i][j] = Maths::max(price[i][j], putValue);
            }
        }
    }
    
    // adjust pricing if contingent conversion, not done at start of function because
    // it may use conversion info from the previous EOD and currently we can not place
    // time line points at EOD to avoid conflict with SOD time points from conversion schedule
    // (closeby points may be erased during timeline construction and causing problem)
    if( cvb->hasContConversion() )
    {
        cvb->contConversion->adjustPrices(model1F->TimePts.StepDates[step], top + bot, 
            &s[-bot], &price[layerCoCoNCV][-bot], &price[layerCoCoCV][-bot], callBool[step]?callLevel[step]:(-1), cvb);
    }
    
    
    //////////////////////////////////////////////////////////////
    //                                                          //
    //  IV: OCB                                                 //
    //                                                          //
    //////////////////////////////////////////////////////////////
    
    if (isOCB && pEnd > pStart) { // we're an OCB
        for (i=numCVBPriceArrays ; i<= pEnd ; ++i) {
            for (j=-bot; j<=top; j++) {
                price[i][j] += ocbRiskyAdj[step];
            }
            
            if (ocbExerBool[step] == true) {
                for (j=-bot; j<=top; j++)
                {    
                    price[i][j] = Maths::max((price[i-numCVBPriceArrays][j] - ocbStrike[step]), price[i][j]);
                }
                
                // be careful -- we don't want to add value for bonds that have already been called so min with call 
                // unless we're going to convert just before. Also be careful not to subtract value when the strike 
                // is above the call level. Avoid this by never MINing with a negative number.
                if (convBeforeCall[step] == false) {
                    double callMinusStrike = Maths::max(callLevel[step] - ocbStrike[step], 0.);
                    double conversionRatio;
                    if (hasReset) {
                        conversionRatio = (*conversionRatioGrid)[step][i-numCVBPriceArrays];
                    } else {
                        conversionRatio = convRatio[step];
                    }
                    
                    for (j=-bot; j<=top; j++) {   
                        if (conversionRatio*s[j]+conversionCash >= callLevel[step]) {
                            price[i][j] = Maths::min(price[i][j], callMinusStrike);
                        }
                    }
                }
            }
        }
    }
    
    //////////////////////////////////////////////////////////////
    //                                                          //
    //  IV: Remove bond floor                                   //
    //                                                          //
    //////////////////////////////////////////////////////////////
    
    for (i=pStart; i<=endIdx; ++i) {
        for (j=-bot; j<=top; ++j)
        {
            if( isBondFlrGeneric ) bondFloorTmp = bondFloorGeneric[j];
            price[i][j] -= bondFloorTmp;
        }
    }

    // for Fwd Starting, addint the coupon before starting date but after issue date.
    if (step ==0 && fwdStarting){
        for (i=pStart; i<=endIdx; ++i) {
            for (j=-bot; j<=top; ++j)
                price[i][j] += pvCFbeforeStart;
        }
    }
}


void ConvBond1fProd::mergePriceForReset(const double* s, int step, int bot, int top, int pStart, int pEnd,
                                        double * const * price)
{
    if ( priceCopy.numCols() == 0 ) {
        priceCopy = CDoubleMatrix(pEnd+1, top-bot+1);
    }

    // back up the input prices 
    int layer, j;
    for(layer= 0;layer<=pEnd;++layer) {
        for(j=-bot; j<=top; ++j) {
            priceCopy[layer][j] = price[layer][j];
        }
    }
    
    double  maxStrike   = (*maxConvPrice)[step];
    double  minStrike   = (*minConvPrice)[step];
    double  resetPremium  = (*resetParity)[step];
    
    for (layer=pStart; layer<=pEnd; ++layer) {
        
        double conversionRatio = (*conversionRatioGrid)[step][layer];
        
        for (j=-bot; j<=top; ++j) {
            
            // need to cross over the grids depending on current spot here
            int lowIdx  = 0;
            int highIdx = 0;
            double targetCR;

            double resetRatio;

            if (step == 0) {
                // valueDate = resetDate so we know what the conversion ratio will reset to.
                // We should just use the reset layer(s) corresponding to this
                // conversion ratio and not those implied by grid spot levels.
                resetRatio = adjustedResetSchedule->getCurrentConversionRatio(
                                        convRatio[step],
                                        model1F->TimePts.StepDates[step],
                                        faceValue);
            } else {
                resetRatio = ResetSchedule::getResetRatio(faceValue,
                                                          maxStrike,
                                                          minStrike,
                                                          resetPremium,
                                                          s[j]);
            }
            
            string resetTypeStr = adjustedResetSchedule->getResetType(model1F->TimePts.StepDates[step]);
            if (  resetTypeStr == "UP_DOWN" ) {
                // UP_DOWN
                targetCR   = resetRatio;
            } else {
                if (adjustedResetSchedule->isFlat(model1F->TimePts.StepDates[step])) {
                    targetCR   = Maths::max(Maths::max(resetRatio,conversionRatio),convRatio[step]);
                } else {
                    targetCR   = Maths::max(resetRatio,conversionRatio);
                }
            }
            
            
            // new way of setting up conversion ratios
            double cr;
            int    idx = 0;
            double lowCR, highCR;
            while (idx < pEnd ) {
                cr = (*conversionRatioGrid)[step][idx];
                if (targetCR >= cr) {
                    lowIdx = idx;
                    lowCR  = cr;
                    if ( targetCR == cr ) {
                        highIdx = idx;
                        highCR  = cr;
                    } else {
                        highIdx = idx+1;
                        highCR  = (*conversionRatioGrid)[step][highIdx];
                    }
                }
                ++idx;
            }
            
            double weight = (highCR != lowCR)?((targetCR - lowCR) / ( highCR - lowCR)):1.0;
            price[layer][j]  = weight * priceCopy[highIdx][j] + (1.-weight) * priceCopy[lowIdx][j];

            if (adjustedResetSchedule->getReceiveIntrinsic(
                                            model1F->TimePts.StepDates[step])) {
                // Adjust option prices for any intrinsic value
                if (step > 0) {
                    // If intrinsic has not yet been paid add intrinsic to the option value
                    // where intrinsic is based on the conversion ratio prior to reset
                    double lastStepCR = (*conversionRatioGrid)[step-1][layer];
                    double parity = s[j] * lastStepCR;
                    double intrinsic = Maths::max(parity-faceValue,0.0);
                    price[layer][j] += intrinsic;
                }
            }
        }
    }
}


// the input price holds (bondfloor + option)
void ConvBond1fProd::calcPutLevelNProb(int step, int bot, int top, const double *s, 
                                       const double *priceLayer, double putValue, double &putStockLevel, double &putProbability)
{
    static const string method = "ConvBond1fProd::calcPutLevelNProb";
    
    // don't fail if this doesn't work
    try {
        putStockLevel  = 0.;
        putProbability = 1.0;
        if ( isStaticSpread() && putValue < bondFloor[step]) {
            putProbability = 0.0;
        } else if (Maths::areEqualWithinTol(callLevel[step],putValue,1.e-11)) {
            if( convBool[step] )
            {
                if( cvb->hasAddOnConvRatios() )
                    putStockLevel = cvb->getImpliedStrike(model1F->TimePts.StepDates[step], putValue);
                else
                    putStockLevel = (putValue-convCash[step])/convRatio[step];
            }
        } else if (isStaticSpread() && putValue > priceLayer[top]) {
            putProbability = 1.0;
        } else if (isStaticSpread() && putValue <= priceLayer[-bot] ) {
            putProbability = 0.0;
        } else {
            for (int idx=top-1; idx >= -bot; --idx) {
                
                if ( priceLayer[idx] <= putValue ) {
                    if (!Maths::isZero((priceLayer[idx+1] - priceLayer[idx]))) {
                        double a;
                        a = (priceLayer[idx+1] - putValue)/(priceLayer[idx+1] - priceLayer[idx]);
                        putStockLevel = (1-a)*s[idx+1] + a*s[idx];
                    } else {
                        putStockLevel = s[idx+1];
                    }
                    break;
                }
            }
        }
        
        if (putStockLevel > 0.) {
            // create vol request and get processed vol
            CVolRequestConstSP volRequest = GetLNRequest();
            CVolProcessedSP  procVol(getEquityAsset()->getProcessedVol(volRequest.get()));
            // cast to the type of vol we're expecting
            CVolProcessedBSSP volBS = CVolProcessedBSSP::dynamicCast(procVol);
            // calculate the variance
            double variance = volBS->CalcVar(cvb->getValueDate(),
                model1F->TimePts.StepDates[step]);
            
            // get the forward at maturity - could possibly pull this one out the the FWD_AT_MAT
            // results if some of the calculations are done in recordOutputRequests
            double fwd = getEquityAsset()->fwdValue(model1F->TimePts.StepDates[step]);
            if ( !Maths::isZero(putStockLevel)) {
                putProbability = 1. - N1((log(fwd/putStockLevel) - (0.5 * variance))/sqrt(variance));
            } else {
                putProbability = 0.0;
            }
        }
        
        calcedFirstPutProb = true;
    } 
    catch (exception& e) {
        putStockLevel = 999.;
        ModelException myE = ModelException(e);
        myE.addMsg(method + ": Failed to compute put stock level. Continuing...\n");
        myE.errorLog();
    }
}

// compute value of put if is put to stock but not already put. 
// code adapted from CorpAct.cpp
double ConvBond1fProd::calcPutToStockValue(int step)
{ 
    DateTime putDate = model1F->TimePts.StepDates[step];
    SampleListSP samples = cvb->getPutToStockSample(putDate);

    // make sure the first date is past value date, otherwise, can not put!
    if( samples->getFirstDate().getDate() < cvb->getValueDate().getDate() ) return 0.0;

    double factor = cvb->calcPutToStockFactor(putDate, samples.get(), GetLNRequest().get());
    return putLevel[step] * factor;
}

double ConvBond1fProd::adjConversionValueDECDelay(int step, double s, double convertValue, double conversionValue)
{
    // price the call
    double convRatio = conversionValue / s;
    DateTime fwdDate   = model1F->TimePts.StepDates[step].rollDate(cvb->delayReset);
    double   discFact  = cvb->discount->pv(model1F->TimePts.StepDates[step], fwdDate);
    double   strike    = convertValue / convRatio;
    double   fwd       = getEquityAsset()->fwdFwd(model1F->TimePts.StepDates[step],s,fwdDate);
    double dummy = s;
    vector<double> fwdVol;
    model1F->GetStepVol(step, fwdVol, &dummy,0,0);
    double variance = fwdVol[0] * fwdVol[0] * (double(cvb->delayReset) / 365.25);
    double convOption = conversionOption(true, model1F->TimePts.StepDates[step], fwdDate, fwd, strike, variance, discFact);
    convOption *= convRatio;
    conversionValue += convOption;
    
    // price the put spread if necessary
    if ( cvb->includePutSpread && convRatio * s > convertValue )
    {
        double highStrike = s;
        double lowStrike   = convertValue / convRatio;
        double put1 = conversionOption(false /*put*/, model1F->TimePts.StepDates[step], fwdDate, fwd, highStrike, variance, discFact);
        double put2 = conversionOption(false /*put*/, model1F->TimePts.StepDates[step], fwdDate, fwd, lowStrike, variance, discFact);
        
        conversionValue += (put2 - put1) * convRatio;
    }
    
    return conversionValue;
}


double ConvBond1fProd::calcAccelerateDECCallValue(int step, double s, double bondFloor)
{
    // accelerate maturity feature
    double callValue = bondFloor;
    
    if ( cvb->dividendPassThrough ) {
        // subtract the pv of the dividends
        CashFlowArraySP dividendCashFlows(new CashFlowArray(0));
        
        DividendListSP divSchedule = DividendCollector::divsBetweenDates(getEquityAsset().get(),
            model1F->TimePts.StepDates[step],
            model1F->TimePts.StepDates[step],
            cvb->bond->getMaturityDate());
        
        const DividendArray& divArray = divSchedule->getArray();
        
        int k;
        for (k=0 ; k<divArray.size() ; ++k) {
            callValue -= divArray[k].getDivAmount() * cvb->initialConvRatio * cvb->dividendPassThroughPct * 
                cvb->discount->pv(model1F->TimePts.StepDates[step], 
                divArray[k ].getPayDate());
            
        }
    }
    
    if (s < cvb->initialPrice) {
        callValue += s * cvb->initialConvRatio;
    } else if (s < cvb->convPrice){
        callValue += cvb->initialPrice * cvb->initialConvRatio;
    } else {
        callValue += cvb->initialPrice * cvb->initialConvRatio + (s - cvb->convPrice)*cvb->minConvRatio;
    }
    
    return callValue;
}


// do some precomputing to speed up later calculations.
void ConvBond1fProd::preCalcBlackOnCall(int step, double midGrid, double &bocLowerStockPrice, double &bocUpperStockPrice, 
                                        DoubleArray &offset, DoubleArray &coefficient)
{
    static const string method = "ConvBond1fProd::preCalcBlackOnCall";
    
    DateTime callRedemptionDate = model1F->TimePts.StepDates[step].rollDate(cvb->callNotification);
    
    FD1FLNGeneric* lnModel;
    lnModel = dynamic_cast<FD1FLNGeneric*>(genericFDModel);
    
    if ( !lnModel ) {
        throw ModelException(method, "Only log-normal supported at the moment");
    }
    
    CVolProcessedBSConstSP processedVol = lnModel->getProcessedVol();
    
    // calculate the variance
    double variance = processedVol->CalcVar(model1F->TimePts.StepDates[step],
        callRedemptionDate);
    
    // double call stock price 
    double stockPrice  = (callLevel[step] - convCash[step])/convRatio[step];
    if( convBool[step] && cvb->hasAddOnConvRatios() )
        stockPrice = cvb->getImpliedStrike(model1F->TimePts.StepDates[step], callLevel[step]);
    
    // get forward price - currently set to stock price for performance
    double forward = stockPrice;
    
    // calculate upper and lower boundary (I am ignoring the vol correction term for simplicity for now)
    double threeSigmaExp = exp(sqrt(variance)*3);
    bocLowerStockPrice = forward / threeSigmaExp;
    bocUpperStockPrice = forward * threeSigmaExp;
    
    // calculate offset and slope for fast linear forwards
    int i;
    DateTimeArray fwdDates(cvb->callNotification);
    for (i=0 ; i<fwdDates.size() ; ++i) {
        fwdDates[i] = model1F->TimePts.StepDates[step].rollDate(i+1);
    }
    DoubleArray baseFwdPrices(fwdDates.size());
    DoubleArray tweakedFwdPrices(fwdDates.size());
    offset      = DoubleArray(fwdDates.size());
    coefficient = DoubleArray(fwdDates.size());
    
    getEquityAsset()->fwdFwd(model1F->TimePts.StepDates[step],midGrid,fwdDates,baseFwdPrices);
    getEquityAsset()->fwdFwd(model1F->TimePts.StepDates[step],midGrid  * 1.01 ,fwdDates,tweakedFwdPrices);
    
    for (i=0 ; i<fwdDates.size() ; ++i) {
        coefficient[i] = (tweakedFwdPrices[i] - baseFwdPrices[i])/(midGrid * .01);
        offset[i]      = baseFwdPrices[i] - midGrid * coefficient[i];
    }
}


/** Penultimate step smoothing. Should do a better job. At a minimum it should check for the call level.
    Could do an up and out call with a rebate. The old model did a call spread which is probably
    good enough for most purposes. Was probably overkill there as it would smooth the calls anyway. */
void ConvBond1fProd::PenultStraight(const double* s, int step, int bot, int top, double *price)
{
    int j;
    double fwdPrice;
    double strike = bondFloor[step+1]-convCash[step+1];
    vector <double> vol_arr, drift_arr;
    double dt = model1F->TimePts.TradeYrFrac[step+1];
    int nn = model1F->GetStepVol(step, vol_arr, s, -bot, top);
    if (nn > 1){
        return;
//        throw ModelException("ConvBond1fProd::PenultStraight", "Penultimate step smoothing not implemented for local vol");
    }
    double variance = vol_arr[0]*vol_arr[0]*dt;
    double df = cvb->discount->pv(model1F->TimePts.StepDates[step],
                                  model1F->TimePts.StepDates[step+1]);

    double drift = getEquityAsset()->fwdValue(model1F->TimePts.StepDates[step+1])/
                   getEquityAsset()->fwdValue(model1F->TimePts.StepDates[step]);

    double conversionRatio = convRatio[step+1], conversionCash = convCash[step+1];

    for (j=-bot; j<=top; j++)
    {
        if (s[j]>0.0) {
            fwdPrice = s[j]*drift;
            
            if( convBool[step] && cvb->hasAddOnConvRatios() )
            {
                bool isConvertible;
                cvb->getConversionInfo(model1F->TimePts.StepDates[step+1], fwdPrice, &isConvertible, &conversionRatio, &conversionCash);
                strike =  bondFloor[step+1]-conversionCash;
            }
            if (fabs(log(fwdPrice*conversionRatio/strike)) < model1F->TruncationStd*sqrt(variance))
            {
                price[j] = conversionRatio*Black::price(true, fwdPrice, 
                    strike/conversionRatio, df, variance);
            }
        }
    }
}

/** Penultimate step smoothing for mandatory convertibles. This is a distinct routine
    from the straight routine even though there is a lot of overlap because I anticipate 
    improving the straight routine. This one should be fine for practical purposes
    as it is. */
void ConvBond1fProd::PenultMandatory(const double* s, int step, int bot, int top, double *price)
{
    int j;
    double fwdPrice;

    vector <double> vol_arr, drift_arr;
    double dt = model1F->TimePts.TradeYrFrac[step+1];
    int nn = model1F->GetStepVol(step, vol_arr, s, -bot, top);
    /*
    if (nn > 1){
        throw ModelException("ConvBond1fProd::PenultMandatory", "Penultimate step smoothing not implemented for local vol");
    }
    */
    double variance = vol_arr[0]*vol_arr[0]*dt;
    double df = cvb->discount->pv(model1F->TimePts.StepDates[step],
                                  model1F->TimePts.StepDates[step+1]);

    double drift = getEquityAsset()->fwdValue(model1F->TimePts.StepDates[step+1])/
                   getEquityAsset()->fwdValue(model1F->TimePts.StepDates[step]);
    for (j=-bot; j<=top; j++)
    {
        if (s[j]>0.0) {
            fwdPrice = s[j]*drift;

            if (nn > 1){
                variance = vol_arr[j+bot]*vol_arr[j+bot]*dt;
            }


            if (cvb->DECS == true) {
                if (fabs(log(fwdPrice/cvb->initialPrice)) < model1F->TruncationStd*sqrt(variance) ||
                    fabs(log(fwdPrice/cvb->convPrice)) < model1F->TruncationStd*sqrt(variance))
                {
                    price[j] = cvb->initialConvRatio*fwdPrice*df;

                    price[j] -= cvb->initialConvRatio*Black::price(true, fwdPrice, 
                        cvb->initialPrice, df, variance);

                    price[j] += cvb->minConvRatio*Black::price(true, fwdPrice, 
                        cvb->convPrice, df, variance);
                }
            } else { // cvb->PERCS == true
                if (fabs(log(fwdPrice/cvb->initialPrice)) < model1F->TruncationStd*sqrt(variance))
                {
                    price[j] = cvb->initialConvRatio*fwdPrice*df;

                    price[j] -= cvb->initialConvRatio*Black::price(true, fwdPrice, 
                        cvb->initialPrice, df, variance);
                }
            }
        }
    }
}

double ConvBond1fProd::getCoupon(int step, const double* s, int start, int end)
{
    // this needs to get implemented for the E2C model
    double accruedToday = (cvb->payAccruedUponDefault)?cvb->getAccruedAtDate(model1F->TimePts.StepDates[step]):0.0;
    double recovery     = cvb->useAssetRecovery?(faceValue * cvb->recoveryPct):0.0;

    return recovery + accruedToday * cvb->recoveryPct;
}

bool ConvBond1fProd::hasEquityLayer()
{
    return true;
}

/** premium scaling */
double ConvBond1fProd::scalePremium(const double& fairValue)
{
    if (isOCB == false) {
        // return fairValue + bondFloor[0];
        double bondFloorTmp;
        if( isStaticSpread() )
            if (hasE2CLayer)
                bondFloorTmp = model1F->PriceEnd[model1F->PriceEnd.size()-1];
            else
                bondFloorTmp = bondFloor[0];
        else
        {
            int numPriceArray = (isOCB?2:1)*numCVBPriceArrays; // # if price array excl of the riskybond
            bondFloorTmp = model1F->PriceEnd[numPriceArray];
        }
        return fairValue + Maths::max(bondFloorTmp, putLevel[0]);
    } else {
        // return model1F->PriceEnd[1];
        return model1F->PriceEnd[numCVBPriceArrays];
    }
}

/** extra output requests */
void ConvBond1fProd::recordOutputRequests(Control* control, Results* results, double fairValue)
{
    static const string method = "ConvBond1fProd::recordOutputRequests";
    OutputRequest* request = NULL;

    if ( control && control->isPricing() ) {
        double nakedBond = 0.0;
        if (hasE2CLayer) {
             nakedBond = model1F->PriceEnd[genericFDModel->numPriceArrays-1];
        } else {
            nakedBond = bondFloor[0]+pvCFbeforeStart;
            if (putBool[0] == true) {
                nakedBond = Maths::max(nakedBond, putLevel[0]);
            }
        }
        cvb->recordOutputRequests(control, results, fairValue, nakedBond, isOCB);

        // record the bond floor as calculated by the no static spread
        if( !isStaticSpread() &&
            control->requestsOutput(OutputRequest::NAKED_BOND_PRICE2, request) )
        {
            int numPriceArray = (isOCB?2:1)*numCVBPriceArrays; // # if price array excl of the riskybond
            double bondFloorTmp = model1F->PriceEnd[numPriceArray];
            results->storeRequestResult(request, Maths::max(bondFloorTmp, putLevel[0]));
        }

        // Indicative vol
        if (control->requestsOutput(OutputRequest::IND_VOL, request)) {
            double             indVol     = 0.0;
            CVolRequestConstSP volRequest   = GetLNRequest();
            CVolRequestSP      ncVolRequest =CVolRequestSP::constCast(volRequest);

            LinearStrikeVolRequest* lsVolRequest = dynamic_cast<LinearStrikeVolRequest*>(ncVolRequest.get());
            // interpolate the vol using our LN request
            CVolProcessedBSSP volBS(getEquityAsset()->getProcessedVol(lsVolRequest));
            // calculate the indicative vol
            try {
                indVol = volBS->CalcVol(cvb->valueDate, cvb->bond->getMaturityDate());
            }
            catch (exception& ) {
                indVol = 0.0;
            }
            results->storeRequestResult(request, indVol);
        }

        // Duration
        if ( control->requestsOutput(OutputRequest::BOND_DURATION, request) ) {
            results->storeRequestResult(request, cvb->bond->duration(cvb->valueDate,risky));
        }
        // Convexity
        if ( control->requestsOutput(OutputRequest::BOND_CONVEXITY, request) ) {
            results->storeRequestResult(request, cvb->bond->convexity(cvb->valueDate,risky));
        }


        if (isOCB == true) {
            ocb->recordOutputRequests(control, results, model1F->PriceEnd[0] + Maths::max(bondFloor[0], putLevel[0]), bondFloor[0], fairValue);
        } else {
            // put probability
            if ( firstPutIndex >= 0 && control->requestsOutput(OutputRequest::PUT_PROBABILITY, request)) {
                if (calcedFirstPutProb == true) {
                    results->storeRequestResult(request, putProbability);
                } else {
                    string errorMsg = "could not calculate put probability";
                    UntweakableSP untweakablePutProb(new Untweakable(errorMsg));
                    results->storeRequestResult(request,untweakablePutProb);
                }
            }

            // put stock level
            if ( firstPutIndex >= 0 && control->requestsOutput(OutputRequest::PUT_STOCK_LEVEL, request)) {
                if (calcedFirstPutProb == true) {
                    results->storeRequestResult(request, putStockLevel);
                } else {
                    string errorMsg = "could not calculate put stock level";
                    UntweakableSP untweakablePutStock(new Untweakable(errorMsg));
                    results->storeRequestResult(request,untweakablePutStock);
                }
            }

            // put strike
            if ( firstPutIndex >= 0 && control->requestsOutput(OutputRequest::PUT_STRIKE, request) ) {
                double firstPutLevel = putLevel[firstPutIndex];
                double putConvRatio  = (convBool[firstPutIndex])?convRatio[firstPutIndex]:0.0;
                double putConvCash   = (convBool[firstPutIndex])?convCash[firstPutIndex]:0.0;
                double spotFX        = 1.0;
                double putStrike;

                if( convBool[firstPutIndex] && cvb->hasAddOnConvRatios() )
                {    
                    bool isConvertible;
                    cvb->getConversionInfo(model1F->TimePts.StepDates[firstPutIndex], firstPutLevel, &isConvertible, &putConvRatio, &putConvCash, false/*is not spot*/);
                }

                if ( Maths::isPositive(firstPutLevel) && Maths::isPositive(putConvRatio) ) {
                    if (StruckEquity::TYPE->isInstance(getEquityAsset().get())) {
                        // cast to struck equity object
                        const IObject* obj           = dynamic_cast<const IObject*>(getEquityAsset().get());
                        const StruckEquity* struckEq = dynamic_cast<const StruckEquity*>(obj);
                        if ( struckEq ) {
                            spotFX = struckEq->getFX()->getSpot();
                        }
                    }
                    putStrike = (firstPutLevel-putConvCash)/(putConvRatio*spotFX);
                    results->storeRequestResult(request,putStrike);
                } else {
                    UntweakableSP untweakablePutStrike(new Untweakable(
                        "firstPutLevel(" + 
                        Format::toString(firstPutLevel) + 
                        ") and conversion ratio(" +
                        Format::toString(putConvRatio) + 
                        ") must be strictly positive at first put date"));

                    results->storeRequestResult(request,untweakablePutStrike);
                }
            }

            // call strike
            if ( control->requestsOutput(OutputRequest::CALL_STRIKE, request) ) {

                if ( !(!cvb->callSchedule) && cvb->callSchedule->length() > 0 ) {
                   DateTime firstCallDate;
                   if ( cvb->callSchedule->getNextDate(cvb->valueDate,
                                                       cvb->bond->getMaturityDate(),
                                                       firstCallDate)) {
                       // AS: should possibly put in the correct fwd vol here
                      bool    isCallable, canConvert, isAdjustedForAccrued;
                      double  firstCallLevel, callConvRatio, callConvCash;
                      double spotFX        = 1.0;
                      double callStrike;
                      bool  isHardCall;

                      cvb->getCallLevel(firstCallDate, 0.0, &isCallable, &isAdjustedForAccrued, &isHardCall, &firstCallLevel);

                      bool isSpot = !cvb->hasAddOnConvRatios();
                      double spotOrLevel = isSpot?0:firstCallLevel;
                      cvb->getConversionInfo(firstCallDate, spotOrLevel, &canConvert, &callConvRatio, &callConvCash, isSpot);

                      if ( Maths::isPositive(firstCallLevel) && Maths::isPositive(callConvRatio) ) {
                         if ( StruckEquity::TYPE->isInstance(getEquityAsset().get())) {

                            // cast to struck equity object
                            const IObject* obj           = dynamic_cast<const IObject*>(getEquityAsset().get());
                            const StruckEquity* struckEq = dynamic_cast<const StruckEquity*>(obj);
                            if ( struckEq ) {
                               spotFX = struckEq->getFX()->getSpot();
                            }
                         }
                         callStrike = (firstCallLevel-callConvCash)/(callConvRatio*spotFX);
                         results->storeRequestResult(request,callStrike);
                      } else {
                         UntweakableSP untweakableCallStrike(new Untweakable(
                            "firstCallLevel(" + 
                            Format::toString(firstCallLevel) + 
                            ") and conversion ratio(" +
                            Format::toString(callConvRatio) + 
                            ") must be strictly positive at first call date"));

                            results->storeRequestResult(request,untweakableCallStrike);
                       }
                   }
               }
            }
        }
    }
}

/** returns a vol request for log-normal vol */
CVolRequestConstSP ConvBond1fProd::GetLNRequest()
{ 
    CVolRequestConstSP volRequest(
        new LinearStrikeVolRequest(cvb->getStrike(), 
                                   cvb->valueDate,
                                   cvb->bond->getMaturityDate(), 
                                   false));
  
    return volRequest;
}

/** Implementation of DDEInitiator interface, built on FD1FGeneric
 */
class ConvBond1fProdDDE: public ConvBond1fProd, virtual public IDDEInitiator
{
public:
    ConvBond1fProdDDE(const ConvBond *cvb) : ConvBond1fProd(cvb)
    {
        if (cvb->riskyGrowth)
            throw ModelException("ConvBond1fProdDDE", "Do not allow riskyGrowth flag under DDE");
    }

    DateTime maxMaturity() const { return cvb->getBondMaturityDate(); }

    void sensitiveDates( DateTimeArray    &dates ) const
    {
        dates.clear();

        // stop inserting extra dates to avoid cases where strike change
        // drastically between 2 very close dates, eg. strike on and right 
        // after a put date with INTERP_NONE. 
        return;

        /************ DEC/PERC need more thinking ***************/
        if ( isOCB || cvb->DECS == true || cvb->PERCS == true || !cvb->conversionRatios ) 
        {
            IDDEInitiator::sensitiveDates(dates);
            return;
        }

        // collect conv schedule dates and put dates
        const DateTimeArray &conversionDates = cvb->conversionRatios->getDateArray();
        dates = conversionDates;
        if ( cvb->putSchedule.get() ) 
        {
            const DateTimeArray& putDates = cvb->putSchedule->getDateArray();
            dates = DateTime::merge(dates, putDates);
        }
        // remove dates outside conversion schedule span
        DateTime::removeOutliers(dates, conversionDates[0], conversionDates[conversionDates.size()-1]);

        // remove not convertible dates
        for (vector<DateTime>::iterator iter(dates.begin()); iter != dates.end(); /* inc in loop */)
        {
            double convRatio, convCash;
            bool isConv;
            cvb->getConversionInfo(*iter, 0, &isConv, &convRatio, &convCash, true);

            if( !isConv )
                iter = dates.erase(iter);
            else
                ++iter;
        }
    }

    void sensitiveStrikes(    const DateTimeArray     dates,
                            DoubleArray             &strikes,    // same dimension as dates
                            bool                    &strikeIsPct) const
    {
        if( dates.size() != strikes.size() )
            throw ModelException("ConvBond1fProdDDE::sensitiveStrikes", "Dates and strikes dimension mismatch");

        /************ DEC/PERC need more thinking ***************/
        if ( isOCB || cvb->DECS == true || cvb->PERCS == true || !cvb->conversionRatios ) 
        {
            IDDEInitiator::sensitiveStrikes(dates, strikes, strikeIsPct);
            return;
        }

        strikeIsPct = false;

        // calculate bond floor (not using call level) and sensitive strike
        CashFlowArrayConstSP myCFs = cvb->bond->getCashFlowRef();
        int n=myCFs->size()-1;
        DateTime prevDate = (*myCFs).back().date;
        DateTimeArray putDates;
        if ( cvb->putSchedule.get() ) putDates = cvb->putSchedule->getDates();
        int nn = putDates.size()-1;
        double floor=0;
        for (int i=dates.size()-1; i>=0; i--) // loop backwards
        {
            while(n>=0 && dates[i]<=(*myCFs)[n].date)
            {
                // we assume put is INTERP_NONE
                while(nn>=0 && (*myCFs)[n].date < putDates[nn])
                {
                    floor *= risky->pv(putDates[nn], prevDate);
                    // get max(redemption, put level)
                    bool isPut;
                    double putLevel;
                    cvb->getPutLevel(dates[i], &isPut, &putLevel);
                    if( isPut && floor < putLevel ) floor = putLevel;
                    prevDate = putDates[nn];
                    nn--;
                }

                floor *= risky->pv((*myCFs)[n].date, prevDate);
                floor += (*myCFs)[n].amount;
                prevDate = (*myCFs)[n].date;
                n--;
            }
            // we assume put is INTERP_NONE
            while(nn>=0 && dates[i] <= putDates[nn])
            {
                floor *= risky->pv(putDates[nn], prevDate);
                // get max(redemption, put level)
                bool isPut;
                double putLevel;
                cvb->getPutLevel(dates[i], &isPut, &putLevel);
                if( isPut && floor < putLevel ) floor = putLevel;
                prevDate = putDates[nn];
                nn--;
            }
            floor *= risky->pv(dates[i], prevDate);
            prevDate = dates[i];

            bool isConv, isSpot = !cvb->hasAddOnConvRatios();
            double convRatio, convCash, spotOrLevel = isSpot?0:floor;
            cvb->getConversionInfo(dates[i], spotOrLevel, &isConv, &convRatio, &convCash, isSpot);

            if( isConv )
            {
                floor -= convCash;
                if( floor <= 0 )
                    throw(ModelException("ConvBond1fProdDDE::sensitiveStrikes","Strike < 0. Check redemption/put/conversion cash level."));

                strikes[i] = floor/convRatio;
            }
            else
                strikes[i] = (i==(dates.size()-1)?asset->getSpot():strikes[i+1]); // use next strike if not convertible

        }

    } // end of sensitiveStrikes

    bool isStaticSpread() const { return false; }
};

/** Implementation of CFDGridPass::IntoProduct interface */
FD1FGeneric::IProduct* ConvBond::createProduct(
    FD1FGeneric* model) const{
    ConvBond1fProd *tmpConvBond1fProd;


    if(FD1FDDE::TYPE->isInstance(model)) {
        tmpConvBond1fProd = new ConvBond1fProdDDE(this);
        tmpConvBond1fProd->hasE2CLayer = false;
    } else {
        tmpConvBond1fProd = new ConvBond1fProd(this);
        tmpConvBond1fProd->hasE2CLayer    = FD1FE2C::TYPE->isInstance(model) && 
                                            FirmAsset::TYPE->isInstance(asset.get());
    }
    tmpConvBond1fProd->genericFDModel = model;
    tmpConvBond1fProd->model1F        = model;

    return tmpConvBond1fProd;
}

/** Implementation of CFDGridPass::IntoProduct interface */
FD1F::IProduct* ConvBond::createProduct(
    FD1F* model) const {

    ConvBond1fProd *tmpConvBond1fProd = new ConvBond1fProd(this);
    tmpConvBond1fProd->fdModel = model;
    tmpConvBond1fProd->model1F = model;

    tmpConvBond1fProd->hasE2CLayer    = FD1FE2C::TYPE->isInstance(model) &&
                                        FirmAsset::TYPE->isInstance(asset.get());

    return tmpConvBond1fProd;
}


/** Implementation of CFDGridPass::IntoProduct interface */
FD1F::IProduct* OptOnConvBond::createProduct(
    FD1F* model) const{
    ConvBond1fProd *tmpConvBond1fProd = new ConvBond1fProd(cvb.get());
    tmpConvBond1fProd->fdModel = model;
    tmpConvBond1fProd->model1F = model;

    tmpConvBond1fProd->isOCB = true;
    tmpConvBond1fProd->ocb = this;

    tmpConvBond1fProd->hasE2CLayer    = FD1FE2C::TYPE->isInstance(model) &&
                                        FirmAsset::TYPE->isInstance(cvb->asset.get());

    return tmpConvBond1fProd;
}

// adjust risky bond floor price layers, assign risky bond floor to bondFloorGeneric
// return the new pEnd for non-risky-bond related price layers.
// notice that the address of bondFloorGeneric is assigned here as well.
int ConvBond1fProd::calcRiskyBondNAdjPriceLayer(int step, int bot, int top, 
                    int pEnd, double * const * price, double *&bondFloorGeneric)
{
    bondFloorGeneric = bondFloorGenericMem + bot;
    if( !isStaticSpread() )
    {
        calcBondRiskyNFloorBeforeMat(step, bot, top, price[pEnd], bondFloorGeneric);
        pEnd--; // make sure to decrease the number of layer to return the remaining layers
    }
    else // isE2CLayer
    {
        // this calculation is inheritted from prior version and is flawed
        // mainly it missed discouting of coupon for 1 time step
        // it does not contain other logic related to riskFree etc
        int layerBond = pEnd + 2;
        for(int layer = pEnd+1; layer <= pEnd+2; layer++)
        {
            for(int j=-bot; j<=top; j++)
            {
                // check coupon paid today ...
                if (addCouponBool[step+1] == true)
                    price[layer][j] += coupon[step+1];

                // check callable ... this should only be done for hard calls
                if (callBool[step] && callIsHard[step])
                    price[layer][j] = Maths::min(price[layer][j], callLevel[step]);
        
                if( layer==layerBond ) bondFloorGeneric[j] = price[layer][j];
    
                // check puttable ....
                if (putBool[step])
                    price[layer][j] = Maths::max(price[layer][j], putLevel[step]);
            }
        }
    }

    return pEnd;
}


// update bondRisky and bondRiskFree
void ConvBond1fProd::calcBondRiskyNFloorBeforeMat(int i/*step*/, int bot, int top, double *bondRiskyArray, double *bondFloorArray)
{
    for(int j=-bot; j<=top; j++)
    {
        // use the non-equity-participation spread to further discount the riskybond
        // this simple treatment requires the bond to have 0 recovery
        // use PV recovery rule for preferreds
        double &bondRisky = bondRiskyArray[j];
        if(  nonEqCreditPvs )  bondRisky *= nonEqCreditPvs[i];

        double bondRiskFree = bondRiskFreeArray[i];

        if (callIsHard[i] && callLevel[i] < (bondRisky + bondRiskFree) )
        {
            bondRisky = Maths::max((double)callLevel[i]-bondRiskFree, 0.);
        }
        if(bondFloorArray) bondFloorArray[j] = bondRisky + bondRiskFree;

        if( i>0 )
        {
            if ( putBool[i] == true &&
                !( model1F->TimePts.StepDates[i-1].getDate() == 
                model1F->TimePts.StepDates[model1F->TimePts.StepDates.size()-1].getDate()) ) 
            {
                bondRisky = Maths::max(bondRisky, putLevel[i]);
            }

            if (addCouponBool[i] == true) {
                if( !(model1F->TimePts.StepDates[i] <= cvb->riskFreeCouponEndDate) )
                    bondRisky += coupon[i];
            }
        }

    }
    return;
}






/**************************************************************************************************************************************/
/**************************************************************************************************************************************/
/**************************************************************************************************************************************/
/**************************************************************************************************************************************/
/**************************************************************************************************************************************/
/**************************************************************************************************************************************/

/** private class */
class ConvBondFDProd: public LatticeProdEDRIns
{   
public:

    ConvBondFDProd(const ConvBond* cvb, FDModel* model, bool isOCB = false);

    virtual ~ConvBondFDProd();
    
    virtual CAssetConstSP GetAssetRef() 
    {       
        return CAssetConstSP::dynamicCast((IObjectConstSP)getEquityAsset());
    }
    
    virtual bool GetFwdStartLV()
    {
        if ( !!inst->resetSchedule && inst->resetSchedule->length() > 0 ) 
        {
            return inst->resetSchedule->isFwdStart();
        }
        else
            return false;                   
    }
    
    virtual DateTime GetFwdStartDateLV()
    {
        if ( !!inst->resetSchedule && inst->resetSchedule->length() > 0 && inst->resetSchedule->isFwdStart()) 
        {
            return inst->resetSchedule->getStartDate(inst->valueDate);
        }
        else
            return inst->valueDate;
    }
    
    virtual YieldCurveConstSP GetDiscCurveRef() const;
    
    virtual CVolRequestConstSP GetLNRequest() const;
    
    virtual YieldCurveConstSP GetRiskFreeCurve() const;
 
    virtual void init(Control* control) const;
    
private:
    double InitProd1();    
    void InitProd2(double);

public:
    virtual void initProd();
    
    virtual bool Positive() 
    {
        // option payoff may < 0 if soft call level is less than bond flr.
        return false;
    }
    
    virtual string getCcyTreatment() const
    { 
        return inst->ccyTreatment;
    }

    double scalePremium(vector<double> & P, 
                        YieldCurveConstSP disc);
    
    /** extra output requests */    
    virtual void recordOutput(Control* control, 
                              YieldCurveConstSP disc, 
                              Results* results); 
      

    virtual void preCalc(int step);                         
    
     /** product payoff method at maturity */
    void prod_BWD_T(
        const TreeSlice & spot,
        int step,
        int bot,
        int top,
        int pStart,
        int pEnd,
        const vector< TreeSliceSP > & price);

    // this is payoff boundary condition, for KO, early exercise etc.
    void prod_BWD(
        const TreeSlice & spot,
        int step,
        int bot,
        int top,
        int pStart,
        int pEnd,
        const vector< TreeSliceSP > & price);

    virtual void update(int& step, FDProduct::UpdateType type);

    /** ignore start date if not forward starting */
    virtual DateTime getStartDate() const
    {
        return inst->valueDate;
    }

    void PenultStraight(const double* s, 
                              int step, 
                              int bot, 
                              int top, 
                              double *price);
    
    void PenultMandatory(const double* s, 
                               int step, 
                               int bot, 
                               int top, 
                               double *price);

    virtual double getCoupon(int step, 
                       const double* s, 
                             int start, 
                             int end);

    virtual bool hasEquityLayer();

    CAssetSP     getEquityAsset() 
    {
        if ( FirmAsset::TYPE->isInstance(asset.get()) ) 
        {
            FirmAsset* firmAsset = 0;
            if ( FirmAsset::TYPE->isInstance(asset.get())) 
            {
                CAsset*    ncAsset   = const_cast<CAsset*>(asset.get());
                firmAsset = dynamic_cast<FirmAsset*>(ncAsset);
            } 
            else 
            {
                throw ModelException("ConvBondProd::getEquityAsset",
                                     "Underlying asset must be a FirmAsset");
            }
            return firmAsset->getEquityAsset().getSP();
        } 
        else 
        {
            return asset.getSP();
        }
    }


    // following functions for spot dependent spread
    virtual bool isStaticSpread() const { return true; }

    void mergePriceForReset(const double* s, 
                                  int step, 
                                  int bot, 
                                  int top, 
                                  int pStart, 
                                  int pEnd,
                                  const vector< double * > & p);

    inline double adjConversionValueDECDelay(int step, 
                                             double s, 
                                             double convertValue, 
                                             double conversionValue);
    
    inline double calcAccelerateDECCallValue(int step,  
                                             double s, 
                                             double bondFloor);

    void calcPutLevelNProb(int step, 
                           int bot, 
                           int top, 
                           const double *s, 
                           const double *priceLayer, 
                           double putValue, 
                           double &putStockLevel, 
                           double &putProbability);

    double calcPutToStockValue(int step);

    int calcRiskyBondNAdjPriceLayer(int step, 
                                    int bot, 
                                    int top, 
                                    int pEnd, 
                                    const vector< double * > & p, 
                                    double *&bondFloorGeneric);

    void calcBondRiskyNFloorBeforeMat(int step, 
                                      int bot, 
                                      int top, 
                                      double *bondRiskyArray, 
                                      double *bondFloorArray);
public:
    const OptOnConvBond *ocb;

protected:

    const ConvBond* inst;

    bool  isOCB;
    int   numCVBPriceArrays;

    double *coupon;
    bool   *addCouponBool;
    double *callLevel;
    bool   *callBool;
    bool   *softPutBool;
    double *softPutLevel;
    double *softPutTrigger;
    bool   *callIsAdjForAccrued;
    bool   *callIsHard;
    double *putLevel;
    bool   *putBool;
    double *convRatio;
    double *convCash;
    double *parityFloor;
    bool   *convBool;
    bool   *convBeforeCall;
    double *bondFloor;

    double *bondFloorGenericMem;
    double *bondRiskFreeArray;
    double *nonEqCreditPvs;    // extra pv between steps for the portion of clean spread curve that eq does not participate 

    bool   *ocbExerBool;
    double *ocbStrike;
    double *ocbRiskyAdj;
    
    YieldCurveSP        risky;
    CAssetWrapper       asset;

    // members for resettable convertible bond - need to get these into a proper class
    // when the dust has settled
    CDoubleArraySP   minConvPrice;
    CDoubleArraySP   maxConvPrice;
    CDoubleArraySP   resetParity;

    CBoolArraySP     resetArray;
    int              resetType;
    int              firstGridStep;
    CDoubleMatrixSP  conversionRatioGrid;
    CDoubleMatrixSP  conversionRatioGridnoReset; // for triggerReset
    double           actualCR;

    // members required for put probability calculation
    int            firstPutIndex;
    double         putStockLevel;
    double         putProbability;
    bool           calcedFirstPutProb;

    // processed vol is required for Black on call
    CVolProcessedBSSP processedVol;

    // some members used to speed up performance;
    double          faceValue;
    bool            hasReset;
    CDoubleMatrix   priceCopy;

    ResetScheduleSP adjustedResetSchedule;

    // for fwdStarting Resettable
    bool            fwdStarting;
    double          pvCFbeforeStart;

};

FDProductSP ConvBond::createProduct(FDModel* model) const
{
    return FDProductSP( new ConvBondFDProd(this, model) );
}

ConvBondFDProd::ConvBondFDProd(const ConvBond* cvb, FDModel* model, bool isOCB) :
    LatticeProdEDRIns(model, 0, 0), 
    inst(cvb), isOCB( isOCB )
{
    // first: set discount curve
    if( tree1f )
        tree1f->setDiscountCurve( inst->discount.getSP() );

    // second: create spot payoff
    payoffIndex = model->createProduct( IProdCreatorSP( new
        IndexSpecEQ( inst->asset.getName(), inst->asset, inst->ccyTreatment ) ) );

    coupon              = 0;
    addCouponBool       = 0;
    callLevel           = 0;
    callBool            = 0;
    softPutBool         = 0;
    softPutLevel        = 0;
    softPutTrigger      = 0;
    callIsHard          = 0;
    callIsAdjForAccrued = 0;
    putLevel            = 0;
    putBool             = 0;
    convRatio           = 0;
    convCash            = 0;
    parityFloor         = 0;
    convBool            = 0;
    convBeforeCall      = 0;
    bondFloor           = 0;
    bondFloorGenericMem = 0;
    bondRiskFreeArray   = 0;
    nonEqCreditPvs        = 0;
    putStockLevel       = 0.0;
    putProbability      = 0.0;
    firstPutIndex       = -1;
    calcedFirstPutProb  = false;
    resetType           = 0;
    numCVBPriceArrays   = 0;
    pvCFbeforeStart     = 0.0; 

    ocbExerBool = 0;
    ocbStrike = 0;
    ocbRiskyAdj = 0;

    // make the risky curve
    // make sure that we have a BootstrappedYieldCurve otherwise we can't add the spreadCurve
    if (!BootstrappedYieldCurve::TYPE->isInstance(inst->discount.get()))
    {
        throw ModelException("ConvBondFDProd::ConvBondFDProd", 
                             "Yield curve must be a cash swap curve");
    }

    DateTime tmp(inst->getEffMaturityDate());
    DateTimeSP matDate(new DateTime(tmp));

    risky = inst->getRiskyCurve(matDate.get());
    asset = CAssetWrapper(copy(inst->asset.get()));
    if (inst->riskyGrowth && ICanBeRisky::TYPE->isInstance(getEquityAsset().get())) 
    {
        // cast ICanBeRisky pointer
        IObject* obj            = dynamic_cast<IObject*>(getEquityAsset().get());
        ICanBeRisky* riskyAssetI = dynamic_cast<ICanBeRisky*>(obj);
        CreditCurveWrapper* cs = 0;
        ICreditCurveSP stockCurve;
        if ( inst->useJointDefaultCorrelation ) 
        {
            ExpiryArraySP stockExp(new ExpiryArray(1));
            DoubleArraySP stockCredit(new DoubleArray(1));

            (*stockExp)[0]    = ExpirySP(new MaturityPeriod("30Y"));
            (*stockCredit)[0] = inst->stockCreditSpread;

            stockCurve = CreditSpreadCurveSP(new CreditSpreadCurve("STOCK_CREDIT", stockExp.get(), *(stockCredit.get())));

            riskyAssetI->makeRisky(stockCurve,matDate.get());
        } 
        else 
        {
            cs = (CreditCurveWrapper*)(&(inst->creditSpreads));
            riskyAssetI->makeRisky(cs->getSP(),matDate.get());
        }
    }

    numIns = 1;

    if (isOCB) 
    {
        if ( cvb->triggerReset ) 
        {
            throw ModelException("ConvBondFDProd::ConvBondFDProd", "Cannot price asset swaps on trigger resettables");
        }
        if ( cvb->hasContConversion() ) 
        {
            throw ModelException("ConvBondFDProd::ConvBondFDProd", "Cannot price asset swaps on cvb with contingenet conversion");
        }

        if ( !!cvb->resetSchedule && cvb->resetSchedule->length() > 0 ) 
        {
            // default number of conversion ratio samples
            numPrices   = 30;
            numCVBPriceArrays = 15;
        } 
        else 
        {
            // four price arrays for trigger reset, 1 otherwise
            // left trigger reset in here in case I ever get around to implement it
            if (cvb->triggerReset) 
            {
                numPrices   =  4;
                numCVBPriceArrays =  2;
            } 
            else 
            {
                numPrices   =  2;
                numCVBPriceArrays =  1;
            }
        }
    } 
    else 
    {
        if ( !!cvb->resetSchedule && cvb->resetSchedule->length() > 0 ) 
        {
            // default number of conversion ratio samples
            numPrices = 15;
        } 
        else 
        {
            // two price arrays for trigger reset, 1 otherwise
            if (cvb->triggerReset || 
                // if contingent conversion and convertible for rest of life, 1 array
                cvb->hasContConversion() )
            {
                numPrices =  2;
            } 
            else 
            {
                numPrices =  1;
            }
        }            
    }

    // add a last price array for bondrisky
    if( !isStaticSpread() ) numPrices++;
}

ConvBondFDProd::~ConvBondFDProd() {
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
    if (softPutBool != 0) {
        delete [] softPutBool;
        softPutBool = 0;}
    if (softPutLevel != 0) {
        delete [] softPutLevel;
        softPutLevel = 0;}
    if (softPutTrigger != 0) {
        delete [] softPutTrigger;
        softPutTrigger = 0;}
    if (callIsHard != 0) {
        delete [] callIsHard;
        callIsHard = 0;}
    if (callIsAdjForAccrued != 0) {
        delete [] callIsAdjForAccrued;
        callIsAdjForAccrued = 0;}
    if (putLevel != 0) {
        delete [] putLevel;
        putLevel = 0;}
    if (putBool != 0) {
        delete [] putBool;
        putBool = 0;}
    if (convRatio != 0) {
        delete [] convRatio;
        convRatio = 0;}
    if (convCash != 0) {
        delete [] convCash;
        convCash = 0;}
    if (parityFloor != 0) {
        delete [] parityFloor;
        parityFloor = 0;}
    if (convBool != 0) {
        delete [] convBool;
        convBool = 0;}
    if (convBeforeCall != 0) {
        delete [] convBeforeCall;
        convBeforeCall = 0;}
    if (bondFloor != 0) {
        delete [] bondFloor;
        bondFloor = 0;}
    if (bondFloorGenericMem != 0) {
        delete [] bondFloorGenericMem;
        bondFloorGenericMem = 0;}
    if (bondRiskFreeArray != 0) {
        delete [] bondRiskFreeArray;
        bondRiskFreeArray = 0;}
    if (nonEqCreditPvs != 0) {
        delete [] nonEqCreditPvs;
        nonEqCreditPvs = 0;}
    if (ocbExerBool != 0) {
        delete [] ocbExerBool;
        ocbExerBool = 0;}
    if (ocbStrike != 0) {
        delete [] ocbStrike;
        ocbStrike = 0;}
    if (ocbRiskyAdj != 0) {
        delete [] ocbRiskyAdj;
        ocbRiskyAdj = 0;}
}

YieldCurveConstSP ConvBondFDProd::GetDiscCurveRef() const
{
    if (inst->convertIntoIssuerShares == true && inst->riskyGrowth == false) 
    {
        return YieldCurveConstSP::dynamicCast((IObjectConstSP)inst->discount.getSP());
    } 
    else 
    {
        // exchangeables will be discounted at the risky rate
        if ( inst->useJointDefaultCorrelation ) 
        {
            DateTime tmp(inst->getEffMaturityDate());
            DateTimeSP matDate(new DateTime(tmp));
            YieldCurveConstSP riskyDiscountCurve = inst->getRiskyCurve(matDate.get(), inst->useJointDefaultCorrelation);
            return riskyDiscountCurve;
        } 
        else 
        {
            return YieldCurveConstSP::dynamicCast((IObjectConstSP)risky);
        }
    }
}

YieldCurveConstSP ConvBondFDProd::GetRiskFreeCurve() const
{
    return YieldCurveConstSP::dynamicCast((IObjectConstSP)inst->discount.getSP());
}

/** initialise model */
void ConvBondFDProd::init(CControl* control) const
{
    static const string method = "ConvBondFDProd::Init";
    try 
    {
        int i;
        
        // critical dates
        DateTimeArray critDates;

        CashFlowArraySP myCFs = inst->bond->getExAdjCashFlows(inst->valueDate, risky); 
        for (i=0; i<myCFs->size(); i++) 
        {
            critDates.push_back((*myCFs)[i].date);
        }
        if (inst->putSchedule.get()) 
        {
            const DateTimeArray& putDates = inst->putSchedule->getDateArray();
            for (i=0; i<putDates.size(); i++)
                critDates.push_back(putDates[i]);
        }
        if (inst->callSchedule.get()) 
        {
            const DateTimeArray& callDates = inst->callSchedule->getDateArray();
            for (i=0; i<callDates.size(); i++)
                critDates.push_back(callDates[i]);
        }
        if ( inst->softCallTriggerSchedule.get() ) 
        {
            const DateTimeArray& softcallDates = inst->softCallTriggerSchedule->getDateArray();
            for (i=0; i<softcallDates.size(); i++)
                critDates.push_back(softcallDates[i]);
        }
        if ( inst->softPutTriggerSchedule.get() ) 
        {
            const DateTimeArray& softPutDates = inst->softPutTriggerSchedule->getDateArray();
            for (i=0; i<softPutDates.size(); i++)
                critDates.push_back(softPutDates[i]);
        }

        if ( inst->conversionRatios.get() ) 
        {
            const DateTimeArray& conversionDates = inst->conversionRatios->getDateArray();
            for (i=0; i<conversionDates.size(); i++)
                critDates.push_back(conversionDates[i]);
        }
        if (isOCB == true) 
        {
            const DateTimeArray& exerDates = ocb->exerSched->getDateArray();
            for (i=0; i<exerDates.size(); i++)
                critDates.push_back(exerDates[i]);
        }

        // push in additional addon conversion dates for embedded warrant
        if( inst->hasAddOnConvRatios() )
            inst->addOnConvRatios->insertCritDates(critDates);

        // push in additional dates for contingenet conversion
        if( inst->hasContConversion() )
            inst->contConversion->insertCritDates(critDates);

        if(tree1f)
        {
            // default to DollarDivTreament = false !!! to be removed after EAS/EIS can handle it correctly !!!
            tree1f->SetDivAmountTreatment(false);    

            // default to NODE_INSERTION smoothing
            
            tree1f->equalTime = true;
            tree1f->NumOfPrice = numPrices;
            tree1f->NumOfInsertNode = numIns;
        }

        // add critical dates
        model->addCritDates( critDates );

        DateTimeArray segDates;
        segDates.resize(2);
        segDates[0] = inst->valueDate;
        segDates[1] = inst->bond->getMaturityDate();
        IntArray density( 1, 1 );
       
        // prepare timeline set up
        model->initSegments( segDates, density );
    }
    catch (exception& e) 
    {
        throw ModelException(e, method);
    }  
}

// just to be able to access protected members of Tree1f
class AccessTree1f : public CTree1f
{
    // this class cannot be instantiated
    AccessTree1f();
    friend class NAMESPACE::ConvBondFDProd;
};

// just to be able to access protected members of FD1D
class AccessFD1DMulti : public FD1DMulti
{
    // this class cannot be instantiated
    AccessFD1DMulti();
    friend class NAMESPACE::ConvBondFDProd;
};

void ConvBondFDProd::initProd() 
{
    if (tree1f){
		static_cast<AccessTree1f*>(tree1f)->Underlier = GetAssetRef();
		static_cast<AccessTree1f*>(tree1f)->discYC = GetDiscCurveRef();
	}else if (fd1dMulti){
		static_cast<AccessFD1DMulti*>(fd1dMulti)->underlying = GetAssetRef();
		static_cast<AccessFD1DMulti*>(fd1dMulti)->discYC = GetDiscCurveRef();
	} else {
		throw ModelException("initProd","Only tree1f or fd1dMulti is supported.");
	}
    
    initSlices( numPrices );
    initInsertNode();
}

double ConvBondFDProd::InitProd1()
{
    static const string method = "ConvBondFDProd::InitProd";

    int numStep           = model->getLastStep();
    coupon                = new double[numStep+1];
    addCouponBool         = new bool[numStep+1];

    // for performance
    faceValue = inst->bond->getFaceValue();
    hasReset  = (!!inst->resetSchedule && inst->resetSchedule->length() > 0 );

    // prepare for fwdStarting case
    fwdStarting = hasReset?inst->resetSchedule->isFwdStart():false;
    DateTime startDate = fwdStarting ? inst->resetSchedule->getStartDate(inst->valueDate):inst->valueDate;
    double fwdAtStart = inst->asset->fwdValue(startDate);
    if (fwdStarting) 
    {
        // :-( Under Solaris/gcc3.2, we dump core here if getStartLevel throws
        // an exception (test resettableinp/fwdRst24.xml), because the
        // compiler tries to be too clever about deconstructing the objects
        // that thereby go out of scope.

        fwdAtStart = inst->resetSchedule->getStartLevel(inst->valueDate, fwdAtStart);   // use Known Value
    }
    return fwdAtStart;
}

void ConvBondFDProd::InitProd2(double fwdAtStart)
{
    static const string method = "ConvBondFDProd::InitProd";
    try
    {   
        int i;
        CashFlowArraySP myCFs = inst->bond->getExAdjCashFlows(inst->valueDate, risky); 
        int numStep           = model->getLastStep();
        int cfIdx             = 0;
        double putLevelAtMat  = 0.0;
        DateTime startDate = fwdStarting ? inst->resetSchedule->getStartDate(inst->valueDate):inst->valueDate;

        int cfIdxFirst = 0;
        if (fwdStarting){
            while (startDate.getDate()>(*myCFs)[cfIdxFirst].date.getDate()){
                cfIdxFirst++;   // finding for the first cfIdx after TimePts.StepDates[0]
            }
        }
        cfIdx = cfIdxFirst;
            
        // generate adjusted reset schedule 
        if ( !!inst->resetSchedule) {
            adjustedResetSchedule = ResetScheduleSP(inst->resetSchedule.clone());

            // convert the schedule to bond currency for currency struck bonds
            if ( StruckEquity::TYPE->isInstance(getEquityAsset().get())) {
                StruckEquity* struckAsset = dynamic_cast<StruckEquity*>(getEquityAsset().get());
                double       fxSpot      = struckAsset->getFXSpot();

                // get all reset dates
                DateTimeArray resetDates = adjustedResetSchedule->getDates();
                DoubleArraySP fxFwds(new DoubleArray(resetDates.size()));
                for(int i=0 ; i<fxFwds->size() ; ++i) {
                    if (resetDates[i] >= inst->valueDate ) {
                        (*fxFwds)[i] = struckAsset->fxFwdValue(resetDates[i]);
                    } else {
                        (*fxFwds)[i] = 0.0;
                    }
                }

                // calculate initial conversion price
                double initialCR        = inst->conversionRatios->firstValue();
                if (fwdStarting){
                    initialCR /= fwdAtStart;
                    // fwdAtStart is already in payoff currency Unit (not Underlying Ccy Unit)
                    // in a few lines later, adjustedResetSchedule->preProcessSchedule(....)
                    // fx rate is multiplying ResetLevels, by assuming that reset level are given by Underlying Ccy.
                    // Thus, the Reset Level at heare should be Ccy Unit of underling.  Not payoff.
                    adjustedResetSchedule->scaleLevels(fwdAtStart/struckAsset->fxFwdValue(startDate));
                }
                if (Maths::isZero(initialCR)) {
                    initialCR = DBL_EPSILON;
                }
                double initialConvPrice = inst->bond->getFaceValue() / initialCR;

                // convert the reset schedule into a currency struck schedule, if necessary
                adjustedResetSchedule->preProcessSchedule(true, fxSpot, fxFwds, initialConvPrice);
            } else {
                double       fxSpot      = 1.0;

                // get all reset dates
                DateTimeArray resetDates = adjustedResetSchedule->getDates();
                DoubleArraySP fxFwds(new DoubleArray(resetDates.size()));
                for(int i=0 ; i<fxFwds->size() ; ++i) {
                    (*fxFwds)[i] = 1.0;
                }

                // calculate initial conversion price
                double initialCR        = inst->conversionRatios->firstValue();
                if (fwdStarting){
                    initialCR /= fwdAtStart;
                    adjustedResetSchedule->scaleLevels(fwdAtStart);
                }
                if (Maths::isZero(initialCR)) {
                    initialCR = DBL_EPSILON;
                }
                double initialConvPrice = inst->bond->getFaceValue() / initialCR;

                // convert the reset schedule into a currency struck schedule, if necessary
                adjustedResetSchedule->preProcessSchedule(false, fxSpot, fxFwds, initialConvPrice);
            }
        }


        // assume the cashFlows are in order
        for (i=0; i<=numStep; i++) {
            if ( cfIdx <  myCFs->size() && 
                model->getDate(i).getDate() == (*myCFs)[cfIdx].date.getDate() ) {
                if ( i == 0 || 
                     model->getDate(i).getDate() != 
                     model->getDate(i-1).getDate() ) {
                    addCouponBool[i] = true;
                    // I might have a problem here with coupons at maturity
                    if (cfIdx == myCFs->size()-1) {
                        if (inst->getCouponAtMat == true) {
                            addCouponBool[i] = true;
                            coupon[i] = (*myCFs)[cfIdx].amount - inst->bond->getRedemption();
                        } else {
                            addCouponBool[i] = false;
                            coupon[i] = 0.;
                        }
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

        // add the dividend schedule 
        if ( inst->dividendPassThrough ) {
            DividendListSP divSchedule = DividendCollector::divsBetweenDates(getEquityAsset().get(),
                                                              inst->valueDate,
                                                              inst->valueDate,
                                                              inst->bond->getMaturityDate());

            if (!!divSchedule) {
                const DividendArray& divArray = divSchedule->getArray();

                double convRatio = 0.0;
                if (inst->DECS == true ) {
                   convRatio = inst->initialConvRatio;
                } else {
                   if ( inst->conversionRatios.get() && inst->conversionRatios->length() > 0) {
                      convRatio = inst->conversionRatios->firstValue();
                   }
                }

                // assume the dividends are in order
                cfIdx = cfIdxFirst;
                for (i=0; i<=numStep; i++) {
                    if ( cfIdx <  divArray.size() && 
                        model->getDate(i).getDate() == divArray[cfIdx].getPayDate().getDate() ) {
                        if ( i == 0 || 
                             model->getDate(i).getDate() != 
                             model->getDate(i-1).getDate() ) {
                            addCouponBool[i] = true;
                            // I might have a problem here with coupons at maturity
                            if (cfIdx == divArray.size()-1) {
                                if (inst->getCouponAtMat == true) {
                                    addCouponBool[i] = true;
                                    coupon[i] += convRatio * inst->dividendPassThroughPct * divArray[cfIdx].getDivAmount();
                                }
                            } else {
                                coupon[i] += convRatio * inst->dividendPassThroughPct * divArray[cfIdx].getDivAmount();
                            }
                            cfIdx++;
                        }
                    }
                }
            }
        }

        callLevel           = new double[numStep+1];
        callBool            = new bool[numStep+1];
        callIsAdjForAccrued = new bool[numStep+1];
        callIsHard          = new bool[numStep+1];

        double dummy = getEquityAsset()->getSpot();
        vector<double> fwdVol;
        for (i=0; i<=numStep; i++) {
            if(tree1f) {
				tree1f->GetStepVol(i, fwdVol, &dummy,0,0);
			} else if (fd1dMulti){
				fd1dMulti->GetStepVol(i, fwdVol, &dummy, 0, 0);
			}
            inst->getCallLevel(model->getDate(i), fwdVol[0], &(callBool[i]), &(callIsAdjForAccrued[i]), 
                &(callIsHard[i]), &(callLevel[i]));
        } 

        putLevel   = new double[numStep+1];
        putBool    = new bool[numStep+1];

        for (i=0; i<=numStep; i++) {
            inst->getPutLevel(model->getDate(i), &(putBool[i]), &(putLevel[i]));
            // first put date is required for put probability calculation
            if ( putBool[i] && firstPutIndex < 0 ) {
                firstPutIndex = i;
            }
        }

        softPutBool         = new bool[numStep+1];
        softPutLevel        = new double[numStep+1];
        softPutTrigger      = new double[numStep+1];

        for (i=0; i<=numStep; i++) {
            inst->getSoftPutLevel(model->getDate(i), &(softPutBool[i]), &(softPutTrigger[i]), &(softPutLevel[i]));
        }

        convRatio   = new double[numStep+1];
        convCash    = new double[numStep+1];
        parityFloor = new double[numStep+1];
        convBool    = new bool[numStep+1];

        for (i=0; i<=numStep; i++) {
            inst->getConversionInfo(model->getDate(i), 0/*dummyspot*/, &(convBool[i]), &(convRatio[i]), &(convCash[i]), true, &(parityFloor[i]));
            if ( Maths::isZero(convRatio[i])) {
                convBool[i] = false;
            }
            if (fwdStarting)
                convRatio[i] /= fwdAtStart;
        }

        convBeforeCall = new bool[numStep+1];

        convBeforeCall[0] = true; // handle the first point specially. no harm in setting it true. OCB will be wrong if set to false.
        // There can be some trouble with setting convBeforeCall[0] = true. The cvb value can be above the call level 
        // when parity is below. Should fix this in way that doesn't screw up OCB.
        for (i=1; i<=numStep; i++) {
            if (callBool[i] == false){
                // can't call now
                convBeforeCall[i] = true; // must be true or OCB will screw up
            } else {
                if (callBool[i-1] == false) {
                    // can call now but not at previous step
                    convBeforeCall[i] = true;
                } else {
                    // can call now and at previous step
                    if (callLevel[i] < callLevel[i-1]*.9975) { // give 25 bp leeway
                        convBeforeCall[i] = true;
                    } else {
                        convBeforeCall[i] = false;
                    } 
                }
            }
        }

        if (firstPutIndex > 0 && model->getDate(firstPutIndex) <= inst->riskFreeCouponEndDate)
            throw ModelException("ConvBondFDProd::InitProd", 
                                 "the model does not handle puts that are at or before the end date for risk free coupons");

        bondFloor = new double[numStep+1];
        if (inst->getCouponAtMat == true){
            if (inst->DECS == true || inst->PERCS == true) {
                bondFloor[numStep] = 0.; // note that getCouponAtMat is validated to be true.
            } else {
                if (putBool[numStep] == true) {
                    bondFloor[numStep] = Maths::max(inst->bond->getRedemption(), putLevel[numStep]);
                } else {
                    bondFloor[numStep] = inst->bond->getRedemption();
                }
            }
        } else {
            if (putBool[numStep] == true) {
                // we always get the coupon if we put at mat - only in the first case is the accrued interest
                // not yet included in the put level
                if ( !inst->putAdjustForAccrued ) {
                    putLevelAtMat = putLevel[numStep] + inst->getAccruedAtDate(model->getDate(numStep));
                } else {
                    putLevelAtMat = putLevel[numStep];
                }
                bondFloor[numStep] = Maths::max((*myCFs)[myCFs->size()-1].amount, putLevelAtMat);
            } else {
                bondFloor[numStep] = (*myCFs)[myCFs->size()-1].amount;
            }
        }

        double  bondRiskFree = 0;
        double  bondRisky = bondFloor[numStep];

        // check for stock being in Escrow account
        if (inst->stockInEscrow) {
            double parity     = getEquityAsset()->getSpot() * convRatio[numStep];
            double difference = parity - bondRiskFree;
            bondRiskFree += difference;
            bondRisky    -= difference;
        }

        if( !isStaticSpread() )
        {
            bondRiskFreeArray = new double[numStep+1];
            bondRiskFreeArray[numStep] = bondRiskFree;
            
            // make sure to get non equity portion of spread
            if( dynamic_cast<IHaveNonEqSpread *>(tree1f) )
            {
                const CleanSpreadCurve *nonEqCleanSprds = dynamic_cast<IHaveNonEqSpread *>(tree1f)->getNonEqSpreads();
                if( nonEqCleanSprds == 0 )
                    throw ModelException(method, "Model has non eq spread but no eq spread is returned");

                nonEqCreditPvs = new double[numStep];
                for(i=0; i<numStep; i++)
                    nonEqCreditPvs[i] = nonEqCleanSprds->getDefaultPV(model->getDate(i), model->getDate(i+1));
            }
        }

        for (i=numStep-1; i>=0; i--) {
            double accruedIR = inst->getAccruedAtDate(model->getDate(i));
            // We assume here that puts cannot be before riskFreeCouponEndDate
            if ( putBool[i+1] == true && 
                 !( model->getDate(i).getDate() == 
                    model->getDate(numStep-1).getDate()) &&
                    !((inst->DECS == true || inst->PERCS == true) && i == numStep-1)) {
                bondRisky = Maths::max(bondRisky, putLevel[i+1]);
            } 

            if (addCouponBool[i+1] == true) {
                if (model->getDate(i+1) <= inst->riskFreeCouponEndDate)
                    bondRiskFree += coupon[i+1];
                else
                    bondRisky += coupon[i+1];
            }

            if ( (model->getDate(i+1) <= inst->riskFreeCouponEndDate && bondRiskFree > 0) ||
                 inst->stockInEscrow ) {
                bondRiskFree *= inst->discount.getSP()->pv(model->getDate(i), 
                                     model->getDate(i+1));
            }

            // use PV recovery rule for preferreds
            if (inst->DECS || inst->PERCS) {
                bondRisky = risky->riskyPV(model->getDate(i),
                                           model->getDate(i+1),
                                           bondRisky,
                                           bondRisky,
                                           inst->useAssetRecovery,
                                           inst->recoveryPct);
            } else {
                bondRisky = risky->riskyPV(model->getDate(i),
                                           model->getDate(i+1),
                                           bondRisky,
                                           inst->bond->getNotional(model->getDate(i)) + accruedIR,
                                           inst->useAssetRecovery,
                                           inst->recoveryPct);
            }

            if (callIsHard[i] && callLevel[i] < bondRisky + bondRiskFree)
            {
                bondRisky = Maths::max((double)callLevel[i]-bondRiskFree, 0.);
                if (bondRisky == 0)
                    bondRiskFree = callLevel[i];
            }
            bondFloor[i] = bondRisky + bondRiskFree;

            if( !isStaticSpread() )
            {
                bondRiskFreeArray[i] = bondRiskFree;
            }
        }

        // add all coupons before fwdStartdate, after value date.
        if (fwdStarting){
            pvCFbeforeStart = inst->bond->couponsPV(inst->valueDate,startDate,risky); // pv of cpns before starting.
            double dfRiskFree = inst->discount->pv(inst->valueDate,startDate);
            double dfRisky    = risky->pv(inst->valueDate,startDate);
            double riskyBondPVAfterStart = inst->bond->presentValue(inst->valueDate, risky) - pvCFbeforeStart;
            riskyBondPVAfterStart *= (1.0-dfRiskFree/dfRisky);
            pvCFbeforeStart += riskyBondPVAfterStart;
            pvCFbeforeStart /= dfRiskFree;
        }

        if (isOCB) {
            try {
                ((OptOnConvBond *)ocb)->strikeAtOptionMat = inst->bond->getNotional(ocb->exerSched->lastDate());
            } catch (exception&) {
                ((OptOnConvBond *)ocb)->strikeAtOptionMat = faceValue;
            }
        }


        if (isOCB == true) {
            ocbExerBool  = new bool[numStep+1];
            ocbStrike    = new double[numStep+1];
            for (i=0; i<=numStep; i++) {
                ocb->getStrike(model->getDate(i), bondFloor[i], true, &(ocbExerBool[i]), &(ocbStrike[i]));  
            }

            ocbRiskyAdj = new double[numStep+1];
            ocbRiskyAdj[numStep] = 0.;
            long lastExerIdx = -1;
            for (i=numStep-1; i>=0; i--) {
                if (ocbExerBool[i+1] == true) {
                    ocbRiskyAdj[i] = Maths::max(bondFloor[i+1] - ocbStrike[i+1], 0.);
                    lastExerIdx = i+1;
                } else if (lastExerIdx < 0) {
                    ocbRiskyAdj[i] = 0;
                } else {
                    ocbRiskyAdj[i] = Maths::max(bondFloor[lastExerIdx] - ocbStrike[lastExerIdx], 0.);
                }
                if (lastExerIdx > i+1) {
                     ocbRiskyAdj[i] *= 
                        risky->pv(model->getDate(i+1), model->getDate(lastExerIdx));
                }
                if (lastExerIdx >= 0) {
                    ocbRiskyAdj[i] *= 
                        (risky->pv(model->getDate(i), model->getDate(i+1)) -
                        inst->discount->pv(model->getDate(i), model->getDate(i+1)));
                }
            }
            // a little validation that's hard to do earlier
            bool isExerAfter = false;
            for (i=numStep-1; i>=0; i--) {
                if (isExerAfter == false && ocbExerBool[i] == true) {
                    isExerAfter = true;
                }
                if (isExerAfter == true && callBool[i] == true && ocbExerBool[i] == false) {
                    throw ModelException("ConvBondFDProd::InitProd", 
                                         "option on convertible must be exercisable whenever it is callable");
                }
            }
        }

        if (hasReset) 
        {
            int numStep             = model->getLastStep();
            minConvPrice            = CDoubleArraySP(new CDoubleArray(numStep + 1));
            maxConvPrice            = CDoubleArraySP(new CDoubleArray(numStep + 1));
            resetParity             = CDoubleArraySP(new CDoubleArray(numStep + 1));
    
            resetArray              = CBoolArraySP(new BoolArray(numStep + 1));
            firstGridStep           = numStep+1;
            int numOfPriceArray     = 0;
    
            numOfPriceArray =  15;

            conversionRatioGrid              = CDoubleMatrixSP(new CDoubleMatrix(numStep + 1, numOfPriceArray));
            conversionRatioGridnoReset       = CDoubleMatrixSP(new CDoubleMatrix(numStep + 1, numOfPriceArray));
    
            double initCVR = inst->conversionRatios->firstValue();
            if (fwdStarting) 
            {
                initCVR /=  fwdAtStart;
            }
    
            double crAtValueDate    = adjustedResetSchedule->getCurrentConversionRatio(
                initCVR,
                model->getDate(0),
                faceValue);
    
            for (i=0; i<=numStep; i++) 
            {
                // check if there is a reset today
                if ( adjustedResetSchedule->hasReset(model->getDate(i))) 
                {
                    (*resetArray)[i]     = true;

                    // maximum conversion price
                    (*maxConvPrice)[i]   = adjustedResetSchedule->getMaxResetStrike(model->getDate(i)); 

                    // minimum conversion price
                    (*minConvPrice)[i]   = adjustedResetSchedule->getMinResetStrike(model->getDate(i)); 
           
                    // minimum conversion price
                    (*minConvPrice)[i]   = Maths::min((*minConvPrice)[i],(*maxConvPrice)[i]);
                    (*resetParity)[i]    = adjustedResetSchedule->getParity(model->getDate(i)); 
            
                    if (i < firstGridStep) 
                    {
                        firstGridStep =  i;
                    }
                } 
                else 
                {
                    if ( i > 0 ) 
                    {
                        (*resetArray)[i]     = false;
                        (*maxConvPrice)[i]   = (*maxConvPrice)[i-1];
                        (*minConvPrice)[i]   = (*minConvPrice)[i-1];
                        (*resetParity)[i]    = (*resetParity)[i-1];
                    } 
                    else 
                    {
                        (*resetArray)[i]     = false;
                        double initialCR;
                        if ( model->getDate(i) < inst->conversionRatios->firstDate()) 
                        {
                            initialCR = inst->conversionRatios->firstValue();
                            if (fwdStarting)
                                initialCR /= fwdAtStart;
                        } 
                        else 
                        {
                            initialCR = convRatio[i];
                        }

                        double currentCR     = adjustedResetSchedule->getCurrentConversionRatio(
                            initialCR,
                            model->getDate(i),
                            faceValue);
                        (*maxConvPrice)[i]   = faceValue / currentCR;
                        (*minConvPrice)[i]   = faceValue / currentCR;
                        (*resetParity)[i]    = 1.0;
                        actualCR=adjustedResetSchedule->getCurrentConversionRatio(
                            initialCR, model->getDate(i), faceValue);
                    }
                }

                // set up the conversion ratio grid
                int j;
                double minCR = faceValue / (*maxConvPrice)[i];
                double maxCR = faceValue / (*minConvPrice)[i];

                string resetTypeStr = adjustedResetSchedule->getResetType(model->getDate(i));

                if ( resetTypeStr == "UP" && adjustedResetSchedule->isFlat(model->getDate(i))) 
                {
                    // up reset only
                    minCR = Maths::max(minCR, crAtValueDate);
                    maxCR = Maths::max(minCR, maxCR);
                }

                for (j=0;j<numOfPriceArray;++j) 
                {
                    (*conversionRatioGrid)[i][j] = exp(log(minCR) + (log(maxCR) - log(minCR)) /( numOfPriceArray-1) * j);
                    //(*conversionRatioGridnoReset)[i][j] = exp(log(minCR) + (log(maxCR) - log(minCR)) /( numOfPriceArray-1) * j);
                }
            }
        } 
        else 
        {
            firstGridStep = 0;
        }

        // set vol for trigger adjustment in contingenet conversion
        if( inst->hasContConversion() ) {
            CVolRequestConstSP volRequest = GetLNRequest();
            CVolProcessedSP  procVol(inst->asset->getProcessedVol(volRequest.get()));
            CVolProcessedBSSP volBS = CVolProcessedBSSP::dynamicCast(procVol);
            inst->contConversion->setVol(volBS);
        }
    }
    catch (exception& e) 
    {
        throw ModelException(e, method);
    }
}

/** calculate barriers and place barriers at inserted node if needed */
void ConvBondFDProd::preCalc(int step)
{
    static const string method("ConvBondDProd::preCalc");

    if (step == model->getLastStep())
        InitProd2(InitProd1());

    try 
    {
        	if (tree1f){
			int idx = tree1f->getSliceIndex(step);
			if (tree1f->GetSmoothMethod() == CTree1f::NODE_INSERTION) 
			{
				if( inst->hasAddOnConvRatios() )
					throw ModelException("ConvBondFDProd::preCalcTree", "Don't know how to insert node yet if hasAddOnConvRatio");
	        
				if (convRatio[step] > 0) 
				{
					tree1f->SetInsertNode(idx, 0, (callLevel[step]-convCash[step])/convRatio[step], 0);
				}
				else 
				{    
					tree1f->SetInsertNode(idx, 0, 0., 0);
				}
			}
		}
    }
    catch (exception& e) 
    {
        throw ModelException(e, method);
    }
}

/** product payoff method at steps at maturity */
void ConvBondFDProd::prod_BWD_T(const TreeSlice & spot,
                                      int step, 
                                      int bot, 
                                      int top, 
                                      int pStart, 
                                      int pEnd,
                                      const vector< TreeSliceSP > & price)   
{
    double * s = spot.getValues();
    const vector< double * > & p = getValues( price );

    int     i,j;
    double  strike = 0.;
    
    // allocate local memory for bondFloorGeneric. assume that tree/fd need max bond floor memory at maturity
    if( bondFloorGenericMem == 0 && (!isStaticSpread()) ) 
        bondFloorGenericMem = new double[-bot + top + 1];
    
    // if not static spread, need to initiate the risky bond prices and get the non-eq related spread curve
    if (!isStaticSpread()) 
    {
        int i = step;
        double  bondRiskFree = 0;
        double  bondRisky = bondFloor[step];
        
        // ********** escrow code is suspicious since bondRiskFree is always 0 prior to this calculation!
        // check for stock being in Escrow account
        if (inst->stockInEscrow) 
        {
            double parity     = getEquityAsset()->getSpot() * convRatio[i];
            double difference = parity - bondRiskFree;
            bondRisky    -= difference;
        }

        // We assume here that puts cannot be before riskFreeCouponEndDate
        if ( putBool[i] == true && 
             !( model->getDate(i-1).getDate() == 
                model->getDate(model->getLastStep()-1).getDate()) &&
             !(inst->DECS == true || inst->PERCS == true) ) 
        {
            bondRisky = Maths::max(bondRisky, putLevel[i]);
        } 
        
        if (addCouponBool[i] == true) 
        {
            if ( !( model->getDate(i) <= inst->riskFreeCouponEndDate ) ) 
            {
                bondRisky += coupon[i];
            }
        }
        for (j=bot; j<=top; j++) 
        {
            (p[pEnd])[j] = bondRisky;
        }
        
        pEnd--;
    }
    
    // set numCVBPriceArrays if we're pricin on a tree - this should possibly be in preCalcTree, which
    // would require touching the interface in many payouts.
    
    if (isOCB) 
    {
        numCVBPriceArrays = pEnd / 2;
    } 
    else 
    {
        numCVBPriceArrays = pEnd;
    }
   

    int endIdx = (isOCB && hasReset)?numCVBPriceArrays-1:pEnd;
    
    // when we put on the maturity date we always get the coupon
    double put = putLevel[step];
    if ( putBool[step] == true && inst->putAdjustForAccrued == false ) 
    {
        put += inst->getAccruedAtDate(inst->bond->getMaturityDate());
    }

    for (i=pStart; i<=endIdx; i++) 
    {
        if (inst->DECS == true) 
        {
            for (j=bot; j<=top; j++) 
            {
                if (s[j] < inst->initialPrice) 
                {
                    (p[i])[j] = s[j]*inst->initialConvRatio;
                } 
                else if (s[j] < inst->convPrice)
                {
                    (p[i])[j] = inst->initialPrice * inst->initialConvRatio;
                } 
                else if ( inst->decsHasCutoff && s[j] > inst->decsCutoffLevel ) 
                {
                    (p[i])[j] = inst->initialPrice * inst->initialConvRatio + (inst->decsCutoffLevel - inst->convPrice)*inst->minConvRatio;
                } 
                else 
                {
                    (p[i])[j] = inst->initialPrice * inst->initialConvRatio + (s[j] - inst->convPrice)*inst->minConvRatio;
                }               
            }        
        } 
        else if (inst->PERCS == true) 
        {
            for (j=bot; j<=top; j++) 
            {
                if (s[j] < inst->initialPrice) 
                {
                    (p[i])[j] = s[j]*inst->initialConvRatio;
                } 
                else 
                {
                    (p[i])[j] = inst->initialPrice*inst->initialConvRatio;
                } 
            }
        } 
        else 
        {
            strike = bondFloor[step];
            double conversionRatio, conversionCash = convCash[step];
            for (j=bot; j<=top; j++) 
            {                
                if (hasReset && (endIdx - pStart > 0)) 
                {
                    // trigger no reset code to has here
                    conversionRatio = (*conversionRatioGrid)[step][i];
                }
                else 
                {
                    if( convBool[step] && inst->hasAddOnConvRatios() ) 
                    {
                        bool isConvertible;
                        inst->getConversionInfo(model->getDate(step), s[j], &isConvertible, &conversionRatio, &conversionCash);
                    }
                    else 
                    {
                        conversionRatio = convRatio[step];
                    }
                }
                
                double convertValue = 0.0;

                if (inst->triggerReset                                                                      &&
                    inst->resetType == "UP_RESET"                                                           &&
                    model->getDate(step).getDate() >= inst->resetObservationStartDate.getDate() &&
                    model->getDate(step).getDate() <= inst->resetObservationEndDate.getDate()   &&
                    inst->resetTrigger <=  s[j] ) 
                {
                    conversionRatio = inst->resetConversionRatio;
                } 
                else if ( 
                    inst->triggerReset                                                                       &&
                    inst->resetType == "DOWN_RESET"                                                          &&
                    model->getDate(step).getDate() >= inst->resetObservationStartDate.getDate()  &&
                    model->getDate(step).getDate() <= inst->resetObservationEndDate.getDate()    &&
                    inst->resetTrigger >=  s[j]  ) 
                {
                    conversionRatio = inst->resetConversionRatio;
                }
                
                if( inst->hasContConversion() && inst->contConversion->isTrigActiveAtMat() &&
                    i==(pStart+inst->contConversion->isHistEnabled(inst->getValueDate())) ) 
                {
                    convertValue = 0.0; // non-convertible price array
                } 
                else if ( inst->contingentConversion && inst->triggerActiveAtMat && 
                          conversionRatio*s[j] < inst->contingentConversionTrigger ) 
                {
                    convertValue = 0.0;
                } 
                else 
                {
                    convertValue = conversionRatio*s[j]+conversionCash;
                }
                double holderValue = Maths::max(convertValue, put);
                
                double issuerValue = strike;
                
                (p[i])[j] = Maths::max((Maths::max(holderValue, issuerValue) - strike), 0.0);
                
                if ( inst->triggerReset ) 
                {
                    convertValue = inst->resetConversionRatio*s[j]+convCash[step];
                    double holderValue = Maths::max(convertValue, put);
                    double issuerValue = bondFloor[step];
                    
                    (p[i+1])[j] = Maths::max((Maths::max(holderValue, issuerValue) - strike), 0.0);
                }
            }    
        }
    }

    // does not handle options on resettables yet
    if (isOCB && pEnd > pStart) 
    { // we're an OCB
        for (i=numCVBPriceArrays ; i<=pEnd ; ++i) {
            for (j=bot; j<=top; j++) 
            {

                if (ocbExerBool[step] == true) 
                {
                    (p[i])[j] = Maths::max((p[i-numCVBPriceArrays])[j] + strike - ocbStrike[step], 0.);
                } 
                else 
                {
                    (p[i])[j] = 0.;
                }
            }
        }
    }
}

/** product payoff method at steps earlier than maturity */
void ConvBondFDProd::prod_BWD(const TreeSlice & spot,
                                    int step, 
                                    int bot, 
                                    int top, 
                                    int pStart, 
                                    int pEnd,
                                    const vector< TreeSliceSP > & price) 
{
    static const string method = "ConvBondFDProd::prod_BWD";

    double * s = spot.getValues();
    const vector< double * > & p = getValues( price );

    int         j;
    int         layer;
    double      convertValue;
    bool        aboveCall;
    double      slope;
    double      conversionValue;
    double      conversionValueForCall; // distinguish between conversion value from convertibility and from right to convert when called
    DoubleArray offset;
    DoubleArray coefficient;
    double      conversionRatioForCall, conversionCash = convCash[step];
    double      conversionRatio = convRatio[step];

    int         layerCoCoCV=-1, layerCoCoNCV=-1, layerTrigResetRatio=-1;
    
    //////////////////////////////////////////////////////////////
    //                                                          //
    //  For debt equity model, need to calc spot dependent      //
    //  bond floor. also adjust pEnd if needed                  //
    //                                                          //
    //////////////////////////////////////////////////////////////
    
    bool        isBondFlrGeneric = ( !isStaticSpread() );
    double      bondFloorTmp = bondFloor[step], *bondFloorGeneric=0;
    if( isBondFlrGeneric ) 
    {
        pEnd = calcRiskyBondNAdjPriceLayer(step, bot, top, pEnd, p, bondFloorGeneric);
    }
    
    int endIdx = pStart;
    if ( hasReset || inst->triggerReset || inst->hasContConversion() ) 
    {
        if ( step >= firstGridStep ) 
        {
            endIdx = (isOCB)?numCVBPriceArrays-1:pEnd;
        } 
        else 
        {
            // only need single grid if there is no reset before the current step
            endIdx = pStart;
        }
        
        if( inst->triggerReset ) 
        {
            layerTrigResetRatio = pStart + 1;
            if( endIdx != pStart + 1 ) 
            {
                throw ModelException(method, "Internal error. 2 cvb pricing grids needed for trigger reset");
            }
        }
        else if ( inst->hasContConversion() ) 
        {   
            if( endIdx != pStart + 1 ) 
            {
                throw ModelException(method, "Internal error. 2 cvb pricing grids needed for contingent conversion");
            }
            
            bool idxFlag = inst->contConversion->isHistEnabled(inst->getValueDate());
            layerCoCoCV = pStart + !idxFlag;
            layerCoCoNCV = pStart + idxFlag;
        }
    }
    
    // first do Penult adjustment if required
    if (step == model->getLastStep()-1) 
    {
        for (layer=pStart; layer<=endIdx; ++layer) 
        {
            // disable smoothing for spot dependent spread
            if( isStaticSpread() && 
                !( inst->hasContConversion() && inst->contConversion->isTrigActiveAtMat() && layer==layerCoCoNCV ) && 
                // disable smoothing for no-convert price slice if contingent conversion and trigger active at maturity
                ( !adjustedResetSchedule || adjustedResetSchedule->length() == 0) ) 
            {
                // no penultimate step smoothing for resettables yet, I'm afraid
                if (inst->DECS == false && inst->PERCS == false)
                    PenultStraight(s, step, bot, top, (p[layer]));
                else
                    PenultMandatory(s, step, bot, top, (p[layer]));
            }
        }
    }

    // calculate relevant stock interval for Black on call and some coefficients
    // initial stock interval to be not-do-Black-on-call-adjustments
    
    //////////////////////////////////////////////////////////////
    //                                                          //
    //  Merge price layers of diff conv ratios to get cvb       //
    //  option values for conversion/call/put decision making   //
    //  Needed for: reset, triggerReset                         //
    //  Ideally should do CoCo (done at end of function)        //
    //                                                          //
    //////////////////////////////////////////////////////////////
    if ( hasReset && (pEnd - pStart) > 0 && (*resetArray)[step] ) 
    {
        mergePriceForReset(s, step, bot, top, pStart, endIdx, p);
    }
    
    if (inst->triggerReset &&
        model->getDate(step).getDate() >= inst->resetObservationStartDate.getDate() &&
        model->getDate(step).getDate() <= inst->resetObservationEndDate.getDate() ) {
        bool isUp = (inst->resetType == "UP_RESET"); 
        for (j=bot; j<=top; ++j) 
        {
            if( isUp == (inst->resetTrigger <=  s[j]) ) (p[pStart])[j] = (p[layerTrigResetRatio])[j];
        }
    }
     
    
    //////////////////////////////////////////////////////////////
    //                                                          //
    //  Pricing for each layer                                  //
    //                                                          //
    //////////////////////////////////////////////////////////////

    // normally, FD calculation is up to callLevel and need to linear extend beyond
    // for certain cases, no upBarrier is used and should not use slope due to potential inaccuracy/instability
    // therefore, ideally should be in-synch with boundary treatment in preCalcFD()
    bool useSlopeAboveCall = true;
    if ( inst->hasContConversion() ) useSlopeAboveCall = false;

    for (layer=pStart; layer<=endIdx; ++layer) 
    {
        // initialize for each loop over spots. used be addOnConvRatios
        double nextAddOnConvTrigger = inst->hasAddOnConvRatios()?(-1):(s[top]*1.1); 

        for (j=bot; j<=top; ++j) 
        {         
            aboveCall = false;
            slope   = 0.0;

            // if not static spread, we replace the bond floor by the appropriate bond floor per node
            if ( isBondFlrGeneric ) 
            {
                bondFloorTmp = bondFloorGeneric[j];
            }

            convertValue = (p[layer])[j] + bondFloorTmp;
            
            //////////////////////////////////////////////////////////////
            //                                                          //
            //  I: conversionValue and conversionValue/Ratio            //
            //  in case of call                                         //
            //                                                          //
            //////////////////////////////////////////////////////////////
            
            if ( !convBool[step] ) 
            {
                conversionValue = 0.0;                
                conversionRatioForCall = 0.0;
                conversionValueForCall = 0.0;
            } 
            // if DEC but not continuousReset, use convRatio[step] (see following logic)
            else if ( inst->DECS && inst->continuousReset ) 
            {
                if (s[j] < inst->initialPrice)
                    conversionValue = s[j] * inst->initialConvRatio;
                else if (s[j] < inst->convPrice)
                    conversionValue = inst->initialPrice * inst->initialConvRatio;
                else if( !inst->decsHasCutoff || s[j] < inst->decsCutoffLevel)
                    conversionValue = inst->initialPrice * inst->initialConvRatio + (s[j] - inst->convPrice)*inst->minConvRatio;
                else
                    conversionValue = inst->initialPrice * inst->initialConvRatio + (inst->decsCutoffLevel - inst->convPrice)*inst->minConvRatio;
                
                // calculate conversion option
                if ( inst->delayReset > 0 )
                    conversionValue = adjConversionValueDECDelay(step, s[j], convertValue, conversionValue);
                
                conversionValue += conversionCash;
                
                conversionRatioForCall = (conversionValue - conversionCash)/s[j];
                conversionValueForCall = conversionValue;
                
            } 
            else if ( hasReset && (pEnd - pStart > 0) ) 
            {
                // currently using linearly interpolated conversion ratios
                conversionRatio = (*conversionRatioGrid)[step][layer];
                conversionValue = conversionRatio*s[j]+conversionCash;
                
                if ( inst->maxWithParity == true ) 
                {
                    double parityCR, parityValue;
                    double  minStrike   = (*minConvPrice)[step];
                    if ( s[j] * inst->manualResetConvPricePct < minStrike) 
                    {
                        parityCR = faceValue / minStrike;
                    } 
                    else 
                    {
                        parityCR = faceValue / (s[j] * inst->manualResetConvPricePct);
                    }
                    parityValue = parityCR * s[j];
                    conversionValue = Maths::max(parityValue,conversionValue);
                }
                
                conversionRatioForCall = conversionRatio;
                conversionValueForCall = conversionValue;
            } 
            else if( inst->triggerReset && layer==layerTrigResetRatio )  
            {
                conversionValue = inst->resetConversionRatio*s[j]+conversionCash;
                conversionRatioForCall = inst->resetConversionRatio;
                conversionValueForCall = conversionValue;
            }
            else 
            {
                // if has addOnConvRatio, update conversionRatio/Cash across addOnConvRatio
                // trigger level. we know that s[j] >0 and isConvertible
                if( inst->hasAddOnConvRatios() && s[j] >= nextAddOnConvTrigger ) 
                {
                    bool isConvertible;
                    inst->getConversionInfo(model->getDate(step), s[j], &isConvertible, &conversionRatio, &conversionCash);
                    nextAddOnConvTrigger = inst->addOnConvRatios->getNextTrigger(s[j], s[top] * 1.1);
                }
                conversionValue = conversionRatio*s[j]+conversionCash;
                conversionRatioForCall = conversionRatio;
                conversionValueForCall = conversionValue;
                
                // apply contingent conversion trigger
                // does not affect conversionValue/Ratio for call purpose
                if( (inst->hasContConversion() && layer==layerCoCoNCV ) ||
                    (inst->contingentConversion && conversionValue < inst->contingentConversionTrigger) )
                {
                    conversionValue = 0.0; 
                }
            }
            
            //////////////////////////////////////////////////////////////
            //                                                          //
            //  II: call level and DEC's call-like features             //
            //                                                          //
            //////////////////////////////////////////////////////////////
            
            if (callBool[step]) 
            {
                if (convBeforeCall[step]) 
                {
                    convertValue = Maths::min(convertValue, callLevel[step]);
                    // first make sure convertValue is below parityFloor
                    if( convertValue < parityFloor[step] * conversionValue )
                        convertValue = Maths::max(conversionValue, convertValue);
                } 
                else 
                {
                    // convert value now given conversionValue and convertValue-if-not-conversion-now
                    // first make sure convertValue is below parityFloor
                    if( convertValue < parityFloor[step] * conversionValue )
                        convertValue = Maths::max(conversionValue, convertValue);
                    
                    // notice that convertValue may be less than conversionValueForCall, eg for certain ContConversion
                    // notice that (p[layer])[j-1] is already sum of (bondfloor + option)
                    if (useSlopeAboveCall && !aboveCall && !Maths::isNegative( convertValue - conversionValueForCall ) &&
                        (inst->notResetCallTrigger? inst->conversionRatios->firstValue()*
                          s[j]:conversionValueForCall) >= callLevel[step] ) 
                    {
                        aboveCall = true;
                        if (j > bot)
                            slope = (callLevel[step] - (p[layer])[j-1])/((callLevel[step]-conversionCash)/conversionRatioForCall-s[j-1]);
                        else
                            slope = conversionRatioForCall;
                    }
                    
                    if (!aboveCall) 
                    {                        
                        // still need to min with call
                        convertValue = Maths::min(convertValue, callLevel[step]);                        
                    }
                    else 
                    {
                        convertValue = callLevel[step] + slope*(s[j] - (callLevel[step]-conversionCash)/conversionRatioForCall);
                    }
                }
                
            }
            // no calls. but may have call-like feature for DECS
            else if ( inst->DECS && inst->accelerateDECS && model->getDate(step) >= inst->accelerateDECSStartDate && 
                    model->getDate(step) <= inst->accelerateDECSEndDate && s[j] >= inst->accelerateDECSTriggerLevel )
            {
                double callValue = calcAccelerateDECCallValue(step, s[j], bondFloorTmp);
                convertValue = Maths::max(conversionValue, Maths::min(convertValue, callValue));
            }   
            else if ( inst->DECS && inst->decsHasCutoff && !inst->triggerStartDate.empty() && model->getDate(step) >= inst->triggerStartDate &&
                s[j] > inst->decsCutoffLevel * inst->cappedDECSTrigger )
            {
                convertValue = inst->initialPrice * inst->initialConvRatio + (inst->decsCutoffLevel - inst->convPrice)*inst->minConvRatio;
            } 
            else 
            { 
                // first make sure convertValue is below parityFloor
                if( convertValue < parityFloor[step] * conversionValue )
                {
                    convertValue = Maths::max(conversionValue, convertValue);
                }
            }

            // add soft put condition
            if ( softPutBool[step] && conversionValue <= softPutTrigger[step] ) 
            {
                convertValue = Maths::max(convertValue, softPutLevel[step]);
            }
            
            //////////////////////////////////////////////////////////////
            // Notice that price array holds (convert = floor + option) //
            //////////////////////////////////////////////////////////////
            (p[layer])[j] = convertValue;
        }  
    }
    
    
    //////////////////////////////////////////////////////////////
    //                                                          //
    //  III: put strike/prob calc and adj price for put level   //
    //                                                          //
    //////////////////////////////////////////////////////////////
    
    // if put to stock, need to calc put value, which does not affect bond floor calculation.
    double putValue = putLevel[step];
    if( putBool[step] == true && !!inst->putToStock ) putValue = calcPutToStockValue(step);

    // check if we want to calculate the put stock levels and put probs here
    // we should really check the control pointer here to avoid unnecessarily duplicating work, but this requires
    // some changes to callbacks - still need to be done.
    // Don't calc if puts a  re American style i.e. if you can put at next step.
    // if contingent conversion, the put level/prob is not unique. so we just use the cv/ncv price depending on whether today is cv or ncv      
    if ((pStart == pEnd || inst->hasContConversion() ) && 
        firstPutIndex == step && !Maths::isPositive(putStockLevel) && putBool[step+1] == false)
        calcPutLevelNProb(step, bot, top, s, (p[pStart]), putValue, putStockLevel, putProbability);
    
    // max with put - this has to be done after the put probability has been calculated
    // also adjust the bond floor by put level
    int i;
    if (putBool[step] == true ) {
        bondFloorTmp = Maths::max(bondFloorTmp, putLevel[step]);
        if ( isBondFlrGeneric )
        {
            for (j=bot; j<=top; ++j)
                bondFloorGeneric[j] = Maths::max(bondFloorGeneric[j], putLevel[step]);
        }
        
        for (i=endIdx; i>=pStart; --i) 
        {
            for (j=bot; j<=top; ++j)
            {
                (p[i])[j] = Maths::max((p[i])[j], putValue);
            }
        }
    }
    
    // adjust pricing if contingent conversion, not done at start of function because
    // it may use conversion info from the previous EOD and currently we can not place
    // time line points at EOD to avoid conflict with SOD time points from conversion schedule
    // (closeby points may be erased during timeline construction and causing problem)
    if( inst->hasContConversion() )
    {
        inst->contConversion->adjustPrices(model->getDate(step), top - bot, 
            &s[bot], &(p[layerCoCoNCV])[bot], &(p[layerCoCoCV])[bot], callBool[step]?callLevel[step]:(-1), inst);
    }
        
    //////////////////////////////////////////////////////////////
    //                                                          //
    //  IV: OCB                                                 //
    //                                                          //
    //////////////////////////////////////////////////////////////
    
    if (isOCB && pEnd > pStart) 
    { // we're an OCB
        for (i=numCVBPriceArrays ; i<= pEnd ; ++i) 
        {
            for (j=bot; j<=top; j++) 
            {
                (p[i])[j] += ocbRiskyAdj[step];
            }
            
            if (ocbExerBool[step] == true) 
            {
                for (j=bot; j<=top; j++)
                {    
                    (p[i])[j] = Maths::max(((p[i-numCVBPriceArrays])[j] - ocbStrike[step]), (p[i])[j]);
                }               
                // be careful -- we don't want to add value for bonds that have already been called so min with call 
                // unless we're going to convert just before. Also be careful not to subtract value when the strike 
                // is above the call level. Avoid this by never MINing with a negative number.
                if (convBeforeCall[step] == false) 
                {
                    double callMinusStrike = Maths::max(callLevel[step] - ocbStrike[step], 0.);
                    double conversionRatio;
                    if (hasReset) 
                    {
                        conversionRatio = (*conversionRatioGrid)[step][i-numCVBPriceArrays];
                    } 
                    else 
                    {
                        conversionRatio = convRatio[step];
                    }
                    
                    for (j=bot; j<=top; j++) 
                    {   
                        if (conversionRatio*s[j]+conversionCash >= callLevel[step]) 
                        {
                            (p[i])[j] = Maths::min((p[i])[j], callMinusStrike);
                        }
                    }
                }
            }
        }
    }
    
    //////////////////////////////////////////////////////////////
    //                                                          //
    //  IV: Remove bond floor                                   //
    //                                                          //
    //////////////////////////////////////////////////////////////
    
    for (i=pStart; i<=endIdx; ++i) 
    {
        for (j=bot; j<=top; ++j)
        {
            if( isBondFlrGeneric ) bondFloorTmp = bondFloorGeneric[j];
            (p[i])[j] -= bondFloorTmp;
        }
    }

    // for Fwd Starting, addint the coupon before starting date but after issue date.
    if (step ==0 && fwdStarting)
    {
        for (i=pStart; i<=endIdx; ++i) 
        {
            for (j=bot; j<=top; ++j)
                (p[i])[j] += pvCFbeforeStart;
        }
    }   
}

void ConvBondFDProd::mergePriceForReset(const double* s, 
                                              int step, 
                                              int bot, 
                                              int top, 
                                              int pStart, 
                                              int pEnd,
                                              const vector< double * > & p)
{
    if ( priceCopy.numCols() == 0 ) 
    {
        priceCopy = CDoubleMatrix(pEnd+1, top+bot+1);
    }
    
    // back up the input prices 
    int layer, j;
    for(layer= 0;layer<=pEnd;++layer) 
    {
        for(j=bot; j<=top; ++j) 
        {
            priceCopy[layer][j] = (p[layer])[j];
        }
    }
    
    double  maxStrike   = (*maxConvPrice)[step];
    double  minStrike   = (*minConvPrice)[step];
    double  resetPremium  = (*resetParity)[step];
    
    for (layer=pStart; layer<=pEnd; ++layer) 
    {    
        double conversionRatio = (*conversionRatioGrid)[step][layer];
        
        for (j=bot; j<=top; ++j) 
        {           
            // need to cross over the grids depending on current spot here
            int lowIdx  = 0;
            int highIdx = 0;
            double targetCR;
            
            double resetRatio;
            if ( resetPremium * s[j] < minStrike ) 
            {
                resetRatio = faceValue / minStrike;
            } 
            else if ( resetPremium * s[j] > maxStrike ) 
            {
                resetRatio = faceValue / maxStrike;
            } 
            else 
            {
                if ( step == 0 ) 
                {
                    resetRatio = adjustedResetSchedule->getCurrentConversionRatio(
                        convRatio[step],
                        model->getDate(step),
                        faceValue);
                } 
                else 
                {
                    resetRatio = faceValue / (s[j] * resetPremium);
                }
            }
            
            string resetTypeStr = adjustedResetSchedule->getResetType(model->getDate(step));
            if (  resetTypeStr == "UP_DOWN" ) 
            {
                // UP_DOWN
                targetCR   = resetRatio;
            } 
            else 
            {
                if (adjustedResetSchedule->isFlat(model->getDate(step))) 
                {
                    targetCR   = Maths::max(Maths::max(resetRatio,conversionRatio),convRatio[step]);
                } 
                else 
                {
                    targetCR   = Maths::max(resetRatio,conversionRatio);
                }
            }
            
            
            // new way of setting up conversion ratios
            double cr;
            int    idx = 0;
            double lowCR, highCR;
            while (idx < pEnd ) 
            {
                cr = (*conversionRatioGrid)[step][idx];
                if (targetCR >= cr) 
                {
                    lowIdx = idx;
                    lowCR  = cr;
                    if ( targetCR == cr ) 
                    {
                        highIdx = idx;
                        highCR  = cr;
                    } 
                    else 
                    {
                        highIdx = idx+1;
                        highCR  = (*conversionRatioGrid)[step][highIdx];
                    }
                }
                ++idx;
            }
            double weight = (highCR != lowCR)?((targetCR - lowCR) / ( highCR - lowCR)):1.0;
            (p[layer])[j]  = weight * priceCopy[highIdx][j] + (1.-weight) * priceCopy[lowIdx][j];            
        }
    }
}

// the input price holds (bondfloor + option)
void ConvBondFDProd::calcPutLevelNProb(int step, 
                                       int bot, 
                                       int top, 
                                       const double *s, 
                                       const double *priceLayer, 
                                       double putValue, 
                                       double &putStockLevel, 
                                       double &putProbability)
{
    static const string method = "ConvBondFDProd::calcPutLevelNProb";
    
    // don't fail if this doesn't work
    try 
    {
        putStockLevel  = 0.;
        putProbability = 1.0;
        if ( isStaticSpread() && putValue < bondFloor[step]) 
        {
            putProbability = 0.0;
        } 
        else if (Maths::areEqualWithinTol(callLevel[step],putValue,1.e-11)) 
        {
            if( convBool[step] )
            {
                if( inst->hasAddOnConvRatios() )
                    putStockLevel = inst->getImpliedStrike(model->getDate(step), putValue);
                else
                    putStockLevel = (putValue-convCash[step])/convRatio[step];
            }
        } 
        else if (isStaticSpread() && putValue > priceLayer[top]) 
        {
            putProbability = 1.0;
        } 
        else if (isStaticSpread() && putValue <= priceLayer[bot] ) 
        {
            putProbability = 0.0;
        } 
        else 
        {
            for (int idx=top-1; idx >= bot; --idx) 
            {              
                if ( priceLayer[idx] <= putValue ) 
                {
                    if (!Maths::isZero((priceLayer[idx+1] - priceLayer[idx]))) 
                    {
                        double a;
                        a = (priceLayer[idx+1] - putValue)/(priceLayer[idx+1] - priceLayer[idx]);
                        putStockLevel = (1-a)*s[idx+1] + a*s[idx];
                    } 
                    else 
                    {
                        putStockLevel = s[idx+1];
                    }
                    break;
                }
            }
        }
        
        if (putStockLevel > 0.) 
        {
            // create vol request and get processed vol
            CVolRequestConstSP volRequest = GetLNRequest();
            CVolProcessedSP  procVol(getEquityAsset()->getProcessedVol(volRequest.get()));
            // cast to the type of vol we're expecting
            CVolProcessedBSSP volBS = CVolProcessedBSSP::dynamicCast(procVol);
            // calculate the variance
            double variance = volBS->CalcVar(inst->getValueDate(),
                model->getDate(step));
            
            // get the forward at maturity - could possibly pull this one out the the FWD_AT_MAT
            // results if some of the calculations are done in recordOutputRequests
            double fwd = getEquityAsset()->fwdValue(model->getDate(step));
            if ( !Maths::isZero(putStockLevel)) 
            {
                putProbability = 1. - N1((log(fwd/putStockLevel) - (0.5 * variance))/sqrt(variance));
            } 
            else 
            {
                putProbability = 0.0;
            }
        }        
        calcedFirstPutProb = true;
    } 
    catch (exception& e) 
    {
        putStockLevel = 999.;
        ModelException myE = ModelException(e);
        myE.addMsg(method + ": Failed to compute put stock level. Continuing...\n");
        myE.errorLog();
    }
}

// compute value of put if is put to stock but not already put. 
// code adapted from CorpAct.cpp
double ConvBondFDProd::calcPutToStockValue(int step)
{ 
    DateTime putDate = model->getDate(step);
    SampleListSP samples = inst->getPutToStockSample(putDate);

    // make sure the first date is past value date, otherwise, can not put!
    if( samples->getFirstDate().getDate() < inst->getValueDate().getDate() ) return 0.0;

    double factor = inst->calcPutToStockFactor(putDate, samples.get(), GetLNRequest().get());
    return putLevel[step] * factor;
}

double ConvBondFDProd::adjConversionValueDECDelay(int step, double s, double convertValue, double conversionValue)
{
    // price the call
    double convRatio = conversionValue / s;
    DateTime fwdDate   = model->getDate(step).rollDate(inst->delayReset);
    double   discFact  = inst->discount->pv(model->getDate(step), fwdDate);
    double   strike    = convertValue / convRatio;
    double   fwd       = getEquityAsset()->fwdFwd(model->getDate(step),s,fwdDate);
    double dummy = s;
    vector<double> fwdVol;
    	if (tree1f){
		tree1f->GetStepVol(step, fwdVol, &dummy,0,0);
	} else if (fd1dMulti){
		fd1dMulti->GetStepVol(step, fwdVol, &dummy, 0, 0);
	}
    double variance = fwdVol[0] * fwdVol[0] * (double(inst->delayReset) / 365.25);
    double convOption = conversionOption(true, model->getDate(step), fwdDate, fwd, strike, variance, discFact);
    convOption *= convRatio;
    conversionValue += convOption;
    
    // price the put spread if necessary
    if ( inst->includePutSpread && convRatio * s > convertValue )
    {
        double highStrike = s;
        double lowStrike   = convertValue / convRatio;
        double put1 = conversionOption(false /*put*/, model->getDate(step), fwdDate, fwd, highStrike, variance, discFact);
        double put2 = conversionOption(false /*put*/, model->getDate(step), fwdDate, fwd, lowStrike, variance, discFact);
        
        conversionValue += (put2 - put1) * convRatio;
    }
    
    return conversionValue;
}

double ConvBondFDProd::calcAccelerateDECCallValue(int step, 
                                                  double s, 
                                                  double bondFloor)
{
    // accelerate maturity feature
    double callValue = bondFloor;
    
    if ( inst->dividendPassThrough ) 
    {
        // subtract the pv of the dividends
        CashFlowArraySP dividendCashFlows(new CashFlowArray(0));
        
        DividendListSP divSchedule = DividendCollector::divsBetweenDates(getEquityAsset().get(),
            model->getDate(step),
            model->getDate(step),
            inst->bond->getMaturityDate());
        
        const DividendArray& divArray = divSchedule->getArray();
        
        int k;
        for (k=0 ; k<divArray.size() ; ++k) 
        {
            callValue -= divArray[k].getDivAmount() * inst->initialConvRatio * inst->dividendPassThroughPct * 
                inst->discount->pv(model->getDate(step), divArray[k ].getPayDate());
        }
    }
    
    if (s < inst->initialPrice) 
    {
        callValue += s * inst->initialConvRatio;
    } 
    else if (s < inst->convPrice)
    {
        callValue += inst->initialPrice * inst->initialConvRatio;
    } 
    else 
    {
        callValue += inst->initialPrice * inst->initialConvRatio + (s - inst->convPrice)*inst->minConvRatio;
    }
    
    return callValue;
}

/** Penultimate step smoothing. Should do a better job. At a minimum it should check for the call level.
    Could do an up and out call with a rebate. The old model did a call spread which is probably
    good enough for most purposes. Was probably overkill there as it would smooth the calls anyway. */
void ConvBondFDProd::PenultStraight(const double* s, 
                                          int step, 
                                          int bot, 
                                          int top, 
                                          double *price)
{
    int j;
    double fwdPrice;
    double strike = bondFloor[step+1]-convCash[step+1];
    vector <double> vol_arr, drift_arr;
    int nn;
	double dt, truncationStd;
	if (tree1f) {
		dt = tree1f->getTrdYrFrac(step+1);
		nn = tree1f->GetStepVol(step, vol_arr, s, bot, top);
		truncationStd = tree1f->TruncationStd;
	} else if (fd1dMulti){
		dt = fd1dMulti->getTrdYrFrac(step+1);
		nn = fd1dMulti->GetStepVol(step, vol_arr, s, bot, top);
		truncationStd = fd1dMulti->getTruncationStd();
	}
    if (nn > 1){
        return;
//        throw ModelException("ConvBondFDProd::PenultStraight", "Penultimate step smoothing not implemented for local vol");
    }
    double variance = vol_arr[0]*vol_arr[0]*dt;
    double df = inst->discount->pv(model->getDate(step), model->getDate(step+1));

    double drift = getEquityAsset()->fwdValue(model->getDate(step+1))/
                   getEquityAsset()->fwdValue(model->getDate(step));

    double conversionRatio = convRatio[step+1], conversionCash = convCash[step+1];

    for (j=bot; j<=top; j++)
    {
        if (s[j]>0.0) 
        {
            fwdPrice = s[j]*drift;
            
            if( convBool[step] && inst->hasAddOnConvRatios() )
            {
                bool isConvertible;
                inst->getConversionInfo(model->getDate(step+1), fwdPrice, &isConvertible, &conversionRatio, &conversionCash);
                strike =  bondFloor[step+1]-conversionCash;
            }
            if (fabs(log(fwdPrice*conversionRatio/strike)) < truncationStd*sqrt(variance))
            {
                price[j] = conversionRatio*Black::price(true, fwdPrice, 
                    strike/conversionRatio, df, variance);
            }
        }
    }
}

/** Penultimate step smoothing for mandatory convertibles. This is a distinct routine
    from the straight routine even though there is a lot of overlap because I anticipate 
    improving the straight routine. This one should be fine for practical purposes
    as it is. */
void ConvBondFDProd::PenultMandatory(const double* s, 
                                           int step, 
                                           int bot, 
                                           int top, 
                                           double *price)
{
    int j;
    double fwdPrice;

    vector <double> vol_arr, drift_arr;
	double dt, truncationStd;
	int nn;
	if (tree1f){
		dt = tree1f->getTrdYrFrac(step+1);
		nn = tree1f->GetStepVol(step, vol_arr, s, bot, top);
		truncationStd = tree1f->TruncationStd;
	} else if (fd1dMulti){
		dt = fd1dMulti->getTrdYrFrac(step+1);
		nn = fd1dMulti->GetStepVol(step, vol_arr, s, bot, top);
		truncationStd = fd1dMulti->getTruncationStd();
	}
    /*
    if (nn > 1){
        throw ModelException("ConvBondFDProd::PenultMandatory", "Penultimate step smoothing not implemented for local vol");
    }
    */
    double variance = vol_arr[0]*vol_arr[0]*dt;
    double df = inst->discount->pv(model->getDate(step), model->getDate(step+1));

    double drift = getEquityAsset()->fwdValue(model->getDate(step+1))/
                   getEquityAsset()->fwdValue(model->getDate(step));
    for (j=bot; j<=top; j++)
    {
        if (s[j]>0.0) {
            fwdPrice = s[j]*drift;

            if (nn > 1)
            {
                variance = vol_arr[j-bot]*vol_arr[j-bot]*dt;
            }


            if (inst->DECS == true) 
            {
				if (fabs(log(fwdPrice/inst->initialPrice)) < truncationStd*sqrt(variance) ||
					fabs(log(fwdPrice/inst->convPrice)) < truncationStd*sqrt(variance))
                {
                    price[j] = inst->initialConvRatio*fwdPrice*df;

                    price[j] -= inst->initialConvRatio*Black::price(true, fwdPrice, 
                        inst->initialPrice, df, variance);

                    price[j] += inst->minConvRatio*Black::price(true, fwdPrice, 
                        inst->convPrice, df, variance);
                }
            } 
            else 
            { // inst->PERCS == true
                if (fabs(log(fwdPrice/inst->initialPrice)) < truncationStd*sqrt(variance))
                {
                    price[j] = inst->initialConvRatio*fwdPrice*df;

                    price[j] -= inst->initialConvRatio*Black::price(true, fwdPrice, 
                        inst->initialPrice, df, variance);
                }
            }
        }
    }
}

double ConvBondFDProd::getCoupon(int step, 
                                 const double* s,  
                                 int start, 
                                 int end)
{
    // this needs to get implemented for the E2C model
    double accruedToday = (inst->payAccruedUponDefault)?inst->getAccruedAtDate(model->getDate(step)):0.0;
    double recovery     = inst->useAssetRecovery?(faceValue * inst->recoveryPct):0.0;

    return recovery + accruedToday * inst->recoveryPct;
}

bool ConvBondFDProd::hasEquityLayer()
{
    return true;
}

/** premium scaling */
double ConvBondFDProd::scalePremium(vector<double> & P, 
                                     YieldCurveConstSP disc)
{    
    if (isOCB == false) 
    {
        // return fairValue + bondFloor[0];
        double bondFloorTmp;
        if( isStaticSpread() )         
            bondFloorTmp = bondFloor[0];
        else
        {
            int numPriceArray = (isOCB?2:1)*numCVBPriceArrays; // # if price array excl of the riskybond
            bondFloorTmp = P[numPriceArray];
        }
        return (P[0] + Maths::max(bondFloorTmp, putLevel[0])) ;
    } 
    else 
    {        
        return ( P[numCVBPriceArrays] ) ;
    }
}

/** extra output requests */
void ConvBondFDProd::recordOutput(Control*            control, 
                                  YieldCurveConstSP   disc, 
                                  Results*            results)
{
    static const string method = "ConvBondFDProd::recordOutput";
    OutputRequest* request = NULL;

    // get prices at t=0
    int size = slices.size();
    vector< double > price0( size );
    for( int i = 0; i < size; ++i )
        price0[i] = model->getPrice0( *slices[i] );

    // save price
    double price = scalePremium(price0, disc);
    results->storePrice(price, disc->getCcy());

    if ( control && control->isPricing() ) 
    {
        double nakedBond = 0.0;        
        nakedBond = bondFloor[0]+pvCFbeforeStart;
        if (putBool[0] == true) 
        {
            nakedBond = Maths::max(nakedBond, putLevel[0]);
        }
        
        inst->recordOutputRequests(control, results, price, nakedBond, isOCB);

        // record the bond floor as calculated by the no static spread
        if( !isStaticSpread() &&
            control->requestsOutput(OutputRequest::NAKED_BOND_PRICE2, request) )
        {
            int numPriceArray = (isOCB?2:1)*numCVBPriceArrays; // # if price array excl of the riskybond
            double bondFloorTmp = model->getPrice0( *slices[numPriceArray] );
            results->storeRequestResult(request, Maths::max(bondFloorTmp, putLevel[0]));
        }

        // Indicative vol
        if (control->requestsOutput(OutputRequest::IND_VOL, request)) 
        {
            double             indVol     = 0.0;
            CVolRequestConstSP volRequest   = GetLNRequest();
            CVolRequestSP      ncVolRequest =CVolRequestSP::constCast(volRequest);

            LinearStrikeVolRequest* lsVolRequest = dynamic_cast<LinearStrikeVolRequest*>(ncVolRequest.get());
            // interpolate the vol using our LN request
            CVolProcessedBSSP volBS(getEquityAsset()->getProcessedVol(lsVolRequest));
            // calculate the indicative vol
            try 
            {
                indVol = volBS->CalcVol(inst->valueDate, inst->bond->getMaturityDate());
            }
            catch (exception& ) 
            {
                indVol = 0.0;
            }
            results->storeRequestResult(request, indVol);
        }

        // Duration
        if ( control->requestsOutput(OutputRequest::BOND_DURATION, request) ) 
        {
            results->storeRequestResult(request, inst->bond->duration(inst->valueDate,risky));
        }
        // Convexity
        if ( control->requestsOutput(OutputRequest::BOND_CONVEXITY, request) ) 
        {
            results->storeRequestResult(request, inst->bond->convexity(inst->valueDate,risky));
        }


        if (isOCB == true) 
        {
            ocb->recordOutputRequests(control, results, price0[0] + Maths::max(bondFloor[0], putLevel[0]), bondFloor[0], price);
        } 
        else 
        {
            // put probability
            if ( firstPutIndex >= 0 && control->requestsOutput(OutputRequest::PUT_PROBABILITY, request)) 
            {
                if (calcedFirstPutProb == true) 
                {
                    results->storeRequestResult(request, putProbability);
                } 
                else 
                {
                    string errorMsg = "could not calculate put probability";
                    UntweakableSP untweakablePutProb(new Untweakable(errorMsg));
                    results->storeRequestResult(request,untweakablePutProb);
                }
            }

            // put stock level
            if ( firstPutIndex >= 0 && control->requestsOutput(OutputRequest::PUT_STOCK_LEVEL, request)) 
            {
                if (calcedFirstPutProb == true) 
                {
                    results->storeRequestResult(request, putStockLevel);
                } 
                else 
                {
                    string errorMsg = "could not calculate put stock level";
                    UntweakableSP untweakablePutStock(new Untweakable(errorMsg));
                    results->storeRequestResult(request,untweakablePutStock);
                }
            }

            // put strike
            if ( firstPutIndex >= 0 && control->requestsOutput(OutputRequest::PUT_STRIKE, request) ) 
            {
                double firstPutLevel = putLevel[firstPutIndex];
                double putConvRatio  = (convBool[firstPutIndex])?convRatio[firstPutIndex]:0.0;
                double putConvCash   = (convBool[firstPutIndex])?convCash[firstPutIndex]:0.0;
                double spotFX        = 1.0;
                double putStrike;

                if( convBool[firstPutIndex] && inst->hasAddOnConvRatios() )
                {    
                    bool isConvertible;
                    inst->getConversionInfo(model->getDate(firstPutIndex), firstPutLevel, &isConvertible, &putConvRatio, &putConvCash, false/*is not spot*/);
                }

                if ( Maths::isPositive(firstPutLevel) && Maths::isPositive(putConvRatio) ) 
                {
                    if (StruckEquity::TYPE->isInstance(getEquityAsset().get())) 
                    {
                        // cast to struck equity object
                        const IObject* obj           = dynamic_cast<const IObject*>(getEquityAsset().get());
                        const StruckEquity* struckEq = dynamic_cast<const StruckEquity*>(obj);
                        if ( struckEq ) 
                        {
                            spotFX = struckEq->getFX()->getSpot();
                        }
                    }
                    putStrike = (firstPutLevel-putConvCash)/(putConvRatio*spotFX);
                    results->storeRequestResult(request,putStrike);
                } 
                else 
                {
                    UntweakableSP untweakablePutStrike(new Untweakable(
                        "firstPutLevel(" + 
                        Format::toString(firstPutLevel) + 
                        ") and conversion ratio(" +
                        Format::toString(putConvRatio) + 
                        ") must be strictly positive at first put date"));

                    results->storeRequestResult(request,untweakablePutStrike);
                }
            }

            // call strike
            if ( control->requestsOutput(OutputRequest::CALL_STRIKE, request) ) 
            {
                if ( !(!inst->callSchedule) && inst->callSchedule->length() > 0 ) 
                {
                   DateTime firstCallDate;
                   if ( inst->callSchedule->getNextDate(inst->valueDate,
                                                       inst->bond->getMaturityDate(),
                                                       firstCallDate)) 
                   {
                       // AS: should possibly put in the correct fwd vol here
                      bool    isCallable, canConvert, isAdjustedForAccrued;
                      double  firstCallLevel, callConvRatio, callConvCash;
                      double spotFX        = 1.0;
                      double callStrike;
                      bool  isHardCall;

                      inst->getCallLevel(firstCallDate, 0.0, &isCallable, &isAdjustedForAccrued, &isHardCall, &firstCallLevel);

                      bool isSpot = !inst->hasAddOnConvRatios();
                      double spotOrLevel = isSpot?0:firstCallLevel;
                      inst->getConversionInfo(firstCallDate, spotOrLevel, &canConvert, &callConvRatio, &callConvCash, isSpot);

                      if ( Maths::isPositive(firstCallLevel) && Maths::isPositive(callConvRatio) ) 
                      {
                         if ( StruckEquity::TYPE->isInstance(getEquityAsset().get())) 
                         {
                            // cast to struck equity object
                            const IObject* obj           = dynamic_cast<const IObject*>(getEquityAsset().get());
                            const StruckEquity* struckEq = dynamic_cast<const StruckEquity*>(obj);
                            if ( struckEq ) 
                            {
                               spotFX = struckEq->getFX()->getSpot();
                            }
                         }
                         callStrike = (firstCallLevel-callConvCash)/(callConvRatio*spotFX);
                         results->storeRequestResult(request,callStrike);
                      } 
                      else 
                      {
                         UntweakableSP untweakableCallStrike(new Untweakable(
                            "firstCallLevel(" + 
                            Format::toString(firstCallLevel) + 
                            ") and conversion ratio(" +
                            Format::toString(callConvRatio) + 
                            ") must be strictly positive at first call date"));

                            results->storeRequestResult(request,untweakableCallStrike);
                       }
                   }
               }
            }
        }
    }
}

/** returns a vol request for log-normal vol */
CVolRequestConstSP ConvBondFDProd::GetLNRequest() const
{ 
    CVolRequestConstSP volRequest(
        new LinearStrikeVolRequest(inst->getStrike(), 
                                   inst->valueDate,
                                   inst->bond->getMaturityDate(), 
                                   false));
  
    return volRequest;
}

// adjust risky bond floor price layers, assign risky bond floor to bondFloorGeneric
// return the new pEnd for non-risky-bond related price layers.
// notice that the address of bondFloorGeneric is assigned here as well.
int ConvBondFDProd::calcRiskyBondNAdjPriceLayer(int step, 
                                                int bot, 
                                                int top, 
                                                int pEnd, 
                                                const vector< double * > & p,
                                                double *&bondFloorGeneric)
{
    bondFloorGeneric = bondFloorGenericMem - bot;
    if( !isStaticSpread() )
    {
        calcBondRiskyNFloorBeforeMat(step, bot, top, (p[pEnd]), bondFloorGeneric);
        pEnd--; // make sure to decrease the number of layer to return the remaining layers
    }
    else // isE2CLayer
    {
        // this calculation is inheritted from prior version and is flawed
        // mainly it missed discouting of coupon for 1 time step
        // it does not contain other logic related to riskFree etc
        int layerBond = pEnd + 2;
        for(int layer = pEnd+1; layer <= pEnd+2; layer++)
        {
            for(int j=bot; j<=top; j++)
            {
                // check coupon paid today ...
                if (addCouponBool[step+1] == true)
                    (p[layer])[j] += coupon[step+1];

                // check callable ... this should only be done for hard calls
                if (callBool[step] && callIsHard[step])
                    (p[layer])[j] = Maths::min((p[layer])[j], callLevel[step]);
        
                if( layer==layerBond ) bondFloorGeneric[j] = (p[layer])[j];
    
                // check puttable ....
                if (putBool[step])
                    (p[layer])[j] = Maths::max((p[layer])[j], putLevel[step]);
            }
        }
    }
    return pEnd;
}


// update bondRisky and bondRiskFree
void ConvBondFDProd::calcBondRiskyNFloorBeforeMat(int i/*step*/, 
                                                  int bot, 
                                                  int top, 
                                                  double *bondRiskyArray, 
                                                  double *bondFloorArray)
{
    for(int j=bot; j<=top; j++)
    {
        // use the non-equity-participation spread to further discount the riskybond
        // this simple treatment requires the bond to have 0 recovery
        // use PV recovery rule for preferreds
        double &bondRisky = bondRiskyArray[j];
        if(  nonEqCreditPvs )  bondRisky *= nonEqCreditPvs[i];

        double bondRiskFree = bondRiskFreeArray[i];

        if (callIsHard[i] && callLevel[i] < (bondRisky + bondRiskFree) )
        {
            bondRisky = Maths::max((double)callLevel[i]-bondRiskFree, 0.);
        }
        if(bondFloorArray) bondFloorArray[j] = bondRisky + bondRiskFree;

        if( i>0 )
        {
            if ( putBool[i] == true &&
                !( model->getDate(i-1).getDate() == 
                model->getDate(model->getLastStep()-1).getDate()) ) 
            {
                bondRisky = Maths::max(bondRisky, putLevel[i]);
            }

            if (addCouponBool[i] == true) {
                if( !(model->getDate(i) <= inst->riskFreeCouponEndDate) )
                    bondRisky += coupon[i];
            }
        }

    }
    return;
}

void ConvBondFDProd::update(int& step, 
                            FDProduct::UpdateType type)
{
    // we assume just need one und level for spot here
    const TreeSlice & s = payoffIndex->getValue( step );
    int bot, top;
    s.getCalcRange( bot, top );

    const vector< TreeSliceSP > & price = slices;
    int pStart = 0, pEnd = price.size() - 1;

    if (type == FDProduct::BWD_T)
    {
        prod_BWD_T(s,
                   step,
                   bot,
                   top,
                   pStart, 
                   pEnd,
                   price);

        //insert nodes
        if (tree1f && tree1f->NumOfInsertNode>0)
        {
            prod_BWD_T(*insNodes,
                        step,
                        0,
                        tree1f->NumOfInsertNode-1,
                        pStart, 
                        pEnd,
                       *insPrices);
        }
    }
    else if(type == FDProduct::BWD)
    {
        prod_BWD(s,
                 step,
                 bot,
                 top,
                 pStart, 
                  pEnd,
                 price);

        //insert nodes
        if (tree1f && tree1f->NumOfInsertNode>0)
        {
            prod_BWD(*insNodes,
                      step,
                      0,
                       tree1f->NumOfInsertNode-1,
                      pStart, 
                      pEnd,
                     *insPrices);
        }
    }   
};



FDProductSP OptOnConvBond::createProduct(FDModel* model) const
{
    ConvBondFDProd * treeProd = new ConvBondFDProd( cvb.get(), model, true );

    treeProd->ocb = this;

    return FDProductSP( treeProd );
}

DRLIB_END_NAMESPACE
