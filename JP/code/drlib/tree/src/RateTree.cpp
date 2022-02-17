//----------------------------------------------------------------------------
//
//   Group       : Global QR&D
//
//   Filename    : RateTree.cpp
//
//   Description : top level model class for rates trees
//
//----------------------------------------------------------------------------

#include "edginc/config.hpp"
#include "edginc/RateTree.hpp"
#include "edginc/RatesUtils.hpp"
#include "edginc/IndexSpecIR.hpp"
#include "edginc/IndexSpecFX.hpp"
#include "edginc/IndexSpecEQ.hpp"
#include "edginc/ModelFilter.hpp"
#include "edginc/IrConverter.hpp"
#include "edginc/IRSmileMQ.hpp"
#include "edginc/VolRequestRaw.hpp"
#include "edginc/VolSurface.hpp"
#include "edginc/SRMEQVol.hpp"


DRLIB_BEGIN_NAMESPACE

START_PUBLIC_ENUM_DEFINITION(ZeroInterpStyle::Enum, "");
ENUM_VALUE_AND_NAME(ZeroInterpStyle::LINEAR, "LINEAR", "Linear interpolation");
ENUM_VALUE_AND_NAME(ZeroInterpStyle::FLAT_FWD, "FLAT_FWD", "Flat forward interpolation");
END_ENUM_DEFINITION(ZeroInterpStyle::Enum);

START_PUBLIC_ENUM_DEFINITION(ZeroBankMode::Enum, "");
ENUM_VALUE_AND_NAME(ZeroBankMode::ZEROBANK, "ZEROBANK", "");
ENUM_VALUE_AND_NAME(ZeroBankMode::CLAIMBANK, "CLAIMBANK", "");
END_ENUM_DEFINITION(ZeroBankMode::Enum);


/****************************** IndexSpecIRTree ********************************/

class IndexSpecIRTree : public FDProduct {
public:
    IndexSpecIRConstSP  indexSpecIR;
    RateTree*           rateTree;
    DateTimeArray       modelResetDates;
    DateTimeArray       lateResetDates;

    /******************** methods ********************/
    IndexSpecIRTree(const IndexSpecIRConstSP &indexSpecIR, RateTree* rateTree);
    virtual DateTime getStartDate(void) const { return model->getDate(0);}

    // FDProduct virtual methods - begin
    virtual void addModelResetDates(const DateTimeArray &modelResetDates, const DateTimeArray &lateResetDates);
    virtual void init(Control*) const;
    virtual void initProd(void);
    virtual void update(int& step, UpdateType type) {}
    virtual const TreeSlice& getValue(int step, DateTime eventDate) const;
    virtual void recordOutput(Control* ctrl, YieldCurveConstSP, Results* results) {}
    // FDProduct virtual methods - end

    void calculateNamedZeroBank1(RateTree::NamedZeroBank& zeroBank) const;
    void calculateNamedZeroBank2(RateTree::NamedZeroBank& zeroBank) const;
    void calculateNamedZeroBank2bis(RateTree::NamedZeroBank& zeroBank) const;
    void calculateNamedZeroBank3(RateTree::NamedZeroBank& zeroBank) const;
    void calculateNamedZeroBankExplicit(RateTree::NamedZeroBank& zeroBank) const;
    void calculateNamedZeroBankCB(RateTree::NamedZeroBank& zeroBank) const;
private:
    TreeSliceSP slice;
};

IndexSpecIRTree::IndexSpecIRTree(const IndexSpecIRConstSP &indexSpecIR, RateTree* rateTree) :
    FDProduct(rateTree), indexSpecIR(indexSpecIR), rateTree(rateTree)
{}

void IndexSpecIRTree::addModelResetDates(
    const DateTimeArray &modelResetDatesP, const DateTimeArray &lateResetDatesP)
{
    modelResetDates.insert(modelResetDates.end(), modelResetDatesP.begin(), modelResetDatesP.end());
    lateResetDates.insert(lateResetDates.end(), lateResetDatesP.begin(), lateResetDatesP.end());
}

void IndexSpecIRTree::init(Control*) const{

    if (rateTree->getZeroBankMode() == ZeroBankMode::ZEROBANK) {

        RateTree::NamedZeroBank zeroBank(indexSpecIR->getName());
        if (modelResetDates.size()>=1) {
            switch (indexSpecIR->zeroBankMethod) {
                case IndexSpecIR::NamedZeroBankType::STANDARD:
                case IndexSpecIR::NamedZeroBankType::SWAPTION_T: 
                    calculateNamedZeroBank1(zeroBank); 
                    break;
                case IndexSpecIR::NamedZeroBankType::CALLTURFLOWS_T:
                    calculateNamedZeroBank2(zeroBank); 
                    break;
                case IndexSpecIR::NamedZeroBankType::RIBOBS:
                    calculateNamedZeroBank2bis(zeroBank); 
                    break;
                case IndexSpecIR::NamedZeroBankType::CALLFXRIB_T:
                    calculateNamedZeroBank3(zeroBank); 
                    break;
                case IndexSpecIR::NamedZeroBankType::EXPLICIT:
                    calculateNamedZeroBankExplicit(zeroBank);
                    break;
                case IndexSpecIR::NamedZeroBankType::CLAIMBANK:
                    calculateNamedZeroBankCB(zeroBank);
                    break;
                default: 
                    throw ModelException("Unsupported zeroBankMethod for IndexSpec "+
                                        indexSpecIR->getName());
            }
        }
        YieldCurveConstSP yc = indexSpecIR->factor.getSP();
        zeroBank.currency = yc->getCcy();
        zeroBank.curveName = yc->getName();

        rateTree->insertNamedZeroBank(zeroBank);
        model->addCritDates(zeroBank.zeroMatDateList);
    }
    else {
        int i;
        string curveName = indexSpecIR->factor.getSP()->getName();

        // ??? for now, hardwiring claimBank zero dates as critical
        // somehow product has to be able to set this flag so it will probably end up
        // as a field in indexSpecIR, but it is a shame to to export that directly as it
        // is a model related parameter rather than instrument related.  if we put it in 
        // the model it is one size fits all - can't specify some components to make zeros
        // critical and other components to make non critical without sending the info down
        // the instrument dependency tree
        for (i = 0; i < modelResetDates.size(); i++) {
            rateTree->insertIRIndex(*indexSpecIR, modelResetDates[i], lateResetDates[i], true);
        }
    }
}

void IndexSpecIRTree::initProd(void){
    string curveName = indexSpecIR->getFactor()->getName();
    slice = model->createSlice();
    slice->name = "IndexSpecIR";
}

const TreeSlice & IndexSpecIRTree::getValue(int step, DateTime eventDate) const {
    try {
        if (historicalValueAvailable(*indexSpecIR, *slice, eventDate))
            return *slice;

        DateTime currentDate = model->getDate(step);

        if (rateTree->getZeroBankMode() == ZeroBankMode::CLAIMBANK) {
            rateTree->getIRIndex(*indexSpecIR, currentDate, eventDate, *slice);
        }
        else {
            if (currentDate != eventDate) {
                throw ModelException("Tree model zeroBankMode ZEROBANK does not support offset "
                                     "reset dates (resetDate = " + eventDate.toString() + ", "
                                     "currentDate = " + currentDate.toString() + ")");
            }
            rateTree->getIRIndex(*indexSpecIR, model->getDate(step), *slice);
        }
        return *slice;
    }
    catch (exception& e) {
        throw ModelException(e, __FUNCTION__, "(indexSpecIR=" + indexSpecIR->getName() + ")");
    }
}

// must use esl_zerobank.c: ZbkDLFromIdx()


void IndexSpecIRTree::calculateNamedZeroBank1(RateTree::NamedZeroBank& zeroBank) const {
    try {
        DateTime resetFirstDate = modelResetDates.front();
        DateTime resetLastDate = modelResetDates.back();
        int resetFrequency;
        if (modelResetDates.size() > 1) {
            MaturityPeriodSP mp(MaturityPeriod::dateSubtract(modelResetDates[1], modelResetDates[0]));
            resetFrequency = mp->approxAnnualFrequency();
        }
        else resetFrequency=1;

        int tenorInMonths = indexSpecIR->tenor->toMonths();

        /* Last date in tree is maturity of the swap plus maturity
         * of index used in the floating reset. For reset-in-arrears, this
         * is exactly the required maturity (as last reset is on the swap
         * maturity), and for reset-in-advance this gives and extra period
         * beyond the last reset maturity.
         * One extra zero maturity one month after the 'real' last one. */
        DateTime lastZeroMatDate = RatesUtils::nextMonth(resetLastDate, tenorInMonths+1, true);

        /* First zero maturity is first reset date, and if it does not exits,
         * it is the value date. (There is a need for zero bank if there is
         * an exercise date in the last period of swap reseting-in-advance) */
        DateTime firstZeroMatDate = resetFirstDate.max(model->getValueDate());

        /* Zero dates are generated based on
         * Libor index: the minimum of frequency of the floating index and
         *              the payment frequency of the floating leg;
         * Cms index:   subject to a maximum interval of 6 months between
         *              zeros
         * The zero dates start from the first reset date (thereby ensuring
         * that the first reset is not interpolated). */

        int rateFreq  = indexSpecIR->frequency->annualFrequency();

        /* The only allowable freq exceeding semi is annual */
        int zeroFreq;
        if (tenorInMonths <= 12)
            zeroFreq = max(rateFreq, resetFrequency);  /* Libor index, 1yCms */
        else if (rateFreq==1)
            zeroFreq = 2;
        else
            zeroFreq = rateFreq;         /* XyCms index, X = 2, 3, ... */

        /* Create a date list with the ZeroFreq */
        RatesUtils::genDateList(firstZeroMatDate, lastZeroMatDate,
            zeroFreq, true, zeroBank.zeroMatDateList);

        /* Finally determine number of zeros needed at any one time.
         * Force an extra zero to cover cases where
         * we might end up extrapolating */
         zeroBank.nbZeros = (tenorInMonths*zeroFreq)/12 + 1;
    }
    catch (exception& e) {
        throw ModelException(e, __FUNCTION__, "calculateNamedZeroBank, name = " + zeroBank.name);
    }
}

void IndexSpecIRTree::calculateNamedZeroBankCB(RateTree::NamedZeroBank& zeroBank) const {
    try {
        DateTimeArray swapStart = modelResetDates;
        if (indexSpecIR->fwdRateOffset.get()) {
            for (int i=0; i<swapStart.size(); ++i) {
                swapStart[i] = indexSpecIR->fwdRateOffset->toDate(swapStart[i]);
            }
        }

        DateTimeArray idxZUse;
        RatesUtils::zbkDLFromIdx(modelResetDates, swapStart, indexSpecIR->tenor, 
            indexSpecIR->frequency, zeroBank.zeroMatDateList, idxZUse);
//        zeroBank.nbZeros = zeroBank.zeroMatDateList.size();
        int tenorInMonths = indexSpecIR->tenor->toMonths();
        int zeroFreq = indexSpecIR->frequency->annualFrequency();
        zeroBank.nbZeros = (tenorInMonths*zeroFreq)/12+1;

        ASSERT(zeroBank.zeroMatDateList.size() == idxZUse.size());
    }
    catch (exception& e) {
        throw ModelException(e, __FUNCTION__);
    }
}

void IndexSpecIRTree::calculateNamedZeroBank2(RateTree::NamedZeroBank& zeroBank) const {
    /* Last date in tree is max of last coupon payment and maturity  */
    /* of  index  used  in the last floating reset.                  */
    /* So, set last zero date = "last reset" + rate mat + 1          */
    /* "last reset" is the max of Last reset date from flt schedule  */
    /* and last exercise date (partway exercise)                     */

    int tenorInMonths = indexSpecIR->tenor->toMonths();
    DateTime resetFirstDate = modelResetDates.front();
    DateTime resetLastDate = modelResetDates.back();
    DateTime lastZeroMatDate = RatesUtils::nextMonth(resetLastDate, tenorInMonths+1, true);
    DateTime firstZeroMatDate = resetFirstDate.max(model->getValueDate());

    int zeroFreq;
    if (tenorInMonths <= 12)
        zeroFreq = 12; /* for Libor index */
    else if (tenorInMonths <= 36)
        zeroFreq = 4;  /* for 1,2,3 CMS */
    else
        zeroFreq = 2;  /* for longer CMS index */

    /* Create a date list with the ZeroFreq */
    RatesUtils::genDateList(firstZeroMatDate, lastZeroMatDate, zeroFreq, true, zeroBank.zeroMatDateList);

    zeroBank.nbZeros = (tenorInMonths*zeroFreq)/12 + 1;
}

void IndexSpecIRTree::calculateNamedZeroBank2bis(RateTree::NamedZeroBank& zeroBank) const {
    /* Last date in tree is max of last coupon payment and maturity  */
    /* of  index  used  in the last floating reset.                  */
    /* So, set last zero date = "last reset" + rate mat + 1          */
    /* "last reset" is the max of Last reset date from flt schedule  */
    /* and last exercise date (partway exercise)                     */

    int tenorInMonths = indexSpecIR->tenor->toMonths();
    DateTime resetFirstDate = modelResetDates.front();
    DateTime resetLastDate = modelResetDates.back();
    DateTime lastZeroMatDate = RatesUtils::nextMonth(resetLastDate, tenorInMonths, true);
    DateTime firstZeroMatDate = resetFirstDate.max(model->getValueDate());

    int zeroFreq;
    if (tenorInMonths <= 12)
        zeroFreq = 12; /* for Libor index */
    else if (tenorInMonths <= 36)
        zeroFreq = 4;  /* for 1,2,3 CMS */
    else
        zeroFreq = 2;  /* for longer CMS index */

    /* Create a date list with the ZeroFreq */
    RatesUtils::genDateList(firstZeroMatDate, lastZeroMatDate, zeroFreq, true, zeroBank.zeroMatDateList);

    zeroBank.nbZeros = (tenorInMonths*zeroFreq)/12 + 1;
}

void IndexSpecIRTree::calculateNamedZeroBank3(RateTree::NamedZeroBank& zeroBank) const {
    try {
        int tenorInMonths = indexSpecIR->tenor->toMonths();
        int frontIdx=0;
        while (frontIdx<modelResetDates.size() && modelResetDates[frontIdx].getDate()<model->getToday().getDate()) ++frontIdx;
        if (frontIdx==modelResetDates.size())
            return; // no reset in the future
        DateTime resetFirstDate = modelResetDates[frontIdx];
        DateTime resetLastDate = modelResetDates.back();
        DateTime lastZeroMatDate = RatesUtils::nextMonth(resetLastDate, tenorInMonths+1, true);
        // lastZeroMatDate = lastZeroMatDate.rollDate(1);
        DateTime firstZeroMatDate = resetFirstDate.max(model->getValueDate());

        int resetFrequency;
        if (modelResetDates.size() > 1) {
            MaturityPeriodSP mp(MaturityPeriod::dateSubtract(modelResetDates[1], modelResetDates[0]));
            resetFrequency = mp->approxAnnualFrequency();
        }
        else 
            resetFrequency=1;

        int rateFreq = indexSpecIR->frequency->annualFrequency();
        int zeroFreq;

        if (tenorInMonths <= 12) {  // observing on libor rates
            if (resetFrequency >= 52) {  // if weekly or daily
                if (rateFreq <= 4)  // rate frequency A, S or Q
                    zeroFreq = 12;
                else
                    zeroFreq = rateFreq;
            }
            else
                zeroFreq = max(rateFreq, resetFrequency);
        }
        else {   // observing on CMS rates
            if (rateFreq == 1)
                zeroFreq = 2;
            else
                zeroFreq = rateFreq;
        }

        /* Create a date list with the ZeroFreq */
        RatesUtils::genDateList(firstZeroMatDate, lastZeroMatDate, zeroFreq, true, zeroBank.zeroMatDateList);

        zeroBank.nbZeros = (tenorInMonths*zeroFreq)/12 + 1;
    }
    catch (exception& e) {
        throw ModelException(e, "KFloatLegTree::calculateNamedZeroBank3");
    }
}

void IndexSpecIRTree::calculateNamedZeroBankExplicit(RateTree::NamedZeroBank& zeroBank) const {
    try {
        if (indexSpecIR->nbZerosAtOneTime == -1)
            throw ModelException("if using explicit zeroDates mode, must supply nbZerosAtOneTime field");

        zeroBank.zeroMatDateList = indexSpecIR->explicitZeroDates;
        zeroBank.nbZeros = indexSpecIR->nbZerosAtOneTime;
    }
    catch (exception& e) {
        throw ModelException(e, "KFloatLegTree::calculateNamedZeroBankExplicit");
    }
}

/****************************** IndexSpecFXTree ********************************/

class IndexSpecFXTree : public FDProduct {
public:
    IndexSpecFXConstSP  indexSpecFX;
    RateTree*           rateTree;

    /******************** methods ********************/
    IndexSpecFXTree(const IndexSpecFXConstSP &indexSpecFX, RateTree* rateTree);
    DateTime getStartDate(void) const {return model->getValueDate();}

    // FDProduct virtual methods - begin
    virtual void init(Control*) const {}
    virtual void initProd(void);
    virtual void update(int& step, UpdateType type) {}
    virtual const TreeSlice& getValue(int step, DateTime eventDate) const;
    virtual void recordOutput(Control* ctrl, YieldCurveConstSP, Results* results) {}
    virtual bool isElementary() const { return true; }
    virtual FDDirection getDirectionPref() const { return FDProduct::FWD_BACK_INDUCTION; }
    // FDProduct virtual methods - end

private:
    TreeSliceSP slice;
};

IndexSpecFXTree::IndexSpecFXTree(const IndexSpecFXConstSP &indexSpecFX, RateTree* rateTree) :
    FDProduct(rateTree), indexSpecFX(indexSpecFX), rateTree(rateTree)
{}

void IndexSpecFXTree::initProd(void) {
    slice = model->createSlice(); 
    slice->name = "IndexSpecFX";
}

const TreeSlice& IndexSpecFXTree::getValue(int step, DateTime eventDate) const {
    try {
        // no need to bother with FX adjustments if single currency
        if (indexSpecFX->sameCurve) {
            *slice = 1.;
            return *slice;
        }
      
        // retrive past fixing from assest history associated with indexSpec itself
        if (historicalValueAvailable(*indexSpecFX, *slice, eventDate))
            return *slice;

        DateTime currentDate = rateTree->getDate(step);
        if (rateTree->getZeroBankMode() == ZeroBankMode::CLAIMBANK) {

            rateTree->getFXIndex(*indexSpecFX, currentDate, eventDate, *slice);
        }
        else {
            if (currentDate != eventDate) {
                throw ModelException("Current Date " + currentDate.toString() 
                + " must equal eventDate " + eventDate.toString() 
                + " - date offset logic not yet supported");
            }

            rateTree->getFXIndex(*slice);
        }
        return *slice;
    }
    catch (exception& e) {
        throw makeException(e, __FUNCTION__);     
    }
}

/****************************** IndexSpecFXTree ********************************/

class IndexSpecEQTree : public FDProduct {
public:
    IndexSpecEQConstSP  indexSpecEQ;
    RateTree*           rateTree;

    /******************** methods ********************/
    IndexSpecEQTree(const IndexSpecEQConstSP &indexSpecEQ, RateTree* rateTree);
    DateTime getStartDate(void) const {return model->getValueDate();}

    // FDProduct virtual methods - begin
    virtual void init(Control*) const {}
    virtual void initProd(void);
    virtual void update(int& step, UpdateType type) {}
    virtual const TreeSlice & getValue(int step, DateTime eventDate) const;
    virtual void recordOutput(Control* ctrl, YieldCurveConstSP, Results* results) {}
    virtual bool isElementary() const { return true; }
    virtual FDDirection getDirectionPref() const { return FDProduct::FWD_BACK_INDUCTION; }
    // FDProduct virtual methods - end

private:
    TreeSliceSP slice;
};

IndexSpecEQTree::IndexSpecEQTree(const IndexSpecEQConstSP &indexSpecEQ, RateTree* rateTree) :
    FDProduct(rateTree), indexSpecEQ(indexSpecEQ), rateTree(rateTree)
{}

void IndexSpecEQTree::initProd(void) {
    slice = model->createSlice();
    slice->name = "IndexSpecEQ";
}

const TreeSlice & IndexSpecEQTree::getValue(int step, DateTime eventDate) const {
    try {
        if (historicalValueAvailable(*indexSpecEQ, *slice, eventDate))
            return *slice;

        DateTime currentDate = rateTree->getDate(step);
        if (rateTree->getZeroBankMode() == ZeroBankMode::CLAIMBANK) {

            rateTree->getEQIndex(*indexSpecEQ, currentDate, eventDate, *slice);
        }
        else {
            if (currentDate != eventDate) {
                throw ModelException("Current Date " + currentDate.toString() 
                + " must equal eventDate " + eventDate.toString() 
                + " - date offset logic not yet supported");
            }

            rateTree->getEQIndex(*slice);
        }
        return *slice;
    }
    catch (exception& e) {
        throw makeException(e, __FUNCTION__);     
    }
}

/********************************** RateTree::NamedZeroBank ********************************/

RateTree::NamedZeroBank::~NamedZeroBank() {
    if (!zeroBank.empty()) // cannot throw an exception from destructor
        fprintf(stderr,"%s: zeroBank not empty\n",__FUNCTION__);
}

/********************************** RateTree ********************************/

RateTree::RateTree(const CClassConstSP &type) : FDModel(type),
    mTpIdx(-1), CBLegacyPricingMode(false) 
{}


IRVolRawSP RateTree::getIRVolRaw(RateTree::CcyIRParamsSP modelIRParams) {
    try
    {
        VolRequestRaw tmpVolRequest;  // not used except to construct vol
        IRVolRawSP volSP(dynamic_cast<IRVolRaw*>(
            modelIRParams->curveToDiffuse->getProcessedVol(&tmpVolRequest)));
        if (!volSP)
            throw ModelException("IRVolPair required in market for diffusion curve");

        return volSP;
    }
    catch (exception& e) {
        throw ModelException(e, __FUNCTION__);
    }
}

// register required IR index/fixings with the engine in order for it to
// configure the required claim banks
void RateTree::insertIRIndex(const IndexSpecIR& rateSpec, DateTime resetDate, DateTime lateResetDate,  bool isCrit) {
    try {

        string curveName = rateSpec.getFactor()->getName();
        DateTime todayDate = getToday();

        map<string, NamedClaimBank>::iterator cb;

        DateTimeArray resetDates(1);
        DateTimeArray indexStartDates(1);
        DateTimeArray idxZUse;
        DateTimeArray idxZMat;

        resetDates[0] = resetDate;
        if (rateSpec.fwdRateOffset.get())
            // ??? indexStartDate should be holiday adjusted
            indexStartDates[0] = rateSpec.fwdRateOffset->toDate(resetDate);
        else
            indexStartDates[0] = resetDate;  //  no offset

        // ??? rateStartDate should be holiday adjusted.  The question will be do we want 
        // to regenerate these dates when calling getIRIndex or cache them somewhere - 
        // maybe indexSpecIR can calculate all the index dates itself then just insert 
        // them as zeroDates here
        if (indexStartDates[0] < resetDate) {
            throw ModelException("Negative fwdRateOffset not allowed "
                "(IndexSpecIR \""+rateSpec.name+"\").");
        }

        // only insert future resets in tree - this also covers AM/PM pricing if resetDate = today
        if (resetDate < todayDate)
            return;

        RatesUtils::zbkDLFromIdx(resetDates, indexStartDates, rateSpec.tenor, rateSpec.frequency,
            idxZMat, idxZUse);

        // if claimBank doesn't exist for this curveName, create a new one, otherwise
        // add the dates to the existing bank
        cb = claimBank.find(curveName);

        // if claimBank doesn't exist, create a new one and insert into map
        if (cb == claimBank.end()) {
            NamedClaimBank newClaimBank;
            newClaimBank.curveName = curveName;
            claimBank[curveName] = newClaimBank;
            cb = claimBank.find(curveName);
        }

        // compute latest zero date
        DateTime latestDate = lateResetDate;
        MaturityPeriodSP mp = RatesUtils::getMaturityPeriod(rateSpec.fwdRateOffset, resetDate);
        if (mp.get())
            latestDate = mp->toDate(latestDate);
        latestDate = rateSpec.tenor->toDate(latestDate);

        // Add to zero bank dates       
        if (isCrit) {
            for (int i=0; i<idxZUse.size(); i++) {
                cb->second.critZeroUseDates.push_back(idxZUse[i].toIrDate());  // reset/observation date
                cb->second.critZeroMatDates.push_back(idxZMat[i].toIrDate());  // actual zero dates
                cb->second.critZeroLabel.push_back(string("RATE - IndexSpecName:") + 
                    rateSpec.getName());
            }
            // if in legacy claim bank mode, don't include this final critical date
            // as it makes it difficult to compare pricing with legacy models
            if (!CBLegacyPricingMode) {
                cb->second.critZeroUseDates.push_back(resetDate.toIrDate());  // reset/observation date
                cb->second.critZeroMatDates.push_back(latestDate.toIrDate());  // actual zero dates
                cb->second.critZeroLabel.push_back(string("RATE - Extra late zero "
                    "for IndexSpecName:") + rateSpec.getName());
            }
        }
        else {
            for (int i=0; i<idxZUse.size(); i++) {
                cb->second.optZeroUseDates.push_back(idxZUse[i].toIrDate());
                cb->second.optZeroMatDates.push_back(idxZMat[i].toIrDate());
                // ??? implement optionZeroLabel
            }
            cb->second.optZeroUseDates.push_back(resetDate.toIrDate());
            if (!CBLegacyPricingMode) {
                cb->second.optZeroMatDates.push_back(latestDate.toIrDate());
            }
        }
    }
    catch (exception &e) {
        throw ModelException(e,__FUNCTION__);
    }
}


FDProductSP RateTree::makeProduct(const IProdCreatorSP & creator) {
    try {
        const IndexSpec *indexSpec = dynamic_cast<const IndexSpec *>(creator.get());
        if (indexSpec && !indexSpec->setupCalled)
            throw ModelException(indexSpec->getClass()->getName()+" \""+indexSpec->getName()+"\" not setup");

        if (IndexSpecIR::TYPE->isInstance(creator)) {
            return FDProductSP(new IndexSpecIRTree(IndexSpecIRConstSP::dynamicCast(creator), this));
        }
        if (IndexSpecFX::TYPE->isInstance(creator)) {
            return FDProductSP(new IndexSpecFXTree(IndexSpecFXConstSP::dynamicCast(creator), this));
        }
        if (IndexSpecEQ::TYPE->isInstance(creator)) {
            return FDProductSP(new IndexSpecEQTree(IndexSpecEQConstSP::dynamicCast(creator), this));
        }
        return FDModel::makeProduct(creator);
    }
    catch (exception& e) {
        throw ModelException(e, "RateTree::makeProduct");
    }
}

void RateTree::load(CClassSP& clazz)
{
    clazz->setPublic();
    REGISTER(RateTree, clazz);
    SUPERCLASS(FDModel);

    // these fields are transient to allow tweaking to copy
    // the model after ::getMarket() has been called
    FIELD(indexSpecs,"")
    FIELD_MAKE_TRANSIENT_BUT_TWEAKABLE(indexSpecs)
    FIELD(fxFactors,"");
    FIELD_MAKE_TRANSIENT_BUT_TWEAKABLE(fxFactors)

    FIELD(engineSet, "(opt) IR Engine key for the engine numerics (eg. tree, MC, FD etc).");
    FIELD_MAKE_OPTIONAL(engineSet);
    FIELD(engineTable, "(opt) MarketTable map containing the collection of IR engine objects with corresponding keys.");
    FIELD_MAKE_OPTIONAL(engineTable);
}

CClassConstSP const RateTree::TYPE = CClass::registerClassLoadMethod(
    "RateTree", typeid(RateTree), RateTree::load );

// for class loading
bool RateTreeLoad(void) {
    return (RateTree::TYPE != 0);
}

TreeSliceSP RateTree::createSimpleSlice(
    const string & curveToDEV ) const
{
    TreeSliceRatesSP newSlice( new TreeSliceRates(*range, curveToDEV, getCurveIdx(curveToDEV)));
    newSlice->treeStep = range->treeStep; //???? mTpIdx;
    return newSlice;
}

/** creates new slice */
TreeSliceSP RateTree::createSlice(
    const string & curveToDEV,
    const string & factorName1,
    const string & factorName2,
    const string & factorName3 ) const
{
    TreeSliceSP newSlice = createSimpleSlice(curveToDEV);

    // state variable support
    return createLayer(newSlice);
}

/** returns an ExposureHighlighter - a "model" that does everything except
actually price, so you get to see what market data it uses */
ExposureHighlighter* RateTree::exposureHighlighter() {
    return new IRVegaPointwiseExposureHighlighter(this);
}

IModel::WantsRiskMapping RateTree::wantsRiskMapping() const {
    return riskMappingIrrelevant;
}


///  solver loop : decide if forward or backward induction
void RateTree::roll()
{
    InductionType::Enum solverType = getInductionType();
    switch ( solverType )
    {
        case InductionType::FWD : rollFront();
                break;
        case InductionType::EXPRESSBACK: produceStatePrices(); // roll forward and store the state prices
        case InductionType::BACK : rollBack();
                break;
        default: 
            throw ModelException("RateTree::roll: only FWD, BACK and EXPRESSBACK induction implemented");
            break;

    }
}
///  solver loop 
// helper to loop slice layers
void RateTree::loopDev(const vector< TreeSliceSP >& slices, int curveIdx)
{
    for (int j=0; j<(int)slices.size(); j++) {
        // check slice treeStep only if current treeStep is below previos to last
        if (slices[j]->treeStep != range->treeStep+1
            && range->treeStep < getLastStep()-1) 
        {
            throw ModelException("Cannot DEV slice " + slices[j]->name  + 
            " from step "+Format::toString(slices[j]->treeStep)+" to step "
            +Format::toString(range->treeStep));
        }

        // do DEV
        TreeSliceLayer* layer = dynamic_cast<TreeSliceLayer*>(slices[j].get());
        if (layer) {
            loopDev(layer->getSlices(), layer->getSlices()[0]->getCurveToDEVIdx());
        }
        else if (!slices[j]->isZero()) {
            if( curveIdx < 0 ) {
                curveIdx = (*slices[j]).getCurveToDEVIdx();
            }
            if (curveIdx<0) {
                throw ModelException("Invalid curve to DEV index ("
                +Format::toString(curveIdx)+") for slice "+slices[j]->name);
            }
            sliceDev(*slices[j], curveIdx);
        }
        // update treeStep
        slices[j]->treeStep = range->treeStep;
    }
}

void RateTree::loopExpand(TreeSlice &slice, const string &curveToDev) {

    if (slice.isZero()) {
        return;
    }
    TreeSliceLayer *layer = dynamic_cast<TreeSliceLayer *>(&slice);
    if (!layer) {
        if (!slice.isZero()) {
            sliceExpand(slice, curveToDev);
        }
        return;
    }
    const vector< TreeSliceSP >& slices = layer->getSlices();
    for (int i=0; i<(int)slices.size(); ++i) {
        loopExpand(*slices[i], curveToDev);
    }
}

///  solver loop
// Update tree by sweeping timeline
void RateTree::rollBack()
{
    int t;
    int T = getLastStep();
    FDProduct *curProd=0; // product being used currently (for debug messages)
    string prodType, prodOutputName;
    
    try {

        const FDProductArray & products = getProducts();  //for easy reference
        FDProduct::UpdateType type = FDProduct::BWD_T;

        for (t = T; t >= 0; t--) {

            // pre-DEV expantion
            if (t!=T) {
                for(size_t k=0; k < products.size(); k++) 
                {
                    curProd = products[k].get(); // to improve error messages
                    if (!products[k]->isCalcOff()) 
                    {
                        const vector< TreeSliceSP > & slices 
                            = products[k]->getSlicesToDEV();

                        for (int j=0; j<(int)slices.size(); j++) {
                            string sliceName = slices[j]->name;  // to help when debugging
                            loopExpand(*slices[j],slices[j]->getCurveToDEV());
                        }
                    }
                }
            }
            curProd = 0; // to improve error messages and run time debugging

            // update the engine prob's and internal variables
            update(t, type);

            // preCalc() and DEV
            for(size_t k=0; k < products.size(); ++k) 
            {
                curProd = products[k].get(); // to improve error messages

                if (!products[k]->isCalcOff()) {

                    // product with no slices still calls preCalc() and update()
                    products[k]->preCalc(t);

                    // handles both cases state var and not
                    loopDev( products[k]->getSlicesToDEV() );
                }
            }

            for(size_t k=0; k < products.size(); ++k) 
            {
                curProd = products[k].get(); // to improve error messages

                if (!products[k]->isCalcOff()) {
                    products[k]->update(t, type);
                }
            }
            curProd = 0; // to improve error messages and run time debugging

            updateFinalize(t, type);
            type = FDProduct::BWD;
        }
    }
    catch (exception& e) {
        DateTime treeStepDate = getDate(t);
        string prodDescr;

        if (curProd) {
            prodDescr = "while processing product (type: \"" 
                + string(typeid(*curProd).name())
                + "\" product name: \""+ curProd->getOutputName()+"\") ";
        }

        throw ModelException(e, __FUNCTION__, "failed " + prodDescr
            + "on timestep " + Format::toString(t) + " (" 
            + treeStepDate.toString() + ") of " + Format::toString(T)
            + ",  model: \"" + getClass()->getName()+"\"");
    }
}

///  solver loop 
// Update tree by sweeping timeline
void RateTree::rollFront()
{
    int t;
    int T = getLastStep();
    string prodType, prodOutputName;
    
    try {
    
        const FDProductArray&  products = getProducts();  //for easy reference
        FDProduct::UpdateType type = FDProduct::FWD_0; // ??

        for (t = 0; t <= T; ++t) {

            // update engine prob's and internal variables
            updateStatePrice(t, type); // update state prices

            // expand slices to the correct dimension for dev 
            for(size_t k=0; k < products.size(); k++) {
                if (products[k]->isCalcOff())
                    continue;

                // update the product (i.e. multiply statePrices with payout if needed)
                products[k]->update(t, type); 

                // to assist debugging print out type and name if available of current product
                prodType = typeid((*products[k])).name();
                prodOutputName = products[k]->getOutputName();
            }

            // updateFinalize(t, type); // ?? probably not needed
            type = FDProduct::FWD;
        }
    }
    catch (exception& e) {
        throw ModelException(e, "RateTree::Solver::rollFront failed for product type \"" + 
                             prodType + "\", product outputName \"" +
                             prodOutputName + "\" on timestep " +
                             Format::toString(t) + " of " +
                             Format::toString(T)+ " (model is "+getClass()->getName()+")");
    }
}

///  solver loop 
// Update tree by sweeping timeline
void RateTree::produceStatePrices()
{
    int t;
    int T = getLastStep();
    string prodType, prodOutputName;

    try {
        DateTimeArray::const_iterator date_it = statePriceDates.begin();

        while ( (date_it != statePriceDates.end()) &&  (*(date_it) <= getDate(0) ) )
            date_it++ ; 

        // only do this part if we need any state prices
        const FDProductArray&  products = getProducts();  //for easy reference
        FDProduct::UpdateType type = FDProduct::FWD_0; // ??

        for (t = 0; (t <= T) && (date_it!=statePriceDates.end()) ; ++t) {

            // update engine prob's and internal variables
            updateStatePrice(t, type); // update state prices

            // do we need to store the stateprices 
            DateTime currentDate = getDate(t);

            if (currentDate == (*date_it) )
            {
                TreeSliceSP currStatePriceSP = TreeSliceSP( getStatePriceSlice(t).clone() ); 
                statePrices.insert( make_pair(t, currStatePriceSP) );

                date_it++; // go to the next StatePrice Date
            }
            type = FDProduct::FWD;

            if ( (date_it != statePriceDates.end()) && ( currentDate > (*date_it) ) )
            {
                throw ModelException("Cannot produce cached state-prices at step " 
                    + Format::toString(t) + " for date point " + (*date_it).toString()  + 
                    ", state-price date point not in timeline!!!" );
            }
        }
    }
    catch (exception& e) {
        throw ModelException(e, "RateTree::Solver::produceStatePrices failed for product type \"" + 
            prodType + "\", product outputName \"" +
            prodOutputName + "\" on timestep " +
            Format::toString(t) + " of " +
            Format::toString(T)+ " (model is "+getClass()->getName()+")");
    }
}


/** interpolate/integrate for prices at t=0 with asset=asset(t=0) etc. */
double RateTree::getPrice0( const TreeSlice & price ) const
{
    try {
        const TreeSliceLayer * layer = dynamic_cast< const TreeSliceLayer * >( &price );
        if( layer )
            return layer->getPrice0();

        return price.getCentre();
    }
    catch (exception& e) {
        throw ModelException(e, "RateTree::getPrice0");
    }
}

void RateTree::ExtendTreeCurve(T_CURVE *tc, const DateTime &newStartDate)
{
    // Convert dates and delegate to ESL function.  The first date is
    // used to extend the start of the curve (i.e. before base/value date) while
    // the second date is used to extend the end of the curve.  We pass the
    // same date twice to effectively disable extension of end of the curve.
    // NOTE: This function is not (yet) implemented for IRX curves. 
    ExtendTCurve(tc, newStartDate.toIrDate(), newStartDate.toIrDate());
}

/*************************** ZeroBond ************************/

ZeroBond::ZeroBond(const DateTimeArray& startDatesL,
                   const DateTimeArray& matDatesL,
                   const string& curveName)
    : CObject(TYPE), curveName(curveName),
    // for security
    zeroSlice(0), currencyEnum(RateTree::DOMESTIC), sliceDim(0), devMode(0)
{
    static int uniqueId=0;
    name = "ZeroBond-"+Format::toString(++uniqueId);

    if (startDatesL.size() != matDatesL.size()) throw ModelException(
        "ZeroBond", "startDatesL.size() != matDatesL.size()");

    // copy startDatesL, matDatesL to startDates, matDates removing null-lenght bonds
    int nb = startDatesL.size();
    startDates.reserve(nb);
    matDates.reserve(nb);
    for (int i=0; i<nb; ++i) {
        if (startDatesL[i]!=matDatesL[i]) {
            startDates.push_back(startDatesL[i]);
            matDates.push_back(matDatesL[i]);
        }
    }
}

ZeroBond::~ZeroBond() {
    if (zeroSlice!=NULL) // cannot throw an exception from destructor
        fprintf(stderr,"%s: zeroSlice not NULL\n",__FUNCTION__);
}

/** returns critical dates */
DateTimeArraySP ZeroBond::getCritDates(void) const {
    DateTimeArraySP critDates(new DateTimeArray(startDates));
    critDates->insert(critDates->end(), matDates.begin(), matDates.end());
    return critDates;
}

FDProductSP ZeroBond::createProduct( FDModel * model ) const {
    return FDProductSP(new ZeroBondProd(ZeroBondConstSP(this), model));
}

void ZeroBond::load(CClassSP& clazz)
{
    REGISTER(ZeroBond, clazz);
    SUPERCLASS(CObject);
    IMPLEMENTS(FDModel::IIntoProduct);
    EMPTY_SHELL_METHOD(defaultConstructor);
}

CClassConstSP const ZeroBond::TYPE = CClass::registerClassLoadMethod(
    "ZeroBond", typeid(ZeroBond), load);

/****************************** ZeroBondProd *************************/

ZeroBondProd::ZeroBondProd(const ZeroBondConstSP &inst, FDModel* model)
    : FDProduct(model), inst(inst)
{
    rateTree = dynamic_cast<RateTree*>(model);
    if (!rateTree) throw ModelException("ZeroBondProd", "Model must derive from RateTree");
}

void ZeroBondProd::init(Control*) const {
    try {
        rateTree->insertZero(*inst);
    }
    catch (exception& e) {
        throw ModelException(e, "ZeroBondProd::init, Error for instrument type " +
                             inst->getClass()->getName() + " and name " +
                             inst->getName());
    }
}

void ZeroBondProd::initProd(){
    // create a single slice to store the tree Par_Yield function results
    // this slice is never DEV'd and is reused for each reset date
    mainSlice = model->createSlice();
    mainSlice->name = "ZeroBond";
}


const TreeSlice & ZeroBondProd::getValue(int step, DateTime matDate) const {
    DateTime stepDate = model->getDate(step);
    if (stepDate==matDate) {
        *mainSlice = 1.;
    } else {
        rateTree->getZero(*inst, step, *mainSlice);
    }
    return *mainSlice;
}

/****************************** CCyIRParams *************************/

void RateTree::CcyIRParams::validate(string const&ccy, string const& discName) const {
    try {
        // check curves are same currency
        if (curveToDiffuse->getCcy() != ccy) {
            throw ModelException("curveToDiffuse is of currency "
            + curveToDiffuse->getCcy() + " != instrument ccy " + ccy);
        }
        if (curveToDiscount->getCcy() != ccy) {
            throw ModelException("curveToDiscount is of currency "
            + curveToDiffuse->getCcy() + " != instrument ccy " + ccy);
        }
        if (discName.size() && curveToDiscount->getName() != discName) {
            throw ModelException("curveToDiscount is " + curveToDiscount->getName()
            + " != instrument discount " + discName);
        }
    }
    catch (exception& e) {
        throw ModelException(e, __FUNCTION__);
    }
}

void RateTree::CcyIRParams::getData(CModel *model, const MarketData*  market)
{
    static const string method = "RateTree::CcyIRParams::getData";
    if (!smileTable.isEmpty())
    {
        smileTable.getData(model, market);

        // Lookup Smile
        MarketObjectSP ptrSmile = smileTable->getMarketObject(smileSet);
        MarketObject*  pSmile = ptrSmile.get();
        if (!pSmile)
            throw ModelException(method, " - Smile object " + smileSet + " not found in Smile Table " + smileTable.getName() + "!");

        // Initialisation for IRSmileMQ
        IRSmileMQ* pIRSmileMQ = dynamic_cast<IRSmileMQ*>(pSmile);
        if (pIRSmileMQ)
            pIRSmileMQ->initialiseData(market, model);
    }

    if (!modelTable.isEmpty())
        modelTable.getData(model, market);

    // Phasing out of irCalib - do not use irCalib if using smileTable or modelTable.
    if (!irCalib.isEmpty() && smileTable.isEmpty() && modelTable.isEmpty())
        irCalib.getData(model, market);

    curveToDiffuse.getData(model, market);
    if(!curveToDiscount.isEmpty())
        curveToDiscount.getData(model, market);
}

void RateTree::CcyIRParams::load(CClassSP& clazz){
    clazz->setPublic();                     /* can be seen from outside the library */
    REGISTER(RateTree::CcyIRParams, clazz); /* always the same, mandatory */
    EMPTY_SHELL_METHOD(defaultConstructor); /* always the same, mandatory */
    SUPERCLASS(CObject);

    string volCalibIndexDescr("swaption/caplet volatility calibration index eg. 5yFix=5y column of swpn matrix");
    volCalibIndexDescr += "10yCMS=diagonal of swpn matrix, 3m=3m caplet column";
    FIELD(volCalibIndex, volCalibIndexDescr);
    FIELD_MAKE_OPTIONAL(volCalibIndex);
    FIELD(irCalib,"market cache collection of IR model parameters");
    FIELD_MAKE_OPTIONAL(irCalib);
    FIELD(nbFactors,"number of factors for this currency (default: 1)");
    FIELD_MAKE_OPTIONAL(nbFactors);
    FIELD(smileSet,"name of smile set required to be used the irCalib market object");
    FIELD(modelSet,"name of model parameters set to be used within the irCalib market object");
    FIELD(curveToDiffuse, "curve name to diffuse in tree model");
    FIELD(curveToDiscount, "curve name to discount from the model's perspective");
    FIELD_MAKE_OPTIONAL(curveToDiscount);
    FIELD(smileTable,"A MarketTable object containing collection of IR smile objects, accessed using the smile key given above.");
    FIELD_MAKE_OPTIONAL(smileTable);
    FIELD(modelTable,"A MarketTable object containing collection of IR model objects, accessed using the model key given above.");
    FIELD_MAKE_OPTIONAL(modelTable);
    FIELD(calibVolType, "Calibration vol surface type eg. BASE_VOL or SWAP_VOL.");
    FIELD_MAKE_OPTIONAL(calibVolType);
    FIELD(calibType,    "Calibration type eg. FIXED_TENOR, CMS or FIXED_MATURITY etc.");
    FIELD_MAKE_OPTIONAL(calibType);
    FIELD(calibTenor, "Calibration Tenor (or explicit maturity date when using FIXED_MATURITY.")
    FIELD_MAKE_OPTIONAL(calibTenor);
    FIELD(calibVolOverride, "Calibration spot vol override value - vol surfaces are ignored when used.")
    FIELD_MAKE_OPTIONAL(calibVolOverride);
}

CClassConstSP const RateTree::CcyIRParams::TYPE =
    CClass::registerClassLoadMethod("RateTree::CcyIRParams", typeid(RateTree::CcyIRParams),
    RateTree::CcyIRParams::load);

DEFINE_TEMPLATE_TYPE(RateTree::CcyIRParamsArray);

/****************************** RateTreeMDF *************************/

RateTreeMDF::RateTreeMDF(int num, const MarketObjectQualifiersSP q): qualifiers(q){
    ASSERT(num==1 || num==2 || num==3);
    // get 1,2 or 3 factor ir parameters
    setRetrievalMode(IRCalib::Model::TYPE, true, 
                     IRCalib::getModelType(num));

    // needs this for IRVol (base or swaption vol)
    setRetrievalMode(IRVol::TYPE, IYieldCurve::TYPE, true, IRVol::TYPE);

    // get 'IRVolPair' vols when asked for ir vols if appropriate
    // This will let through IRCalibs unchanged (since they are not
    // derived from IRVolPair)
    setRetrievalMode(IRVolBase::TYPE, IYieldCurve::TYPE, true, IRVolPair::TYPE);

    setRetrievalMode(CVolBase::TYPE, CAsset::TYPE, true, SRMEQVol::TYPE);
}

/** Resolve from multiple market data set of same name and type */
MarketObjectConstSP RateTreeMDF::resolveMultiData(const MarketObjectArray& marketObjs, CClassConstSP) const{
    if (marketObjs.size() ==1)
        return marketObjs[0]; // no need to further select
    else if (!qualifiers)
        throw ModelException("RateTreeMDF::resolveMultiData", 
            "no qualifiers supplied for resolving multiple data sets");

    MarketObjectArray::const_iterator i;
    for (i = marketObjs.begin();i != marketObjs.end(); ++i){
        if (qualifiers->equals((*i)->getMarketObjectQualifier().get()))
            return *i;
    }
    return MarketObjectConstSP();
}

// *** IRVolSelector stuff ***

// Predicate for finding IRGridPointAbs objects by value (note that this class
// is not derived from IObject).
struct IsSameGridPointAs : public std::unary_function<IRGridPointAbsSP,bool>
{
    IsSameGridPointAs(const IRGridPointAbsSP &b) : b(b) {}
    bool operator()(const IRGridPointAbsSP &a) const
    {
        return
            a->getExpiry()->equalTo(b->getExpiry().get()) &&
            a->getTenor()->equalTo(b->getTenor().get());
    }
private:
    const IRGridPointAbsSP& b;
};

// Add grid point to array if its not already there (by value).
inline void addUniqueVolExposure(
    IRGridPointAbsArray &volExposures,
    IRGridPointAbsSP volExposureSP)
{
    if (std::find_if(
            volExposures.begin(),
            volExposures.end(),
            IsSameGridPointAs(volExposureSP)) == volExposures.end())
        volExposures.push_back(volExposureSP);
}

/******************************************/

RateTree::IRVolSelector::IRVolSelector(
    IRVolRawSP     volSP,
    const T_CURVE& diffusionCurve,
    const string&  volCalibIndex,
    bool           smoothing) 
    :
    mVolSP(volSP),
    mDiffusionCurve(diffusionCurve),
    mVolCalibIndex(volCalibIndex),
    mSmoothing(smoothing) 
{ 
    validate(); 
}


RateTree::IRVolSelector::IRVolSelector(
    IRVolRawSP     volSP,
    const T_CURVE& diffusionCurve,
    const string&  volCalibIndex,
    bool           smoothing,
    CcyIRParamsSP  ccyIRParams) 
    :
    mVolSP(volSP),
    mDiffusionCurve(diffusionCurve),
    mVolCalibIndex(volCalibIndex),
    mSmoothing(smoothing) 
{
    if (ccyIRParams.get())
    {
        mCalibVolType = ccyIRParams->calibVolType;
        mCalibType = ccyIRParams->calibType;
        mCalibTenor = ccyIRParams->calibTenor;
        mCalibVolOverride = ccyIRParams->calibVolOverride;
    }
    validate(); 
}

// Helper method to process information from last vol selection and generate
// swaption vol exposures.
void RateTree::IRVolSelector::addSwapVolExposures(
    IRGridPointAbsArraySP volExposuresSP)
{
    for (int i = 0; i < mSelectedSV.NbPoints; ++i)
    {
        IRGridPointAbsSP gridPointSP(new IRGridPointAbs(
            mSwapVolExpiries[mSelectedSV.SwaptionExpiryIndices[i]],
            DateTime::fromIrDate(mSelectedSV.SwaptionExpiryDates[i]),
            mSwapVolTenors[mSelectedSV.SwapTenorIndices[i]],
            DateTime::fromIrDate(mSelectedSV.SwapMatDates[i])));

        addUniqueVolExposure(*volExposuresSP, gridPointSP);
    }
}

// Helper method to process information from last vol selection and generate
// swaption vol exposures.
void RateTree::IRVolSelector::addBaseVolExposures(
    IRGridPointAbsArraySP volExposuresSP)
{
    if (mSelectedBV.NbPoints > 0)
    {
        // Assumption is that we select a single column from base vol
        // matrix (i.e. that is why we use mBaseVolTenors[0]).
        DateTime baseVolMatDate(mBaseVolTenors[0]->toDate(mBaseVolBaseDate));

        for (int i = 0; i < mSelectedBV.NbPoints; ++i)
        {
            IRGridPointAbsSP gridPointSP(new IRGridPointAbs(
                mBaseVolExpiries[mSelectedBV.VolIndices[i]],
                DateTime::fromIrDate(mSelectedBV.VolDates[i]),
                mBaseVolTenors[0],
                baseVolMatDate));

            addUniqueVolExposure(*volExposuresSP, gridPointSP);
        }
    }
}

void RateTree::IRVolSelector::validate()
{
    if (!mVolSP)
        throw ModelException("IRVolSelector::validate()",
                             "No volatility object was supplied!");
}

// Helper function to select volatilties from IRVolRaw into MKTVOL_DATA
// and, potentially, record info for volatility exposure recording.
void RateTree::IRVolSelector::selectVols(
    MKTVOL_DATA&          mktVolData,
    IRGridPointAbsArraySP volExposuresSP,
    OutputNameConstSP     irVolNameSP)
{
    try
    {
        BASEVOL_DATA bv;
        SWAPVOL_DATA sv;

        const IRVol* baseVol = mVolSP->getBaseVol().get();
        const IRVol* swapVol = mVolSP->getSwapVol().get();

        if (!mVolCalibIndex.empty())
        {
            // USING VolCalibindex string
            if (baseVol)
            {
                IrConverter::to_BASEVOL_DATA(bv, mBaseVolTenors, mBaseVolExpiries,
                                            mBaseVolBaseDate, baseVol, mVolCalibIndex);
            }
            if (swapVol)
            {
                IrConverter::to_SWAPVOL_DATA(sv, mSwapVolTenors, mSwapVolExpiries, swapVol);
            }
            IrConverter::to_MKTVOL_DATA(mktVolData, &mSelectedBV, &mSelectedSV,
                mVolCalibIndex, mDiffusionCurve, (baseVol ? &bv : NULL), (swapVol ? &sv : NULL), mSmoothing);
        }
        else
        {
            // USING calibration enums
            const IRVol* useVol = (mCalibVolType == IRVol::BASE_VOL) ? baseVol : swapVol;
            IrConverter::to_MKTVOL_DATA(mktVolData, &mSelectedBV, &mSelectedSV,
                                        mDiffusionCurve, mSmoothing, useVol, mCalibVolType, mCalibType, mCalibTenor,
                                        mCalibVolOverride);
        }

        // Sensitivities
        if (volExposuresSP.get())
        {
            if (!irVolNameSP || irVolNameSP->equals(baseVol->getName()))
                addBaseVolExposures(volExposuresSP);

            if (!irVolNameSP || irVolNameSP->equals(swapVol->getName()))
                addSwapVolExposures(volExposuresSP);
        }
    }
    catch (exception& e)
    {
        throw ModelException(e, __FUNCTION__, "with IRVolPair \"" + mVolSP->getName() + "\"");
    }
}

// Populates volExposuresSP with relevant grid points for the named IRVol.
void RateTree::IRVolSelector::getVolExposures(
    IRGridPointAbsArraySP volExposuresSP,
    OutputNameConstSP irVolNameSP)
{
    if (!volExposuresSP)
        throw ModelException("RateTree::IRVolSelector::getVolExposures",
                             "NULL volExposure container encountered!");

    MKTVOL_DATA mktVolData;  // "Dummy"
    selectVols(mktVolData, volExposuresSP, irVolNameSP);
}

bool RateTree::happensNow(const DateTimeArray& array, DateTime stepDate, int *arrayPos)
{
    for (int i=0; i < array.size(); ++i) {
        if (array[i] == stepDate) {
            if (arrayPos)  // if requested, store array position index
                *arrayPos = i;
            return true;
        }
    }
    return false;
}

RateTree::IRModelParams::IRModelParams() {
    memset(this, 0, sizeof(*this));
    for (int i = 0; i < 3; ++i)
    {
        FactorVol[i]=-999;
        FactorMR[i]=-999;
        FactorCorr[i]=-999;
    }
}

// take copy of registered critical dates to use when printing out
// the dates in the print function, along with the description labels
void RateTree::NamedClaimBank::saveUnsorted() {
    critZeroMatDatesUnsorted = critZeroMatDates;
    critZeroUseDatesUnsorted = critZeroUseDates;
} 

YieldCurveConstSP RateTree::collectYC(string const& ycName, IObjectConstSP instruments) {
    try {
        // recurse the supplied instruments to retrieve the first YieldCurve object
        // type encountered.  This should be the top level instrument discount
        // curve object.

        class Action : public ObjectIteration::IActionConst {
            const string & name;
        public:
            YieldCurveConstSP yc;
            Action(const string & name) : name(name) {}

            virtual bool invoke(const ObjectIteration::State& state, IObjectConstSP obj) {
                const YieldCurveConstSP & ryc = YieldCurveConstSP::dynamicCast(obj);
                if (name == ryc->getName()) {
                    yc = ryc;
                    state.quitRecursion(true); // Found. Quit now.
                }
                // don't recurse inside the YieldCurve
                return false;
            }
        };

        Action action(ycName);
        ObjectIteration iteration(YieldCurve::TYPE);
        iteration.recurse(action, instruments);

        if (!action.yc) {
            throw ModelException("Object does not contain a YieldCurve field.");
        }
        return action.yc;
    }
    catch (exception& e) {
        throw ModelException(e, __FUNCTION__);
    }

}

void RateTree::collectFactors(FDModel * model,                // IN
                              IObjectSP instruments,          // IN
                              IMarketFactorArray & factors)   // OUT
{
    try {
        // recurse the instrument components to find and store all the IMarketFactor types found
        // in the exported fields - these are yield curve types for fix3
        class Action : public ObjectIteration::IAction {
            FDModel * model;
            IMarketFactorArray & factors;

        public:
            Action(FDModel* model, IMarketFactorArray& factors) :
                model(model), factors(factors) {factors.clear();}

            virtual bool invoke(ObjectIteration::State& state, IObjectConstSP obj) {
                IMarketFactorSP factor = IMarketFactorSP::dynamicCast(state.getObject());

                string name = factor->getName();
                string type = factor->getClass()->getName();
                // check if already added
                for (int i = 0; i < factors.size(); ++i) {
                    if (name == factors[i]->getName() 
                    && type == factors[i]->getClass()->getName())
                        return true;
                }
                // add the factor to the list
                if (model->acceptFactor(factor.get())) {
                    factors.push_back(factor);
                }
                return false;
            }
        };
        Action action(model, factors);
        ObjectIteration iteration(IMarketFactor::TYPE);
        iteration.recurse(action, instruments);
    }
    catch (exception& e) {
        throw ModelException(e, __FUNCTION__);
    }
}

DRLIB_END_NAMESPACE
