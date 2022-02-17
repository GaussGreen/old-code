//----------------------------------------------------------------------------
//
//   Group       : Global QR&D
//
//   Filename    : KRib2.cpp
//
//   Description : Rib component
//
//----------------------------------------------------------------------------

#include "edginc/config.hpp"
#include "edginc/KRib2.hpp"
#include "edginc/Control.hpp"
#include "edginc/Results.hpp"
#include "edginc/Addin.hpp"

DRLIB_BEGIN_NAMESPACE

KRib2Tree::KRib2Tree(const KRib2ConstSP &inst, FDModel* model) 
: FDProduct(model), inst(inst), nbCoupons(inst->payDates.size()) {
    try {
        KRib2Tree::CompressorSP compressor = model->getProdModifier<KRib2Tree::Compressor>();
        if (compressor.get()) {
            compressor->lop(model->getToday(), inst->obsDates, modelObsDates, inst->obsToCpn);
            inst->assignObsDatesToCoupons(modelObsToCpn, modelObsDates);
        }
        else {
            modelObsDates = inst->obsDates;
            modelObsToCpn = inst->obsToCpn;
        }

        // create products obsProd[] for obsUnds[]
        int nbObs = modelObsDates.size();
        obsProd.resize(nbObs);
        if (inst->obsUnds.size()==1) {
            // same underlying for all observations
            FDProductSP prod = model->createProduct(inst->obsUnds[0]);
            for (int i=0; i<nbObs; ++i)
                obsProd[i] = prod;
        }
        else {
            // different underlyings for each observation
            for (int i=0; i<nbObs; ++i)
                obsProd[i] = model->createProduct(inst->obsUnds[i]);
        }
        // send modelResetDates to obsProd
        for (size_t obsIdx=0; obsIdx<obsProd.size(); ++obsIdx) {
            int cpnIdx = modelObsToCpn[obsIdx];
            obsProd[obsIdx]->addModelResetDate(modelObsDates[obsIdx], inst->accEnd[cpnIdx]);
        }
    } 
    catch (exception& e) {
        throw makeException(e, __FUNCTION__);     
    }
}

void KRib2Tree::init(Control*) const {
    try {
        DateTime today = model->getToday();
        model->addCritDates(modelObsDates);

        // register zeros from pay dates to observation dates
        for (int obsIdx = modelObsDates.size()-1; obsIdx>=0; --obsIdx) {
            DateTime obsDate = modelObsDates[obsIdx];
            int cpnIdx = modelObsToCpn[obsIdx];

            if (inst->resetDates[cpnIdx] < obsDate) {
                model->registerZero(obsDate, inst->payDates[cpnIdx], inst->discountCurveName);
                // otherwise we use couponFixings[] instead of a zero slice
            }
        }
    } 
    catch (exception& e) {
        throw makeException(e, __FUNCTION__);     
    }
}

void KRib2Tree::initProd(void) {
    try {
        nbObsPerCoupon.resize(nbCoupons, 0);
        couponFixings.resize(nbCoupons);

        // create ribCoupons[] slices
        ribCoupons.resize(nbCoupons);
        for (int i=0; i<nbCoupons; ++i) {
            ribCoupons[i] = model->createSlice(inst->discountCurveName);
            ribCoupons[i]->name = inst->outputName + "_ribFractions_"+Format::toString(i);
            *ribCoupons[i] = 0.;
            startDEV(ribCoupons[i]);
        }
    } 
    catch (exception& e) {
        throw makeException(e, __FUNCTION__);     
    }
}


void KRib2Tree::couponHasReset(TreeSliceSP couponFixing, int cpnIdx) {
    try {
        (*ribCoupons[cpnIdx]) *= (*couponFixing);
        couponFixings[cpnIdx] = couponFixing;
    } 
    catch (exception& e) {
        throw makeException(e, __FUNCTION__);     
    }
}


void KRib2Tree::update(int &step, FDProduct::UpdateType updateType) {
    try {
        DateTime currDate = model->getDate(step);
        int nbObs = modelObsDates.size();

        for (int obsIdx=0; obsIdx<nbObs; ++obsIdx) {
            // look for an observation date
            if (modelObsDates[obsIdx] !=currDate)
                continue;

            int cpnIdx = modelObsToCpn[obsIdx];
            addObservation(step, *ribCoupons[cpnIdx], nbObsPerCoupon[cpnIdx], cpnIdx, obsIdx);
        }
    } 
    catch (exception& e) {
        throw makeException(e, __FUNCTION__);     
    }
}

void KRib2Tree::deleteCoupon(int couponIdx) {
    ASSERT(couponIdx < (int)ribCoupons.size());

    stopDEV(ribCoupons[couponIdx]);
    ribCoupons[couponIdx] = TreeSliceSP();
    couponFixings[couponIdx] = TreeSliceSP();
}

void KRib2Tree::addObservation(int step, TreeSlice &accuSlice, int &nbObs, int couponIdx, int obsIdx) {
    try {
        DateTime currentDate = model->getDate(step);
        TreeSliceSP couponFixing;
        if (couponFixings[couponIdx].get()) {
            couponFixing = couponFixings[couponIdx];
        }
        else {
            // get a zero bond from the payment date
            model->getZero(currentDate, inst->payDates[couponIdx], inst->discountCurveName, couponFixing);
        }
        // update the rib fraction of the coupon
        const TreeSlice &obsValue = obsProd[obsIdx]->getValue(step, modelObsDates[obsIdx]);
        accuSlice += obsValue * (*couponFixing);
        ++nbObs;
    } 
    catch (exception& e) {
        throw makeException(e, __FUNCTION__);     
    }
}

void KRib2Tree::getRibCoupon(int step, int couponIdx, TreeSliceSP &ribCoupon) {
    try {
        DateTime currentDate = model->getDate(step);
        if (!ribCoupons[couponIdx]) {
            throw ModelException("Slice ribCoupon["+Format::toString(couponIdx)
            +"] is not available anymore (deleteCoupon() was called)");
        }
        ribCoupon = ribCoupons[couponIdx]->clone(true);
        int nbObs = nbObsPerCoupon[couponIdx];

        if (inst->accStart[couponIdx] <= currentDate) {
            // If the RIB observation is not complete for the coupon, 
            // add the missing ones now. The inderlying observation
            // index may take its values from the past fixings
            // or estimate them

            for (int obsIdx = modelObsDates.size()-1; obsIdx>=0; --obsIdx) {

                // skip observations not related to this coupon
                if (modelObsToCpn[obsIdx] != couponIdx)
                    continue;

                // skip the observations that have been added in update()
                DateTime obsDate = modelObsDates[obsIdx];
                if (currentDate <= obsDate)
                    continue;

                addObservation(step, *ribCoupon, nbObs, couponIdx, obsIdx);
            }
        }
        if (nbObs==0) {
            throw ModelException("No RIB observation between ["
            +inst->accStart[couponIdx].toString()+" and "
            +inst->accEnd[couponIdx].toString());
        }
        if (nbObs > 1) {
            *ribCoupon /= (double)nbObs;
        }
    } 
    catch (exception& e) {
        throw makeException(e, __FUNCTION__);     
    }
}

void KRib2Tree::recordOutput(Control* ctrl, YieldCurveConstSP, Results* results) {
    try {
        if (!ctrl->requestsOutput(OutputRequest::DBG))
            return;

        if (model->getProdModifier<KRib2Tree::Compressor>().get()) {
            DateTimeArraySP dates(new DateTimeArray(modelObsDates));
            results->storeGreek(dates, Results::DEBUG_PACKET, 
                OutputNameConstSP(new OutputName("RIB_COMPRESSOR_DATES_"+inst->outputName)));
        }
    } 
    catch (exception& e) {
        throw makeException(e, __FUNCTION__);     
    }
}

/********************************** KRib2 **********************************/

void KRib2::assignObsDatesToCoupons(
        IntArray &convTable,
        DateTimeArray const &dates) const
{
    {   
        convTable.resize(dates.size());
        int cpnIdx = 0;
        DateTime accEndL = accEnd[0];
        for (int obsIdx = 0; obsIdx < dates.size(); ++obsIdx) {
            DateTime obsDate = dates[obsIdx];

            // move to the next cpnIdx if needed
            bool hasMoved = false;
            while (accEndL < obsDate
            || (includeAccStartObs && accEndL == obsDate)) {
                if (hasMoved) {
                    // two increments in a row -->error
                    throw ModelException("No rib observation for the "
                        +Format::toString(cpnIdx)+"-th coupon");
                }
                hasMoved = true;
                if (cpnIdx == accStart.size()-1) { // has reached the last coupon already
                    ASSERT(accEndL == obsDate);
                    break; // allowed, do not re-check, for this observation is attached to the last coupon anyway
                }
                // increment cpnIdx
                ++cpnIdx;
                accEndL = accEnd[cpnIdx];
            }
            convTable[obsIdx] = cpnIdx;
        }
    }
}

/** implement KComponent virtual function */
void KRib2::setup(const IModel* model, const MarketData* market) {
    try {
        KComponent::setup(model, market);

        if (discount.getName().size()) {
            throw ModelException("\"discount\" is not used by KRib2. "
            "It should not be provided.");
        }
        if (discountCurveName.empty()) {
            throw ModelException("\"discountCurveName\" not provided by parent component");
        }

        int nbCoupons = payDates.size();
        if (accStart.size() != nbCoupons)
            throw ModelException("accStart.size() != payDates.size()");
        if (accEnd.size() != nbCoupons)
            throw ModelException("accEnd.size() != payDates.size()");
        if (resetDates.size() != nbCoupons)
            throw ModelException("resetDates.size() != payDates.size()");

        // check obsDates
        if (obsDates.size() < 1) 
            throw ModelException("At least one  observation date is required");

        DateTime::ensureIncreasing(obsDates,"obsDates",true);

        if (obsDates[0] < accStart[0]) {
            throw ModelException("first observation date "+obsDates[0].toString()
            + " < first accrual start date " + accStart[0].toString());
        }
        if (accEnd.back() < obsDates.back()) {
            throw ModelException("last observation date "+obsDates.back().toString()
            + " > last accrual end date " + accEnd.back().toString());
        }

        assignObsDatesToCoupons(obsToCpn, obsDates);

        // check that payment date is after all corresponding observations
        DateTime::ensureIncreasing(payDates,"payDates",true);
        for (int obsIdx = obsDates.size()-1; obsIdx>=0; --obsIdx) {
            DateTime obsDate = obsDates[obsIdx];
            int cpnIdx = obsToCpn[obsIdx];
            DateTime payDate = payDates[cpnIdx];

            if (payDate < obsDate) {
                throw ModelException("Coupon #"+Format::toString(cpnIdx)
                    +" is attached to a rib observation "+obsDate.toString()
                    +" occuring before the payment date "
                    +payDate.toString());
            }
        }

        // initialize underlyings
        if (obsUnds.size()==1) {
            obsUnds[0]->addResetDates(obsDates);
            obsUnds[0]->setup(model, market);
        }
        else if (obsUnds.size() == obsDates.size()) {
            for (int obsIdx = obsDates.size()-1; obsIdx>=0; --obsIdx) {
                obsUnds[obsIdx]->addResetDate(obsDates[obsIdx]);
                obsUnds[obsIdx]->setup(model, market);
            }
        }
        else throw ModelException("obsUnds.size() must be 1 or nb pay dates ("
            +Format::toString(nbCoupons)+")");
    }
    catch (exception& e) {
        throw makeException(e, __FUNCTION__);
    }
}


FDProductSP KRib2::createProduct(FDModel * model) const {
    try {
        if (!setupCalled) 
            throw ModelException("steup() has not been called by parent component");
        return FDProductSP(new KRib2Tree(KRib2ConstSP(this), model));
    } 
    catch (exception& e) {
        throw makeException(e, __FUNCTION__);     
    }
}

double KRib2::getRibFraction(int cpnIdx, CashflowInfo &cfi ) const {
    try {
        double cumulObs=0.;
        int nbObs = 0;
        CashflowInfo cfiLoc;
        for (int obsIdx=0; obsIdx<obsDates.size(); ++obsIdx) {
            if (obsToCpn[obsIdx]!=cpnIdx)
                continue;
            int obsUndsIdx = (obsUnds.size() == 1 ? 0 : obsIdx);
            double observation;
            {
                CashflowInfo dummyCfi;
                observation = obsUnds[obsUndsIdx]->getValue(obsDates[obsIdx], dummyCfi);
                cfiLoc.updateAmountType(dummyCfi.amountType);
            }
            cumulObs += observation;
            ++nbObs;
            cfiLoc.push("observation["+Format::toString(obsIdx)+"] on "+obsDates[obsIdx].toString(), 0.5 < observation);
        }
        if (nbObs==0)
            throw ModelException("No observation for coupon "+Format::toString(cpnIdx));

        double ribFraction = cumulObs/nbObs;
        cfiLoc.push("ribFraction", Format::toString((int)(cumulObs+.5))+"/"+Format::toString(nbObs));
        cfi.merge(cfiLoc, outputName);
        return ribFraction;
    } 
    catch (exception& e) {
        throw makeException(e, __FUNCTION__);     
    }
}


void KRib2::load(CClassSP& clazz)
{
    clazz->setPublic();
    clazz->setDescription("Rib component");
    REGISTER(KRib2, clazz)
    SUPERCLASS(KComponent);
    EMPTY_SHELL_METHOD(KRib2::defaultConstructor);
    IMPLEMENTS(FDModel::IIntoProduct);

    FIELD(obsUnds, "RIB events observation composite index");
    FIELD(obsDates, "RIB actual observation date");
    FIELD(includeAccStartObs, "If true, an observation falling on an accrual "
        "start date belongs to the following period. Otherwise it belongs to "
        "the previous period.");
    FIELD_MAKE_OPTIONAL(includeAccStartObs);

    FIELD(accStart,"");           FIELD_MAKE_TRANSIENT(accStart);
    FIELD(accEnd,"");             FIELD_MAKE_TRANSIENT(accEnd);
    FIELD(resetDates,"");         FIELD_MAKE_TRANSIENT(resetDates);
    FIELD(payDates,"");           FIELD_MAKE_TRANSIENT(payDates);
    FIELD(discountCurveName,"");  FIELD_MAKE_TRANSIENT(discountCurveName);
    FIELD(obsToCpn, "");    FIELD_MAKE_TRANSIENT(obsToCpn);

    Addin::registerConstructor(Addin::UTILITIES, KRib2::TYPE);
}

CClassConstSP const KRib2::TYPE = CClass::registerClassLoadMethod(
    "KRib2", typeid(KRib2), KRib2::load);

/************************* KRib2Tree::Compressor ********************/

void KRib2Tree::Compressor::lop(DateTime today,
    DateTimeArray const &inDates, DateTimeArray &outDates, 
    IntArray const &obsToCpn) 
{
    outDates.clear();
    if (!inDates.size()) {
        return;
    }
    outDates.push_back(inDates[0]);
    if (inDates.size()==1) {
        return;
    }
    ASSERT(inDates.size()==obsToCpn.size());
    int i=0;
    // date before which we keep all dates similar
    DateTime verbatimDate = verbatimPeriod->toDate(today);
    // curPeriod will not go higher than targetPeriod
    double targetPeriod = 365.0 * targetFreq->toYears();
    // growth rate of curPeriod after exactDate
    double periodGrowth = ::log(2.)/(365. * decayHalfLife->toYears());
    int lastCouponId=-1;
    int lastDateIdx=0;
    while (i < inDates.size()-1) {
        ++i;
        DateTime curDate = inDates[i];

        if (verbatimDate < curDate && i!=(inDates.size()-1)) {
            // not the last date, and not protected by verbatimPeriod: we might want to drop this date
            int nb1 = curDate.daysDiff(outDates.back());
            int nb2 = (i<(inDates.size()-1) ? inDates[i+1].daysDiff(outDates.back()) : 0);
            if (nb1 < targetPeriod && (nb2-targetPeriod < targetPeriod-nb1)) {
                // not protected by targetPeriod and next date is closer to targetPeriod
                double curPeriod = ::exp(periodGrowth * curDate.daysDiff(verbatimDate));
                if ((i-lastDateIdx)<(curPeriod-0.5)) {
                    // not protected by curPeriod: we might want to drop this date

                    if (obsToCpn[i] == lastCouponId) {
                        continue; // we already have an observation for this coupon --> drop
                    }
                    if (i < inDates.size()-1
                        && (obsToCpn[i+1] == obsToCpn[i]))
                    {
                        // in the worst case we can still use the next observation for this coupon
                        // --> drop this one
                        continue; 
                    }
                }
            }
        }
        // [end] decide whether we want to keep the current date
        outDates.push_back(curDate); // save this date
        lastCouponId = obsToCpn[i];
        lastDateIdx = i;
    }
}

void KRib2Tree::Compressor::load(CClassSP& clazz)
{
    clazz->setPublic();
    clazz->setDescription("Rib Compressor component");
    REGISTER(KRib2Tree::Compressor, clazz)
    SUPERCLASS(CObject);
    IMPLEMENTS(FDModel::IProdModifier);
    EMPTY_SHELL_METHOD(KRib2Tree::Compressor::defaultConstructor);

    FIELD(verbatimPeriod, "Date (can be expressed as a maturity period from today) before which all observations must be kept.");
    FIELD(targetFreq, "");
    FIELD(decayHalfLife, "");

    Addin::registerConstructor(Addin::UTILITIES, KRib2Tree::Compressor::TYPE);
}

typedef KRib2Tree::Compressor KRib2TreeCompressor;
CClassConstSP const KRib2TreeCompressor::TYPE = CClass::registerClassLoadMethod(
    "KRib2Tree::Compressor", typeid(KRib2Tree::Compressor), KRib2TreeCompressor::load);

/*************************** KRib2Tree::Compressor::Tester ******************/

IObjectSP KRib2Tree::Compressor::Tester::run() 
{
    DateTimeArraySP outDates(new DateTimeArray);
    compressor->lop(today, inDates, *outDates, obsToCpn);
    return IObjectSP(outDates);
}

void KRib2Tree::Compressor::Tester::load(CClassSP& clazz)
{
    clazz->setPublic();
    clazz->setDescription("Rib Compressor tester component");
    REGISTER(KRib2Tree::Compressor::Tester, clazz)
    SUPERCLASS(CObject);
    EMPTY_SHELL_METHOD(KRib2Tree::Compressor::Tester::defaultConstructor);
    FIELD(today, "");
    FIELD(obsToCpn, "");
    FIELD(compressor, "");
    FIELD(inDates, "");

    Addin::registerObjectMethod(
            "KRib2Tree_Compressor_Tester", Addin::UTILITIES,
            "", false, Addin::expandSimple,
			&KRib2Tree::Compressor::Tester::run);
}

typedef KRib2Tree::Compressor::Tester KRib2TreeCompressorTester;
CClassConstSP const KRib2TreeCompressorTester::TYPE = CClass::registerClassLoadMethod(
    "KRib2Tree::Compressor::Tester", typeid(KRib2Tree::Compressor::Tester), KRib2TreeCompressorTester::load);

/**********************************/
bool KRib2Load()
{
    return KRib2::TYPE !=0;
}

DRLIB_END_NAMESPACE
