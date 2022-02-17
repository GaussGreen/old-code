//----------------------------------------------------------------------------
//
//   Group       : QR&D Interest Rates
//
//   Filename    : CallableRib.cpp
//
//   Description : instrument shell that makes use of components 
//                 with customized internal logic to fill optional data fields
//                 for each components. 
//
//----------------------------------------------------------------------------

#include "edginc/config.hpp"
#include "edginc/KRibTurbo.hpp"
#include "edginc/KRib2.hpp"
#include "edginc/KKnockOut.hpp"
#include "edginc/QuasiMQ.hpp"
#include "edginc/Business252.hpp"
#include "edginc/Results.hpp"
#include "edginc/IrConverter.hpp"
#include "edginc/IndexSpecIR.hpp"
#include "edginc/PhysicalSettlement.hpp"
#include "edginc/KKnockOut.hpp"
#include "edginc/KSum.hpp"

DRLIB_BEGIN_NAMESPACE


/********************************** Declarations *******************************/

// NEW version of the shell instrument for CallableRib
// It looks like we could factor out a new base class for "shell instruments"
// that are based on the KComponent architecture and move common functionality
// into the new base.
class SimpleRib : public CInstrument,
                  virtual public FDModel::IIntoProduct,
                  virtual public QuasiMQ::IIntoProduct,
                  virtual public LastSensDate
{
public:
    static CClassConstSP const TYPE;

    /* CInstrument:: */
    virtual void Validate(void) {}

    virtual DateTime getValueDate() const
    {
        return compRoot->getValueDate();
    }
    
    virtual string discountYieldCurveName() const
    {
        return discount.getName();
    }
    
    virtual bool priceDeadInstrument(CControl* control, CResults* results) const
    {
        return compRoot->priceDeadInstrument(control, results);
    }
    
    virtual void GetMarket(const IModel *model, const CMarketDataSP market);

    /* FDModel::IIntoProduct:: */
    virtual FDProductSP createProduct(FDModel* model) const;

    /* Quasi::IIntoProduct:: */
    virtual QuasiMQ::ProductSP createProduct(QuasiMQ* model) const;

    /* LastSensDate */
    virtual DateTime endDate(const Sensitivity* sensControl) const
    {
        return compRoot->endDate(sensControl);
    }

private:
    friend class SimpleRibMQProd;

    SimpleRib() : CInstrument(TYPE) {}
    static void load(CClassSP& clazz);
    static IObject* defaultConstructor(void) { return new SimpleRib(); }

    /************** exported fields ************/
    KRibTurboSP         paymentLeg;
    KRib2SP             rib;
    YieldCurveWrapper   discount;

    /************* transiend fields ************/
    KComponentSP        compRoot;  // for FDModel
};

/*********************/

class SimpleRibMQProd : public QuasiMQ::Product
{
    /*********** variables **********/

    SimpleRib const * inst;

    /************ methods ************/

    void getPaySpecs(std::vector<IndexSpecIRSP>& paySpecs) const;
    void getRIBSpecs(std::vector<IndexSpecIRSP>& ribSpecs) const;

    // setups the inputs required for MQ pricing. return argument is number of variates
    int setupRIBPricingInputs(int                               cpnIndex,
                              const DateTime&                   obsDate,
                              const std::vector<IndexSpecIRSP>& paySpecs,
                              const std::vector<IndexSpecIRSP>& ribSpecs,
                              IndexSpecIRSP&                    rateSpec1,
                              DateTime&                         rateExpiry1,
                              IndexSpecIRSP&                    rateSpec2,
                              DateTime&                         rateExpiry2,
                              DoubleArray&                      payoffParams,
                              long&                             optionType) const;
                              

public:
    SimpleRibMQProd(QuasiMQ* model, SimpleRib const * inst);
    virtual void price(Control* control, CResults*  results);

};

/******************************* SimpleRibMQProd ***************************/

SimpleRibMQProd::SimpleRibMQProd(QuasiMQ* model, SimpleRib const* inst) : QuasiMQ::Product(model), inst(inst)
{}


void SimpleRibMQProd::price(Control* control, CResults*  results)
{
    try
    {
        //////////////////////////////////////////////////////
        //                                                  //
        // This data does not change inside the loop        //
        //                                                  //
        //////////////////////////////////////////////////////
        DateTime                localValDt = model->spotDate;
        DateTime                localToday = model->today;
        long optType;

        //////////////////////////////////////////////////////
        //                                                  //
        // This data is coupons dependent                   //
        //                                                  //
        //////////////////////////////////////////////////////
        // ??? make more sensible later
        KRibTurboConstSP paymentLeg = KRibTurboConstSP::dynamicCast(inst->compRoot);
        if (!paymentLeg.get())
            throw ModelException("Internal error - unable to cast compRoot to KRibTurbo type");

        const CouponSchedDates &sched = *paymentLeg->sched;

        //These are all output arrays (i.e. for each coupon)
        DoubleArraySP undiscPriceSP(new CDoubleArray(sched.nbCoupons));
        DoubleArraySP couponStrikesSP(new CDoubleArray(sched.nbCoupons));
        DoubleArraySP dscFactorsSP(new CDoubleArray(sched.nbCoupons));
        DoubleArraySP yearFractionSP(new CDoubleArray(sched.nbCoupons));
//        IntArraySP          couponNbSP(       new CIntArray(sched.nbCoupons));

        KRibTurbo::IndexLeg& leg = *(*paymentLeg->legs)[0];
        
        QuasiMQ::VariateDebugData debugData1; // IndexSpecIR1 debug data
        QuasiMQ::VariateDebugData debugData2; // IndexSpecIR2 debug data - bivariate only

        // Payment IndexSpecIRs
        std::vector<IndexSpecIRSP> paymentSpecs;
        getPaySpecs(paymentSpecs);

        // RIB IndexSpecIRs
        std::vector<IndexSpecIRSP> ribSpecs;
        getRIBSpecs(ribSpecs);

        DayCountConventionSP busDCC(new Business252);
        InstrumentSettlementSP settleType(new PhysicalSettlement);

        double dscPrice = 0;
        const DayCountConventionSP dcc = leg.dcc;
        
        ///////////////////////////////////////////////////////////////////////////////////////////////////////////
        // Price each coupon
        ///////////////////////////////////////////////////////////////////////////////////////////////////////////
        for (int i = 0; i < sched.nbCoupons; ++i)
        {
            double price, fwdRate;
            double currCoupon = 0;

            if (!inst->rib.get())
            {
                // Despite this being a "SimpleRib" class, we are allowing this class to
                // be able to price non-rib products - i.e. standard puts, calls, collars for
                // both uni-variate and bi-variate cases - actually bi-variate doesn't work due to
                // restrictions of KRibTurbo on number of IndexSpecIRs
                int numPaySpecs = Maths::isZero((*leg.weights)[i]) ? 0 : paymentSpecs.size();

                switch (numPaySpecs)
                {
                    case 1:
                    {
                        if (paymentLeg->floorRates.get())
                        {
                            model->univarPricer(*paymentSpecs[0],
                                                sched.reset[i],
                                                sched.pay[i],
                                                (*paymentLeg->floorRates)[i],
                                                busDCC,
                                                2,
                                                DoubleArray(1),
                                                settleType,
                                                price,
                                                fwdRate,
                                                debugData1);
                            currCoupon += price;
                        }
                        if (paymentLeg->capRates.get())
                        {
                            model->univarPricer(*paymentSpecs[0],
                                                sched.reset[i],
                                                sched.pay[i],
                                                (*paymentLeg->capRates)[i],
                                                busDCC,
                                                1,
                                                DoubleArray(1),
                                                settleType,
                                                price,
                                                fwdRate,
                                                debugData1);
                            currCoupon += price;
                        }
                        break;
                    }
                    case 2:
                    {
                        DoubleArray payoffParams(2);
                        payoffParams[1] = (*leg.weights)[i];
                        if (paymentLeg->floorRates.get())
                        {
                            payoffParams[0] = (*paymentLeg->floorRates)[i];
                            model->bivarPricer(*paymentSpecs[0],
                                               sched.reset[i],
                                               *paymentSpecs[1],
                                               sched.reset[i],
                                               sched.pay[i],
                                               busDCC,
                                               302,
                                               payoffParams,
                                               settleType,
                                               price,
                                               fwdRate,
                                               debugData1,
                                               debugData2);
                            currCoupon += price;
                        }
                        if (paymentLeg->capRates.get())
                        {
                            payoffParams[0] = (*paymentLeg->capRates)[i];
                            model->bivarPricer(*paymentSpecs[0],
                                               sched.reset[i],
                                               *paymentSpecs[1],
                                               sched.reset[i],
                                               sched.pay[i],
                                               busDCC,
                                               301,
                                               payoffParams,
                                               settleType,
                                               price,
                                               fwdRate,
                                               debugData1,
                                               debugData2);
                            currCoupon += price;
                        }
                        break;
                    }
                    default:
                    {
                        throw ModelException("Unable to price non-RIB product with " + Format::toString(numPaySpecs) + " payment IndexSpecs");
                        break;
                    }
                }
            }
            else
            {
                // count number of rib observations in current accrual period and loop over them
                const FlexDates& obsDates = inst->rib->obsDates;

                // if includeAccStartObs, then [accStart, accEnd), else (accStart, accEnd]
                int obsStart = (inst->rib->includeAccStartObs ? std::lower_bound(obsDates.begin(), obsDates.end(), sched.accStart[i])
                                                              : std::upper_bound(obsDates.begin(), obsDates.end(), sched.accStart[i])) - obsDates.begin();

                int obsEnd   = (inst->rib->includeAccStartObs ? std::lower_bound(obsDates.begin(), obsDates.end(), sched.accEnd[i])
                                                              : std::upper_bound(obsDates.begin(), obsDates.end(), sched.accEnd[i])) - obsDates.begin();

                if (obsStart == obsDates.size() || obsStart == obsEnd || obsDates[obsStart] > sched.accEnd[i])
                {
                    throw ModelException("No observation dates for accrual period " + sched.accStart[i].toString() + " to " + sched.accEnd[i].toString());
                }

                for (int j = obsStart; j < obsEnd; ++j)
                {
                    IndexSpecIRSP rateSpec1;
                    DateTime rateExpiry1;

                    IndexSpecIRSP rateSpec2;
                    DateTime rateExpiry2;

                    DoubleArray payoffParams;

                    int numVariates = setupRIBPricingInputs(i,
                                                            obsDates[j],
                                                            paymentSpecs,
                                                            ribSpecs,
                                                            rateSpec1,
                                                            rateExpiry1,
                                                            rateSpec2,
                                                            rateExpiry2,
                                                            payoffParams,
                                                            optType);

                    switch (numVariates)
                    {
                        case 1:
                            model->univarPricer(*rateSpec1,
                                                rateExpiry1,
                                                sched.pay[i],
                                                (*couponStrikesSP)[i],
                                                busDCC,
                                                optType,
                                                payoffParams,
                                                settleType,
                                                price,
                                                fwdRate,
                                                debugData1);
                            break;
                        case 2:
                            model->bivarPricer(*rateSpec1,
                                               rateExpiry1,
                                               *rateSpec2,
                                               rateExpiry2,
                                               sched.pay[i],
                                               busDCC,
                                               optType,
                                               payoffParams,
                                               settleType,
                                               price,
                                               fwdRate,
                                               debugData1,
                                               debugData2);
                            break;
                        default:
                            throw ModelException("Unable to price " + Format::toString(numVariates) + "-variate RIB product");
                            break;
                    }

                    currCoupon += price / (obsEnd - obsStart);

                }// for j, riblets
            } // end if rib

            (*yearFractionSP)[i] = dcc->years(sched.accStart[i], sched.accEnd[i]);
            currCoupon *= (*yearFractionSP)[i] * (*leg.notionals)[i];

            (*dscFactorsSP)[i] = model->getZero(localValDt, sched.pay[i], inst->discountYieldCurveName());
            (*undiscPriceSP)[i] = currCoupon;
            dscPrice += currCoupon * (*dscFactorsSP)[i];

        }//for i, coupons
        
        //generate output
        results->storePrice(dscPrice, model->domesticYC->getCcy());
        results->storeScalarGreek(dscPrice, Results::DEBUG_PACKET, OutputNameSP(new OutputName("DISCOUNDED_PRICE")));
        results->storeGreek(undiscPriceSP, Results::DEBUG_PACKET, OutputNameSP(new OutputName("UNDISCOUNTED_PREMIUM")));
        results->storeGreek(dscFactorsSP, Results::DEBUG_PACKET, OutputNameSP(new OutputName("DISCOUNT_FACTOR")));

        results->storeGreek(DateTimeArraySP(new DateTimeArray(sched.pay)),      Results::DEBUG_PACKET, OutputNameSP(new OutputName("PAYMENT_DATE")));
        results->storeGreek(DateTimeArraySP(new DateTimeArray(sched.accStart)), Results::DEBUG_PACKET, OutputNameSP(new OutputName("ACC_START")));
        results->storeGreek(DateTimeArraySP(new DateTimeArray(sched.accEnd)),   Results::DEBUG_PACKET, OutputNameSP(new OutputName("ACC_END")));
        results->storeGreek(yearFractionSP,                                     Results::DEBUG_PACKET, OutputNameSP(new OutputName("ACC_DCF")));
        results->storeGreek(DoubleArraySP(new DoubleArray(*leg.notionals)),     Results::DEBUG_PACKET, OutputNameSP(new OutputName("NOTIONAL")));

        results->storeGreek(couponStrikesSP, Results::DEBUG_PACKET, OutputNameSP(new OutputName("STRIKES")));
//        results->storeGreek(couponNbSP, Results::DEBUG_PACKET, OutputNameSP(new OutputName("COUPON_NBS")));

        // riblet Debug data
        results->storeGreek(debugData1.expiryDates, Results::DEBUG_PACKET, OutputNameSP(new OutputName("RIBLET_SPEC1_EXP")));
        results->storeGreek(debugData1.startDates,  Results::DEBUG_PACKET, OutputNameSP(new OutputName("RIBLET_SPEC1_START")));
        results->storeGreek(debugData1.endDates,    Results::DEBUG_PACKET, OutputNameSP(new OutputName("RIBLET_SPEC1_END")));
        results->storeGreek(debugData1.payDates,    Results::DEBUG_PACKET, OutputNameSP(new OutputName("RIBLET_SPEC1_PAY")));
        results->storeGreek(debugData1.atmVols,     Results::DEBUG_PACKET, OutputNameSP(new OutputName("RIBLET_SPEC1_ATMVOL")));
        results->storeGreek(debugData1.prices,      Results::DEBUG_PACKET, OutputNameSP(new OutputName("RIBLET_SPEC1_PX")));
        results->storeGreek(debugData1.fwdRates,    Results::DEBUG_PACKET, OutputNameSP(new OutputName("RIBLET_SPEC1_FWDRATE")));

        results->storeGreek(debugData2.expiryDates, Results::DEBUG_PACKET, OutputNameSP(new OutputName("RIBLET_SPEC2_EXP")));
        results->storeGreek(debugData2.startDates,  Results::DEBUG_PACKET, OutputNameSP(new OutputName("RIBLET_SPEC2_START")));
        results->storeGreek(debugData2.endDates,    Results::DEBUG_PACKET, OutputNameSP(new OutputName("RIBLET_SPEC2_END")));
        results->storeGreek(debugData2.payDates,    Results::DEBUG_PACKET, OutputNameSP(new OutputName("RIBLET_SPEC2_PAY")));
        results->storeGreek(debugData2.atmVols,     Results::DEBUG_PACKET, OutputNameSP(new OutputName("RIBLET_SPEC2_ATMVOL")));
        results->storeGreek(debugData2.prices,      Results::DEBUG_PACKET, OutputNameSP(new OutputName("RIBLET_SPEC2_PX")));
        results->storeGreek(debugData2.fwdRates,    Results::DEBUG_PACKET, OutputNameSP(new OutputName("RIBLET_SPEC2_FWDRATE")));
    }
    catch (exception& e)
    {
        throw ModelException(e, __FUNCTION__);
    }
}

void SimpleRibMQProd::getPaySpecs(std::vector<IndexSpecIRSP>& paySpecs) const
{
    paySpecs.clear();

    KRibTurboConstSP paymentLeg = KRibTurboConstSP::dynamicCast(inst->compRoot); // guaranteed to work - see SimpleRib::GetMarket
    KRibTurbo::IndexLeg& leg = *(*paymentLeg->legs)[0];

    // plonk this somewhere else - e.g. in an object validation function
    KSumConstSP ksum = KSumConstSP::dynamicCast(leg.indexes[0]);
    if (!ksum.get()) throw ModelException("KRibTurbo::IndexLeg.indexes components should be KSum objects");

    CModel::IProdCreatorArray indexes = ksum->getListK();
    
    for (int j = 0; j < indexes.size(); ++j)
    {
        IndexSpecIRSP paySpec = IndexSpecIRSP::dynamicCast(indexes[j]);
        if (paySpec.get())
        {
            paySpecs.push_back(paySpec);
        }
        else
        {
            throw ModelException("Floating index on rib payment leg must be indexSpecIR, type supplied = " +
                                 indexes[j]->getClass()->getName());
        }
    }
}

void SimpleRibMQProd::getRIBSpecs(std::vector<IndexSpecIRSP>& ribSpecs) const
{
    ribSpecs.clear();

    if (inst->rib.get())
    {
        KKnockOutSP ko = KKnockOutSP::dynamicCast(inst->rib->obsUnds[0]);
        for (int j = 0; j < ko->knockouts.size(); ++j)
        {
            IndexSpecIRSP rSpec = IndexSpecIRSP::dynamicCast(ko->knockouts[j]->und);
            if (rSpec.get())
            {
                ribSpecs.push_back(rSpec);
            }
            else
            {
                throw ModelException("Rib observation index must be indexSpecIR, type supplied = " +
                                     ko->knockouts[j]->und->getClass()->getName());
            }
        }
    }
}

int SimpleRibMQProd::setupRIBPricingInputs(int                               cpnIndex,
                                           const DateTime&                   obsDate,
                                           const std::vector<IndexSpecIRSP>& paySpecs,
                                           const std::vector<IndexSpecIRSP>& ribSpecs,
                                           IndexSpecIRSP&                    rateSpec1,
                                           DateTime&                         rateExpiry1,
                                           IndexSpecIRSP&                    rateSpec2,
                                           DateTime&                         rateExpiry2,
                                           DoubleArray&                      payoffParams,
                                           long&                             optionType) const
{
    payoffParams.clear();
    int numVariates = 0;

    KKnockOutSP ko = KKnockOutSP::dynamicCast(inst->rib->obsUnds[0]);

    KRibTurboConstSP payLeg = KRibTurboConstSP::dynamicCast(inst->compRoot);  // guaranteed to work - see SimpleRib::GetMarket
    KRibTurbo::IndexLegConstSP indexLeg = (*payLeg->legs)[0];

    int numPaySpecs = Maths::isZero((*indexLeg->weights)[cpnIndex]) ? 0 : paySpecs.size();
    
    rateExpiry1 = obsDate; // initial guess - overridden further down if wrong
    rateExpiry2 = obsDate; // initial guess - overridden further down if wrong

    payoffParams.push_back(ko->knockouts[0]->loBarrier->interpolate(obsDate)); // LB1
    payoffParams.push_back(ko->knockouts[0]->hiBarrier->interpolate(obsDate)); // HB1

    switch (ribSpecs.size())
    {
        case 1:  // ribSepcs.size() == 1
        {
            switch (numPaySpecs)
            {
                case 0:   // ribSpecs.size() == 1, paymentSpecs.size() == 0
                {
                    // univariate Fixed Rib
                    numVariates = 1;
                    rateSpec1 = ribSpecs[0];

                    payoffParams.push_back(0.0); // LE1
                    payoffParams.push_back(0.0); // HE1

                    switch (ko->barrierType)
                    {
                        case KKnockOut::Barrier::ONE_IN:
                        case KKnockOut::Barrier::IN:
                            optionType = 13;
                            break;
                        case KKnockOut::Barrier::ONE_OUT:
                        case KKnockOut::Barrier::OUT:
                            optionType = 14;
                            break;
                        default:
                            throw ModelException("Invalid barrierType set");
                            break;
                    }

                    break;
                } 
                case 1:   // ribSpecs.size() == 1, paymentSpecs.size() == 1
                {
                    // inside/outside collared rib
                    numVariates = 2;
                    rateSpec1 = paySpecs[0];
                    rateSpec2 = ribSpecs[0];
                    rateExpiry1 = payLeg->sched->reset[cpnIndex];
                    payoffParams.push_back((*indexLeg->weights)[cpnIndex]);   // (L)everage
                    payoffParams.push_back((*indexLeg->spreads)[cpnIndex]);   // (S)pread
                    payoffParams.push_back((*payLeg->floorRates)[cpnIndex]);  // (F)loor
                    payoffParams.push_back((*payLeg->capRates)[cpnIndex]);    // (C)ap

                    switch (ko->barrierType)
                    {
                        case KKnockOut::Barrier::ONE_IN:
                        case KKnockOut::Barrier::IN:
                            optionType = 313;
                            break;
                        case KKnockOut::Barrier::ONE_OUT:
                        case KKnockOut::Barrier::OUT:
                            optionType = 314;
                            break;
                        default:
                            throw ModelException("Invalid barrierType set");
                            break;
                    }

                    break;
                }
                default:  // ribSpecs.size() == 1, paymentSpecs.size() == ???
                {
                    throw ModelException("Number of payment IndexSpecs must be 0 or 1, supplied : " +
                                        Format::toString(ribSpecs.size()));
                    break;
                }
            }
            break;
        }

        case 2: // ribSpecs.size() == 2
        {
            switch (numPaySpecs)
            {
                case 0:   // ribSpecs.size() == 2, paymentSpecs.size() == 0
                {
                    // Dual Rib
                    numVariates = 2;
                    rateSpec1 = ribSpecs[0];
                    rateSpec2 = ribSpecs[1];
                    payoffParams.push_back(ko->knockouts[1]->loBarrier->interpolate(obsDate)); // LB2
                    payoffParams.push_back(ko->knockouts[1]->hiBarrier->interpolate(obsDate)); // HB2

                    switch (ko->barrierType)
                    {
                        case KKnockOut::Barrier::IN_IN:
                            optionType = 361;
                            break;
                        case KKnockOut::Barrier::ONE_IN:
                            optionType = 362;
                            break;
                        case KKnockOut::Barrier::IN_OUT: // interpreted as "in and out", as opposed to "in or out" (optionType 364)
                            optionType = 363;
                            break;
                        case KKnockOut::Barrier::OUT_IN: // interpreted as "out and in", as opposed to "out or in" (optionType 366)
                            optionType = 365;
                            break;
                        case KKnockOut::Barrier::OUT_OUT:
                            optionType = 367;
                            break;
                        case KKnockOut::Barrier::ONE_OUT:
                            optionType = 368;
                            break;
                        default:
                            throw ModelException("Invalid barrierType set");
                            break;
                    }

                    break;
                } 
                case 1:   // ribSpecs.size() == 2, paymentSpecs.size() == 1
                {
                    numVariates = 2;
                    rateSpec1 = paySpecs[0];

                    // for 3 IndexSpecIRs - resetDate for payment must match RIB date for a bi-variate product
                    // otherwise it is a different underlying and the product is tri-variate
                    if (payLeg->sched->reset[cpnIndex] == obsDate)
                    {
                        // the IndexSpecIR for payment index must match one of the RIB indexes for bi-variate
                        if (paySpecs[0]->tenor->equalTo(ribSpecs[0]->tenor.get()) && paySpecs[0]->frequency->equalTo(ribSpecs[0]->frequency.get()))
                        {
                            rateSpec2 = ribSpecs[1];
                        }
                        else if (paySpecs[0]->tenor->equalTo(ribSpecs[1]->tenor.get()) && paySpecs[0]->frequency->equalTo(ribSpecs[1]->frequency.get()))
                        {
                            rateSpec2 = ribSpecs[0];
                        }
                        else
                        {
                            throw ModelException("Unable to price trivariate product");
                        }
                    }
                    else
                    {
                        throw ModelException("Unable to price trivariate product");
                    }

                    payoffParams.push_back((*indexLeg->weights)[cpnIndex]);   // (L)everage
                    payoffParams.push_back((*indexLeg->spreads)[cpnIndex]);   // (S)pread
                    payoffParams.push_back((*payLeg->floorRates)[cpnIndex]);  // (F)loor
                    payoffParams.push_back((*payLeg->capRates)[cpnIndex]);    // (C)ap

                    switch (ko->barrierType)
                    {
                        case KKnockOut::Barrier::ONE_IN:
                        case KKnockOut::Barrier::IN:
                            optionType = 323;
                            break;
                        case KKnockOut::Barrier::ONE_OUT:
                        case KKnockOut::Barrier::OUT:
                            optionType = 324;
                            break;
                        default:
                            throw ModelException("Invalid barrierType set");
                            break;
                    }

                    break;
                }
                default:  // ribSpecs.size() == 2, paymentSpecs.size() == ???
                {
                    throw ModelException("Number of payment IndexSpecs must be 0 or 1, supplied : " +
                                         Format::toString(ribSpecs.size()));
                    break;
                }
            }
            break;
        }

        default: // ribSpecs.size() == ???
        {
            throw ModelException("Number of rib IndexSpecs, i, must be 0 < i < 3, supplied : " +
                                 Format::toString(ribSpecs.size()));
            break;
        }
    }

    return numVariates;
}

/******************************* SimpleRibMQProd ***************************/

FDProductSP SimpleRib::createProduct(FDModel* model) const
{
    return model->createProduct(compRoot);
}

QuasiMQ::ProductSP SimpleRib::createProduct(QuasiMQ* model) const
{
    return QuasiMQ::ProductSP(new SimpleRibMQProd(model, this));
}

// export the members of the class through the library interface
void SimpleRib::load(CClassSP& clazz)
{
    clazz->setPublic();
    REGISTER(SimpleRib, clazz);
    SUPERCLASS(CInstrument);
    IMPLEMENTS(FDModel::IIntoProduct);
    EMPTY_SHELL_METHOD(defaultConstructor);

    // exported fields
    FIELD(paymentLeg, "payment leg for rib - may be either fixed or floating rate");
    FIELD(rib, "rib component");
    FIELD_MAKE_OPTIONAL(rib);
    FIELD(discount,"pricing discount curve for product");
    
    // transient fields - for FDModel
    FIELD(compRoot,""); 
    FIELD_MAKE_TRANSIENT(compRoot);
}

// construct and link the rib component dependencies in instrumet::getMarket
void SimpleRib::GetMarket(const IModel* model, const CMarketDataSP market)
{
    try
    {
        discount.getData(model, market.get());

        if (paymentLeg->outputName.empty())
            paymentLeg->outputName = "PAYMENT_LEG";
        if (rib.get() && rib->outputName.empty())
            rib->outputName = "RIB";
        //if (observationEvent->outputName.empty())
        //    observationEvent->outputName = "OBSERVATION_EVENT";

        if (paymentLeg->discountYieldCurveName().empty())
            paymentLeg->discount = discount;

        // setup component heirachy used by tree pricer
        paymentLeg->rib = rib;

        compRoot = paymentLeg;

        if (compRoot->discountYieldCurveName().empty())
            compRoot->discount = discount;

        // additional checks on the components
        if (paymentLeg->legs->size()!=1)
        {
            throw ModelException("paymentLeg.legs.size() != 1. Is " + Format::toString(paymentLeg->legs->size()));
        }

        KRibTurbo::IndexLeg &leg = *(*paymentLeg->legs)[0];
        int nbCoupons = paymentLeg->sched->nbCoupons;

        //for (int idxIdx=0; idxIdx < leg.indexes.size(); ++idxIdx)
        //{
        //    IndexSpecIR *spec = dynamic_cast<IndexSpecIR*>(leg.indexes[idxIdx].get());
        //    if (!spec)
        //    {
        //        throw ModelException("Index for indexes[" + Format::toString(idxIdx) + "] of paymentLeg->legs[0] is not IndexSpecIR");
        //    }
        //}

        // cleanup fields so that iterator do not get lost
        paymentLeg.reset(0);
        //rib.reset(0);
        //observationEvent.reset(0);

        compRoot->GetMarket(model, market);

        // Implement cross component checks
        // ??? toDo - need access to private members to do any sensible
        // checks, but mad to make all components friends of shell instruments
    }
    catch (exception& e)
    {
        throw ModelException(e, __FUNCTION__);
    }
}

// register the SimpleRib in the library framework at startup
CClassConstSP const SimpleRib::TYPE = CClass::registerClassLoadMethod("SimpleRib", typeid(SimpleRib), SimpleRib::load);

// to ensure linker doesn't optimise out the class
bool SimpleRibLoad(void)
{
    return (SimpleRib::TYPE != 0);
}

DRLIB_END_NAMESPACE
