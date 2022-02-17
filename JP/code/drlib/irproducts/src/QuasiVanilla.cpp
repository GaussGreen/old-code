//----------------------------------------------------------------------------
//
//   Group       : QR Rates
//
//   Filename    : QuasiVanilla.cpp
//
//   Description : First attempt at implementing quasi-vanilla in qlib
//                 subject to change in the future
//
//----------------------------------------------------------------------------

#include "edginc/config.hpp"
#include "edginc/QuasiVanilla.hpp"
#include "edginc/IrConverter.hpp"
#include "edginc/MarketData.hpp"
#include "edginc/IInstrumentCollection.hpp"
#include "edginc/Results.hpp"
#include "edginc/VolRequestRaw.hpp"
#include "edginc/Business252.hpp"
#include "edginc/PhysicalSettlement.hpp"


DRLIB_BEGIN_NAMESPACE

/********************************* UnivarIRCpn *********************************/

void UnivarIRCpn::UnivarIRCpnProdMQ::price(
                  Control*   control,
                  CResults*  results) {
    try {
        //////////////////////////////////////////////////////
        //                                                  //
        // This data does not change inside the loop        //
        //                                                  //
        //////////////////////////////////////////////////////
        double                  dscPrice = 0;
        DateTimeArray           zeroDatesIdx;
        DoubleArray             zeroRatesIdx;
        DateTimeArray           zeroDatesDisc;
        DoubleArray             zeroRatesDisc;
        DoubleArray             outputs(2);
        long                    rateFreq;
        string                  rateDCC;
        int                     numPayoffParams = inst->payoffParams.getLength();
        DayCountConventionSP    localExpiryDCC = inst->expiryDCC;
        DateTime                localValDt = model->getValueDate();
        DateTime                localToday = model->today;
        int                     idx=0, i = 0;
        CouponSchedDatesSP      localSched = inst->sched;

        //      if (cms)
        {
            rateFreq = inst->rate->frequency->approxAnnualFrequency();
        }
//      else
//      {
//          rateFreq = 0;
//      }

        rateDCC = inst->rate->dcc->toString();

        //Get zero curves
//SAMY      IrConverter::q3GetZerosData(zeroDatesIdx, zeroRatesIdx, model->idxCurve.get());
//SAMY      IrConverter::q3GetZerosData(zeroDatesDisc, zeroRatesDisc, model->discCurve.get());


        //////////////////////////////////////////////////////
        //                                                  //
        // This data is coupons dependent                   //
        //                                                  //
        //////////////////////////////////////////////////////
        int nbCoupons = localSched->nbCoupons;
        //These are all output arrays (i.e. for each coupon)
        DoubleArraySP       premiumsSP(new CDoubleArray(nbCoupons));
        DoubleArraySP       fwdRatesSP(new CDoubleArray(nbCoupons));
        DoubleArraySP       scaledSigATMsSP(new CDoubleArray(nbCoupons));
        DoubleArraySP       couponNotionalsSP(new CDoubleArray(nbCoupons));
        DoubleArraySP       couponStrikesSP(new CDoubleArray(nbCoupons));
        DoubleArraySP       dscFactorsSP(new CDoubleArray(nbCoupons));// = (model->discCurve).pv(localValDt, localPaymentDate);
        DateTimeArraySP     expiryDatesSP(new DateTimeArray(nbCoupons));
        DateTimeArraySP     startDatesSP(new DateTimeArray(nbCoupons));
        DateTimeArraySP     matDatesSP(new DateTimeArray(nbCoupons));
        DateTimeArraySP     payDatesSP(new DateTimeArray(nbCoupons));
        DateTimeArraySP     accrualStartSP(new DateTimeArray(nbCoupons));
        DateTimeArraySP     accrualEndSP(new DateTimeArray(nbCoupons));
        DoubleArray         smile(4);
		double				dcf;

        QuasiMQ::VariateDebugData debug;

		for (i=0; i< nbCoupons; ++i)
        {
            (*accrualStartSP)[i]    = localSched->accStart[i];
            (*accrualEndSP)[i]      = localSched->accEnd[i];
            (*expiryDatesSP)[i]     = localSched->reset[i];
            (*startDatesSP)[i]      = localSched->resetEff[i];
			(*matDatesSP)[i]		= (inst->rate->tenor)->toDate( (*startDatesSP)[i] );
            (*payDatesSP)[i]        = localSched->pay[i];
            (*couponNotionalsSP)[i] = inst->notionals[i];
            (*couponStrikesSP)[i]   = inst->strikes[i];
			dcf						= (inst->accrualDCC)->years((*accrualStartSP)[i], (*accrualEndSP)[i]);

            YieldCurveWrapper factor;
            IndexSpecIR rateSpec("UnivarIRCpn rateSpec", factor); // rateFreq, &(rateDCC[0]), ((*matDatesSP)[i]),
            // to do: populate: rateSpec.frequency..

            model->univarPricer(
                *inst->rate,
                ((*expiryDatesSP)[i]),
                ((*payDatesSP)[i]),
                (*couponStrikesSP)[i],
                DayCountConventionSP(new Business252),
                inst->payoffType,
                inst->payoffParams,
                InstrumentSettlementSP(new PhysicalSettlement),//              setlType,
                outputs[0],
                outputs[1],
                debug);

			(*premiumsSP)[i] = outputs[0] * (*couponNotionalsSP)[i] * dcf;
            (*fwdRatesSP)[i] = outputs[1];
            (*dscFactorsSP)[i] = ( (model->domesticYC).get() )->pv(localValDt, (*payDatesSP)[i]);
            dscPrice += outputs[0] * (*dscFactorsSP)[i];
        }

        // Create Caplet output
        results->storePrice(dscPrice, "USD");
        results->storeScalarGreek(dscPrice, Results::DEBUG_PACKET,
                                  OutputNameSP(new OutputName("DISCOUNDED_PRICE")));
        results->storeGreek(premiumsSP, Results::DEBUG_PACKET,
                                  OutputNameSP(new OutputName("UNDISCOUNTED_PREMIUM")));
        results->storeGreek(fwdRatesSP, Results::DEBUG_PACKET,
                                  OutputNameSP(new OutputName("UNADJUSTED_FWD_RATE")));
        results->storeGreek(dscFactorsSP, Results::DEBUG_PACKET,
                                  OutputNameSP(new OutputName("DISCOUNT_FACTOR")));
        results->storeGreek(expiryDatesSP, Results::DEBUG_PACKET,
                                  OutputNameSP(new OutputName("EXPIRY")));
        results->storeGreek(startDatesSP, Results::DEBUG_PACKET,
                                  OutputNameSP(new OutputName("RATE_START")));
        results->storeGreek(matDatesSP, Results::DEBUG_PACKET,
                                  OutputNameSP(new OutputName("RATE_MAT")));
        results->storeGreek(payDatesSP, Results::DEBUG_PACKET,
                                  OutputNameSP(new OutputName("PAYMENT_DATE")));
        results->storeGreek(accrualStartSP, Results::DEBUG_PACKET,
                                  OutputNameSP(new OutputName("ACC_START")));
        results->storeGreek(accrualEndSP, Results::DEBUG_PACKET,
                                  OutputNameSP(new OutputName("ACC_END")));
        results->storeGreek(couponNotionalsSP, Results::DEBUG_PACKET,
                                  OutputNameSP(new OutputName("NOTIONAL")));
        results->storeGreek(couponStrikesSP, Results::DEBUG_PACKET,
                                  OutputNameSP(new OutputName("STRIKES")));
    }
    catch (exception& e) {
        throw ModelException(e, __FUNCTION__);
    }
}



void UnivarIRCpn::validatePop2Object() {
    try {
    /*  if (start > end) {
            start = start;
            throw ModelException("Start must be before end");
        }

        if (notionals.getLength() != (amortisationDates).getLength())
            throw ModelException("Nb of notional step up dates must equal nb Notionals");

        if (   amortisationDates[ (amortisationDates).getLength()-1 ] < end   )
            throw ModelException("Final Amortisation Date must be on/after cap end date otherwise some Notionals undefined!");

        if (strikes.getLength() != (strikeStepUpDates).getLength())
            throw ModelException("Nb of strike step up dates must equal nb strikes");

        if (   strikeStepUpDates[ (strikeStepUpDates).getLength()-1 ] < end   )
            throw ModelException("Final strike date must be on/after cap end date otherwise some strikes undefined!");
        */

        if (notionals.getLength() != sched->nbCoupons)
            throw ModelException("Nb of notionals must equal nb coupons");

        if (strikes.getLength() != sched->nbCoupons)
            throw ModelException("Nb of strike  must equal nb coupons");
    }
    catch (exception& e) {
        throw ModelException(e, __FUNCTION__);
    }
}


QuasiMQ::ProductSP UnivarIRCpn::createProduct(QuasiMQ* model) const {
    try {
        return QuasiMQ::ProductSP(new UnivarIRCpnProdMQ(model, this));
    }
    catch (exception& e) {
        throw ModelException(e, __FUNCTION__);
    }
}


void UnivarIRCpn::load(CClassSP& clazz){
    REGISTER(UnivarIRCpn, clazz);
    clazz->setPublic(); // make externally visible
    clazz->setDescription("Q3 quasi-vanilla style instrument definition");
    SUPERCLASS(KComponent);
    EMPTY_SHELL_METHOD(defaultConstructor);

    FIELD(start,            "start date of instrument");FIELD_MAKE_OPTIONAL(start);
    FIELD(end,              "end date of instrument");  FIELD_MAKE_OPTIONAL(end);
    FIELD(payoffType,       "option payoff number");
    FIELD(payoffParamsLabel,"label");                   FIELD_MAKE_OPTIONAL(payoffParamsLabel);
    FIELD(settlementType,   "cash or physical");
    FIELD(stub,             "'F'ront or 'B'ack stub");  FIELD_MAKE_OPTIONAL(stub);
    FIELD(strikes,          "Strikes");                 FIELD_MAKE_OPTIONAL(strikes);
    FIELD(strikeStepUpDates,"Strike Step up dates");    FIELD_MAKE_OPTIONAL(strikeStepUpDates);
    FIELD(payoffParams,     "payoff params");           FIELD_MAKE_OPTIONAL(payoffParams);
    FIELD(arrears,          "coupon in arrears");       FIELD_MAKE_OPTIONAL(arrears);
    FIELD(paymentDelay,     "Delay (from reset)");      FIELD_MAKE_OPTIONAL(paymentDelay);
    FIELD(notionals,        "Amortisation Schedule");
    FIELD(amortisationDates,"Amortisation Dates");      FIELD_MAKE_OPTIONAL(amortisationDates);

    FIELD(couponFreq,       "freq of coupon payments"); FIELD_MAKE_OPTIONAL(couponFreq);
    FIELD(accrualDCC,       "daycount convention for coupon accrual");
    FIELD(expiryDCC,        "daycount convention for time to expiry");
    FIELD(resetBadDayConv,  "reset BDC");               FIELD_MAKE_OPTIONAL(resetBadDayConv);
    FIELD(accrualBadDayConv,"accrual BDC");             FIELD_MAKE_OPTIONAL(accrualBadDayConv);
    FIELD(paymentBadDayConv,"payment BDC");             FIELD_MAKE_OPTIONAL(paymentBadDayConv);
    FIELD(rate,             "observation rate");
    FIELD(sched,            "Coupon Schedule");
}

CClassConstSP const UnivarIRCpn::TYPE = CClass::registerClassLoadMethod(
    "UnivarIRCpn", typeid(UnivarIRCpn), load);

START_PUBLIC_ENUM_DEFINITION(UnivarIRCpn::PayoffType::Enum, "");
ENUM_VALUE_AND_NAME(UnivarIRCpn::PayoffType::CALL, "CALL", "");
ENUM_VALUE_AND_NAME(UnivarIRCpn::PayoffType::PUT, "PUT", "");
ENUM_VALUE_AND_NAME(UnivarIRCpn::PayoffType::BINARY_CALL, "BINARY_CALL", "");
ENUM_VALUE_AND_NAME(UnivarIRCpn::PayoffType::BINARY_PUT, "BINARY_PUT", "");
ENUM_VALUE_AND_NAME(UnivarIRCpn::PayoffType::INSIDE_RIB, "INSIDE_RIB", "");
ENUM_VALUE_AND_NAME(UnivarIRCpn::PayoffType::OUTSIDE_RIB, "OUTSIDE_RIB", "");
END_ENUM_DEFINITION(UnivarIRCpn::PayoffType::Enum);


/***************************************************************************************************************
*
* Define Univariate Interest Rate RIB Coupon object.
*
* This defines the payment details of an instrument with fixed RIB coupon payment (observe on
* floating rate, pay fixed rate).
*
* Example : Cap on CMS rate.
*
* Samy A Mohammed (22 August 2006)
**************************************************************************************************************

// type loading

void UnivarIRRibCpn::validatePop2Object() {
    try {
        if (start > end) {
            start = start;
            throw ModelException("Start must be before end");
        }

        if (notionals.getLength() != (amortisationDates).getLength())
            throw ModelException("Nb of notional step up dates must equal nb Notionals");

        if (   amortisationDates[ (amortisationDates).getLength()-1 ] < end   )
            throw ModelException("Final Amortisation Date must be on/after cap end date otherwise some Notionals undefined!");

        if (strikes.getLength() != (strikeStepUpDates).getLength())
            throw ModelException("Nb of strike step up dates must equal nb strikes");

        if (   strikeStepUpDates[ (strikeStepUpDates).getLength()-1 ] < end   )
            throw ModelException("Final strike date must be on/after cap end date otherwise some strikes undefined!");
    }
    catch (exception& e) {
        throw ModelException(e, __FUNCTION__);
    }
}


QuasiMQ::ProductSP UnivarIRRibCpn::createProduct(QuasiMQ* model) const {
    try {
        return QuasiMQ::ProductSP(new UnivarIRRibCpnProd2Q(model, this));
    }
    catch (exception& e) {
        throw ModelException(e, __FUNCTION__);
    }
}


void UnivarIRRibCpn::UnivarIRRibCpnProd2Q::price(
                           Control*   control,
                           CResults*  results) {
    try {

        //////////////////////////////////////////////////////
        //                                                  //
        // This data does not change inside the loop        //
        //                                                  //
        //////////////////////////////////////////////////////
        double                  dscPrice = 0;
        vector<long>            zeroDatesIdx;
        vector<double>          zeroRatesIdx;
        vector<long>            zeroDatesDisc;
        vector<double>          zeroRatesDisc;
        double                  output[50];
        int                     numVnfmParams;
        double                  vnfmParams[9];
        long                    rateFreq;
        string                  rateDCC;
        int                     numPayoffParams = inst->payoffParams.getLength();
        DayCountConventionSP    localExpiryDCC = inst->expiryDCC;
        DateTime                localValDt = model->getValueDate();
        DateTime                localToday = model->today;
        int                     j=0, i = 0, idx=0;

        //      if (cms)
        {
            rateFreq = inst->rate->frequency->approxAnnualFrequency();
        }
        //      else
        //      {
        //          rateFreq = 0;
        //      }

        rateDCC = inst->rate->dcc->toString();

        //Get zero curves
    //SAMY  IrConverter::q3GetZerosData(zeroDatesIdx, zeroRatesIdx, model->idxCurve.get());

    //SAMY  IrConverter::q3GetZerosData(zeroDatesDisc, zeroRatesDisc, model->discCurve.get());

        //Assume that this data comes from somewhere
        numVnfmParams = 9;
        vnfmParams[0] = 0.005;
        vnfmParams[1] = 0.6215;
        vnfmParams[2] = 1.5119;
        vnfmParams[3] = 0.0964;
        vnfmParams[4] = 0.45;
        vnfmParams[5] = 0.48;
        vnfmParams[6] = -0.0602;
        vnfmParams[7] = -0.1036;
        vnfmParams[8] =-0.95;


        //////////////////////////////////////////////////////
        //                                                  //
        // This data is coupons dependent                   //
        //                                                  //
        //////////////////////////////////////////////////////
        int nbCoupons = ( (inst->end).daysDiff(inst->start)/365 + 1 ) *
            (inst->couponFreq)->approxAnnualFrequency()  *
                (inst->obsFreq)->approxAnnualFrequency();
        //These are all output arrays (i.e. for each coupon)
        DoubleArraySP       premiumsSP(new CDoubleArray(nbCoupons));
        DoubleArraySP       fwdRatesSP(new CDoubleArray(nbCoupons));
        DoubleArraySP       scaledSigATMsSP(new CDoubleArray(nbCoupons));
        DoubleArraySP       couponNotionalsSP(new CDoubleArray(nbCoupons));
        DoubleArraySP       couponStrikesSP(new CDoubleArray(nbCoupons));
        DoubleArraySP       dscFactorsSP(new CDoubleArray(nbCoupons));// = (model->discCurve).pv(localValDt, localPaymentDate);
        IntArraySP          couponNbSP(new CIntArray(nbCoupons));
        DateTimeArraySP     expiryDatesSP(new DateTimeArray(nbCoupons));
        DateTimeArraySP     startDatesSP(new DateTimeArray(nbCoupons));
        DateTimeArraySP     matDatesSP(new DateTimeArray(nbCoupons));
        DateTimeArraySP     ribAccrualStartSP(new DateTimeArray(nbCoupons));
        DateTimeArraySP     ribAccrualEndSP(new DateTimeArray(nbCoupons));
        DateTimeArraySP     payDatesSP(new DateTimeArray(nbCoupons));
        DateTimeArraySP     accrualStartSP(new DateTimeArray(nbCoupons));
        DateTimeArraySP     accrualEndSP(new DateTimeArray(nbCoupons));
        double              smile[4];
        DateTime            localAccrualSt, localAccrualEnd, localPayDate;

        //Organize the caplet dates using some stubbing convention
        if (inst->stub == "B")
        {
            i = 0;
            nbCoupons = 0;
            do
            {
                localAccrualSt = i > 0 ? localAccrualEnd:  inst->start;
                localAccrualEnd = MIN(
                                        inst->couponFreq->toDate( localAccrualSt),
                                        inst->end
                                    );
                localPayDate =  ( localAccrualEnd ).rollDateInMonths(inst->paymentDelay->toMonths());;

                j = 0;
                do
                {
                    (*accrualStartSP)[nbCoupons]        = localAccrualSt;
                    (*accrualEndSP)[nbCoupons]          = localAccrualEnd;
                    (*ribAccrualStartSP)[nbCoupons]     = j > 0 ? (*ribAccrualEndSP)[nbCoupons-1] :  localAccrualSt;
                    (*ribAccrualEndSP)[nbCoupons]       = MIN(
                                                                inst->obsFreq->toDate( (*ribAccrualStartSP)[nbCoupons] ),
                                                                localAccrualEnd
                                                            );
                    if (inst->arrears)
                    {
                        (*expiryDatesSP)[nbCoupons]     = (*startDatesSP)[nbCoupons] = (*ribAccrualEndSP)[nbCoupons];
                    }
                    else
                    {
                        //advance reset
                        (*expiryDatesSP)[nbCoupons]     = (*startDatesSP)[nbCoupons] = (*ribAccrualStartSP)[nbCoupons];
                    }

                    (*matDatesSP)[nbCoupons]            = (inst->rate->tenor)->toDate( (*startDatesSP)[nbCoupons]    ); //inst->rate->tenor defines the tnor of the underlying rate e.g. 5Y in 5Y CMS
                    (*payDatesSP)[nbCoupons]            = localPayDate;
                    (*couponNbSP)[nbCoupons]            = i+1;

                    //Call pricer
                    //calculate notional (staircase interpolation)
                    if ( (*ribAccrualEndSP)[nbCoupons] <= (inst->amortisationDates)[idx])
                    {
                        (*couponNotionalsSP)[nbCoupons] = inst->notionals[idx];
                    }
                    else
                    {
                        while ( (*ribAccrualEndSP)[nbCoupons]  >  (inst->amortisationDates)[idx]  &&
                                (*ribAccrualEndSP)[nbCoupons]  >  (inst->amortisationDates)[idx+1] )
                        {
                            ++idx;
                        }
                        (*couponNotionalsSP)[nbCoupons] = inst->notionals[idx+1];
                    }

                    idx = 0;
                    //calculate strike (staircase interpolation)
                    if ( (*ribAccrualEndSP)[nbCoupons] <= (inst->strikeStepUpDates)[idx])
                    {
                        (*couponStrikesSP)[nbCoupons] = inst->strikes[idx];
                    }
                    else
                    {
                        while ( (*ribAccrualEndSP)[nbCoupons]  >  (inst->strikeStepUpDates)[idx]  &&
                                (*ribAccrualEndSP)[nbCoupons]  >  (inst->strikeStepUpDates)[idx+1] )
                        {
                            ++idx;
                        }
                        (*couponStrikesSP)[nbCoupons] = inst->strikes[idx+1];
                    }

                    //Assume that this data comes from somewhere
                    (*scaledSigATMsSP)[nbCoupons] = 0.12;

                    //Assume that this data comes from somewhere
                    smile[0] = 30.0;
                    smile[1] = model->qLeft;
                    smile[2] = model->qRight;
                    smile[3] = model->fwdShift;


                    Q3MQQuasiPricer(
                        localToday.getDate(),
                        localToday.getDate(),
                        localValDt.getDate(),
                        zeroDatesIdx.size(),
                        &(zeroDatesIdx[0]),
                        &(zeroRatesIdx[0]),
                        zeroDatesDisc.size(),
                        &(zeroDatesDisc[0]),
                        &(zeroRatesDisc[0]),
                        numVnfmParams,
                        vnfmParams,
                        rateFreq,
                        &(rateDCC[0]),
                        ((*expiryDatesSP)[nbCoupons]).getDate(),
                        ((*startDatesSP)[nbCoupons]).getDate(),
                        ((*matDatesSP)[nbCoupons]).getDate(),
                        (*scaledSigATMsSP)[nbCoupons],
                        0.,
                        smile,
                        inst->payoffType,
                        &( (*couponStrikesSP)[nbCoupons] ) ,
                        numPayoffParams,
                        (double*)&(inst->payoffParams[0]),
                        ((*payDatesSP)[nbCoupons]).getDate(),
                        1,//              setlType,
                        NULL,//   char              *holidayFile,
                        NULL,           //char              *BusVolDCC,
                        output) ;

                    if (nbCoupons< nbCoupons-1)
                    {
                        (*expiryDatesSP)[nbCoupons+1] = inst->couponFreq->toDate( (*expiryDatesSP)[nbCoupons] );
                    }
                    (*premiumsSP)[nbCoupons] = output[0];
                    (*fwdRatesSP)[nbCoupons] = output[1];
                    (*dscFactorsSP)[nbCoupons] = ( (model->domesticYC).get() )->pv(localValDt, (*payDatesSP)[nbCoupons]);
                    dscPrice += output[0] * (*dscFactorsSP)[nbCoupons];

                //
                    ++nbCoupons;
                    ++j;
                }
                while (   (*ribAccrualEndSP)[nbCoupons-1] < localAccrualEnd   );
                ++i;

            }
            while ( localAccrualEnd < inst->end  );
        }// if stub = F
        else
        {
            throw ModelException("Not currently supporting Front Stub");
        }


        //generate output
    //  results->storePrice(dscPrice, "USD");
        results->storeScalarGreek(dscPrice, Results::DEBUG_PACKET,
            OutputNameSP(new OutputName("DISCOUNDED_PRICE")));
        results->storeGreek(premiumsSP, Results::DEBUG_PACKET,
            OutputNameSP(new OutputName("UNDISCOUNTED_PREMIUM")));
        results->storeGreek(fwdRatesSP, Results::DEBUG_PACKET,
            OutputNameSP(new OutputName("UNADJUSTED_FWD_RATE")));
        results->storeGreek(dscFactorsSP, Results::DEBUG_PACKET,
            OutputNameSP(new OutputName("DISCOUNT_FACTOR")));
        results->storeGreek(expiryDatesSP, Results::DEBUG_PACKET,
            OutputNameSP(new OutputName("EXPIRY")));
        results->storeGreek(startDatesSP, Results::DEBUG_PACKET,
            OutputNameSP(new OutputName("RATE_START")));
        results->storeGreek(matDatesSP, Results::DEBUG_PACKET,
            OutputNameSP(new OutputName("RATE_MAT")));
        results->storeGreek(ribAccrualStartSP, Results::DEBUG_PACKET,
            OutputNameSP(new OutputName("RIB_ACC_START")));
        results->storeGreek(ribAccrualEndSP, Results::DEBUG_PACKET,
            OutputNameSP(new OutputName("RIB_ACC_END")));
        results->storeGreek(payDatesSP, Results::DEBUG_PACKET,
            OutputNameSP(new OutputName("PAYMENT_DATE")));
        results->storeGreek(accrualStartSP, Results::DEBUG_PACKET,
            OutputNameSP(new OutputName("ACC_START")));
        results->storeGreek(accrualEndSP, Results::DEBUG_PACKET,
            OutputNameSP(new OutputName("ACC_END")));
        results->storeGreek(couponNotionalsSP, Results::DEBUG_PACKET,
            OutputNameSP(new OutputName("NOTIONAL")));
        results->storeGreek(couponStrikesSP, Results::DEBUG_PACKET,
            OutputNameSP(new OutputName("STRIKES")));
        results->storeGreek(couponNbSP, Results::DEBUG_PACKET,
            OutputNameSP(new OutputName("COUPON_NBS")));
    }
    catch (exception& e) {
        throw ModelException(e, __FUNCTION__);
    }
}


void UnivarIRRibCpn::load(CClassSP& clazz){
    REGISTER(UnivarIRRibCpn, clazz);
    clazz->setPublic(); // make externally visible
    clazz->setDescription("Q3 quasi-vanilla fixed RIB style instrument definition");
    SUPERCLASS(KComponent);
    EMPTY_SHELL_METHOD(defaultConstructor);

    FIELD(start,             "start date of instrument");
    FIELD(end,               "end date of instrument");
    FIELD(payoffType,        PayoffTypeEnum::comment() + "option payoff number");
    FIELD(payoffParamsLabel, "label");
    FIELD(settlementType,    "cash or physical");
    FIELD(stub,              "'F'ront or 'B'ack stub");
    FIELD(strikes,           "Strikes");
    FIELD(strikeStepUpDates, "Strike Step up dates");
    FIELD(payoffParams,      "payoff params");
    FIELD(arrears,           "coupon in arrears");
    FIELD(paymentDelay,      "Delay (from reset)");
    FIELD(notionals,         "Amortisation Schedule");
    FIELD(amortisationDates,     "Amortisation Dates");

    FIELD(couponFreq,        "freq of coupon payments");
    FIELD(obsFreq,           "freq of idx obervation");
    FIELD(accrualDCC,        "daycount convention for coupon accrual");
    FIELD(expiryDCC,         "daycount convention for time to expiry");
    FIELD(resetBadDayConv,   "reset BDC");
    FIELD(accrualBadDayConv, "accrual BDC");
    FIELD(paymentBadDayConv, "payment BDC");
    FIELD(rate,              "observation rate");


    FIELD_MAKE_OPTIONAL(payoffParamsLabel);
    FIELD_MAKE_OPTIONAL(strikes);
    FIELD_MAKE_OPTIONAL(strikeStepUpDates);
    FIELD_MAKE_OPTIONAL(payoffParams);
}

CClassConstSP const UnivarIRRibCpn::TYPE = CClass::registerClassLoadMethod(
    "UnivarIRRibCpn", typeid(UnivarIRRibCpn), load);

*/


bool QuasiVanillaLoad(){
    return (UnivarIRCpn::TYPE != 0);
}


DRLIB_END_NAMESPACE
