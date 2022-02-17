#include "edginc/config.hpp"
#include "edginc/SVGenKFloat.hpp"
#include "edginc/Maths.hpp"
#include <fstream>

DRLIB_BEGIN_NAMESPACE

SVIProdCreatorSP SVGenKFloat::getSVIProdCreator(IStateVariableGen::IStateGen* pathGen) const {
    SVDiscFactorSP dfReset(dfResetGen->getSVDiscFactor(pathGen));
    SVDiscFactorSP dfPay(dfPayGen->getSVDiscFactor(pathGen));
    //SVIndexSpecIRSP irSV = SVIndexSpecIRSP::dynamicCast(irSVGen->getSVIProdCreator(pathGen));
	SVIndexSpecIRSP irSV = DYNAMIC_POINTER_CAST<SVIndexSpecIR>(irSVGen->getSVIProdCreator(pathGen));

    return SVIProdCreatorSP(new SVKFloat(sched, notionals, weights, spreads, dcfs, rateType, principalDates, principalPayments, irSV, dfReset, dfPay));
}

void SVGenKFloat::collectStateVars(IStateVariableCollectorSP svCollector) const {
    svCollector->append(irSVGen.get());
    svCollector->append(dfResetGen.get());
    svCollector->append(dfPayGen.get());
}

/*This function takes an index of a coupon reset date and returns the future cashflows (coupon
payments and principal payments) discounted to that reset date. */
double SVKFloat::elem(size_t index) {
    size_t nbCoupons = sched->pay.size();
    double value=0;
    double dfPayVal, dfResetVal;
    for (size_t i=index; i<nbCoupons; i++) {
        /*FIXME: Cannot handle two coupons with the same reset/payment dates because
        dfReset/dfPay would have array size less than the no. of coupons in that case */
        dfPayVal = dfPay->getDF(i);
        dfResetVal = dfReset->getDF(i);
        double discFactor = dfPayVal/dfResetVal;
        double weight = (*weights)[i];          
        double spread = (*spreads)[i];
        double dcf = (*dcfs)[i];
        double notional = (*notionals)[i];

        double annuity = 0;
        //for (size_t j=i+1; )

        if (!Maths::isZero(weight)) {            
             if (rateType==RateType::SIMPLE) {
             //if (isRateSimple) {
                 value += (notional * dcf * weight) * (irSV->elem(i) + spread/weight) * dfPayVal;
             } else {
                 value += notional * (::pow (1. + (weight * irSV->elem(i) + spread), dcf) - 1.) * dfPayVal;
             }
        } else {
            if (rateType==RateType::SIMPLE) {
               value += (notional * dcf * spread) * dfPayVal;
            } else {
               value += notional * (::pow (1. + spread, dcf) - 1.) * dfPayVal;
            }
        }
    }
    /*ofstream out;
    out.open("C:\\temp\\MC_MTM_"+index, ios_base::app);
    out << value << endl;
    out.close();
    out.open("C:\\temp\\MC_index_"+index, ios_base::app);
    out << irSV->elem(index) << endl;
    out.close();
    out.open("C:\\temp\\MC_dfPay_"+index, ios_base::app);
    out << dfPayVal << endl;
    out.close();
    out.open("C:\\temp\\MC_dfReset_"+index, ios_base::app);
    out << dfResetVal << endl;
    out.close();*/

    return value;
}

bool SVKFloat::doingPast() const {return dfPay->doingPast();}

DRLIB_END_NAMESPACE
