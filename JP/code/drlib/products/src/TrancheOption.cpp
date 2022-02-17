//----------------------------------------------------------------------------
//
//   Group       : Credit QR
//
//   Filename    : TrancheOption.cpp
//
//   Description : a pilot instrument for SCID model 
//
//   Author      : Adrian Bozdog
//
//----------------------------------------------------------------------------


#include "edginc/config.hpp"
#include "edginc/TrancheOption.hpp"
#include "edginc/Maths.hpp"
#include "edginc/Black.hpp"
#include "edginc/SwapTool.hpp"
#include "edginc/EffectiveCurve.hpp"
#include "edginc/SCIDconvolution.hpp"
#include "edginc/SCIDparameters.hpp"
#include "edginc/Maths.hpp"
#include <fstream>

DRLIB_BEGIN_NAMESPACE

TrancheOption::~TrancheOption(){}

TrancheOption::TrancheOption(CClassConstSP clazz) : CDOIndexOption(clazz){
}

/** Do some asset specific validation */
void TrancheOption::Validate() {
    static const string routine("TrancheOption::Validate");
	try 
    {
        QLIB_VERIFY(lowStrike<highStrike,routine + ": lowStrike is not less than highStrike"); 
        CDOIndexOption::Validate();
	}
	catch (exception& e)
	{
		throw ModelException(e, routine);
	}

}

void TrancheOption::setFastMCParameters(SCIDparametersSP &sCIDparamSP) {
    DoubleArray lowStrikes(1, lowStrike), highStrikes(1, highStrike);
    sCIDparamSP->setConvolution(lowStrikes, highStrikes);
}

//returns the loss at expiry
double TrancheOption::getLossAndLegs(long expiryIndexInSim, DateTimeArray & simDates, DateTimeArray & maturities, DateTime& expiry, DoubleArray& RA, DoubleArray& DL, SCID * model) 
{
	SCIDparametersSP sCIDparamSP = model->getSCIDparameters();
    double timeSteps = model->getCFtimeSteps();
    DoubleArray etlTimes(simDates.size());
    for (int i=0; i<etlTimes.size(); i++) etlTimes[i] = valueDate.yearFrac(simDates[i]);
    sCIDparamSP->setFastMC(model->getSeed(), timeSteps, int(etlTimes.back()/timeSteps), model->getNbPathsFastNoJump(), model->getNbPathsFastAtLeastOneJump());
    DoubleMatrix ETL(etlTimes.size(),1); 
	ETL.fill(0);
	double loss = sCIDparamSP->getFutureETL(expiryIndexInSim, etlTimes, 0, ETL, true, model->getConvolutionNoJump(), model->getConvolutionAtLeastOneJump());    
	double div = 1.0 / ( highStrike - lowStrike );
    
    DateTimeArray curveDates(1, valueDate);
    curveDates.insert(curveDates.end(), simDates.begin(), simDates.end());
	DoubleArray riskyDiscount(1,1);
   
    for(int i=0; i < etlTimes.size(); ++i) 
        riskyDiscount.push_back(Maths::max(1e-12, 1-ETL[i][0]*div));           
    EffectiveCurve trancheCurve(valueDate, discount.getSP(), curveDates, riskyDiscount, EffectiveCurve::FLAT_FORWARD);
    long matSize = maturities.size();
	double DLT, DLt, RAT, RAt;
    for (long i=0; i < matSize; ++i){		
        /*DL[i] = trancheCurve.protectionPV(expiry, expiry ,maturities[i],IDiscountCurveRisky::RECOVER_1);
        CashFlowArray cashflows = SwapTool::cashflows(expiry, maturities[i], false, 1.0, coupon, "M", &(*dcc));
        cashflows[cashflows.size()-1].amount -= 1.0; // ugly
        RA[i] = trancheCurve.annuityPV(cashflows,expiry,IDiscountCurveRisky::RECOVER_1);
        */ 
        DLT = trancheCurve.protectionPV(valueDate,valueDate, maturities[i],IDiscountCurveRisky::RECOVER_1);
        DLt = trancheCurve.protectionPV(valueDate,valueDate,expiry,IDiscountCurveRisky::RECOVER_1);
        DL[i] = DLT - DLt; 
        CashFlowArray cashflows = SwapTool::cashflows(valueDate, maturities[i], false, 1.0, coupon, "M", &(*dcc));
        cashflows[cashflows.size()-1].amount -= 1.0; // ugly
        RAT = trancheCurve.annuityPV(cashflows,valueDate,IDiscountCurveRisky::RECOVER_1);
        cashflows = SwapTool::cashflows(valueDate, expiry, false, 1.0, coupon, "M", &(*dcc));
        cashflows[cashflows.size()-1].amount -= 1.0; // ugly
        RAt = trancheCurve.annuityPV(cashflows,valueDate,IDiscountCurveRisky::RECOVER_1);
        RA[i] = RAT - RAt;        
    }
        
    return PayoffTranche(loss, lowStrike, highStrike)/(highStrike - lowStrike);
}

/*=============================================================================
 * Reflection, loading, etc.
 *===========================================================================*/
class TrancheOptionHelper {
public:
    static IObject* defaultTrancheOption();
    /** Invoked when Class is 'loaded' */
    static void load(CClassSP& clazz){
        clazz->setPublic(); // make visible to EAS/spreadsheet
        REGISTER(TrancheOption, clazz);
        SUPERCLASS(CDOIndexOption);
        IMPLEMENTS(SCID::IIntoProduct);
        EMPTY_SHELL_METHOD(defaultTrancheOption);
        FIELD(lowStrike,       "Low strike (absolute value, NOT a percentage)");
        FIELD(highStrike,      "High strike (absolute value, NOT a percentage)");
	}
};

IObject* TrancheOptionHelper::defaultTrancheOption() {
    return new TrancheOption();
}

CClassConstSP const TrancheOption::TYPE = 
    CClass::registerClassLoadMethod("TrancheOption", typeid(TrancheOption),TrancheOptionHelper::load);

/** Included in ProductsLib to force the linker to include this file */
bool TrancheOptionLoad() {
    return (TrancheOption::TYPE != 0);
}

DRLIB_END_NAMESPACE

