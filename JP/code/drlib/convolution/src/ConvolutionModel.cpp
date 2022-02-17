//----------------------------------------------------------------------------
//
//   Group       : QR
//
//   Description : View onto instrument as required by ConvolutionEngine
//
//   Date        : 18th Nov 2005
//
//   Author      : Mark Robson
//
//----------------------------------------------------------------------------

#include "edginc/config.hpp"
#include "edginc/Class.hpp"
#include "edginc/Atomic.hpp"
#include "edginc/CreditAsset.hpp"
#include "edginc/CounterPartyCredit.hpp"
#include "edginc/SingleCreditAsset.hpp"
#include "edginc/CreditTrancheLossConfig.hpp"
#include "edginc/FlatCDO2LossConfig.hpp"
#include "edginc/ConvolutionModel.hpp"


DRLIB_BEGIN_NAMESPACE

IConvolutionModel::IConvolutionModel(){}

IConvolutionModel::~IConvolutionModel(){}

void IConvolutionModel::load(CClassSP& clazz){
    clazz->setPublic(); // make visible to EAS/spreadsheet
    REGISTER_INTERFACE(IConvolutionModel, clazz);
    EXTENDS(IObject);
}

CClassConstSP const IConvolutionModel::TYPE = 
CClass::registerInterfaceLoadMethod(
    "IConvolutionModel", typeid(IConvolutionModel), load);



void IConvolutionModel::createLossCalculatorsFIX(
        const DateTimeArray&                timeline,           /* (I) */
        CounterPartyCreditConstSP           cpty,               /* (I) */
        const DateTime&                     maturity,           /* (I) */
        Control*                            control,            /* (I) */
        Results*						     results,            /* (I) */
        bool                                recoverNotional,    /* (I) */
        YieldCurveConstSP                   discount,           /* (I) */
        IFixedTrancheLossCalculatorConstSP& lossCalculator,                   /* (O) */
        IFixedTrancheLossCalculatorConstSP& recoveredNotionalCalculator,      /* (O) */
        IFixedTrancheLossCalculatorConstSP& conditionalLossCalc,              /* (O) */
		IFixedTrancheLossCalculatorConstSP& conditionalRecNtnlCalc,
		FlatCDO2LossConfigConstSP      ts            /* (I) */) const
{
	throw ModelException("Not to be used");	
};

DRLIB_END_NAMESPACE

