//----------------------------------------------------------------------------
//
//   Group       : EDG Derivatives Research
//
//   Filename    : AddinLib.cpp
//
//   Description : Ensures all classes get linked in (important if there are
//                 no other dependencies upon a class)
//
//   Author      : Mark A Robson
//
//   Date        : 28 Feb 2001
//
//
//----------------------------------------------------------------------------

#include "edginc/config.hpp"
#include "edginc/AddinLib.hpp"
#include "edginc/RiskMgrLib.hpp"
#include "edginc/MarketLib.hpp"
#include "edginc/ProductsLib.hpp"
#include "edginc/IRProductsLib.hpp"
#include "edginc/RadarLib.hpp"
#include "edginc/ToolkitLib.hpp"
#include "edginc/TreeLib.hpp"
#include "edginc/XLTest.hpp"
#include "edginc/DRWrapper.hpp"
#include "edginc/DRWrapperInterface.hpp"
#include "edginc/CompositeInstrument.hpp"  
#include "edginc/UtilLib.hpp"  
#include "edginc/EDRTester.hpp"  
#include "edginc/MonteCarloLib.hpp"
#include "edginc/BadBoy.hpp"
#include "edginc/CreditLib.hpp"
#include "edginc/ConvolutionLib.hpp"
#include "edginc/FourierLib.hpp"
#include "edginc/EventsManager.hpp"
#include "edginc/ToolkitDebug.hpp"

DRLIB_BEGIN_NAMESPACE

extern bool XLHelpExists;
extern bool XLTestsRegistered;
extern bool XLObjectViewRegistered;
extern bool XLObjViewAddinsRegistered;
extern bool XLvbRegistrationCodeRegistered;
extern bool XLGetOutputPacketAddinsRegistered;
extern bool UtilsLink;
extern bool AggregateInstrumentLinked;
extern bool PricingCounterLoad();
extern bool DRInterfaceLink;
extern bool EDGServicesLink;
extern bool CalibratedRiskMappingMatrixLoad();
extern bool CalibratorAddinLoad();
extern bool SRMInitialGuessLoad();
extern bool SRMEquityCalibLoad();
extern bool XMLReadWriteLoad();
extern bool InternalInterfaceLoad();
extern bool SetRiskyGrowthLoad();

bool AddinLib::linkInClasses(){
    void* funcPtr = (void*) XLTest::genericXL;
    return (funcPtr != 0 && XLHelpExists &&
            XLTestsRegistered &&
            XLObjectViewRegistered &&
            UtilsLink &&
            DRInterfaceLink &&
            XLObjViewAddinsRegistered &&
            XLvbRegistrationCodeRegistered &&
            AggregateInstrumentLinked  &&
            XLGetOutputPacketAddinsRegistered  &&
            DRWrapper::load()            &&
            PricingCounterLoad()         &&
            DRWrapperInterface::TYPE     &&
            CompositeInstrument::TYPE    &&
            EventsManager::TYPE          &&
            BadBoy::TYPE                 &&
            EDRTester::TYPE              &&
            EDGServicesLink              &&
            XMLReadWriteLoad()           &&
            CalibratedRiskMappingMatrixLoad() &&
			SRMInitialGuessLoad() &&
			SRMEquityCalibLoad() &&
			CalibratorAddinLoad() &&
			InternalInterfaceLoad() &&
			SetRiskyGrowthLoad() &&
            true);
}

void AddinLib::linkInAllLibs(){
    // ordering not important
    UtilLib::linkInClasses();
    CToolkitLib::linkInClasses();
    CRiskMgrLib::linkInClasses();
    CMarketLib::linkInClasses();
    FourierLib::linkInClasses();
    CProductsLib::linkInClasses();
    IRProductsLib::linkInClasses();
    CTreeLib::linkInClasses();
    MonteCarloLib::linkInClasses();
    CCreditLib::linkInClasses();
    ConvolutionLib::linkInClasses();
	RadarLib::linkInClasses();
    linkInClasses();
}

DRLIB_END_NAMESPACE
