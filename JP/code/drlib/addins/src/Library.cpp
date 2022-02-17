//----------------------------------------------------------------------------
//
//   Group       : EDG Derivatives Research
//
//   Filename    : Library.cpp
//
//   Description : Starts up & shuts down the DR library
//
//   Author      : Andrew J Swain
//
//   Date        : 29 January 2001
//
//
//----------------------------------------------------------------------------

#include "edginc/config.hpp"
#include "edginc/Library.hpp"
#include "edginc/DRAnalyticsInterface.h"
#include "edginc/SetNewHandler.hpp"
#include "edginc/ObjectIteration.hpp"
#include "edginc/Instrument.hpp"
#include "edginc/VolParallelShift.hpp"
#include "edginc/VolLevel.hpp"
#include "edginc/VolBenchmarkShift.hpp"
#include "edginc/VegaSkewPointwise.hpp"
#include "edginc/VegaSkewParallel.hpp"
#include "edginc/VolParallel.hpp"
#include "edginc/ThetaFwdSpot.hpp"
#include "edginc/SpotLevel.hpp"
#include "edginc/RootTimeVega.hpp"
#include "edginc/RhoBorrowParallel.hpp"
#include "edginc/RhoBorrowPointwise.hpp"
#include "edginc/Correl.hpp"
#include "edginc/MuSpecial.hpp"
#include "edginc/MuPointwise.hpp"
#include "edginc/MuParallel.hpp"
#include "edginc/FXVega.hpp"
#include "edginc/FXDelta.hpp"
#include "edginc/Delta.hpp"
#include "edginc/DeltaProxy.hpp"
#include "edginc/VegaMatrix.hpp"
#include "edginc/Asset.hpp"
#include "edginc/Equity.hpp"
#include "edginc/VolBase.hpp"
#include "edginc/FXAsset.hpp"
#include "edginc/JumpRateParallel.hpp"
#include "edginc/JumpRatePointwise.hpp"
#include "edginc/CEVPowerParallel.hpp"
#include "edginc/CEVPowerPointwise.hpp"
#include "edginc/XMLReader.hpp"
#include "edginc/SpreadSheetMode.hpp"
#include "edginc/Version.hpp"
#include "edginc/ErrorHandler.hpp"
#include "edginc/DeltaDDE.hpp"
#include "edginc/IRestorableWithRespectTo.hpp"
#include "edginc/Spot.hpp"
#include "edginc/VolPointwise.hpp"
#include "edginc/IDynamicsParameter.hpp"
#include "edginc/RateParallel.hpp"
#include "edginc/RatePointwise.hpp"
#include "edginc/DerivativeAsset.hpp"
#include "edginc/Smile2QElementwise.hpp"
#include "edginc/IRMeanReversionParallelProperty.hpp"
#include "edginc/IRSmile2QTweak.hpp"


#if defined(WIN32) && !defined(__CYGWIN__)
#include "winsock.h"
#endif

DRLIB_BEGIN_NAMESPACE
// declarations to avoid lots of tedious header files
extern void XLConvertRegister();
extern void XLConvertDoubleRegister();
extern void XLConvertIntRegister();
extern void XLConvertStringRegister();
extern void XLConvertEnumRegister();
extern void XLConvertBoolRegister();
extern void XLConvertDateTimeRegister();
extern void XLConvertDoubleMatrixRegister();
extern void XLConvertExpiryRegister();
extern void XLConvertDoubleArrayRegister();
extern void XLConvertIntArrayRegister();
extern void XLConvertBoolArrayRegister();
extern void XLConvertStringArrayRegister();
extern void XLConvertEnumArrayRegister();
extern void XLConvertObjArrayRegister();
extern void XLConvertDateTimeArrayRegister();
extern void XLConvertExpiryArrayRegister();
extern void XLConvertMarketWrapperRegister();
extern void XLConvertMarketWrapperArrayRegister();

static void registerXLConvertTypes(){
    // this has to be invoked after all the "TYPE"'s have been initialised
    XLConvertRegister();
    XLConvertDoubleRegister();
    XLConvertIntRegister();
    XLConvertStringRegister();
    XLConvertEnumRegister();
    XLConvertBoolRegister();
    XLConvertDateTimeRegister();
    XLConvertDoubleMatrixRegister();
    XLConvertDoubleArrayRegister();
    XLConvertIntArrayRegister();
    XLConvertBoolArrayRegister();
    XLConvertStringArrayRegister();
    XLConvertEnumArrayRegister();
    XLConvertObjArrayRegister();
    XLConvertExpiryRegister();
    XLConvertDateTimeArrayRegister();
    XLConvertExpiryArrayRegister();
    XLConvertMarketWrapperRegister();
    XLConvertMarketWrapperArrayRegister();
}

static int startupCount = 0;

// call before using any library functions
void Library::startup() {
    startupCount++;
    if (startupCount != 1){
        return;
    }
    // default behaviour of 'assert' is to print msg to stderr.
    // we want a dialog, to give us the chance to debug.
#if defined(WIN32) && defined(DEBUG) && !defined (DOING_OFF_BUILD)
    _set_error_mode(_OUT_TO_MSGBOX);
#endif

    // ensure that new throws an exception on failure
    CSetNewHandler::handler();
    // startup the XML library
    XMLReader::initialize();
    // startup NT socket API
#if defined(WIN32) && !defined(__CYGWIN__)
    WSADATA WsaData;
    if (WSAStartup(0x0101, &WsaData)) {
        throw ModelException("Library::startup",
                             "NT socket initialisation failed");
    }
#endif

    // XL stuff - needs to be before loadAllClasses()
    registerXLConvertTypes();
    // and then complete our run time type information
    CClass::loadAllClasses();
    
    // cache information about the classes
    // to do - make this more automatic
    CClassVec  tweaks;
    tweaks.push_back(ITweakableWithRespectTo<Spot>::TYPE);
    tweaks.push_back(IRestorableWithRespectTo<Spot>::TYPE);
    tweaks.push_back(DeltaProxy::Shift::TYPE);
    tweaks.push_back(DeltaDDE::Shift::TYPE);
    tweaks.push_back(DeltaDDE::RestorableShift::TYPE);
    tweaks.push_back(FXDelta::Shift::TYPE);
    tweaks.push_back(FXDelta::RestorableShift::TYPE);
    tweaks.push_back(ITweakableWithRespectTo<Correl>::TYPE);
    tweaks.push_back(IRestorableWithRespectTo<Correl>::TYPE);
    tweaks.push_back(FXVega::IShift::TYPE);
    tweaks.push_back(FXVega::IRestorableShift::TYPE);
    tweaks.push_back(MuParallel::Shift::TYPE);
    tweaks.push_back(MuParallel::RestorableShift::TYPE);
    tweaks.push_back(MuPointwise::IShift::TYPE);
    tweaks.push_back(MuPointwise::IRestorableShift::TYPE);
    tweaks.push_back(MuSpecial::IShift::TYPE);
    tweaks.push_back(MuSpecial::IRestorableShift::TYPE);
    tweaks.push_back(RhoBorrowParallel::Shift::TYPE);
    tweaks.push_back(RhoBorrowParallel::RestorableShift::TYPE);
    tweaks.push_back(RhoBorrowPointwise::IShift::TYPE);
    tweaks.push_back(RhoBorrowPointwise::IRestorableShift::TYPE);
    tweaks.push_back(ITweakableWithRespectTo<RateParallel>::TYPE);
    tweaks.push_back(IRestorableWithRespectTo<RateParallel>::TYPE);
    tweaks.push_back(ITweakableWithRespectTo<RatePointwise>::TYPE);
    tweaks.push_back(IRestorableWithRespectTo<RatePointwise>::TYPE);
    tweaks.push_back(RootTimeVega::IShift::TYPE);
    tweaks.push_back(RootTimeVega::IRestorableShift::TYPE);
    tweaks.push_back(SpotLevel::Shift::TYPE);
    tweaks.push_back(Theta::Shift::TYPE);
    tweaks.push_back(Theta::RestorableShift::TYPE);
    tweaks.push_back(ThetaFwdSpot::Shift::TYPE);
    tweaks.push_back(ThetaFwdSpot::RestorableShift::TYPE);
    tweaks.push_back(VegaMatrix::IShift::TYPE);
    tweaks.push_back(VegaMatrix::IRestorableShift::TYPE);
    tweaks.push_back(ITweakableWithRespectTo<VolParallel>::TYPE);
    tweaks.push_back(IRestorableWithRespectTo<VolParallel>::TYPE);
    tweaks.push_back(ITweakableWithRespectTo<VolPointwise>::TYPE);
    tweaks.push_back(IRestorableWithRespectTo<VolPointwise>::TYPE);
    tweaks.push_back(VegaSkewParallel::IShift::TYPE);
    tweaks.push_back(VegaSkewParallel::IRestorableShift::TYPE);
    tweaks.push_back(VegaSkewPointwise::IShift::TYPE);
    tweaks.push_back(VegaSkewPointwise::IRestorableShift::TYPE);  
    tweaks.push_back(VolBenchmarkShift::Shift::TYPE);
    tweaks.push_back(VolLevel::Shift::TYPE);
    tweaks.push_back(VolParallelShift::Shift::TYPE);
    tweaks.push_back(JumpRateParallel::Shift::TYPE);
    tweaks.push_back(JumpRateParallel::RestorableShift::TYPE);
    tweaks.push_back(JumpRatePointwise::IShift::TYPE);
    tweaks.push_back(JumpRatePointwise::IRestorableShift::TYPE);
    tweaks.push_back(CEVPowerParallel::Shift::TYPE);
    tweaks.push_back(CEVPowerParallel::RestorableShift::TYPE);
    tweaks.push_back(CEVPowerPointwise::IShift::TYPE);
    tweaks.push_back(CEVPowerPointwise::IRestorableShift::TYPE);
    tweaks.push_back(IDynamicsParameter::TYPE);
    tweaks.push_back(DerivativeAsset::TYPE);
    tweaks.push_back(ITweakableWithRespectTo<IRMeanReversionParallelProperty>::TYPE);
    tweaks.push_back(ITweakableWithRespectTo<Smile2QElementwise>::TYPE);
    tweaks.push_back(ITweakableWithRespectTo<IRSmile2QTweak::Property>::TYPE);
                     
// 'tweaks' as viewed by the collectors ie they look, for example,
    // for objects derived from an asset
    tweaks.push_back(CAsset::TYPE);
    tweaks.push_back(Equity::TYPE);
    tweaks.push_back(CVolBase::TYPE);
    tweaks.push_back(FXAsset::TYPE);
    tweaks.push_back(DividendList::TYPE);
    // finally ask to immediately cache information for all
    // instruments, all assets, all vols, and all yield curves
    ObjectIteration::cacheClassInformation(Instrument::TYPE, tweaks);
    ObjectIteration::cacheClassInformation(CAsset::TYPE, tweaks);
    ObjectIteration::cacheClassInformation(Equity::TYPE, tweaks);
    ObjectIteration::cacheClassInformation(CVolBase::TYPE, tweaks);
    ObjectIteration::cacheClassInformation(YieldCurve::TYPE, tweaks);
}


// call before program exit
void Library::shutdown() {
    startupCount--;
    if (startupCount == 0){
        // call the XML lib termination method
      // To be reviewed:        XMLReader::terminate();
#if defined(WIN32) && !defined(__CYGWIN__)
        // tidy up sockets stuff on NT 
        WSACleanup();
#endif
    }
}
//// The prefix to use for Excel functions. Required by 
//// XLAddin and xlregister. Returns "COGS_"
const char* Library::XL_PREFIX ="QLIB_";

static string ERROR_LOG = "c:/edrerror.log";

/** This is invoked at 'start up' immediately after SpreadSheetMode is
    configured */
void Library::configureErrorLog(){
    if (SpreadSheetMode::isOn()){
        // write errors to output file (clear file at start)
        ErrorHandler::logFile(ERROR_LOG, true);
        ErrorHandler::writeMsg("Quantitative Research Library, Version " +
                               CVersion::DRLibVersion() + "\n");
    }
}

const string Library::SERVICE_NAME = "QLIB";
const string Library::THREADSAFE_SERVICE_NAME = "TS_QLIB";
// service version is always the same for both services
const int Library::SERVICE_VERSION = DRLIB_VERSION(1,0,0,0);

DRLIB_END_NAMESPACE

