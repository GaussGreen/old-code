//----------------------------------------------------------------------------
//
//   Group       : EDG Derivatives Research
//
//   Filename    : MonteCarloLib.cpp
//
//   Description : Ensures all classes get linked in (important if there are
//                 no other dependencies upon a class)
//
//   Date        : June 2001
//
//
//----------------------------------------------------------------------------
#include "edginc/config.hpp"
#include "edginc/MonteCarloLib.hpp"
#include "edginc/ModelException.hpp"
#include "edginc/MonteCarlo.hpp" // FILTER
#include "edginc/MCPathConfig.hpp" // FILTER
#include "edginc/ImpliedSample.hpp" // FILTER
#include "edginc/HitSample.hpp" // FILTER

DRLIB_BEGIN_NAMESPACE
extern bool MCPathConfigLNLoad();      // avoid need for header file
extern bool MCPathConfigImpliedLoad(); // avoid need for header file
extern bool MonteCarloTurboLoad();
extern bool MCPathConfigLVLoad();
extern bool MCPathConfigMertonLoad();
extern bool MertonDependenceLoad();
extern bool MCPathConfigLVMertonLoad();
extern bool MCPathConfigSVLoad();
extern bool MCPathConfigSVJLoad();
extern bool MCPathConfigSVCJLoad();
extern bool MCPathConfigDemoLoad();
extern bool MCPathConfigSRMLoad();
extern bool MCPathConfigCCMLoad();
extern bool MCPathConfigE3FLoad();
extern bool SRMIRIVAddinLoad();
extern bool SimpathILoad();
extern bool AdaptiveSamplingLoad();
extern bool IMFSamplingLoad();
extern bool ArbSamplingLoad();
extern bool GaussSamplingLoad();
extern bool ConditionalLossModelLoad();
extern bool ITimelineSpecload();

void MonteCarloLib::linkInClasses(){

    // Please don't add any new ...::TYPE entries here ---
    // use FFFLoad() where FFF.cpp is the source filename.

    bool success =
         MonteCarlo::TYPE != 0 &&
         MCPathConfig::TYPE != 0 &&
         LinearImpliedSample::TYPE != 0 &&
         DegenerateImpliedSample::TYPE != 0 &&
         LinearImpliedSampler::TYPE != 0 &&
         MonteCarloTurboLoad() &&
         MCPathConfigLNLoad() &&
         MCPathConfigImpliedLoad() &&
        MCPathConfigLVLoad() &&
        MCPathConfigMertonLoad() &&
        MertonDependenceLoad() &&
        MCPathConfigLVMertonLoad() &&
        MCPathConfigSVLoad() &&
        MCPathConfigSVJLoad() &&
        MCPathConfigSVCJLoad() &&
        MCPathConfigDemoLoad() &&
		MCPathConfigE3FLoad() &&
		MCPathConfigCCMLoad() &&
        HitSample::TYPE != 0 &&
        MCPathConfigSRMLoad() &&
        SRMIRIVAddinLoad() &&
        SimpathILoad() &&
		SimpathILoad() &&
		AdaptiveSamplingLoad() &&
		GaussSamplingLoad() &&
		IMFSamplingLoad() &&
		ArbSamplingLoad() &&
		ConditionalLossModelLoad() &&
		ITimelineSpecload() &&
        true;

    if (!success){
        throw ModelException("MonteCarloLib::registerClasses",
                             "Registration error");
    }
}

DRLIB_END_NAMESPACE


