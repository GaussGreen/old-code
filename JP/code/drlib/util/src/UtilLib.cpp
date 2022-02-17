//----------------------------------------------------------------------------
//
//   Group       : EDG Derivatives Research
//
//   Filename    : UtilLib.cpp
//
//   Description : ensure all classes are loaded
//
//   Author      : Andrew J Swain
//
//   Date        : 14 June 2001
//
//
//----------------------------------------------------------------------------

#include "edginc/config.hpp"
#include "edginc/UtilLib.hpp"
#include "edginc/MuDates.hpp"
#include "edginc/MultiQDistribution.hpp"
#include "edginc/Spline.hpp"
#include "edginc/Spline2D.hpp"
#include "edginc/OLS.hpp"
#include "edginc/TimePoint2D.hpp"
#include "edginc/ToolkitDebug.hpp"
#include "edginc/Surface.hpp"

DRLIB_BEGIN_NAMESPACE

extern bool IntegratorLoad();
extern bool InterpolatorLoad();
extern bool LinearInterpolatorLoad();
extern bool LinearInterpolatorMDLoad();
extern bool SplineLoad();
extern bool RandomLinkIn();
extern bool PadeCoeffsAddinLinkIn();
extern bool PadeAddinLinkIn();
extern bool MultiQDistributionLinkIn();
extern bool FAMultiQDistributionLinkIn();
extern bool IDistribution1DLoad();
extern bool FunctionIntegrator1DLoad();
extern bool FunctionIntegratorNDLoad();
extern bool VNFMCMSLinkIn();
extern bool IIntegratorLoad();
extern bool DiscreteDistributionLoad();
extern bool RectangleIntegrator1DLoad();
extern bool IntFunc2DGeneralLoad();
extern bool QSpellucciLoad();
extern bool IConvolutorLoad();
extern bool IntegratorNDLoad();
extern bool VNFMCMCDSLinkIn();
extern bool SurfaceLoad();

void UtilLib::linkInClasses() {

    // Please don't add any new ...::TYPE entries here ---
    // use FFFLoad() where FFF.cpp is the source filename.

    bool success = MuDates::load() && 
                   IntegratorLoad() &&
                   InterpolatorLoad() &&
                   LinearInterpolatorLoad() &&
                   LinearInterpolatorMDLoad() &&
                   SplineLoad() &&
                   Spline2D::TYPE && 
                   OLS::TYPE &&
                   MultiQDistribution::TYPE &&
                   RandomLinkIn() &&
                   PadeCoeffsAddinLinkIn() &&
                   PadeAddinLinkIn() &&
                   MultiQDistributionLinkIn() &&
                   FAMultiQDistributionLinkIn() &&
                   IDistribution1DLoad() &&
                   VNFMCMSLinkIn() &&
                   FunctionIntegrator1DLoad() &&
                   FunctionIntegratorNDLoad() &&
                   TimePoint2D::TYPE &&
                   IIntegratorLoad() &&
                   DiscreteDistributionLoad() &&
                   RectangleIntegrator1DLoad() &&
                   IntFunc2DGeneralLoad() &&
		           QSpellucciLoad() &&
		           IConvolutorLoad() &&
				   IntegratorNDLoad() &&
				   VNFMCMCDSLinkIn() &&
				   SurfaceLoad() &&
                   true;

    if (!success){
        throw ModelException("UtilLib::registerClasses",
                             "Registration error");
    }
}

DRLIB_END_NAMESPACE
