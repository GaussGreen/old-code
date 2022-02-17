//----------------------------------------------------------------------------
//
//   Group       : EDG Derivatives Research
//
//   Filename    : MCPathConfigSRM.hpp
//
//   Description : A generator of paths using stochastic rates
//                 SRM = stochastic rate model
//
//   Date        : 27 May 2004
//
//
//----------------------------------------------------------------------------

#ifndef ISUPPORTPATHCONFIGSRM_HPP
#define ISUPPORTPATHCONFIGSRM_HPP

#include "edginc/FORWARD_DECLARE.hpp"
#include "edginc/DateTime.hpp"

DRLIB_BEGIN_NAMESPACE

class MCPathConfigSRM;
class IPastValues;
FORWARD_DECLARE(Dependence);

class ISupportPathConfigSRM {
public:
    virtual DateTimeArray getSRMCriticalDates(
        const MCPathConfigSRM* mcPathConfig,
        const DateTime& start,       // likely to be "today"
        const DateTime& finish) = 0; // the latest requested date

    virtual void createSRMUtil(
        MCPathConfigSRM*          mcPathConfig,
        const DateTime&        today,
        DateTimeArrayConstSP   simDates) = 0;

    virtual void setSRMDiffusibleAsset(
        const DateTime&     today,              // base date of the run
        const MCPathConfigSRM* mcPathConfig,       // to get max/min boundaries, etc
        const IPastValues*  pastValues,         // historic values
        DependenceSP        dependence) = 0;    // factorized correlations

};
DRLIB_END_NAMESPACE
#endif

