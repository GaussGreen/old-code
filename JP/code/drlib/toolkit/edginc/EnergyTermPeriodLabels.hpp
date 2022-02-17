//----------------------------------------------------------------------------
//
//   Group       : GCCG Derivatives Research
//
//   Filename    : EnergyTermPeriodLabels.hpp
//
//   Description : Defines a deal period for energy. 
//
//   Author      : Sean Chen
//
//   Date        : 26 Sept. 2005
//
//----------------------------------------------------------------------------

#ifndef _EnergyTermPeriodLabels_H_
#define _EnergyTermPeriodLabels_H_

#include "edginc/config.hpp"
#include "edginc/Class.hpp"
#include "edginc/EnergyTermPeriod.hpp"

#include <string>
using namespace std;

DRLIB_BEGIN_NAMESPACE


class TOOLKIT_DLL EnergyTermPeriodLabels : public EnergyTermPeriod
{

public:

	static CClassConstSP const TYPE;

    friend class EnergyTermPeriodLabelsHelper;

    //EnergyTermPeriodLabels(const EnergyTermPeriodLabels& thePeriod);

    void validatePop2Object();

    virtual ~EnergyTermPeriodLabels();

protected:

    EnergyTermPeriodLabels();

	string startLabel;
	string endLabel;

};

typedef smartPtr<EnergyTermPeriodLabels> EnergyTermPeriodLabelsSP;
typedef smartConstPtr<EnergyTermPeriodLabels> EnergyTermPeriodLabelsConstSP;

DRLIB_END_NAMESPACE

#endif
