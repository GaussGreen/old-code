//----------------------------------------------------------------------------
//
//   Group       : GCCG Derivatives Research
//
//   Filename    : EnergyTERMPeriod.hpp
//
//   Description : Defines a deal period for energy. 
//
//   Author      : Sean Chen
//
//   Date        : 26 Sept. 2005
//
//----------------------------------------------------------------------------

#ifndef _ENERGYTERMPERIOD_H_
#define _ENERGYTERMPERIOD_H_

#include "edginc/config.hpp"
#include "edginc/Class.hpp"
#include "edginc/DateTime.hpp"
#include "edginc/Object.hpp"
#include "edginc/EnergyContractLabel.hpp"

#include <string>
using namespace std;

DRLIB_BEGIN_NAMESPACE


class TOOLKIT_DLL EnergyTermPeriod : public CObject
{

public:

	static CClassConstSP const TYPE;

    friend class EnergyTermPeriodHelper;

    EnergyTermPeriod(const EnergyTermPeriod& thePeriod);

    void validatePop2Object();

    virtual ~EnergyTermPeriod();

    DateTime getStartDate() const;
    DateTime getEndDate() const;

	EnergyContractLabel getStartContractLabel() const;

    EnergyContractLabel getEndContractLabel() const;

	bool isPrompt() const;

protected:

    EnergyTermPeriod(CClassConstSP clazz);
	EnergyTermPeriod() : CObject(TYPE) {}

    DateTime startDate;
    DateTime endDate;
	EnergyContractLabel startContract;
	EnergyContractLabel endContract;
	bool prompt; // for PROMPT $unregistered

};

typedef smartPtr<EnergyTermPeriod> EnergyTermPeriodSP;
typedef smartConstPtr<EnergyTermPeriod> EnergyTermPeriodConstSP;

DRLIB_END_NAMESPACE

#endif
