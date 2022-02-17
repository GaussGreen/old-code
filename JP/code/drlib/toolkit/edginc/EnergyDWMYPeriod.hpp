//----------------------------------------------------------------------------
//
//   Group       : GCCG Derivatives Research
//
//   Filename    : EnergyDWMYPeriod.hpp
//
//   Description : Defines and decompose nD(WMY) periods used in energy. 
//
//   Author      : Sean Chen
//
//   Date        : 26 Sept. 2005
//
//----------------------------------------------------------------------------

#ifndef _ENERGYDWMYPERIOD_H_
#define _ENERGYDWMYPERIOD_H_

#include "edginc/config.hpp"
#include "edginc/Class.hpp"
#include "edginc/Object.hpp"

#include <string>
using namespace std;

DRLIB_BEGIN_NAMESPACE


class TOOLKIT_DLL EnergyDWMYPeriod : public CObject
{
public:

	static CClassConstSP const TYPE;

    friend class EnergyDWMYPeriodHelper;

    EnergyDWMYPeriod(const string& period);
    virtual ~EnergyDWMYPeriod();
    EnergyDWMYPeriod(const EnergyDWMYPeriod& thePeriod);
	EnergyDWMYPeriod();
	EnergyDWMYPeriod& operator=(const EnergyDWMYPeriod& thePeriod);
//	IObject* clone() const;

    virtual void validatePop2Object();

    string getPeriod() const;
    int getCount() const;
    string getInterval() const;
    int getSize() const { return period.size(); }

	// Overwrite the period to that between contracts. Hidden from constructing interface.
    void setToContractPeriod(int aContractType);

private:

    
    string period;
    int count;
    string interval;

};

typedef smartPtr<EnergyDWMYPeriod> EnergyDWMYPeriodSP;
typedef smartConstPtr<EnergyDWMYPeriod> EnergyDWMYPeriodConstSP;

DRLIB_END_NAMESPACE

#endif
