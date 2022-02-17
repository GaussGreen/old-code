//----------------------------------------------------------------------------
//
//   Group       : GCCG Derivatives Research
//
//   Filename    : EnergyTermPeriodBenchmark.hpp
//
//   Description : Defines a deal period for energy. 
//
//   Author      : Sean Chen
//
//   Date        : 26 Sept. 2005
//
//----------------------------------------------------------------------------

#ifndef _EnergyTermPeriodBenchmark_H_
#define _EnergyTermPeriodBenchmark_H_

#include "edginc/config.hpp"
#include "edginc/Class.hpp"
#include "edginc/EnergyTermPeriod.hpp"

#include <string>
using namespace std;

DRLIB_BEGIN_NAMESPACE


class TOOLKIT_DLL EnergyTermPeriodBenchmark : public EnergyTermPeriod
{

public:

	static CClassConstSP const TYPE;

	enum PeriodFormat {MONTH,QUARTER,SEMIANNUAL,ANNUAL,INJECTION,WITHDRAWL};
	
    friend class EnergyTermPeriodBenchmarkHelper;

    //EnergyTermPeriodBenchmark(const EnergyTermPeriodBenchmark& thePeriod);

    void validatePop2Object();

    virtual ~EnergyTermPeriodBenchmark();

   static bool adjustYear(int* year);

protected:

    EnergyTermPeriodBenchmark();

	string benchmark;

};

typedef smartPtr<EnergyTermPeriodBenchmark> EnergyTermPeriodBenchmarkSP;
typedef smartConstPtr<EnergyTermPeriodBenchmark> EnergyTermPeriodBenchmarkConstSP;

DRLIB_END_NAMESPACE

#endif
