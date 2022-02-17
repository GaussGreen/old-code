/*!
 *
 * Copyright (c) CDC IXIS CM July 2003 Paris
  */



#ifndef _GP_CF_TARNPROXY_H
#define _GP_CF_TARNPROXY_H

/// this header has to come first
/// as env.h defines pre-processor constants for checking and validation
/// namely for std::vector<double> and ARM_Matrix

#include "gpbase/env.h"
#include "gpbase/port.h"
#include "gpbase/rootobject.h"
#include "gpbase/gpvector.h"
#include "gpbase/vectormanip.h"
#include "gpbase/typedef.h"

CC_BEGIN_NAMESPACE( ARM )

class ARM_DensityFunctor;

class ARM_TarnProxy : public ARM_RootObject
{
private:
	std::vector<double>				itsResetDates;
	std::vector<double>				itsDiscountFactors;
	std::vector<double>				itsLevPrec;
	std::vector<double>				itsLev;
	std::vector<double>				itsAdd;
	std::vector<double>				itsCap;
	std::vector<double>				itsFloor;
	std::vector<double>				itsFees;
	std::vector<double>				itsDcf;
	double						itsTarget;
	bool						itsGlobalCap;
	bool						itsGlobalFloor;
	ARM_GP_VectorPtr			itsGrid;
	ARM_VectorPtrVector			itsFunc;

	std::vector<double>				itsStdCpn;
	std::vector<double>				itsStdFund;
	std::vector<double>				itsRbCpn;
	std::vector<double>				itsRbFund;
	std::vector<double>				itsGC;
	std::vector<double>				itsGF;
	std::vector<double>				itsProba;

	double						itsStdCpnLeg;
	double						itsStdFundLeg;
	double						itsRbCpnLeg;
	double						itsRbFundLeg;
	double						itsGCLeg;
	double						itsGFLeg;
	
public:
	ARM_TarnProxy(	const std::vector<double>* resetDates,
					const std::vector<double>* df,
					const std::vector<double>* levprec,
					const std::vector<double>* lev,
					const std::vector<double>* fix,
					const std::vector<double>* cap,
					const std::vector<double>* floor,
					const std::vector<double>* fees,
					const std::vector<double>* dcf,
					double C_Target,
					bool C_Globalcap,
					bool C_Globalfloor);

	ARM_TarnProxy( const ARM_TarnProxy& rhs):itsTarget(rhs.itsTarget){};
	~ARM_TarnProxy(){};

	void Build(const std::vector<double>* fwds,const vector< ARM_DensityFunctor* >& densityFunctors);
	void Price(double C_CorrelInput,int C_NbSimul);
	double Index(int i, double sim) const;
	double GetPrice(int i) const;

	/// Standard ARM object support
	virtual ARM_Object*		Clone() const { return new ARM_TarnProxy(*this); };
	virtual string			ExportShortName() const { return "LPROX";}
	virtual	string			toString(const string& indent="",const string& nextIndent="") const;
};

CC_END_NAMESPACE()

#endif
//-----------------------------------------------------------------------------
/*---- End of file ----*/

