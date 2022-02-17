/*!
 *
 * Copyright (c) CDC IXIS CM July 2003 Paris
  */



#ifndef _GP_CF_TARNPROXY_H
#define _GP_CF_TARNPROXY_H

/// this header has to come first
/// as env.h defines pre-processor constants for checking and validation
/// namely for ARM_GP_Vector and ARM_Matrix

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
	ARM_GP_Vector				itsResetDates;
	ARM_GP_Vector				itsDiscountFactors;
	ARM_GP_Vector				itsLevPrec;
	ARM_GP_Vector				itsLev;
	ARM_GP_Vector				itsAdd;
	ARM_GP_Vector				itsCap;
	ARM_GP_Vector				itsFloor;
	ARM_GP_Vector				itsFees;
	ARM_GP_Vector				itsDcf;
	double						itsTarget;
	bool						itsGlobalCap;
	bool						itsGlobalFloor;
	ARM_GP_VectorPtr			itsGrid;
	ARM_VectorPtrVector			itsFunc;

	ARM_GP_Vector				itsStdCpn;
	ARM_GP_Vector				itsStdFund;
	ARM_GP_Vector				itsRbCpn;
	ARM_GP_Vector				itsRbFund;
	ARM_GP_Vector				itsGC;
	ARM_GP_Vector				itsGF;
	ARM_GP_Vector				itsProba;

	double						itsStdCpnLeg;
	double						itsStdFundLeg;
	double						itsRbCpnLeg;
	double						itsRbFundLeg;
	double						itsGCLeg;
	double						itsGFLeg;
	
public:
	ARM_TarnProxy(	const ARM_GP_Vector& resetDates,
					const ARM_GP_Vector& df,
					const ARM_GP_Vector& levprec,
					const ARM_GP_Vector& lev,
					const ARM_GP_Vector& fix,
					const ARM_GP_Vector& cap,
					const ARM_GP_Vector& floor,
					const ARM_GP_Vector& fees,
					const ARM_GP_Vector& dcf,
					double C_Target,
					bool C_Globalcap,
					bool C_Globalfloor);

	ARM_TarnProxy( const ARM_TarnProxy& rhs):itsTarget(rhs.itsTarget){};
	~ARM_TarnProxy(){};

	void Build(const ARM_GP_Vector& fwds,const vector< ARM_DensityFunctor* >& densityFunctors);
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

