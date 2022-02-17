/*!
 *
 * Copyright (c) CDC IXIS CM November 2003 Paris
 *
 *	\file basket.h
 *
 *  \brief base class for basket calibration
 *
 *	\author  J. Messines
 *	\version 1.0
 *	\date August 2006
 */


#ifndef _INGPCALIB_BASKET_H
#define _INGPCALIB_BASKET_H

/// this header has to come first
/// as env.h defines pre-processor constants for checking and validation
/// namely for ARM_GP_Vector and ARM_Matrix

/// gpbase
#include "gpbase/env.h"
#include "gpbase/port.h"
#include "gpbase/typedef.h"
#include "gpbase/countedptr.h"

#include "gpbase/rootobject.h"
#include "gpbase/assignop.h"
#include "gpbase/gpvector.h"
#include "gpbase/gpmatrix.h"
#include "gpbase/gplinalgtypedef.h"

/// gpcalib
#include "gpcalib/typedef.h"

/// kernel
#include <util/refvalue.h>

class ARM_Security;
class ARM_Model;
class ARM_SpreadOption;
class ARM_Swaption;
class ARM_ZeroCurve;

CC_BEGIN_NAMESPACE( ARM )

class ARM_BasketSensi;
class ARM_PricingModel;
class ARM_DateStrip;

struct ARM_VanillaSpreadOptionArg;


class ARM_BasketCalib : public ARM_RootObject
{
public:
	enum BasketType
	{
		FULL = 0,
		SIMPLE
	};
	enum BasketStrike
	{
		ATM = 0,
		EQUIVALENT
	};
public:
    ARM_BasketCalib(ARM_DateStrip* ds,const ARM_ReferenceValue& notionalProfile,const ARM_ReferenceValue& feesProfile,double side,BasketType bt = FULL, BasketStrike strike = EQUIVALENT);
    ARM_BasketCalib(const ARM_BasketCalib& rhs);
	ASSIGN_OPERATOR(ARM_BasketCalib)
    virtual ~ARM_BasketCalib();

	void Compute(vector<ARM_Security*> securities, vector<ARM_Model*> models, const ARM_GP_Vector& weights);
	void Price(ARM_PricingModel* mkmo);

	void CreateAndSetPortfolio(ARM_ZeroCurve* curve);

	ARM_Swaption* CreateVarNotionalSwaptionAtExer (ARM_ZeroCurve* curve,int exerIndex, double fee);

	bool FillNotionals(ARM_Vector& fixAbs,ARM_Vector& fixNotios,const ARM_GP_Vector& endDates,const ARM_GP_Vector& payDates,double& fixNotional);

	inline  const ARM_Date&		GetRealEndDate() const	{ return itsRealEndDate;}
	inline  const ARM_Date&		GetEndDate() const		{ return itsEndDate;} 
	inline ARM_StdPortfolioPtr	GetPortfolio()			{ return itsPortfolio;};
	inline const ARM_GP_Matrix&	GetBasketWeights()		{ return itsBasketCoefs;};

	ARM_DateStrip* GetExerDateStrip() const { return itsExerDateStrip; }


	
    /// Standard ARM object support
	virtual string ExportShortName() const { return "LBSKT";}
	virtual ARM_Object* Clone() const;
    virtual string toString(const string& indent="", const string& nextIndent="") const;
    virtual void View(char* id = NULL, FILE* ficOut = NULL) const;

private:
	ARM_DateStrip*		itsExerDateStrip;			// updated by constructor (clone)
	
	ARM_ReferenceValue itsNotional;					// ...
	ARM_ReferenceValue itsCallFees;					// ...

	double				itsSide;

	BasketType			itsBasketType;				// ...
	BasketStrike		itsBasketStrike;			// ...

	vector<ARM_BasketSensi*>	itsStructLeg;		// updated in Compute

	ARM_StdPortfolioPtr	itsPortfolio;				// updated in CreateAndSetPortfolio ... (called in Compute)

	ARM_Date itsRealEndDate;						// updated in Compute
	ARM_Date itsEndDate;							// updated in Compute

	ARM_GP_Matrix		itsBasketCoefs;
	
};

CC_END_NAMESPACE()

#endif
//-----------------------------------------------------------------------------
/*---- End of file ----*/
