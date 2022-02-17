/*!
 *
 * Copyright (c) IXIS-CIB March 2006
 *
 *	\file hedge.h
 *  \brief file for the abstract pricer
 *	\author  M. Bernardo
 *	\version 1.0
 *	\date Agust 2006
 */

#ifndef INXXXPROJECT_HEDGE
#define INXXXPROJECT_HEDGE

#include <gpbase/removeidentifiedwarning.h>
#include <gpbase/port.h>
#include <gpbase/rootobject.h>
#include <gpbase/gpvector.h>
#include <gpbase/typedef.h>
#include <gpbase/assignop.h>

#include "xxxproject/mktdatas.h"

#include<vector>
#include<string>
using namespace std;

class ARM_Security;

CC_BEGIN_NAMESPACE( ARM )

class	ARM_MktData;
class	ARM_Scenario;
class	ARM_Pricer;
class	ARM_GramFctorArgDict;


typedef ARM_StringVector::iterator It_Str;

class ARM_Hedge: public ARM_RootObject
{

public:
	ARM_Hedge(ARM_Object* theSecurity, ARM_MktData* theMktData){ Init(	theSecurity,  theMktData);	}
	virtual ~ARM_Hedge(){}
	ARM_Hedge(const ARM_Hedge& rhs);
	ARM_Object*		Clone					(	)	const	{ return new ARM_Hedge(*this); }
	ASSIGN_OPERATOR(ARM_Hedge)
	

	string					toString		(	const string& indent="", const string& nextIndent="") const;
	virtual string			ExportShortName	(	)	const	{ return "LHEDG";}
	void					ComputeHedge	( ARM_Scenario* );
	void					ApplyScenario	( ARM_Scenario*, int nbShift );
	void					Init			(	ARM_Object* ,ARM_MktData* );
	ARM_GramFctorArgDict	GetDictionnary	(	)			{	return itsFunctor; }

private:

	void					Apply_0D_Scenario ( ARM_Scenario*, bool computeHedge);
	void					Apply_1D_Scenario ( ARM_Scenario*, int nbShift, bool computeHedge);
	void					Apply_2D_Scenario ( ARM_Scenario*, int nbShift, bool computeHedge);
	void					Apply_ND_Scenario ( ARM_Scenario*, int nbShift, bool computeHedge);


protected:

	double					itsInitPrice;
	ARM_Pricer*				itsPricer;
	ARM_MktData*			itsMktData;
	ARM_GramFctorArgDict	itsFunctor;

};

CC_END_NAMESPACE()

#endif

/*----------------------------------------------------------------------------*/
/*---- End of file ----*/