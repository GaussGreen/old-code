
/*	*	*	*	*	*	*	*	*	*	*	*	*	*	*	*	*	*	*	*	*	*	*	 /
 *																							 *
 *			Class Hybrid Inf Ir Leg															 *
 *																							 *
 *			This class builds a hybrid inf ir leg from swap legs and  Inf Swap Leg			 *									 *
 *																							 *
 *			Author	: Mathieu Bernardo														 *
 *			Date	: April, 2nd 2007														 *																											 *
 *			Copyright (c) NatIxis															 *
 *			Version	: 1.0																	 *
 *																							 *
 *	*	*	*	*	*	*	*	*	*	*	*	*	*	*	*	*	*	*	*	*	*	*	*/

#ifndef _INGPINFLATION_HYBRIDINFIRLEG_H
#define _INGPINFLATION_HYBRIDINFIRLEG_H

#define LAG			13

#define RateBase	CC_NS( ARM_Constants, rateBase )
#define VolBase		CC_NS( ARM_Constants, volBase )

#include <inst/swapleg.h>

#include <gpbase/gplinalgtypedef.h>
#include <gpbase/gpvector.h>
#include <gpbase/rootobject.h>
#include <gpbase/typedef.h>
#include <gpbase/stringconvert.h>

#include <gpinfra/gramfunctorargdict.h>
#include <gpinfra/gramfunctorarg.h>
#include <gpinfra/mktdatamanagerrep.h>
#include <gpinfra/mktdatamanager.h>

#include "gpinflation/hybridinfindex.h"
#include "gpinflation/typedef.h"

#include <map>

CC_BEGIN_NAMESPACE( ARM )

typedef ARM_GP_VectorPtr								ARM_GP_VecPtr;
typedef map <IndexType,ARM_GP_VecPtr>					ARM_GP_MapPtr;
typedef map <pair<IndexType,IndexType>,ARM_GP_VecPtr >	ARM_GP_CorPtr;
typedef map <IndexType,ARM_InfIrIndex>					ARM_GP_MapIdx;
typedef ARM_GP_MapIdx::iterator							idxIter;
typedef ARM_GP_CorPtr::iterator							corIter; 

class ARM_HybridInfIrLeg: public ARM_Security {

public:
	
	typedef ARM_GP_MapIdx::iterator							idxIter;
	typedef ARM_GP_CorPtr::iterator							corIter; 


	ARM_HybridInfIrLeg( ){ }
	ARM_HybridInfIrLeg(	map<string,ARM_Date> &, map<string,string>	&,	map<string,int>	& );
	ARM_HybridInfIrLeg(   const ARM_HybridInfIrLeg &	);

	ASSIGN_OPERATOR	( ARM_HybridInfIrLeg );

	virtual ARM_Object*	Clone	( ) const	{ return new ARM_HybridInfIrLeg (*this);	}
	virtual ~ARM_HybridInfIrLeg	( )			{ 	};

	virtual void	BuildSchedule	( );
	virtual void	Init			( ARM_MarketData_ManagerRep* );
	virtual void	Validate		( ARM_MarketData_ManagerRep* );
	virtual void	View			( char* id = NULL, FILE* ficOut = NULL);
	virtual string	ViewGlob		( );
	virtual string	ViewInit		( );

	ARM_GP_MapIdx			GetIdx	( ) const { return itsIdx;}
	ARM_GP_MapPtr			GetRes	( ) const { return itsRes;}
	ARM_GP_MapPtr			GetTen	( ) const { return itsTen;}
	ARM_GP_MapPtr			GetBeg	( ) const { return itsBeg;}
	ARM_GP_MapPtr			GetEnd	( ) const { return itsEnd;}
	ARM_GP_MapPtr			GetPay	( ) const { return itsPay;}

	ARM_GP_VecPtr			GetIntTerm( ) const { return itsIntTerms;}
	ARM_GP_VecPtr			GetDisFact( ) const { return itsDiscFact;}
	ARM_GP_VecPtr			GetResDate( ) const { return itsResDates;}

	ARM_Date				GetAsOfDate( ) const { return itsAsOfDate;}
	ARM_InfIrCorrel			GetCorObj  ( ) const { return itsCorObj;  }
	ARM_GramFctorArgDict	GetFunctor ( ) const { return itsFunctor; } 

	template< typename A = string, typename B = string>
	struct Ostream{
		string	ToStream( const A & str1, const B & str2)	{
			CC_Ostringstream	os;
			os << CC_NS(std,setfill(' ')) << CC_NS(std,fixed) << CC_NS(std,setw)(30);
			os <<CC_NS(std,right)<<str1;
			os <<"  :  "<< str2<<"\n";
			return os.str();	
		};
	};

protected:

	ARM_GP_VecPtr	itsBegDates;	
	ARM_GP_VecPtr	itsEndDates;
	ARM_GP_VecPtr	itsResDates;		 
	ARM_GP_VecPtr	itsPayDates;
	ARM_GP_VecPtr	itsNumDates;		
	ARM_GP_VecPtr	itsDemDates;
	
	ARM_GP_VecPtr	itsIntTerms;
	ARM_GP_VecPtr	itsDiscFact;
	ARM_GP_VecPtr	itsIntrDays;

	ARM_GP_MapIdx	itsIdx;

	ARM_GP_MapPtr	itsRes;
	ARM_GP_MapPtr	itsTen;
	ARM_GP_MapPtr	itsBeg;
	ARM_GP_MapPtr	itsEnd;
	ARM_GP_MapPtr	itsPay;

	ARM_GP_MapPtr	itsFwd;
	ARM_GP_MapPtr	itsAdj;
	ARM_GP_MapPtr	itsVol;
	ARM_GP_CorPtr	itsCor;

	ARM_Date		itsAsOfDate;
	ARM_Date		itsStartDate;
	ARM_Date		itsEndDate;

	int						itsNbFlows;
	bool					isHybrid;
	bool					isInit;
	string					itsCurrency;
	ARM_InfIrCorrel			itsCorObj;
	ARM_GramFctorArgDict	itsFunctor;
};



CC_END_NAMESPACE()

#endif
/*--------------------------------------------------------------------------*/
/*---- End of file ----*/


