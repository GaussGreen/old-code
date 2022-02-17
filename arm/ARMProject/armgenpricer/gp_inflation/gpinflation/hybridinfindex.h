
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

#ifndef _INGPINFLATION_HYBRIDINFINDEX_H
#define _INGPINFLATION_HYBRIDINFINDEX_H

#define RateBase	CC_NS( ARM_Constants, rateBase )
#define VolBase		CC_NS( ARM_Constants, volBase )

#include <inst/swapleg.h>

#include <gpbase/gplinalgtypedef.h>
#include <gpbase/gpvector.h>
#include <crv/volcurv.h>

#include <gpbase/rootobject.h>
#include <gpbase/typedef.h>
#include <gpbase/stringconvert.h>

#include <gpinfra/mktdatamanagerrep.h>
#include <gpinfra/mktdatamanager.h>

#include <map>


typedef enum { 	IN,	IR, NO }	IndexType;

CC_BEGIN_NAMESPACE( ARM )

class ARM_InfLeg;

struct ARM_InfIrIndex{

	ARM_InfIrIndex( ){}
	ARM_InfIrIndex( const string& indexName,
					map<string,ARM_Date>	&	mDate, 
					map<string,string>		&	mString,
					map<string,int>			&	mInt);

	ASSIGN_OPERATOR	( ARM_InfIrIndex );
	ARM_InfIrIndex( const ARM_InfIrIndex & rhs );
	virtual ARM_InfIrIndex*	Clone	( ) const	{ return new ARM_InfIrIndex (*this);	}
	virtual ~ARM_InfIrIndex	( )	{ };

	ARM_GP_VectorPtr	GetNumDates ( ) const;
	ARM_GP_VectorPtr	GetDemDates ( ) const;
	ARM_GP_VectorPtr	GetBegDates ( ) const;
	ARM_GP_VectorPtr	GetEndDates ( ) const;

	ARM_GP_VectorPtr	GetResLags  ( ) const;
	ARM_GP_VectorPtr	GetTenLags  ( ) const;

	string				GetIndexName( ) const;
	string				GetName		( ) const;

	void				SetMod( ARM_BSSmiledModel  * mod);

	double				CptFwd( const double & res,
								const double & mat,
								const double & pay) const;

	double				CptVol(	const double & exp,
								const double & ten,
								const double & fwd, 
								const double & str) const;

	double				CptAdj( const double & res,
								const double & mat,
								const double & pay,
								const double & fwd) const;

	double				CptQua(	const double & exp,
								const double & ten,
								const double & fwd, 
								const double & gue,
								const double & str) const;   // compute quantile

	ARM_CountedPtr<ARM_IRIndex>		GetIndex	( ) const;

	static ARM_GP_VectorPtr ConvertToARM_GP_Vector ( const ARM_Vector * vec);

	static 	void		SmoothQua( const ARM_GP_Vector&, ARM_GP_Vector& );

private:
	double	CptIrQuantile(	const double & exp,
							const double & ten,
							const double & fwd, 
							const double & str) const;

	double	CptInQuantile(	const double & exp,
							const double & ten,
							const double & fwd, 
							const double & gue,
							const double & str) const;

private:
	ARM_CountedPtr<ARM_BSSmiledModel>		itsMod;
	ARM_Date								itsAsOfDate;
	
public:
	ARM_CountedPtr<ARM_SwapLeg>		itsIndex;
	IndexType						itsType;
};

typedef ARM_CountedPtr<ARM_InfIrIndex > ARM_InfIrIndexPtr;
typedef ARM_CountedPtr<ARM_VolCurve> ARM_VolCurvePtr;

struct ARM_InfIrCorrel{

	ARM_InfIrCorrel( ){ }
	
	ASSIGN_OPERATOR	( ARM_InfIrCorrel );
	ARM_InfIrCorrel( const ARM_InfIrCorrel & rhs );
	virtual ARM_InfIrCorrel*	Clone	( ) const	{ return new ARM_InfIrCorrel (*this);	}
	virtual ~ARM_InfIrCorrel	( )	{ };

	ARM_VolCurvePtr	GetCorrel	( const string &, const string &) const;
	double			GetCorrel	( const string &, const string &, const double &, const double &) const;
	void			SetCorrel	( const string &, const string &, ARM_VolCurve *);

private:
	map< pair<string,string>,  ARM_VolCurvePtr>  itsCorrel;

};

struct ARM_VolParam{

	ARM_VolParam( const ARM_InfIrIndex & idx, const double & res, const double & ten):itsIdx(idx),itsRes(res),itsTen(ten){}
	ARM_VolParam( const ARM_VolParam & rhs):itsIdx(rhs.itsIdx),itsRes(rhs.itsRes),itsTen(rhs.itsTen){}

	virtual ~ARM_VolParam(){};

	ARM_InfIrIndex	itsIdx;
	double			itsRes;
	double			itsTen;
};


CC_END_NAMESPACE()

#endif
/*--------------------------------------------------------------------------*/
/*---- End of file ----*/


