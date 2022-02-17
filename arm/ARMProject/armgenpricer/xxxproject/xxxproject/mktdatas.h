/*!
 *
 * Copyright (c) IXIS-CIB March 2006
 *
 *	\file mktdatas.h
 *  \brief file for market datas
 *	\author  R. Guillemot
 *	\version 1.0
 *	\date March 2006
 */

#ifndef INXXXPROJECT_MKTDATAS
#define INXXXPROJECT_MKTDATAS



#include <gpbase/removeidentifiedwarning.h>
#include "xxxproject/observator.h"

#include <gpbase/port.h>
#include <gpbase/rootobject.h>
#include <gpbase/assignop.h>
#include <gpbase/gpvector.h>
#include <map>
#include <vector>
#include <string>

CC_USING_NS( std, map		)
CC_USING_NS( std, vector	)
CC_USING_NS( std, string	)

class ARM_ZeroCurve;
class ARM_BasisCurve;
class ARM_VolCurve;
class ARM_BSModel;
class ARM_Forex;
class ARM_Date;

class ARM_Subject;
class ARM_Observator;


CC_BEGIN_NAMESPACE( ARM )


class ARM_MktData : public ARM_RootObject
{

public:

	ARM_MktData							( )										{ }
	ARM_MktData							(	const ARM_Date &, const V_Str	&,	const V_Obj	&	);
	ARM_MktData							(	const V_Sub	& ,	const M_Obs	&	);
	ARM_MktData							(	const ARM_MktData & mkt			);
	ASSIGN_OPERATOR						(	ARM_MktData	);

	virtual ARM_Object*		Clone		( ) const								{ return new ARM_MktData(*this);	}
	// This function create a deep copy of the object
	// The assign operator is not doing the same thing ....
	ARM_MktData*			CreateCopy	( ) const;
	virtual ~ARM_MktData				( )										{ Reset( );	}	
	virtual	ARM_Observator*	GetObs		(	const string &					)	{ return NULL; }
	virtual string			toString	(	const string & indent="", 	const string & nextIndent="") const;

	virtual string			ExportShortName	( )	const { return "LMKTD";}	

private:
	void					InitMktData(const V_Str	&,	const V_Obj	&);
	void					Init			(	const string &, ARM_Object* );
	void					Reset			( );
	static string			GetInstrument	(	const string & );
	static string			GetDomCcy		(	const string & );
	static string			GetLabel		(	const string & );

public:
	M_Obs					itsObs;
	V_Sub					itsSub;

};

template< class T  = ARM_Subject >
class ARM_MktSubject: public ARM_MktData{

public:

	ARM_MktSubject< T >						( const ARM_MktData & mkt ): ARM_MktData(mkt)	{ }

	ASSIGN_OPERATOR							(	ARM_MktSubject< T >	);

	virtual ARM_Object*		Clone			(	) const										{ return new ARM_MktSubject< T >(*this);	}

	
	~ARM_MktSubject							(	)											{ }

public:

	ARM_Object*					Get			(	const string &  ) ;

	void						Set			(	const string & ,  const ARM_Object*	);

	void						Finalize	(	const string &	);

	virtual	ARM_Observator*		GetObs		(	const string &	);

	virtual	V_Str				GetLabel	(	const string &  , const string & ) const;

};

template< class T >
ARM_Object*	ARM_MktSubject< T >::Get(	const string &  ccy) {

	for ( It_Sub its = itsSub.begin(); its != itsSub.end(); its++)
		if ( dynamic_cast<T* > (*its) )  	return (*its)->Get( ccy );	

	return NULL;
}

template< class T >
V_Str	ARM_MktSubject< T >::GetLabel(	const string &  ccy, const string & lab) const{
	V_Str tmp;

	for ( Ct_Sub its = itsSub.begin(); its != itsSub.end(); its++){
		if ( dynamic_cast<T* > (*its) )  {
			tmp =(*its)->GetLabel( ccy, lab );
			break;
		}
	}
	return tmp;
}

template< class T >
void ARM_MktSubject< T >::Set(	const string & ccy, const ARM_Object* obj 	){

	for ( It_Sub its = itsSub.begin(); its != itsSub.end(); its++){
		if ( dynamic_cast<T* > (*its) ) { 
			(*its)->Set( ccy, obj );	
			break;	
		}
	}
}

template<class T>
void			ARM_MktSubject< T >::Finalize	(	const string &	ccy ){

	for ( It_Sub its = itsSub.begin(); its != itsSub.end(); its++)
		if ( dynamic_cast<T* > (*its) ) { 		(*its)->Equalize(ccy);		break;		}
}

template<class T>
ARM_Observator*		ARM_MktSubject< T >::GetObs		(	const string &  ccy){

	for ( It_Sub its = itsSub.begin(); its != itsSub.end(); its++)
		if ( dynamic_cast<T* > (*its) )  	return  (*its)->GetModel( ccy );
	return NULL;
}


CC_END_NAMESPACE()



#endif

