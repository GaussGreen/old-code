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

#ifndef INXXXPROJECT_OBSERVATOR
#define INXXXPROJECT_OBSERVATOR

#define		NBTERM 4

#include <gpbase/removeidentifiedwarning.h>
#include <gpbase/port.h>
#include <gpbase/assignop.h>
#include <gpbase/gpvector.h>
#include <glob\linalg.h>

#include <map>
#include <vector>
#include <string>


CC_USING_NS( std, map )
CC_USING_NS( std, vector )
CC_USING_NS( std, string )

class ARM_ZeroCurve;
class ARM_BasisCurve;
class ARM_VolCurve;
class ARM_BSModel;
class ARM_Forex;
class ARM_Object;
class ARM_BSSmiledModel;
class ARM_PricingModel;

#define CAPCONST	0
#define VFXCONST	4
#define MIXCONST	7
#define IBSCONST	12

CC_BEGIN_NAMESPACE( ARM )


class ARM_Subject;
class ARM_Observator;

typedef map<string,int>				M_Str;
typedef vector<string >				V_Str;

typedef vector<ARM_Object*  >		V_Obj;
typedef V_Obj::iterator				It_Obj;
typedef V_Obj::const_iterator		Ct_Obj;

typedef vector<ARM_Subject* >		V_Sub;
typedef V_Sub::iterator				It_Sub;
typedef V_Sub::const_iterator		Ct_Sub;

typedef vector< ARM_Observator* >	V_Obs;
typedef V_Obs::iterator				It_Obs;
typedef V_Obs::const_iterator		Ct_Obs;

typedef map< string, V_Obs >		M_Inv;
typedef M_Inv::iterator				It_Inv;
typedef M_Inv::const_iterator		Ct_Inv;

typedef map<string, ARM_Object* >	M_Mkt;
typedef M_Mkt::iterator				It_Mkt;
typedef M_Mkt::const_iterator		Ct_Mkt;

typedef map< ARM_Object* ,M_Mkt >	M_Mod;
typedef M_Mod::iterator				It_Mod;
typedef M_Mod::const_iterator		Ct_Mod;

typedef map<string, ARM_Observator*>M_Obs;
typedef M_Obs::iterator				It_Mbs;
typedef M_Obs::const_iterator		Ct_Mbs;

	
class ARM_Observator
{

public:
	ARM_Observator							( )						{ }
	ARM_Observator							( ARM_Object* 	, V_Sub, string & );
	virtual					~ARM_Observator	( );
	void					SetObject		( ARM_Object* obj)		{	itsObject = obj;	}
	ARM_Object*				GetObject		( )						{	return itsObject;	}
	V_Sub					GetSubject		( )						{	return itsSubject;	}
	ARM_Object*				ReduceObject	( );
	M_Str					GetCurrency		( );
	void					SetDomCcy		( string );
		
private:
	M_Str					itsCurrency;					
	ARM_Object*				itsObject;
	V_Sub					itsSubject;

};

class ARM_Subject: public ARM_RootObject{

public:
	ARM_Subject								( ):ARM_RootObject(){	Reset();	}
	ARM_Subject								( const ARM_Subject &	);
	ASSIGN_OPERATOR							( ARM_Subject			);
	virtual ARM_Object*		Clone			( )						const	{	return new ARM_Subject(*this); }
	virtual					~ARM_Subject	( )	{	Reset();	}

	void					Attach			( ARM_Observator*	, string & );
	void					Detach			( ARM_Observator*	);
	virtual void			BuildModel		( ARM_Observator*   , const string & ) { } 

	void					Equalize		( const string  &	);	
	void					Set				( const string  &,		const	ARM_Object* );
	ARM_Object*				Get				( const string  &	)	const;
	ARM_Observator*			GetModel		( const string  &	)	const;
	V_Str					GetCurrency		( )						const;
	V_Str					GetCurrency		( ARM_Object	*	)	const;
	virtual V_Str			GetLabel		( const string  &	,	const string & ) const	{ V_Str tmp; return tmp; };

	static  V_Sub			CreateSubject	( );
 	static  V_Sub			AttributeSubject( const string &, vector<ARM_Subject* >* , ARM_Object*);

	virtual string			toString		( const string &  indent="", const string & nextIndent="") const{return "";}
	virtual string			toString		( const	string & , const string  & indent="", const string & nextIndent="") const{return "";}
	virtual string			GetSubjectInfo  ( ) const { return string (""); }
	static  string			GetInstrument	( )	{ return string (""); }

protected:

	void					Notify			( const string &						);
	virtual void			From			( ARM_Object * ,	ARM_Object *,	int )	{ } 
	virtual ARM_Object*		To				( ARM_Object * ,	int dim=1 )	const	{ return NULL;	}
	void					Reset			( );
	void					Map				( const string&,	const	ARM_Object* );
	virtual string			GetCcy			( ARM_Object *						)	const	{ 	return string(""); }

protected:
	static	int				itsCcyDim;
	M_Inv					itsModelInv;
	M_Mod					itsModel;
	M_Mkt					itsCurrentMkt;
	M_Mkt					itsInitialMkt;

};


class ARM_FX_Subject: public ARM_Subject{

public:
	ARM_FX_Subject():ARM_Subject			( )	{ }
	ARM_FX_Subject							( const ARM_FX_Subject & sub	):	ARM_Subject(sub) { }
	ASSIGN_OPERATOR							( ARM_FX_Subject				);		
	virtual ARM_Object*		Clone			( ) const						{	return new ARM_FX_Subject(*this); }	
	virtual ~ARM_FX_Subject					( )	{ }

	static string			GetInstrument	( )								{	return "DELTA_FX";	}
	virtual string			toString		( const	string & , const string  & indent="", const string & nextIndent="") const;
	virtual string			GetSubjectInfo  ( )	const						{	return "FOREX SPOT";}

	static void				Bump			( ARM_Object *, const ARM_Matrix & , const int &);


private:
	virtual string			GetCcy			( ARM_Object *	obj	)	const;
	virtual void			From			( ARM_Object *,	ARM_Object *,	int dim = 1	);
	virtual ARM_Object*		To				( ARM_Object *,	int dim = 1	)	const;		
};


class ARM_1D_Subject: public ARM_Subject{

public: 
	ARM_1D_Subject():ARM_Subject			( )	{ }
	ARM_1D_Subject							( const ARM_1D_Subject & sub	):	ARM_Subject(sub)	{ }
	ASSIGN_OPERATOR							( ARM_1D_Subject				)
	virtual ARM_Object*		Clone			( ) const						{	return new ARM_1D_Subject(*this); }	
	virtual ~ARM_1D_Subject					( )	{ }

	static  void			Bump			( ARM_Object *, const ARM_Matrix & , const int &);
	virtual string			toString		( const	string & , const string  & indent="", const string & nextIndent="") const;

protected:
	virtual string			GetCcy			( ARM_Object* obj				)	const;
};

class ARM_YC_Subject: public virtual ARM_1D_Subject{

public:
	ARM_YC_Subject():ARM_1D_Subject			( )	{ }
	ARM_YC_Subject							( const ARM_YC_Subject & sub	):ARM_1D_Subject	(sub)	{ }
	ASSIGN_OPERATOR							( ARM_YC_Subject				);
	virtual ARM_Object*		Clone			( ) const						{	return new ARM_YC_Subject(*this); }
	virtual ~ARM_YC_Subject					( ){ }

	static  string			GetInstrument	( )								{	return "DELTA_YC";	}
	virtual string			GetSubjectInfo  ( ) const						{	return "ZERO CURVE";}

protected:
	virtual void			From			( ARM_Object *, ARM_Object *,	int dim = 1);
	ARM_Object*				To				( ARM_Object *,	int dim =1 ) 	const;	
};

class ARM_INF_Subject: public virtual ARM_1D_Subject{

public:
	ARM_INF_Subject():ARM_1D_Subject		( )	{ }
	ARM_INF_Subject							( const ARM_INF_Subject & sub	):ARM_1D_Subject	(sub)	{ }
	ASSIGN_OPERATOR							( ARM_INF_Subject				);
	virtual ARM_Object*		Clone			( ) const						{	return new ARM_INF_Subject(*this); }
	virtual ~ARM_INF_Subject					( ){ }

	static  string			GetInstrument	( )								{	return "DELTA_INF";	}
	virtual string			GetSubjectInfo  ( ) const						{	return "INF CURVE";}

	static void				Bump			( ARM_Object *, const ARM_Matrix & , const int &);

protected:

	virtual void			From			( ARM_Object *, ARM_Object *,	int dim = 1);
	ARM_Object*				To				( ARM_Object *,	int dim =1 ) 	const;
};

class ARM_BS_Subject: public virtual ARM_1D_Subject{

public:
	ARM_BS_Subject():ARM_1D_Subject			( ){ }
	ARM_BS_Subject							(	const ARM_BS_Subject & sub	):ARM_1D_Subject(sub) { }
	ASSIGN_OPERATOR							(	ARM_BS_Subject				);
	virtual ARM_Object*		Clone			( ) const						{	return new ARM_BS_Subject(*this); }	
	virtual ~ARM_BS_Subject					( )								{ }

	static string			GetInstrument	( )								{	return "DELTA_BS";	}
	virtual string			GetSubjectInfo  ( )	const						{	return "BASIS CURVE";}

protected:
	virtual void			From			( ARM_Object * ,	ARM_Object *, int dim =1	);
	virtual ARM_Object*		To				( ARM_Object * ,	int dim =1	)	const;	
};

class ARM_DAT_Subject: public ARM_YC_Subject, ARM_BS_Subject, ARM_INF_Subject{

public:
	ARM_DAT_Subject():ARM_YC_Subject( ),ARM_BS_Subject( ),ARM_INF_Subject ( )	{ }
	ARM_DAT_Subject							( const ARM_DAT_Subject & sub	):	ARM_YC_Subject(sub), ARM_BS_Subject(sub),ARM_INF_Subject (sub){ }
	ASSIGN_OPERATOR							( ARM_DAT_Subject				);		
	virtual ARM_Object*		Clone			( ) const						{	return new ARM_DAT_Subject(*this); }	
	virtual ~ARM_DAT_Subject				( )	{ }

	static  string			GetInstrument	( )								{	return "THETA";	}
	virtual string			toString		( const	string & , const string  & indent="", const string & nextIndent="") const;
	virtual string			GetSubjectInfo  ( )	const						{	return "AS OF DATE";}

	static void				Bump			( ARM_Object *, const ARM_Matrix & , const int &);

private:
	virtual string			GetCcy			( ARM_Object *	obj	)	const;
	virtual void			From			( ARM_Object *,	ARM_Object *,	int dim = 1	);
	virtual ARM_Object*		To				( ARM_Object *,	int dim = 1	)	const;	

};

template< int n> 
class ARM_2D_Subject: public ARM_Subject{

public: 

	ARM_2D_Subject():ARM_Subject			( )	{ }
	ARM_2D_Subject							(	const ARM_2D_Subject & sub	):ARM_Subject(sub) { }
	ASSIGN_OPERATOR							(	ARM_2D_Subject<n>			)
	virtual ARM_Object*		Clone			( ) const						{	return new ARM_2D_Subject(*this); }	
	virtual ~ARM_2D_Subject					( )	{ }

	virtual V_Str			GetLabel		( const string &,const string & )	const;
	static  void			Bump			( ARM_Object *, const ARM_Matrix & , const int &);
	virtual string			toString		( const	string & , const string  & indent="", const string & nextIndent="") const;

protected:

	virtual string			GetCcy			( ARM_Object* obj								)	const;
	virtual ARM_VolCurve*	GetVol			( ARM_Object*									)	const;
	virtual void			SetVol			( ARM_Object*, ARM_VolCurve*					);
	virtual void			From			( ARM_Object *, ARM_Object * , int dim =1		);
	virtual ARM_Object*		To				( ARM_Object *, int dim =1 ) const;	
};


template< int n=CAPCONST> 
class ARM_CAP_Subject: public ARM_2D_Subject<n>{

public: 

	ARM_CAP_Subject(): ARM_2D_Subject<n>	( )									{ }
	ARM_CAP_Subject<n>						( const ARM_CAP_Subject<n> & sub)	: ARM_2D_Subject<n>(sub){ }
	ASSIGN_OPERATOR							( ARM_CAP_Subject<n>				)	
	virtual ARM_Object*		Clone			( ) const							{	return new ARM_CAP_Subject<n>(*this); } 
	virtual ~ARM_CAP_Subject				( )									{ }

	virtual string			GetSubjectInfo  ( )	const;
	static  string			GetInstrument	( );	

protected:
	virtual ARM_VolCurve*	GetVol			( ARM_Object *						)	const;
	virtual	void			SetVol			( ARM_Object *, ARM_VolCurve *		);
};

template< int n=CAPCONST> 
class ARM_OSW_Subject: public ARM_2D_Subject<n>{

public: 

	ARM_OSW_Subject():ARM_2D_Subject<n>		( )									{ }
	ARM_OSW_Subject<n>						( const ARM_OSW_Subject<n> & sub)	: ARM_2D_Subject<n>(sub){ }
	ASSIGN_OPERATOR							( ARM_OSW_Subject<n>				)
	virtual ARM_Object*	Clone				( ) const							{	return new ARM_OSW_Subject<n>(*this); } 
	virtual ~ARM_OSW_Subject				( )									{ }

	virtual string			GetSubjectInfo  ( )	const;
	static string			GetInstrument	( );

protected:
	virtual ARM_VolCurve*	GetVol			( ARM_Object *						)	const;
	virtual	void			SetVol			( ARM_Object *, ARM_VolCurve *		);
};

template< int n=CAPCONST> 
class ARM_SO_Subject: public ARM_2D_Subject<n>{

public: 

	ARM_SO_Subject():ARM_2D_Subject<n>		( )									{ }
	ARM_SO_Subject<n>						( const ARM_SO_Subject<n> & sub)	: ARM_2D_Subject<n>(sub){ }
	ASSIGN_OPERATOR							( ARM_SO_Subject<n>				)
	virtual ARM_Object*	Clone				( ) const							{	return new ARM_SO_Subject<n>(*this); } 
	virtual ~ARM_SO_Subject					( )									{ }

	virtual string			GetSubjectInfo  ( )	const;
	static string			GetInstrument	( );

protected:
	virtual ARM_VolCurve*	GetVol			( ARM_Object *						)		const;
	virtual void			SetVol			( ARM_Object *, ARM_VolCurve*		);

};


template<  int n=VFXCONST > 
class ARM_VFX_Subject: public ARM_2D_Subject<n>{

public: 

	ARM_VFX_Subject():ARM_2D_Subject<n>		( )										{ }
	ARM_VFX_Subject<n>						( const ARM_VFX_Subject<n> & sub	)	: ARM_2D_Subject<n>(sub){ }
	ASSIGN_OPERATOR							( ARM_VFX_Subject<n>				);
	virtual ARM_Object*	Clone				( ) const								{	return new ARM_VFX_Subject<n>(*this); } 
	virtual ~ARM_VFX_Subject				( )										{ }

	static  string			GetInstrument	( );
	virtual string			GetSubjectInfo  ( )	const;

	virtual V_Str			GetLabel		( const string&, const string &		) const;
	static  void			Bump			( ARM_Object *,  const ARM_Matrix & , const int &);

protected:

	virtual ARM_VolCurve*	GetVol			( ARM_Object *						)		const;
	virtual void			SetVol			( ARM_Object *, ARM_VolCurve*		);
	virtual string			GetCcy			( ARM_Object *						)		const;
};


template< int n=MIXCONST > 
class ARM_MIX_Subject: public ARM_2D_Subject<n>{

public:

	ARM_MIX_Subject():ARM_2D_Subject<n>		( ){ }
	ARM_MIX_Subject<n>						( const ARM_MIX_Subject<n> & sub	)	:ARM_2D_Subject<n>(sub) { }
	virtual ARM_Object*		Clone			( ) const								{	return new ARM_MIX_Subject<n> (*this); }	
	ASSIGN_OPERATOR							( ARM_MIX_Subject<n>				);
	virtual ~ARM_MIX_Subject				( )										{ }

	static string			GetInstrument	( );
	virtual string			GetSubjectInfo  ( )	const;
	
	virtual ARM_VolCurve*	GetVol			( ARM_Object *						)	const;
	virtual void			SetVol			( ARM_Object *, ARM_VolCurve*		);
	virtual string			GetCcy			( ARM_Object *						)	const;
	static void				Bump			( ARM_Object *, const ARM_Matrix & , const int& );
};

template< int n=IBSCONST> 
class ARM_IBS_Subject: public ARM_2D_Subject<n>{

public: 

	ARM_IBS_Subject():ARM_2D_Subject<n>		( )	{ }
	ARM_IBS_Subject<n>						( const ARM_IBS_Subject<n> & sub)	: ARM_2D_Subject<n>(sub){ }
	virtual ARM_Object*		Clone			( ) const							{	return new ARM_IBS_Subject<n> (*this); }	
	ASSIGN_OPERATOR							( ARM_IBS_Subject<n>				)
	virtual ~ARM_IBS_Subject				( )									{ }

	virtual string			GetSubjectInfo  ( )	const;
	static string			GetInstrument	( );

protected:

	virtual ARM_VolCurve*	GetVol			( ARM_Object *						)		const;
	virtual void			SetVol			( ARM_Object *, ARM_VolCurve*		);
};

class ARM_COR_Subject: public ARM_2D_Subject<0>{

public:
	ARM_COR_Subject():ARM_2D_Subject<0>		( ){ }
	ARM_COR_Subject							(	const ARM_COR_Subject & sub	):ARM_2D_Subject<0>(sub) { }
	ASSIGN_OPERATOR							(	ARM_COR_Subject				);
	virtual ARM_Object*		Clone			( ) const						{	return new ARM_COR_Subject(*this); }	
	virtual ~ARM_COR_Subject				( )								{ }

	static string			GetInstrument	( )								{	return "CEGA";	}
	virtual string			GetSubjectInfo  ( )	const						{	return "CORREL CURVE";}
	static  void			Bump			( ARM_Object *, const ARM_Matrix & , const int &);
	static ARM_Object*		Convert			( ARM_VolCurve*						);
	virtual void			BuildModel		( ARM_Observator* , const string &  ); 

protected:
	virtual ARM_VolCurve*	GetVol			( ARM_Object *							)	const;
	void					SetVol			( ARM_Object *	obj, ARM_VolCurve* vol	);
	virtual string			GetCcy			( ARM_Object *							)	const;
};

			/*	*	*	*	*	*	*	*	*	*	*	*	*	*	*/
			/*														*/
			/*				Template ARM 2D Subject					*/
			/*														*/
			/*	*	*	*	*	*	*	*	*	*	*	*	*	*	*/

template< int n>
void	ARM_2D_Subject<n>::Bump( ARM_Object* vol , const ARM_Matrix & shift, const int & isRelative){
	double		  val;
	ARM_VolCurve *tmp	= dynamic_cast<ARM_VolCurve * > ( vol );

	for ( int i=0; i<shift.GetNumLines(); i++){
		for( int j=0; j< shift.GetNumCols(); j++){
			if ( isRelative ) 
				val = tmp->GetVolatilities()->Elt(i,j)*shift.Elt(i,j);
			else
				val =  shift.Elt(i,j);
				tmp -> BumpVolatility( val,i+1,j+1,K_NO,K_YES);
		}
	}
}

			/*	*	*	*	*	*	*	*	*	*	*	*	*	*	*/
			/*														*/
			/*				Template ARM CAP Subject				*/
			/*														*/
			/*	*	*	*	*	*	*	*	*	*	*	*	*	*	*/


template< int n>
string			ARM_CAP_Subject<n>::GetInstrument( ){

	string	tmp;

	switch(n){
		case	ATM:
			tmp = "VEGA_CAP_ATM";
			break;
		case	RHO:
			tmp = "VEGA_CAP_RHO";
			break;
		case	NU:
			tmp ="VEGA_CAP_NU";
			break;
		case	BETA:
			tmp ="VEGA_CAP_BETA";
			break;
		default:
			throw Exception(__LINE__, __FILE__, ERR_DATE_NOT_VALID, "The SABR volatility surfaces are missing." );
			break;
		}

	return  tmp;
}

template< int n>
string			ARM_CAP_Subject<n>::GetSubjectInfo  ( )	const{

	string	tmp;

	switch(n){
		case	ATM:
			tmp = "VOL CAP ATM";
			break;
		case	RHO:
			tmp = "VOL CAP RHO";
			break;
		case	NU:
			tmp ="VOL CAP NU";
			break;
		case	BETA:
			tmp ="VOL CAP BETA";
			break;
		default:
			throw Exception(__LINE__, __FILE__, ERR_DATE_NOT_VALID, "The SABR volatility surfaces are missing." );
			break;
		}

	return  tmp;
}


			/*	*	*	*	*	*	*	*	*	*	*	*	*	*	*/
			/*														*/
			/*				Template ARM OSW Subject				*/
			/*														*/
			/*	*	*	*	*	*	*	*	*	*	*	*	*	*	*/

template< int n>
string			ARM_OSW_Subject<n>::GetInstrument( ){

	string	tmp;

	switch(n){
		case	ATM:
			tmp = "VEGA_OSW_ATM";
			break;
		case	RHO:
			tmp = "VEGA_OSW_RHO";
			break;
		case	NU:
			tmp ="VEGA_OSW_NU";
			break;
		case	BETA:
			tmp ="VEGA_OSW_BETA";
			break;
		default:
			throw Exception(__LINE__, __FILE__, ERR_DATE_NOT_VALID, "The SABR volatility surfaces are missing." );
			break;
		}

	return  tmp;
}

template< int n>
string			ARM_OSW_Subject<n>::GetSubjectInfo  ( )	const{

	string	tmp;

	switch(n){
		case	ATM:
			tmp = "VOL OSW ATM";
			break;
		case	RHO:
			tmp = "VOL OSW RHO";
			break;
		case	NU:
			tmp ="VOL OSW NU";
			break;
		case	BETA:
			tmp ="VOL OSW BETA";
			break;
		default:
			throw Exception(__LINE__, __FILE__, ERR_DATE_NOT_VALID, "The SABR volatility surfaces are missing." );
			break;
		}

	return  tmp;
}			/*	*	*	*	*	*	*	*	*	*	*	*	*	*	*/
			/*														*/
			/*				Template ARM VFX Subject				*/
			/*														*/
			/*	*	*	*	*	*	*	*	*	*	*	*	*	*	*/


template< int n>
string			ARM_VFX_Subject<n>::GetInstrument( ){

	string	tmp;

	switch(n){

		case	PIV:
			tmp = "VEGA_FX_PIV";
			break;
		case	RR:
			tmp = "VEGA_FX_RR";
			break;
		case	STR:
			tmp ="VEGA_FX_STR";
			break;
		default:
			throw Exception(__LINE__, __FILE__, ERR_DATE_NOT_VALID, "The SABR volatility surfaces are missing." );
			break;
		}

	return  tmp;
}

template< int n>
string		ARM_VFX_Subject<n>::GetSubjectInfo  ( )	const{
	string	tmp;

	switch(n){

		case	PIV:
			tmp	=	"VOL FXBS PIV";
			break;

		case	RR:
			tmp	=	"VOL FXBS RR";
			break;
		
		case	STR:
			tmp	=	"VOL FXBS STR";
			break;

	}
	return tmp;
}

template< int n>
void		ARM_VFX_Subject<n>::Bump( ARM_Object* curve, const ARM_Matrix & shift, const int & isRelative){
		
	ARM_FXVolCurve*	tmp		= dynamic_cast<ARM_FXVolCurve *> (curve);
	ARM_Matrix		tmpM;
	ARM_Vector*		tmpV	= NULL;
	int				row		= shift.GetNumLines();
	int				col		= shift.GetNumCols();

	switch(n){

		case	PIV:{
			tmpV=  const_cast<ARM_Vector* > (tmp-> GetPivotVols());
			for ( int i= 0 ; i< row; i++){
				if ( isRelative )
					tmpV->Elt( i ) *= ( 1+shift.Elt(i,0) );
				else
					tmpV->Elt( i ) += shift.Elt(i,0) ;
				}
			break;
		}
		case	RR:{
			tmpM = tmp -> GetRR();
			for ( int i= 0 ; i < row; i++){
					for( int j = 0; j < col; j++){
					if ( isRelative )
						tmpM.Elt(i,j) *= ( 1+ shift.Elt(i,j) );
					else
						tmpM.Elt(i,j) +=  shift.Elt(i,j) ;
				}
			}
			tmp -> FXBumpRRorSTR( tmpM, K_YES);
			break;
		}
		case	STR:{
			tmpM = tmp -> GetSTR();
			for ( int i= 0 ; i < row; i++){
				for( int j = 0; j < col; j++){
					if ( isRelative )
						tmpM.Elt(i,j) *= ( 1+ shift.Elt(i,j) );
					else
						tmpM.Elt(i,j) +=  shift.Elt(i,j) ;
				}
			}
			tmp -> FXBumpRRorSTR( tmpM, K_NO);
			break;
		}
	}
}

			/*	*	*	*	*	*	*	*	*	*	*	*	*	*	*/
			/*														*/
			/*				Template ARM MIX Subject				*/
			/*														*/
			/*	*	*	*	*	*	*	*	*	*	*	*	*	*	*/


template< int n> 
string		ARM_MIX_Subject<n>::GetInstrument	( ){
	string	tmp;

	switch(n){

		case	VOL:
			tmp =	"VEGA_MIX_VOL";
			break;
		case	SMILE:
			tmp =	"VEGA_MIX_SMILE";
			break;
		case	SHIFT:
			tmp =	"VEGA_MIX_SHIFT";
			break;
		case	Q:
			tmp =	"VEGA_MIX_Q";
			break;
		default:
			throw Exception(__LINE__, __FILE__, ERR_DATE_NOT_VALID, "The Mixture volatility curve are missing." );
			break;
		}

	return  tmp;
}


template< int n>
string		ARM_MIX_Subject<n>::GetSubjectInfo  ( )	const{
	string	tmp;

	switch(n){

		case	VOL:
			tmp	=	"VOL FXMIX VOL";
			break;

		case	SMILE:
			tmp	=	"VOL FXMIX SMILE";
			break;
		
		case	SHIFT:
			tmp	=	"VOL FXMIX SHIFT";
			break;

		case	Q:
			tmp =	"VOL FXMIX Q";
			break;
	}
	return tmp;
}


template< int n>
void	ARM_MIX_Subject<n>::Bump( ARM_Object* vol , const ARM_Matrix & shift, const int & isRelative){
	
	ARM_Matrix*		tmp	= dynamic_cast<ARM_VolCurve * > ( vol )->GetVolatilities();
	ARM_Matrix*		val = new ARM_Matrix ( shift.GetNumLines(), shift.GetNumCols()  );

	for ( int i=0; i<shift.GetNumLines(); i++){
		for( int j=0; j< shift.GetNumCols(); j++){
			if ( isRelative ) 
				val->Elt(i,j) = tmp->Elt(i,j)*( 1 + shift.Elt(i,j) );
			else
				val->Elt(i,j) = tmp->Elt(i,j) + shift.Elt(i,j);
		}
	}
	dynamic_cast<ARM_VolCurve * > ( vol )->SetVolatilities( val );
}

			/*	*	*	*	*	*	*	*	*	*	*	*	*	*	*/
			/*														*/
			/*				Template ARM IBS Subject				*/
			/*														*/
			/*	*	*	*	*	*	*	*	*	*	*	*	*	*	*/


template< int n> 
string		ARM_IBS_Subject<n>::GetInstrument	( ){
	string	tmp;

	switch(n){

		case	CPI:
			tmp =	"VEGA_INF_CPI";
			break;
		case	YOY:
			tmp =	"VEGA_INF_YOY";
			break;
		default:
			throw Exception(__LINE__, __FILE__, ERR_DATE_NOT_VALID,  "The Inf BS volatility curve are missing." );
			break;
		}

	return  tmp;
}


template< int n>
string		ARM_IBS_Subject<n>::GetSubjectInfo  ( )	const{
	string	tmp;

	switch(n){

		case	CPI:
			tmp	=	"VOL INF CPI";
			break;

		case	YOY:
			tmp	=	"VOL INF YOY";
			break;
	}
	return tmp;
}

			/*	*	*	*	*	*	*	*	*	*	*	*	*	*	*/
			/*														*/
			/*				Template ARM SO Subject					*/
			/*														*/
			/*	*	*	*	*	*	*	*	*	*	*	*	*	*	*/


template< int n>
string			ARM_SO_Subject<n>::GetInstrument( ){
	string	tmp;

	switch(n){
		case	ATM:
			tmp = "VEGA_SO_ATM";
			break;
		case	ADJ:
			tmp = "VEGA_SO_ADJ";
			break;
		default:
			throw Exception(__LINE__, __FILE__, ERR_DATE_NOT_VALID, "The SABR volatility surfaces are missing." );
			break;
		}

	return  tmp;
}

template< int n>
string		ARM_SO_Subject<n>::GetSubjectInfo  ( )	const{
	string	tmp;

	switch(n){
		case	ATM:
			tmp	=	"VOL SO ATM";
			break;
		case	ADJ:
			tmp	=	"VOL SO ADJ";
			break;
	}
	return tmp;
}


CC_END_NAMESPACE()

#endif

