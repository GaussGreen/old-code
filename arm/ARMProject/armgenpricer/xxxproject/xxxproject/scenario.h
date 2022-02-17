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

#ifndef INXXXPROJECT_SCENARIO
#define INXXXPROJECT_SCENARIO

#include <gpbase/removeidentifiedwarning.h>
#include <gpbase/port.h>
#include <gpbase/rootobject.h>
#include <gpbase/gpvector.h>
#include <gpbase/gpmatrix.h>
#include <gpbase/typedef.h>
#include <gpbase/assignop.h>
#include <inst/forex.h>
#include<vector>
#include<string>
using namespace std;

class ARM_BasisCurve;
class ARM_ZeroCurve;
class ARM_VolCurve;
class ARM_ZeroLInterpol;
class ARM_BasisCurve;
class ARM_Forex;




CC_BEGIN_NAMESPACE( ARM )

//using ARM::ARM_ArgConvReverse_MktVolType;
///////////////////////////////////////////////////////////////////////////////
//
// Those class define the scenario to calculate the hedge of any security
// of ARM:
// _ ARM_Scenario			: Base class
// _ ARM_ZCCurveScenario	: for IR delta
// _ ARM_BSScenario			: for FX delta
// _ ARM_FXSpotScenario		: for FX delta
// _ ARM_IRVolCurveScenario : for cap & swaption vol ATM, nu, rho, beta, sensis
// _ ARM_FXVolScenario		: for FX volatility
// _ ARM_SOCorrelScenario	: for spread option correlation
//
////////////////////////////////////////////////////////////////////////////////
# define BIGMAMA 999

typedef enum {	NONE,	REF, VALUE,	REF_VALUE }	ScenarioType;

class ARM_MktData;

class ARM_Scenario : public ARM_RootObject
{

public:
	ARM_Scenario				( ):ARM_RootObject( ){ isInitialized = false; }
	ARM_Scenario				(	const double &,					const string &, 
									const int & isRelative		= K_NO , const int & isCumulInv = K_NO,
									const int & isPerturbative	= K_NO);
	ARM_Scenario				(	const ARM_Scenario & rhs );
	ASSIGN_OPERATOR				(	ARM_Scenario )
	~ARM_Scenario				( ){ };
	virtual ARM_Object*	Clone	( )	const	{ return new ARM_Scenario(*this); }
	
	virtual string		ExportShortName	( )	const { return "LSCEN";}
	virtual string		toString		(	const string& indent="", const string& nextIndent="") const;


public:
	virtual void				InitScenario			(	const string &		)			{ }
	virtual void				InitMktData				(	ARM_MktData *		)			{ isInitialized= true;}
	virtual void				ShiftMktData			(	ARM_MktData *		)			{ }
	virtual void				UnShiftMktData			(	ARM_MktData *		)			{ }
	virtual void				Finalize				(	ARM_MktData *		)	const	{ }
	virtual	void				ShiftNextPlot			( ) ;
	virtual void				InitPosition			( ) ;
	
	int							GetNbShift				( )	const	{	return itsNbBump;	}
	int							GetNbEffShift			( ) const;
	void						SetPosition				( ARM_IntVector  vec ) { itsPosition = vec; }

	ARM_IntVector				GetPosition				( )			{	return itsPosition;		}
	ARM_IntVector				GetPermOrder			( )			{	return itsPermOrder;	}
	vector<ARM_StringVector>	GetPlotsPosition		( )			{	return itsPlotsPos;		}
	vector<ARM_StringVector>	GetPlotsLabel			( )			{	return itsPlotsLab;		}
	map<string,ScenarioType>	GetPriceCritera			( )	const	{	return itsPriceCritera;	}	

	virtual string				GetScenarioKey			( )	const	{	return string("_") ;	}
	virtual string				ViewPosition			( )	const	{	return string("");		}
	string						GetCurrency				( )	const	{	return itsCurrency;		}	

	void						OrderData				( ARM_GP_VectorPtr);
	void						OrderData				( ARM_GP_MatrixPtr);
	
	int							GetDim					( )			{ return itsDim;		}

protected:
	static	vector<string>		ReduceString			(	const string & );
	virtual	void				ConvertString			(   const string & ){ }
	virtual void				InitDimScenario			(	const string & , const int & );

protected:
	double						itsShift;
	string						itsCurrency;
	int							itsNbBump;
	ARM_IntVector				itsPermOrder;
	vector<ARM_StringVector>	itsPlotsLab;
	vector<ARM_StringVector>	itsPlotsPos;
	int							isCumulInv; 
	int							isRelative;
	int							isPerturbative;
	int							itsDim; 
	int							isInitialized;

	mutable	ARM_IntVector		itsPosition;
	map<string,ScenarioType>	itsPriceCritera;
};

template <class T=ARM_Subject> 
class ARM_Post_Scenario : public ARM_Scenario
{

public:
	ARM_Post_Scenario			(	const double & shift,					const string & currency,
									const int & isRelative		= K_NO ,	const int & isCumulInv = K_NO,
									const int & isPerturbative	= K_NO): ARM_Scenario(shift, currency, isRelative , isCumulInv, isPerturbative) {}


	virtual void				ShiftMktData			( ARM_MktData *		) ;
	virtual void				UnShiftMktData			( ARM_MktData *		) ;

protected:
	virtual ARM_Matrix			UpDateShift				( bool )=0;
};

template <class T=ARM_Subject> 
class ARM_0D_Scenario : public ARM_Post_Scenario<T>
{

public:
	ARM_0D_Scenario						(	const double & shift,			const string & currency, 
											const int & isRelative = K_NO , const int & isCumulInv = K_NO,
											const int & isPerturbative = K_NO ) : ARM_Post_Scenario<T>(shift, currency, isRelative , isCumulInv, isPerturbative) { itsDim = 0; }
	ASSIGN_OPERATOR						(	ARM_0D_Scenario)
	virtual ARM_Object*			Clone	( )	const	{ return new ARM_0D_Scenario(*this); }

	virtual void				InitScenario			( const string &	) ;
	virtual void				InitMktData				( ARM_MktData *		) ;
	virtual void				Finalize				( ARM_MktData *		)	const;
	virtual string				GetScenarioKey			( )						const;
	virtual string				ViewPosition			( )						const;

protected:
	virtual ARM_Matrix			UpDateShift				( bool );

};

template <class T=ARM_Subject> 
class ARM_1D_Scenario : public ARM_Post_Scenario<T>
{

public:
	ARM_1D_Scenario(	const double & shift,			const string & currency, 
						const int & isRelative = K_NO , const int & isCumulInv = K_NO,
						const int & isPerturbative = K_NO) : ARM_Post_Scenario<T>(shift, currency, isRelative , isCumulInv, isPerturbative) { itsDim = 1; }
	ASSIGN_OPERATOR(ARM_1D_Scenario)
	virtual ARM_Object*			Clone					( )					const	{ return new ARM_1D_Scenario(*this); }	

	virtual void				InitScenario			( const string &);
	virtual void				InitMktData				( ARM_MktData *	);
	virtual void				Finalize				( ARM_MktData *	)	const;
	virtual string				GetScenarioKey			( )					const;
	virtual string				ViewPosition			( )					const;

protected:
	virtual ARM_Matrix			UpDateShift				( bool );
	void						ConvertString			( const string & );

};

template <class T=ARM_Subject> 
class ARM_2D_Scenario : public ARM_Post_Scenario<T>
{

public:
	ARM_2D_Scenario(	const double & shift,			const string & currency, 
						const int & isRelative = K_NO , const int & isCumulInv = K_NO,
						const int & isPerturbative = K_NO) :ARM_Post_Scenario<T>(shift, currency, isRelative , isCumulInv,  isPerturbative) { itsDim = 2 ;  }
	ASSIGN_OPERATOR(ARM_2D_Scenario)
	virtual ARM_Object*			Clone					( )					const	{ return new ARM_2D_Scenario(*this); }	
	
	virtual void				InitScenario			( const string &);
	virtual void				InitMktData				( ARM_MktData *	);
	virtual void				Finalize				( ARM_MktData *	)	const;
	virtual string				GetScenarioKey			( )					const;
	virtual string				ViewPosition			( )					const;

protected:
	virtual ARM_Matrix			UpDateShift				( bool );
	void						ConvertString			( const string & );

};

 
class ARM_ND_Scenario : public ARM_Scenario
{

public:

	ARM_ND_Scenario( vector<ARM_Scenario*> & );

	ASSIGN_OPERATOR(ARM_ND_Scenario)

	ARM_ND_Scenario(	const ARM_ND_Scenario & rhs);

	virtual ARM_Object*			Clone	( )				const	{ return new ARM_ND_Scenario(*this); }
	
	~ARM_ND_Scenario			();
	
	virtual string ExportShortName() const { return "LSCEN";}
	virtual string toString		(	const string& indent="", const string& nextIndent="") const;
	
	virtual void				InitScenario			( const string & scen="");

	virtual void				InitMktData				( ARM_MktData *	);
	virtual void				ShiftMktData			( ARM_MktData *	);
	virtual void				UnShiftMktData			( ARM_MktData *	);
	virtual void				Finalize				( ARM_MktData *	)	const;
	virtual	void				InitPosition			( );
	
	virtual string				GetScenarioKey			( )					const;
	vector<string>				GetSubScenario			( )					const;
		
protected:

	void						NotifySubScenario		( );
	void						ConvertString			(   const string & );

protected:
	vector<ARM_Scenario *>		itsScenario;

};


			/*	*	*	*	*	*	*	*	*	*	*	*	*	*	*/
			/*														*/
			/*					Post SCENARIO						*/
			/*														*/
			/*	*	*	*	*	*	*	*	*	*	*	*	*	*	*/

template <class T> 
void	ARM_Post_Scenario<T>::ShiftMktData( ARM_MktData * mkt ) {

	ARM_MktSubject< T > tmp	( * mkt );
	ARM_Object*			model	= tmp.Get( itsCurrency );

	T::Bump(model, UpDateShift( true ), isRelative);
	tmp.Set( itsCurrency, model);
}

template <class T> 
void	ARM_Post_Scenario<T>::UnShiftMktData( ARM_MktData * mkt )	{

	ARM_MktSubject< T > tmp	( * mkt );
	ARM_Object*			model	= tmp.Get( itsCurrency );

	T::Bump(model, UpDateShift( false ), isRelative);
	tmp.Set( itsCurrency, model);

}


			/*	*	*	*	*	*	*	*	*	*	*	*	*	*	*/
			/*														*/
			/*						0D SCENARIO						*/
			/*														*/
			/*	*	*	*	*	*	*	*	*	*	*	*	*	*	*/

template <class T> 
string	ARM_0D_Scenario<T>::ViewPosition()const{
	CC_Ostringstream	os;

	os << CC_NS(std,setfill(' ')) << CC_NS(std,fixed) << CC_NS(std,setw)(10)<<itsPlotsLab[0][0] << " : ";
	os << CC_NS(std,setfill(' ')) << CC_NS(std,fixed) << CC_NS(std,setw)(10)<<"("<<itsPlotsLab[0][0]<<")"<< " - "<<"(NO)";

	return os.str( );
}


template <class T> 
void	ARM_0D_Scenario<T>::InitScenario(const string& theScenario ) {	InitDimScenario	(theScenario,1);	}

template <class T> 
void ARM_0D_Scenario<T>::InitMktData( ARM_MktData * mkt){

	ARM_Scenario::InitMktData(mkt);

	itsPlotsLab[0].resize(1);
	itsPlotsPos[0].resize(1);
	itsPlotsLab[0][0]	=	string("Spot");
	itsPlotsPos[0][0]	=	string("Spot");

	itsNbBump = 1;
	itsPriceCritera.insert(pair<string, ScenarioType> ( "Spot", VALUE) );
}


template <class T> 
void	ARM_0D_Scenario<T>::Finalize	(ARM_MktData * mkt ) const	{	
	
	ARM_MktSubject< T > tmp( *mkt );
	tmp.Finalize(itsCurrency) ;
}

template <class T> 
string	ARM_0D_Scenario<T>::GetScenarioKey	()	const	{ 

	return (T::GetInstrument() + string("_") +itsCurrency);	

}

template <class T> 
ARM_Matrix	ARM_0D_Scenario<T>::UpDateShift( bool isUp ) {

	double				val	= isUp*itsShift-(1-isUp)*itsShift ;
	if (isRelative)		val = isUp*itsShift-(1-isUp)*itsShift/(1+itsShift) ;
	ARM_Matrix			vec ( 1 , 1, val);

	return vec;
}

			/*	*	*	*	*	*	*	*	*	*	*	*	*	*	*/
			/*														*/
			/*						1D SCENARIO						*/
			/*														*/
			/*	*	*	*	*	*	*	*	*	*	*	*	*	*	*/

template <class T> 
string	ARM_1D_Scenario<T>::ViewPosition()const {

	string tmp("NO");
	CC_Ostringstream	os;
	map<string,ScenarioType>::const_iterator it;

	for( int i=0; i< itsNbBump; i++){
		it =itsPriceCritera.find(itsPlotsLab[0][i]);
		switch( it->second ){

			case VALUE:{
				os << CC_NS(std,setfill(' ')) << CC_NS(std,fixed) << CC_NS(std,setw)(10)<<itsPlotsLab[0][i] << " : ";
				os << CC_NS(std,setfill(' ')) << CC_NS(std,fixed) << CC_NS(std,setw)(10)<<"("<<itsPlotsLab[0][i]<<")"<< " - "<<"("<<tmp<<")";
				os << "\n";
				break;
			}
			case REF_VALUE:{
				os << CC_NS(std,setfill(' ')) << CC_NS(std,fixed) << CC_NS(std,setw)(10)<<itsPlotsLab[0][i] << " : ";
				os << CC_NS(std,setfill(' ')) << CC_NS(std,fixed) << CC_NS(std,setw)(10)<<"("<<itsPlotsLab[0][i]<<")"<< " - "<<"("<<tmp<<")";
				os << "\n";
				tmp = itsPlotsLab[0][i];
				break;
			}
			case REF:{
				tmp =itsPlotsLab[0][i];
				break;
			}
		}
	}
	return os.str( );
}

template <class T> 
void	ARM_1D_Scenario<T>::InitScenario(const string& theScenario ) {	InitDimScenario	(theScenario,1);	}

template <class T> 
void ARM_1D_Scenario<T>::InitMktData( ARM_MktData * mkt){

	int					NbPlots;
	string				tmp;

	ARM_Scenario::InitMktData(mkt);

	ARM_MktSubject< T > subject( *mkt );
	ARM_ZeroCurve* model = dynamic_cast<ARM_ZeroCurve* > ( subject.Get(itsCurrency) );

	NbPlots					=  model->GetMktData()->itsMktValue->size();
	if (NbPlots==0)
		throw Exception(__LINE__, __FILE__, ERR_DATE_NOT_VALID, "The ZC Curve does not contain any datas." );

	itsPlotsLab[0].resize(NbPlots);


	for(int j=0;j< NbPlots;j++){
		tmp						=	tmp.insert(0,(model->GetMktData()->itsMktTerms)[j]);
		itsPlotsLab[0][j]		=	tmp;
		tmp.erase();
	}

	if ( itsPlotsPos[0][0]=="" ){
		itsPlotsPos[0].resize(NbPlots);
		for(int j=0;j< NbPlots;j++)
			itsPlotsPos[0][j] = itsPlotsLab[0][j];
	}
	else if ( (itsPlotsPos[0].size() == 1) && (itsPlotsPos[0][0] == "ALL") )
		itsPlotsPos[0][0] = itsPlotsLab[0][itsPlotsLab[0].size()-1];
	
	itsNbBump = NbPlots;

	for( int j=0; j< NbPlots; j++){
		It_Str its		=	find( itsPlotsPos[0].begin(), itsPlotsPos[0].end(), itsPlotsLab[0][j]);
		if ( its !=	itsPlotsPos[0].end() ){
			if (j ==0)
				itsPriceCritera.insert(pair<string,ScenarioType>( itsPlotsLab[0][j],REF));
			else
				itsPriceCritera.insert(pair<string,ScenarioType>( itsPlotsLab[0][j],VALUE));
			if ( j>0 && isPerturbative == K_YES ){
				if ( itsPriceCritera[itsPlotsLab[0][j-1]] == VALUE )
					itsPriceCritera[itsPlotsLab[0][j-1] ]= REF_VALUE;
				else
					itsPriceCritera[itsPlotsLab[0][j-1] ]= REF;
				}
			}
		else
			itsPriceCritera.insert(pair<string,ScenarioType>( itsPlotsLab[0][j],NONE) );
	}
}

template <class T> 
ARM_Matrix	ARM_1D_Scenario<T>::UpDateShift( bool isUp ) {

	int					dim = itsPlotsLab[0].size();
	ARM_Matrix			vec ( dim, 1, 0.0 );
	double				tmp = isUp*itsShift-(1-isUp)*itsShift;

	if (isRelative)		tmp = isUp*itsShift-(1-isUp)*itsShift/(1+itsShift);
	
	for ( int i= 0; i<itsPosition[0]+1; i++){
		if (isCumulInv) 
			vec.Elt(dim-i-1,0)	= tmp;
		else
			vec.Elt(i,0)		= tmp;
	}
	return vec;
}


template <class T> 
void	ARM_1D_Scenario<T>::Finalize	(ARM_MktData * mkt ) const	{	
	
	ARM_MktSubject< T > tmp( *mkt );
	tmp.Finalize(itsCurrency) ;
}

template <class T> 
string	ARM_1D_Scenario<T>::GetScenarioKey	()	const	{ 

	return (T::GetInstrument() + string("_") +itsCurrency);	

}

template <class T> 
void	ARM_1D_Scenario<T>::ConvertString( const string & str ){
	
	int							pos;
	string						tmp;
	vector< string >			tmpV_Str;	
	
	tmpV_Str.clear();
	
	if ( str == "" )
		tmpV_Str.push_back("");

	else {

		pos = str.find(";");
		if ( pos>=0 && pos<= str.size() )
			throw Exception(__LINE__, __FILE__, ERR_DATE_NOT_VALID, "the scenario should not contain any ; since we are dealing with a one dimension scenario" );
		
		pos = str.find("T");
		if ( pos<0  || pos > str.size() )
			throw Exception(__LINE__, __FILE__, ERR_DATE_NOT_VALID, "we are expecting T ; since we are dealing with a one dimension scenario" );

		tmpV_Str = ARM_Scenario::ReduceString(str);

	}
	
	itsPlotsPos	[0]	=	tmpV_Str;
	itsPermOrder[0]	=	0;
}


			/*	*	*	*	*	*	*	*	*	*	*	*	*	*	*/
			/*														*/
			/*						2D SCENARIO						*/
			/*														*/
			/*	*	*	*	*	*	*	*	*	*	*	*	*	*	*/

template <class T> 
string	ARM_2D_Scenario<T>::ViewPosition()const{

	string					tmp("NO");
	string					lab;
	bool					goNext;
	CC_Ostringstream		os;

	map<string,ScenarioType>::const_iterator	it;
	vector<string>			::const_iterator	its0;
	vector<string>			::const_iterator	its1;

	ARM_Vector				tmpVec(2);
	ARM_GP_T_Matrix<string> mat(itsPlotsLab[0].size(), itsPlotsLab[1].size());
	ARM_GP_T_Matrix<string> tmpMat(itsPlotsPos[0].size()+1, itsPlotsPos[1].size()+1);

	itsPosition[0]=0;
	itsPosition[1]=0;

	for(int i= 0; i<itsNbBump; i++){

		it	=	itsPriceCritera.find(itsPlotsLab[0][itsPosition[0]]+itsPlotsLab[1][itsPosition[1]]);
		lab =	string("(")+itsPlotsLab[0][itsPosition[0]]+string(",")+itsPlotsLab[1][itsPosition[1]]+string(")");
		switch( it->second ){
			case VALUE:{
				mat.Elt(itsPosition[0],itsPosition[1]) 	= lab+string("-")+tmp;
				break;
			}
			case REF_VALUE:{
				mat.Elt(itsPosition[0],itsPosition[1]) 	= lab+string("-")+tmp;
				tmp = lab;
				break;
			}
			case REF:{
				mat.Elt(itsPosition[0],itsPosition[1])="";
				break;
			}
		}
		const_cast<ARM_2D_Scenario<T>* > (this)->ShiftNextPlot();
	}
	itsPosition[0]=0;
	itsPosition[1]=0;

	tmpMat.Elt(0,0) = string("");
    int i;
	
	for(i=0; i< itsPlotsPos[0].size(); i++)
		tmpMat.Elt(i+1,0) = itsPlotsPos[0][i];

	for(i=0; i< itsPlotsPos[1].size(); i++)
		tmpMat.Elt(0,i+1) = itsPlotsPos[1][i];

	int x=0;
	int y=0;
	for(  i=0; i< itsPlotsLab[0].size(); i++){
		its0 = find( itsPlotsPos[0].begin(), itsPlotsPos[0].end(),itsPlotsLab[0][i] );
		goNext= false;
		for(int  j=0; j< itsPlotsLab[1].size(); j++){
			its1 = find( itsPlotsPos[1].begin(), itsPlotsPos[1].end(),itsPlotsLab[1][j] );

			if( its0!=itsPlotsPos[0].end() && its1!=itsPlotsPos[1].end()){
				tmpMat.Elt(x+1,y+1)=mat.Elt(i,j);
				y++;
				goNext = true;
			}
		}
		if (goNext == true) 	x++;
		y=0;
	}

	for(  i=0; i< itsPlotsPos[0].size()+1; i++){
		for(int  j=0; j< itsPlotsPos[1].size()+1; j++)
			os << CC_NS(std,setfill(' '))<< CC_NS(std,fixed)<< CC_NS(std,setw)(20)<<tmpMat.Elt(i,j);
		os <<"\n";
	}
	return os.str( );
}


template <class T> 
void	ARM_2D_Scenario<T>::InitScenario(const string& theScenario){	InitDimScenario	(theScenario,2);	}

template <class T> 
void ARM_2D_Scenario<T>::InitMktData( ARM_MktData * mkt){

	vector<int>			NbPlots(2);
	vector<It_Str>		its(2);
	string				tmp;
	string				tmpPos;
	string				tmpRef("");



	ARM_Scenario::InitMktData(mkt);
	ARM_MktSubject< T > subject( *mkt );
	ARM_VolCurve* model = dynamic_cast<ARM_VolCurve* > ( subject.Get(itsCurrency) );
	
	itsPlotsLab[0]		= subject.GetLabel( itsCurrency, "X");
	itsPlotsLab[1]		= subject.GetLabel( itsCurrency, "Y");

	NbPlots[0]			= itsPlotsLab[0].size();
	NbPlots[1]			= itsPlotsLab[1].size();
	itsNbBump			= NbPlots[0] *	NbPlots[1];
	
	for (int j=0;j<2;j++){
		if ( NbPlots[j]	==	0 )
			throw Exception(__LINE__, __FILE__, ERR_DATE_NOT_VALID, "The Vol Curve does not contain any datas." );
			
		if ( itsPlotsPos[j].size()==0 ){
			itsPlotsPos[j].resize(	NbPlots[j]	);

			for(int i=0; i< NbPlots[j];	i++)
				itsPlotsPos[j][i] = itsPlotsLab[j][i];
		}
		else if ( (itsPlotsPos[j].size() == 1) && (itsPlotsPos[j][0] == "ALL") )
			itsPlotsPos[j][0]=itsPlotsLab[j][itsPlotsLab[j].size()-1];
	}

	InitPosition();
	ARM_IntVector refPosition(itsPosition);

	for(int i=0; i<itsNbBump; i++){
		ShiftNextPlot();
		for(int j=0; j<2; j++)
			its[j]	=	find( itsPlotsPos[j].begin(), itsPlotsPos[j].end(),	itsPlotsLab[j][itsPosition[j]]);

		tmpPos	= itsPlotsLab[0][itsPosition[0]]+itsPlotsLab[1][itsPosition[1]];
		if ( its[0]	!=	itsPlotsPos[0].end() &&  its[1] !=	itsPlotsPos[1].end() ){

			if (i ==0)
				itsPriceCritera.insert(pair<string,ScenarioType>( tmpPos, REF));
			else
				itsPriceCritera.insert(pair<string,ScenarioType>( tmpPos, VALUE));

			if ( tmpRef.size() !=0 ){
				if( itsPriceCritera[tmpRef] == VALUE )
					itsPriceCritera[tmpRef] = REF_VALUE;
				else
					itsPriceCritera[tmpRef] = REF;
			}
			tmpRef	= tmpPos;

			if (isPerturbative == K_YES && i>0){
				tmp= itsPlotsLab[0][refPosition[0]]+itsPlotsLab[1][refPosition[1]];
				if( itsPriceCritera[tmp] == VALUE )
					itsPriceCritera[tmp] = REF_VALUE;
				else
					itsPriceCritera[tmp] = REF;
			}
		}
		else
			itsPriceCritera.insert(pair<string,ScenarioType>( tmpPos, NONE));

		refPosition=itsPosition;
	}
	InitPosition();
}


template <class T> 
ARM_Matrix	ARM_2D_Scenario<T>::UpDateShift( bool isUp ) {
	
	int			i, j, nb; 
	int			nbAbs, nbOrd;
	int			dimAbs, dimOrd;
	int			row, col;
	double		val = isUp*itsShift-(1-isUp)*itsShift;

	if (isRelative)		val= isUp*itsShift-(1-isUp)*itsShift/(1+itsShift);

	dimAbs		= itsPlotsLab[0].size();
	dimOrd		= itsPlotsLab[1].size();
	
	ARM_Matrix	mat( dimAbs, dimOrd, 0.0) ;

	if ( itsPermOrder[0] == 0 )
		nb		= dimAbs;
	else
		nb		= dimOrd;

	nbAbs		= itsPosition[0]*itsPermOrder[0]+itsPosition[1]*itsPermOrder[1];
	nbOrd		= itsPosition[1]*itsPermOrder[0]+itsPosition[0]*itsPermOrder[1];

	for (	i = 0	; i < nbAbs;	i++){
		for(	j	= 0; j <nb;		j++){
			row		=	(i+1)*itsPermOrder[0]+(j+1)*itsPermOrder[1];
			col		=	(j+1)*itsPermOrder[0]+(i+1)*itsPermOrder[1];
			if ( isCumulInv)
				mat.Elt(dimAbs - row, dimOrd - col) = val;	
			else
				mat.Elt(row-1, col-1) = val;
		}
	}
	for(	j	= 0; j <= nbOrd;	j++)	{
		row		=	(nbAbs+1)*itsPermOrder[0]+(j+1)*itsPermOrder[1];
		col		=	(j+1)*itsPermOrder[0]+(nbAbs+1)*itsPermOrder[1];
		if ( isCumulInv)
			mat.Elt(dimAbs - row, dimOrd - col) = val;	
		else
			mat.Elt(row-1, col-1) = val;	
	}
	return mat;
}

template <class T> 
void	ARM_2D_Scenario<T>::Finalize	(ARM_MktData * mkt ) const	{	
	
	ARM_MktSubject< T > tmp( *mkt );
	tmp.Finalize(itsCurrency) ;
}

template <class T> 
string	ARM_2D_Scenario<T>::GetScenarioKey	()	const	{ 

	return (T::GetInstrument() + string("_") +itsCurrency);	
}

template <class T> 
void	ARM_2D_Scenario<T>::ConvertString( const string & str ){
	
	int							pos;
	int							pos1;
	int							pos2;
	int							test;

	string						tmp;
	string						tmp1;
	string						tmp2;

	vector< string >			tmpV_Str1;	
	vector< string >			tmpV_Str2;
	
	tmpV_Str1.clear();
	tmpV_Str2.clear();
	itsPermOrder[0]	=	0;
	itsPermOrder[1]	=	1;
	
	if ( str == "" ){
		tmpV_Str1.resize(0);
		tmpV_Str2.resize(0);
	}

	else{

		pos = str.find(";");
		if ( pos < 0 || pos > str.size() ){

			test = str.find("E");
			if ( test>=0   && test <= str.size() ){
				tmpV_Str1 = ARM_Scenario::ReduceString(str);
				tmpV_Str2.resize(0);
				}

			test = str.find("T");
			if ( test>=0   && test <= str.size() ){
				tmpV_Str1 .resize(0);
				tmpV_Str2 = ARM_Scenario::ReduceString(str);
				itsPermOrder[0]	=	1;
				itsPermOrder[1]	=	0;
				}
			}
		else {

			pos1	= str.find("T");
			pos2	= str.find("E");
			test	= 1;
			if ( pos1 < 0 || pos1 > str.size() ) test	*= 0;
			if ( pos2 < 0 || pos2 > str.size() ) test	*= 0;
			if ( test ==0 )
					throw Exception(__LINE__, __FILE__, ERR_DATE_NOT_VALID, "E  and T are both mising" );
			else{
				tmp1.erase();
				tmp1.insert(0,str.begin(), str.begin()+pos );
				tmp2.erase();
				tmp2.insert(0,str.begin()+pos+1, str.end() );

				pos1 = tmp1.find("E");
				if ( pos1 >=0	&&	pos1 <= tmp1.size() ){
					tmpV_Str1 = ARM_Scenario::ReduceString(tmp1);
					tmpV_Str2 = ARM_Scenario::ReduceString(tmp2);
					}
				else {
					tmpV_Str1 = ARM_Scenario::ReduceString(tmp2);
					tmpV_Str2 = ARM_Scenario::ReduceString(tmp1);
					itsPermOrder[0]	=	1;
					itsPermOrder[1]	=	0;
					}
				}
			itsPlotsPos[0] = tmpV_Str1;
			itsPlotsPos[1] = tmpV_Str2;
			}
		}
}




CC_END_NAMESPACE()

#endif