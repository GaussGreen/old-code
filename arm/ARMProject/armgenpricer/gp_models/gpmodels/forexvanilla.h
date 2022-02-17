/*!
 *
 * Copyright (c) CDC IXIS CM November 2003 Paris
 * 
 * \file forexvanilla.h
 *  \brief to factorize all vanilla/semi-vanilla fx payoffs
 * 
 *	\author  K. Belkheir
 *	\version 1.0
 *	\date January 2007
 */

#ifndef _INGPCALCULATORS_FOREXVANILLA_H
#define _INGPCALCULATORS_FOREXVANILLA_H

#include "gpbase/gpmatrix.h"
#include "gpbase/assignop.h"
#include "gpbase/rootobject.h"
#include "gpbase/typedef.h"
#include "gpbase/cloneutilityfunc.h"
#include "gpbase/env.h"

#include <math.h>

const long  K_NbLGPoint   =   121;
const int	K_Call		  =     1;

CC_BEGIN_NAMESPACE( ARM )

/// ------------------------------
///	--- Forex Base Class ----
/// ------------------------------
class ARM_FXName : public ARM_RootObject
{
	static const string MktNamesTable [];
	static const string InvMktNamesTable [];
protected:
	bool itsIsInvMkt;
	string itsMktName;
	ARM_GP_StrVector itsMktNamesVect;
	ARM_GP_StrVector itsInvMktNamesVect;
public: 
	/// standard constructor
	ARM_FXName(const string& mktName);
	//copy constructor
	ARM_FXName( const ARM_FXName& rhs );
	ASSIGN_OPERATOR(ARM_FXName)
	//destructor
	virtual ~ARM_FXName() {}
	
	/// what it is for ...
	virtual bool IsInvMkt( const string& fxName );

	//accesors
	virtual bool GetIsInvMkt() const {return itsIsInvMkt;};
	virtual string GetMktName() const {return itsMktName;};

    /// Standard ARM object support
	virtual ARM_Object* Clone()  const { return NULL; };
	virtual string toString(const string& indent="",const string& nextIndent="") const { return string("Forex name ");};
	virtual string ExportShortName() const { return "LFX";}
};


/// ------------------------------
///	--- FX Vanilla Base Class ----
/// ------------------------------
class ARM_FXVanilla : public ARM_RootObject
{
public:

	enum QuantoType
	{
		Quanto1,			/// When PAY/COM has the same model than the first FX  
	    Quanto2,			/// When PAY/COM has the same model than the second FX 
		None,				/// When there is non quanto
	};
protected:
	string itsMktFXname;
	bool itsIsInvG;
	int itsCallPut;
	double itsStrike;
public: 
	/// standard constructor
	ARM_FXVanilla( const string& mktFXname, 
		int callput = K_Call,
		bool isInvG = false)
		:
	//itsFXname(fxName),
	itsMktFXname(mktFXname),
	itsIsInvG(isInvG),
	itsCallPut(callput),
	itsStrike(0){}

	//copy constructor
	ARM_FXVanilla( const ARM_FXVanilla& rhs )
		:
	ARM_RootObject(rhs),
	itsMktFXname(rhs.itsMktFXname),
	itsIsInvG(rhs.itsIsInvG),
	itsCallPut(rhs.itsCallPut),
	itsStrike(rhs.itsStrike){}

	ARM_FXVanilla& operator = ( const ARM_FXVanilla& rhs )
	{
		if( this != &rhs ){
			ARM_RootObject::operator=( rhs );
			itsMktFXname = rhs.itsMktFXname;
			itsIsInvG = rhs.itsIsInvG;
			itsCallPut = rhs.itsCallPut;
			itsStrike = rhs.itsStrike;
		}
		return *this;
	}
	//destructor
	virtual ~ARM_FXVanilla() {}

	/// what it is for ...
	
	virtual double Payoff(double x ) const  = 0;
	virtual double Payoff(double x, double y) const = 0;

	// accessors

	virtual string GetMktFXname() const {return itsMktFXname;};
	virtual void SetMktFXname(  string mktFXname ) {itsMktFXname=mktFXname;};

	virtual bool GetIsInvG() const {return itsIsInvG;};
	virtual void SetIsInvG(  bool isInvG ) {itsIsInvG=isInvG;};

	virtual int GetCallPut() const {return itsCallPut;};
	virtual void SetCallPut(  int callPut ) { itsCallPut=callPut; };

	virtual double GetStrike() const {return itsStrike;};
	virtual void SetStrike(  double strike ) { itsStrike=strike; };

	// pure virtual accessors 

	virtual void SetCoefficients ( const ARM_GP_Vector&  coefficients) {};

	virtual string GetMktFX2name() const { return NULL; };
	virtual void SetMktFX2name(  string mktFx2name ) {};

	virtual string GetComCcy() const { return NULL; };
	virtual void SetComCcy(  string comCcy ) {};

	virtual bool GetIsInvG2() const { return false; };
	virtual void SetIsInvG2(  bool isInvG2 ) {};

	virtual QuantoType GetQuantoType() const { return None; };
	virtual void SetQuantoType(  QuantoType quantoType ) {};


    /// Standard ARM object support
	virtual ARM_Object* Clone()  const { return NULL; };
	virtual string toString(const string& indent="",const string& nextIndent="") const { return string("Forex Vanilla");};
	virtual string ExportShortName() const { return "LFXVA";}
};


/// ------------------------------
///	--- 2FX Vanilla Base Class ----
/// ------------------------------
class ARM_FXVanilla2D : public ARM_FXVanilla
{
protected:
	//ARM_FXName itsFXname2;
	string itsMktFx2name;
	string itsComCcy;
	bool itsIsInvG2;
	QuantoType  itsQuantoType;

public: 
	enum FXVanilla2D
	{
		Fx1_Fx2=0,			/// When we can find a proba where FX1 and FX2 are martingal
	    InvFx1_InvFx2,		/// When we can find a proba where 1/FX1 and 1/FX2 are martingal
		Fx1_InvFx2,			/// When we can find a proba where FX1 and 1/FX2 are martingal
		InvFx1_Fx2,			/// When we can find a proba where 1/FX1 and FX2 are martingal
	};
	/// standard constructor
	ARM_FXVanilla2D( const string& mktFXname1, 
		const string& mktFXname2, 
		const string& ccyComName,
		QuantoType  quantoType = None, 
		int callput = K_Call,
		bool isInvG1 = false,
		bool isInvG2 = false)
	:
ARM_FXVanilla(mktFXname1,callput,isInvG1),
		  itsMktFx2name(mktFXname2),
		  itsComCcy(ccyComName),
		  itsIsInvG2(isInvG2),
		  itsQuantoType(quantoType){}

	//copy constructor
	ARM_FXVanilla2D( const ARM_FXVanilla2D& rhs )
		: ARM_FXVanilla(rhs),
		//itsFXname2(rhs.itsFXname2),
		itsMktFx2name(rhs.itsMktFx2name),
		itsComCcy(rhs.itsComCcy),
		itsIsInvG2(rhs.itsIsInvG2),
		itsQuantoType(rhs.itsQuantoType){}
	ASSIGN_OPERATOR(ARM_FXVanilla2D)

	//destructor
	virtual ~ARM_FXVanilla2D() {}

/// what it is for ...
	virtual double Payoff( double x) const {return 0.0;};
	virtual double Payoff(double x, double y) const {return 0.0;};

	// accessors

	virtual string GetMktFX2name() const {return itsMktFx2name;};
	virtual void SetMktFX2name(  string mktFx2name ) {itsMktFx2name=mktFx2name;};

	virtual string GetComCcy() const {return itsComCcy;};
	virtual void SetComCcy(  string comCcy ) {itsComCcy=comCcy;};

	virtual bool GetIsInvG2() const {return itsIsInvG2;};
	virtual void SetIsInvG2(  bool isInvG2 ) {itsIsInvG2=isInvG2;};

	virtual QuantoType GetQuantoType() const {return itsQuantoType;};
	virtual void SetQuantoType(  QuantoType quantoType ) {itsQuantoType=quantoType;};


    /// Standard ARM object support
	virtual ARM_Object* Clone() const { return new ARM_FXVanilla2D(*this); };
};


/// -------------------------
///	--- ARM_FXCall Payoff --- 
/// -------------------------
class ARM_FXCall : public ARM_FXVanilla2D
{
public: 
	/// standard constructor
	ARM_FXCall( const string& mktFXname1, const string& mktFXname2,
		const string& ccyComName,
		QuantoType  quantoType = None, 
		int callput = K_Call,
		bool isInvG1 = false,
		bool isInvG2 = false)
	:	ARM_FXVanilla2D ( mktFXname1, mktFXname2,  
					  ccyComName,quantoType,
					  callput, isInvG1,isInvG2){}
	
	//copy constructor
	ARM_FXCall( const ARM_FXCall& rhs )
		: ARM_FXVanilla2D(rhs){} 
	ASSIGN_OPERATOR(ARM_FXCall)
	
	//destructor
	virtual ~ARM_FXCall(){}

	/// what it is for ...
	virtual double Payoff(double x ) const;	
	virtual double Payoff(double x, double y) const; 

    /// Standard ARM object support
	virtual ARM_Object* Clone() const { return new ARM_FXCall(*this); };
};


/// ------------------------------
///	--- ARM_FXCallInvFx Payoff --- 
/// ------------------------------
class ARM_FXCallInvFx : public ARM_FXCall
{
public: 
	/// standard constructor
		ARM_FXCallInvFx( const string& mktFXname1, const string& mktFXname2,
		const string& ccyComName,
		QuantoType  quantoType = None, 
		int callput = K_Call,
		bool isInvG1 = false,
		bool isInvG2 = false)
	:	ARM_FXCall( mktFXname1, mktFXname2,  
					  ccyComName,quantoType,
					  callput, isInvG1,isInvG2){}

	//copy constructor
	ARM_FXCallInvFx( const ARM_FXCallInvFx& rhs )
		: ARM_FXCall(rhs){}
	ASSIGN_OPERATOR(ARM_FXCallInvFx)

	//destructor
	virtual ~ARM_FXCallInvFx() {}

/// what it is for ...
	virtual double Payoff( double x) const;	
	virtual double Payoff(double x, double y) const; 

    /// Standard ARM object support
	virtual ARM_Object* Clone() const { return new ARM_FXCallInvFx(*this); };
};


/// ---------------------------------
///	--- ARM_FXSpread Payoff ---
/// ---------------------------------
class ARM_FXSpread : public ARM_FXVanilla2D
{
protected:
	ARM_GP_Vector  itsCoefficients;
public: 
	/// standard constructor
	ARM_FXSpread( const string& mktFXname1, const string& mktFXname2,
		const string& ccyComName,
		QuantoType  quantoType = None, 
		int callput = K_Call,
		bool isInvG1 = false,
		bool isInvG2 = false)
	:	ARM_FXVanilla2D ( mktFXname1, mktFXname2,  
					  ccyComName,quantoType,
					  callput, isInvG1,isInvG2),
		itsCoefficients(2,1.0){}
	
	//copy constructor
	ARM_FXSpread( const ARM_FXSpread& rhs )
		: ARM_FXVanilla2D(rhs), 
		  itsCoefficients(rhs.itsCoefficients){}

	ASSIGN_OPERATOR(ARM_FXSpread)
	
	//destructor
	virtual ~ARM_FXSpread() {}

	/// what it is for ...
	
	virtual double Payoff(double x) const {return 0.0;};
	virtual double Payoff(double x, double y) const;
	
	// accessors
	virtual void SetCoefficients (const ARM_GP_Vector&  coefficients) {itsCoefficients = coefficients;};

    /// Standard ARM object support
	virtual ARM_Object* Clone() const const { return new ARM_FXSpread(*this); };
};


/// ------------------------------------
///	--- ARM_FXSpreadInvFx1Fx2 Payoff --- 
/// ------------------------------------
class ARM_FXSpreadInvFx1Fx2 : public ARM_FXSpread 
{
public: 
	/// standard constructor
	ARM_FXSpreadInvFx1Fx2(const string& mktFXname1, const string& mktFXname2,
		const string& ccyComName,
		QuantoType  quantoType = None, 
		int callput = K_Call,
		bool isInvG1 = false,
		bool isInvG2 = false)
		: ARM_FXSpread( mktFXname1, mktFXname2, ccyComName, quantoType,	callput, isInvG1, isInvG2){}

	//copy constructor
	ARM_FXSpreadInvFx1Fx2( const ARM_FXSpreadInvFx1Fx2& rhs )
		: ARM_FXSpread(rhs){}
	ASSIGN_OPERATOR(ARM_FXSpreadInvFx1Fx2)

	//destructor
	virtual ~ARM_FXSpreadInvFx1Fx2() {}

/// what it is for ...
	//virtual double Payoff(double x) const  {return 0.0;}
	virtual double Payoff(double x, double y) const;	

    /// Standard ARM object support
	virtual ARM_Object* Clone() const { return new ARM_FXSpreadInvFx1Fx2(*this); };
};

/// ------------------------------------
///	--- ARM_FXSpreadFx1InvFx2 Payoff ---
/// ------------------------------------
class ARM_FXSpreadFx1InvFx2 : public ARM_FXSpread
{
public: 
	/// standard constructor
	ARM_FXSpreadFx1InvFx2(const string& mktFXname1, const string& mktFXname2,
		const string& ccyComName,
		QuantoType  quantoType = None, 
		int callput = K_Call,
		bool isInvG1 = false,
		bool isInvG2 = false)
		: ARM_FXSpread( mktFXname1, mktFXname2, ccyComName, quantoType,	callput, isInvG1, isInvG2){}

	//copy constructor
	ARM_FXSpreadFx1InvFx2( const ARM_FXSpreadFx1InvFx2& rhs )
		: ARM_FXSpread(rhs){}
	ASSIGN_OPERATOR(ARM_FXSpreadFx1InvFx2)

	//destructor
	virtual ~ARM_FXSpreadFx1InvFx2() {}

/// what it is for ...
	//virtual double Payoff(double x) const  {return 0.0;}
	virtual double Payoff(double x, double y) const;	

    /// Standard ARM object support
	virtual ARM_Object* Clone() const { return new ARM_FXSpreadFx1InvFx2(*this); };
};

/// ---------------------------------------
///	--- ARM_FXSpreadInvFx1InvFx2 Payoff ---
/// ---------------------------------------
class ARM_FXSpreadInvFx1InvFx2 : public ARM_FXSpread
{
public: 
	/// standard constructor
	ARM_FXSpreadInvFx1InvFx2(const string& mktFXname1, const string& mktFXname2,
		const string& ccyComName,
		QuantoType  quantoType = None, 
		int callput = K_Call,
		bool isInvG1 = false,
		bool isInvG2 = false)
		: ARM_FXSpread( mktFXname1, mktFXname2, ccyComName, quantoType,	callput, isInvG1, isInvG2){}

	//copy constructor
	ARM_FXSpreadInvFx1InvFx2( const ARM_FXSpreadInvFx1InvFx2& rhs )
		: ARM_FXSpread(rhs){}
	ASSIGN_OPERATOR(ARM_FXSpreadInvFx1InvFx2)

	//destructor
	virtual ~ARM_FXSpreadInvFx1InvFx2() {}

/// what it is for ...
	//virtual double Payoff(double x) const  {return 0.0;}
	virtual double Payoff(double x, double y) const;	

    /// Standard ARM object support
	virtual ARM_Object* Clone() const { return new ARM_FXSpreadInvFx1InvFx2(*this); };
};


/// ---------------------------------
///	--- ARM_FXBasket Payoff ---
/// ---------------------------------
class ARM_FXBasket : public ARM_FXVanilla2D
{
protected:
	ARM_GP_Vector  itsCoefficients;
public: 
	/// standard constructor
	ARM_FXBasket( const string& mktFXname1, const string& mktFXname2,
		const string& ccyComName,
		QuantoType  quantoType = None, 
		int callput = K_Call,
		bool isInvG1 = false,
		bool isInvG2 = false)
	:	ARM_FXVanilla2D ( mktFXname1, mktFXname2,  
					  ccyComName,quantoType,
					  callput, isInvG1,isInvG2),
		itsCoefficients(2,1.0){}
	
	//copy constructor
	ARM_FXBasket( const ARM_FXBasket& rhs )
		: ARM_FXVanilla2D(rhs), 
		  itsCoefficients(rhs.itsCoefficients){}

	ASSIGN_OPERATOR(ARM_FXBasket)
	
	//destructor
	virtual ~ARM_FXBasket() {}

	/// what it is for ...
	
	virtual double Payoff(double x) const {return 0.0;};
	virtual double Payoff(double x, double y) const;
	
	// accessors
	virtual void SetCoefficients (const ARM_GP_Vector&  coefficients) {itsCoefficients = coefficients;};

    /// Standard ARM object support
	virtual ARM_Object* Clone() const const { return new ARM_FXBasket(*this); };
};


/// ------------------------------------
///	--- ARM_FXBasketInvFx1Fx2 Payoff --- 
/// ------------------------------------
class ARM_FXBasketInvFx1Fx2 : public ARM_FXBasket 
{
public: 
	/// standard constructor
	ARM_FXBasketInvFx1Fx2(const string& mktFXname1, const string& mktFXname2,
		const string& ccyComName,
		QuantoType  quantoType = None, 
		int callput = K_Call,
		bool isInvG1 = false,
		bool isInvG2 = false)
		: ARM_FXBasket( mktFXname1, mktFXname2, ccyComName, quantoType,	callput, isInvG1, isInvG2){}

	//copy constructor
	ARM_FXBasketInvFx1Fx2( const ARM_FXBasketInvFx1Fx2& rhs )
		: ARM_FXBasket(rhs){}
	ASSIGN_OPERATOR(ARM_FXBasketInvFx1Fx2)

	//destructor
	virtual ~ARM_FXBasketInvFx1Fx2() {}

/// what it is for ...
	//virtual double Payoff(double x) const  {return 0.0;}
	virtual double Payoff(double x, double y) const;	

    /// Standard ARM object support
	virtual ARM_Object* Clone() const { return new ARM_FXBasketInvFx1Fx2(*this); };
};

/// ------------------------------------
///	--- ARM_FXBasketFx1InvFx2 Payoff ---
/// ------------------------------------
class ARM_FXBasketFx1InvFx2 : public ARM_FXBasket
{
public: 
	/// standard constructor
	ARM_FXBasketFx1InvFx2(const string& mktFXname1, const string& mktFXname2,
		const string& ccyComName,
		QuantoType  quantoType = None, 
		int callput = K_Call,
		bool isInvG1 = false,
		bool isInvG2 = false)
		: ARM_FXBasket( mktFXname1, mktFXname2, ccyComName, quantoType,	callput, isInvG1, isInvG2){}

	//copy constructor
	ARM_FXBasketFx1InvFx2( const ARM_FXBasketFx1InvFx2& rhs )
		: ARM_FXBasket(rhs){}
	ASSIGN_OPERATOR(ARM_FXBasketFx1InvFx2)

	//destructor
	virtual ~ARM_FXBasketFx1InvFx2() {}

/// what it is for ...
	//virtual double Payoff(double x) const  {return 0.0;}
	virtual double Payoff(double x, double y) const;	

    /// Standard ARM object support
	virtual ARM_Object* Clone() const { return new ARM_FXBasketFx1InvFx2(*this); };
};

/// ---------------------------------------
///	--- ARM_FXBasketInvFx1InvFx2 Payoff ---
/// ---------------------------------------
class ARM_FXBasketInvFx1InvFx2 : public ARM_FXBasket
{
public: 
	/// standard constructor
	ARM_FXBasketInvFx1InvFx2(const string& mktFXname1, const string& mktFXname2,
		const string& ccyComName,
		QuantoType  quantoType = None, 
		int callput = K_Call,
		bool isInvG1 = false,
		bool isInvG2 = false)
		: ARM_FXBasket( mktFXname1, mktFXname2, ccyComName, quantoType,	callput, isInvG1, isInvG2){}

	//copy constructor
	ARM_FXBasketInvFx1InvFx2( const ARM_FXBasketInvFx1InvFx2& rhs )
		: ARM_FXBasket(rhs){}
	ASSIGN_OPERATOR(ARM_FXBasketInvFx1InvFx2)

	//destructor
	virtual ~ARM_FXBasketInvFx1InvFx2() {}

/// what it is for ...
	//virtual double Payoff(double x) const  {return 0.0;}
	virtual double Payoff(double x, double y) const;	

    /// Standard ARM object support
	virtual ARM_Object* Clone() const { return new ARM_FXBasketInvFx1InvFx2(*this); };
};

/// --------------------------------
///	--- GaussianReplication Base Class
/// --------------------------------
class ARM_GaussReplic : public ARM_RootObject
{
protected:
	int itsNbGLPoints; // number of Gauss Legendre points
	ARM_GP_Matrix itsX_W_Fwd;
	ARM_FXVanilla* itsFXVanilla;

public: 
	/// standard constructor
	ARM_GaussReplic( ARM_FXVanilla* vanilla, int nbGLPoints = K_NbLGPoint )
		:
	itsNbGLPoints(nbGLPoints),
	itsX_W_Fwd(),
	itsFXVanilla(CreateClone(vanilla)){}

	//copy constructor
	ARM_GaussReplic( const ARM_GaussReplic& rhs )
		:itsNbGLPoints(rhs.itsNbGLPoints), 
		 itsX_W_Fwd(rhs.itsX_W_Fwd),
	     itsFXVanilla(CreateClone(rhs.itsFXVanilla)){}
	ARM_GaussReplic& operator = ( const ARM_GaussReplic& rhs )
	{
		if( this != &rhs ){
			ARM_RootObject::operator=( rhs );
			itsNbGLPoints = rhs.itsNbGLPoints;
			itsX_W_Fwd = rhs.itsX_W_Fwd;
			itsFXVanilla = CreateClone(rhs.itsFXVanilla);
		}
		return *this;
	}
	//destructor
	virtual ~ARM_GaussReplic() {delete itsFXVanilla;itsFXVanilla=NULL;}

	/// what it is for ...
	virtual double Price() const  = 0;

	// accessors

	virtual int GetNbGLPoints() const {return itsNbGLPoints;};
	virtual void SetNbGLPoints(  int nbGLPoints ) { itsNbGLPoints=nbGLPoints; };

	virtual ARM_GP_Matrix GetX_W_Fwd() const {return itsX_W_Fwd;};
	virtual void SetX_W_Fwd(  ARM_GP_Matrix x_w_fwd ) { itsX_W_Fwd=x_w_fwd; };

	virtual ARM_FXVanilla* GetFXVanilla() const {return itsFXVanilla;};
	virtual void SetFXVanilla(  ARM_FXVanilla* fxvanilla ) { itsFXVanilla=fxvanilla; };


    /// Standard ARM object support
	virtual ARM_Object* Clone() const { return NULL; }
	virtual string toString(const string& indent="",const string& nextIndent="") const { return string("Gauss Replication");};
	virtual string ExportShortName() const { return "LGREP";}
};

/// --------------------------------
///	--- GaussianReplication1D Class
/// --------------------------------
class ARM_GaussReplic1D : public ARM_GaussReplic
{
public: 
	/// standard constructor
	ARM_GaussReplic1D(ARM_FXVanilla* vanilla,int nbGLPoints = K_NbLGPoint)
		:ARM_GaussReplic( vanilla,nbGLPoints ){}

	//copy constructor
	ARM_GaussReplic1D( const ARM_GaussReplic1D& rhs )
		:ARM_GaussReplic(rhs){}
	ASSIGN_OPERATOR(ARM_GaussReplic1D)
	
	//destructor
	virtual ~ARM_GaussReplic1D() {};

	/// what it is for ...
	virtual double Price() const;	
	
    /// Standard ARM object support
	virtual ARM_Object* Clone() const { return new ARM_GaussReplic1D(*this); };
	//virtual string toString(const string& indent="",const string& nextIndent="") const;
};

/// --------------------------------------
///	--- GaussianReplication1DQuanto Class
/// --------------------------------------
class ARM_GaussReplic1DQuanto : public ARM_GaussReplic1D
{
public: 
	/// standard constructor
	ARM_GaussReplic1DQuanto(ARM_FXVanilla* vanilla,int nbGLPoints = K_NbLGPoint)
		:ARM_GaussReplic1D( vanilla,nbGLPoints ){}

	//copy constructor
	ARM_GaussReplic1DQuanto( const ARM_GaussReplic1DQuanto& rhs )
		:ARM_GaussReplic1D(rhs){}
	
	ASSIGN_OPERATOR(ARM_GaussReplic1DQuanto)
	
	//destructor
	virtual ~ARM_GaussReplic1DQuanto() {}

	/// what it is for ...
	virtual double Price() const;	
	
    /// Standard ARM object support
	virtual ARM_Object* Clone() const { return new ARM_GaussReplic1DQuanto(*this); };
	//virtual string toString(const string& indent="",const string& nextIndent="") const;
};

/// --------------------------------
///	--- GaussianReplication2D Class
/// --------------------------------
class ARM_GaussReplic2D : public ARM_GaussReplic
{
protected:
	int itsNbGLPoints2;
	ARM_GP_Matrix itsX_W_Fwd2;
	double itsRho;
public: 
	/// standard constructor
	ARM_GaussReplic2D(ARM_FXVanilla* vanilla,int nbPoint = K_NbLGPoint)
		:ARM_GaussReplic(vanilla, K_NbLGPoint),
		 itsNbGLPoints2(K_NbLGPoint), 
		 itsX_W_Fwd2(),
		 itsRho(0){}
	//copy constructor
	ARM_GaussReplic2D ( const ARM_GaussReplic2D& rhs )
		:ARM_GaussReplic(rhs),
		 itsNbGLPoints2(rhs.itsNbGLPoints2), 
		 itsX_W_Fwd2(rhs.itsX_W_Fwd2),
		 itsRho(rhs.itsRho){}
	ASSIGN_OPERATOR(ARM_GaussReplic2D)
	
	//destructor
	virtual ~ARM_GaussReplic2D(){}

	/// what it is for ...
	virtual double Price() const;

	// accessors

	virtual int GetNbGLPoints2() const {return itsNbGLPoints2;};
	virtual void SetNbGLPoints2(  int nbGLPoints2 ) { itsNbGLPoints2=nbGLPoints2; };

	virtual ARM_GP_Matrix GetX_W_Fwd2() const {return itsX_W_Fwd2;};
	virtual void SetX_W_Fwd2(  ARM_GP_Matrix x_w_fwd2 ) { itsX_W_Fwd2=x_w_fwd2; };

	virtual double GetRho() const {return itsRho;};
	virtual void SetRho(  double rho ) { itsRho=rho; };

    /// Standard ARM object support
	virtual ARM_Object* Clone() const { return new ARM_GaussReplic2D(*this); };
	//virtual string toString(const string& indent="",const string& nextIndent="") const;
};

/// -----------------------------------------
///	--- GaussianReplication2DQuanto1 Class---
/// -----------------------------------------
class ARM_GaussReplic2DQuanto1 : public ARM_GaussReplic2D
{
public: 
	/// standard constructor
	ARM_GaussReplic2DQuanto1(ARM_FXVanilla* vanilla,int nbGLPoints = K_NbLGPoint)
		:ARM_GaussReplic2D(vanilla,K_NbLGPoint){}

	//copy constructor
	ARM_GaussReplic2DQuanto1 ( const ARM_GaussReplic2DQuanto1& rhs )
		:ARM_GaussReplic2D(rhs){}
	ASSIGN_OPERATOR(ARM_GaussReplic2DQuanto1)
	
	//destructor
	virtual ~ARM_GaussReplic2DQuanto1(){}

	/// what it is for ...
	virtual double Price() const;

    /// Standard ARM object support
	virtual ARM_Object* Clone() const { return new ARM_GaussReplic2DQuanto1(*this); };
};

/// -----------------------------------------
///	--- GaussianReplication2DQuanto2 Class---
/// -----------------------------------------
class ARM_GaussReplic2DQuanto2 : public ARM_GaussReplic2D
{
public: 
	/// standard constructor
	ARM_GaussReplic2DQuanto2(ARM_FXVanilla* vanilla,int nbGLPoints = K_NbLGPoint)
		:ARM_GaussReplic2D(vanilla,K_NbLGPoint){}

	//copy constructor
	ARM_GaussReplic2DQuanto2( const ARM_GaussReplic2DQuanto2& rhs )
		:ARM_GaussReplic2D(rhs){}
	ASSIGN_OPERATOR(ARM_GaussReplic2DQuanto2)
	
	//destructor
	virtual ~ARM_GaussReplic2DQuanto2(){}

	/// what it is for ...
	virtual double Price() const;

    /// Standard ARM object support
	virtual ARM_Object* Clone() const { return new ARM_GaussReplic2DQuanto2(*this); };
};


CC_END_NAMESPACE()

#endif
//-----------------------------------------------------------------------------
/*---- End of file ----*/



