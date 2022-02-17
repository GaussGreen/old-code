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

#include "gpbase/gpvector.h"
#include "gpbase/gpmatrix.h"
#include "gpbase/assignop.h"
#include "gpbase/rootobject.h"
#include "gpbase/typedef.h"
#include "gpbase/cloneutilityfunc.h"
#include "gpbase/env.h"

#include "gpinfra/typedef.h"

#include <math.h>

const long  K_NbLGPoint   =   121;
const int	K_Call		  =     1;
const long  K_Epsilon	  =   10e-4;

CC_BEGIN_NAMESPACE( ARM )
/// ------------------------------
///	--- FX Vanilla Base Class ----
/// ------------------------------
class ARM_FXVanilla : public ARM_RootObject
{

protected:
	bool itsIsInvG;
	int itsCallPut;
	double itsStrike;
public: 
	/// standard constructor
	ARM_FXVanilla(double strike,
		int callput = K_Call,
		bool isInvG = false)
		:
	//itsFXname(fxName),
	itsIsInvG(isInvG),
	itsCallPut(callput),
	itsStrike(strike){}

	//copy constructor
	ARM_FXVanilla( const ARM_FXVanilla& rhs )
		:
	ARM_RootObject(rhs),
	itsIsInvG(rhs.itsIsInvG),
	itsCallPut(rhs.itsCallPut),
	itsStrike(rhs.itsStrike){}

	ARM_FXVanilla& operator = ( const ARM_FXVanilla& rhs )
	{
		if( this != &rhs ){
			ARM_RootObject::operator=( rhs );
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
	virtual bool GetIsInvG() const {return itsIsInvG;};
	virtual void SetIsInvG(  bool isInvG ) {itsIsInvG=isInvG;};

	virtual double GetStrike() const {return itsStrike;};
	virtual void SetStrike(  double strike ) { itsStrike=strike; };

	// pure virtual accessors 
	virtual void SetCoefficients ( const ARM_GP_Vector&  coefficients) {};

	virtual bool GetIsInvG2() const { return false; };
	virtual void SetIsInvG2(  bool isInvG2 ) {};

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
public:
	enum FXVanilla2D
	{
		Fx1_Fx2=0,			/// When we can find a proba where FX1 and FX2 are martingal
		InvFx1_InvFx2,		/// When we can find a proba where 1/FX1 and 1/FX2 are martingal
		Fx1_InvFx2,			/// When we can find a proba where FX1 and 1/FX2 are martingal
		InvFx1_Fx2,			/// When we can find a proba where 1/FX1 and FX2 are martingal
	};

	typedef ARM_FXVanilla2D::FXVanilla2D  ARM_FXVanilla2DType;
protected:
	bool itsIsInvG2;
	ARM_FXVanilla2DType  itsVanilla2dType;

public: 
	
	/// standard constructor
	ARM_FXVanilla2D(
		double strike,
		int callput = K_Call,
		ARM_FXVanilla2DType  fxVanilla2DType = Fx1_Fx2,		
		bool isInvG1 = false,
		bool isInvG2 = false)
	:
ARM_FXVanilla(strike,callput,isInvG1),
		  itsIsInvG2(isInvG2),
		  itsVanilla2dType(fxVanilla2DType){}

	//copy constructor
	ARM_FXVanilla2D( const ARM_FXVanilla2D& rhs )
		: ARM_FXVanilla(rhs),
		itsIsInvG2(rhs.itsIsInvG2),
		itsVanilla2dType(rhs.itsVanilla2dType){}
	ASSIGN_OPERATOR(ARM_FXVanilla2D)

	//destructor
	virtual ~ARM_FXVanilla2D() {}

/// what it is for ...
	virtual double Payoff( double x) const {return 0.0;};
	virtual double Payoff(double x, double y) const {return 0.0;};

	// accessors
	virtual bool GetIsInvG2() const {return itsIsInvG2;};
	virtual void SetIsInvG2(  bool isInvG2 ) {itsIsInvG2=isInvG2;};

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
	ARM_FXCall(double strike,
		int callput = K_Call,
		ARM_FXVanilla2DType  fxVanilla2DType = Fx1_Fx2,		
		bool isInvG1 = false,
		bool isInvG2 = false);
	//copy constructor
	ARM_FXCall( const ARM_FXCall& rhs );
	ASSIGN_OPERATOR(ARM_FXCall)
	
	//destructor
	virtual ~ARM_FXCall(){}

	/// what it is for ...
	virtual double Payoff(double x ) const;	
	virtual double Payoff(double x, double y) const { return Payoff(x);}; 

    /// Standard ARM object support
	virtual ARM_Object* Clone() const { return new ARM_FXCall(*this); };
};

/// ---------------------------------
///	--- ARM_FXSpread Payoff ---
/// ---------------------------------
class ARM_FXSpread : public ARM_FXVanilla2D
{
protected:
	double itsAlpha;
	double itsBeta;
public: 
	/// standard constructor
	ARM_FXSpread(double strike,
		int callput = K_Call,
		ARM_FXVanilla2DType  fxVanilla2Dype = Fx1_Fx2, 
		double alpha = 0.0,
		double beta = 0.0,
		bool isInvG1 = false,
		bool isInvG2 = false);
	//copy constructor
	ARM_FXSpread( const ARM_FXSpread& rhs );
	ASSIGN_OPERATOR(ARM_FXSpread)
	
	//destructor
	virtual ~ARM_FXSpread() {}

	/// what it is for ...	
	virtual double Payoff(double x) const {return 0.0;};
	virtual double Payoff(double x, double y) const;
	
	// accessors
	virtual void SetAlpha (double  alpha) {itsAlpha = alpha;};
	virtual void SetBeta (double  beta) {itsBeta = beta;};

    /// Standard ARM object support
	virtual ARM_Object* Clone() const const { return new ARM_FXSpread(*this); };
};

///	--- ARM_FXBasket Payoff ---
/// ---------------------------------
class ARM_FXBasket : public ARM_FXVanilla2D
{
protected:
	double itsAlpha;
	double itsBeta;
public: 
	/// standard constructor
	ARM_FXBasket(double strike,
		int callput = K_Call,
		ARM_FXVanilla2DType  fxVanilla2Dype = Fx1_Fx2, 
		bool isInvG1 = false,
		bool isInvG2 = false)
	:	
ARM_FXVanilla2D (strike,callput, fxVanilla2Dype,isInvG1, isInvG2),
	itsAlpha(1.0),
	itsBeta(1.0)
	{}
	//copy constructor
	ARM_FXBasket( const ARM_FXBasket& rhs )
		: ARM_FXVanilla2D(rhs), 
	itsAlpha(rhs.itsAlpha),
	itsBeta(rhs.itsBeta){}

	ASSIGN_OPERATOR(ARM_FXBasket)
	
	//destructor
	virtual ~ARM_FXBasket() {}

	/// what it is for ...
	
	virtual double Payoff(double x) const {return 0.0;};
	virtual double Payoff(double x, double y) const;
	
    /// Standard ARM object support
	virtual ARM_Object* Clone() const const { return new ARM_FXBasket(*this); };
};

/// -------------------------
///	--- ARM_FXDigital Payoff --- 
/// -------------------------
class ARM_FXDigital : public ARM_FXVanilla2D
{
protected:
	ARM_DigitType itsDigitType;
	double itsEpsilon;
public:
	/// standard constructor
	ARM_FXDigital(double strike,
		int callput = K_Call,
		ARM_FXVanilla2DType fxVanilla2DType = Fx1_Fx2,		
		bool isInvG1 = false,
		bool isInvG2 = false,
		ARM_DigitType digitType = ARM_FXDigitType::centred,		
		double epsilon = K_Epsilon);
	//copy constructor
	ARM_FXDigital( const ARM_FXDigital& rhs );
	ASSIGN_OPERATOR(ARM_FXDigital)
	
	//destructor
	virtual ~ARM_FXDigital(){}

	/// what it is for ...
	virtual double Payoff(double x ) const;	
	virtual double Payoff(double x, double y) const { return Payoff(x);}; 

    /// Standard ARM object support
	virtual ARM_Object* Clone() const { return new ARM_FXDigital(*this); };
};

/// -------------------------
///	--- ARM_FXQuotient Payoff --- 
/// -------------------------
class ARM_FXQuotient : public ARM_FXVanilla2D
{
protected:
	double itsAlpha;
	double itsBeta;
	double itsStrike2;
public: 
	/// standard constructor
	ARM_FXQuotient(double strike,
		int callput = K_Call,
		ARM_FXVanilla2DType  fxVanilla2DType = Fx1_Fx2,
		double alpha = 0.0,
		double beta = 0.0,
		double strike2 = 1.0,
		bool isInvG1 = false,
		bool isInvG2 = false);
	//copy constructor
	ARM_FXQuotient( const ARM_FXQuotient& rhs );
	ASSIGN_OPERATOR(ARM_FXQuotient)
	
	//destructor
	virtual ~ARM_FXQuotient(){}

	/// what it is for ...
	virtual double Payoff(double x) const {return 0.0;};
	virtual double Payoff(double x, double y) const;

    /// Standard ARM object support
	virtual ARM_Object* Clone() const { return new ARM_FXQuotient(*this); };
};

/*
/// ------------------------------------
///	--- ARM_FXBasketInvFx1Fx2 Payoff --- 
/// ------------------------------------
class ARM_FXBasketInvFx1Fx2 : public ARM_FXBasket 
{
public: 

	ARM_FXBasketInvFx1Fx2( const string& mktFXname1, const string& mktFXname2,
		const string& ccyComName,
		double strike ,
		int callput = K_Call,
		QuantoType  quantoType = None, 
		bool isInvG1 = false,
		bool isInvG2 = false)
	:	ARM_FXBasket ( mktFXname1, mktFXname2, ccyComName,strike, callput, quantoType, isInvG1,isInvG2)
	{}
	
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

	ARM_FXBasketFx1InvFx2( const string& mktFXname1, const string& mktFXname2,
		const string& ccyComName,
		double strike ,
		int callput = K_Call,
		QuantoType  quantoType = None, 
		bool isInvG1 = false,
		bool isInvG2 = false)
		:	ARM_FXBasket ( mktFXname1, mktFXname2, ccyComName,strike, callput, quantoType, isInvG1,isInvG2){}
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
	ARM_FXBasketInvFx1InvFx2( const string& mktFXname1, const string& mktFXname2,
		const string& ccyComName,
		double strike ,
		int callput = K_Call,
		QuantoType  quantoType = None, 
		bool isInvG1 = false,
		bool isInvG2 = false)
		:	ARM_FXBasket ( mktFXname1, mktFXname2, ccyComName,strike, callput, quantoType, isInvG1,isInvG2){}
	

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
};*/

/// --------------------------------
///	--- GaussianReplication Base Class
/// --------------------------------
class ARM_GaussReplic : public ARM_RootObject
{
protected:
	ARM_GP_Matrix itsGLParams;
	ARM_FXVanilla* itsFXVanilla;

public: 
	/// standard constructor
	ARM_GaussReplic( ARM_FXVanilla* vanilla, 
		const ARM_GP_Matrix& glParams )
	:
	ARM_RootObject(),
	itsGLParams(glParams),
	itsFXVanilla(CreateClone(vanilla)){}

	//copy constructor
	ARM_GaussReplic( const ARM_GaussReplic& rhs )
		:itsGLParams(rhs.itsGLParams),
	 itsFXVanilla(CreateClone(rhs.itsFXVanilla))
	{}
	ARM_GaussReplic& operator = ( const ARM_GaussReplic& rhs )
	{
		if( this != &rhs ){
			ARM_RootObject::operator=( rhs );
			itsGLParams = rhs.itsGLParams;
			itsFXVanilla = CreateClone(rhs.itsFXVanilla);
		}
		return *this;
	}
	//destructor
	virtual ~ARM_GaussReplic() {delete itsFXVanilla;itsFXVanilla=NULL;}

	/// what it is for ...
	virtual double Price() const  = 0;

	// accessors
	virtual ARM_GP_Matrix GetGLParams() const {return itsGLParams;};
	virtual void SetGLParams(  const ARM_GP_Matrix& glparams ) { itsGLParams = glparams; };

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
	enum QuantoType
	{
		Quanto,			/// When PAY/COM has the same model than the first FX  
	    InvQuanto,		/// When PAY/COM has the same model than the second FX 
	};
	typedef ARM_GaussReplic1D::QuantoType ARM_1DQuantoType;

	ARM_1DQuantoType  itsQuanto1DType;
public: 
	/// standard constructor
	ARM_GaussReplic1D(ARM_FXVanilla* vanilla,const ARM_GP_Matrix& glParams ,ARM_1DQuantoType type = Quanto);

	//copy constructor
	ARM_GaussReplic1D( const ARM_GaussReplic1D& rhs );
	ASSIGN_OPERATOR(ARM_GaussReplic1D)
	
	//destructor
	virtual ~ARM_GaussReplic1D() {};

	/// what it is for ...
	virtual double Price() const;	
	
    /// Standard ARM object support
	virtual ARM_Object* Clone() const { return new ARM_GaussReplic1D(*this); };
};

/// --------------------------------------
///	--- GaussianReplication1DQuanto Class
/// --------------------------------------
class ARM_GaussReplic1DQuanto : public ARM_GaussReplic1D
{
public: 
	/// standard constructor
	ARM_GaussReplic1DQuanto(ARM_FXVanilla* vanilla, const ARM_GP_Matrix& glParams );
	//copy constructor
	ARM_GaussReplic1DQuanto( const ARM_GaussReplic1DQuanto& rhs );	
	ASSIGN_OPERATOR(ARM_GaussReplic1DQuanto)
	
	//destructor
	virtual ~ARM_GaussReplic1DQuanto() {}

	/// what it is for ...
	virtual double Price() const;	
	
    /// Standard ARM object support
	virtual ARM_Object* Clone() const { return new ARM_GaussReplic1DQuanto(*this); };
};

/// --------------------------------
///	--- GaussianReplication2D Class
/// --------------------------------
class ARM_GaussReplic2D : public ARM_GaussReplic
{
public:
	enum QuantoType
	{
		Quanto1,			/// When PAY/COM has the same model than the first FX  
	    Quanto2,			/// When PAY/COM has the same model than the second FX 
		None
	};
	typedef ARM_GaussReplic2D::QuantoType ARM_QuantoType;
private:
	ARM_GP_Matrix itsGLParams2;
	double itsRho;
	ARM_QuantoType itsQuantoType;
public: 
	/// standard constructor
	ARM_GaussReplic2D(ARM_FXVanilla* vanilla,
		double rho = 0.0,
		ARM_QuantoType  quantoType = None,
		const ARM_GP_Matrix& glmatrix = ARM_GP_Matrix(), 
		const ARM_GP_Matrix& glmatrix2 = ARM_GP_Matrix());
	//copy constructor
	ARM_GaussReplic2D ( const ARM_GaussReplic2D& rhs );
	ASSIGN_OPERATOR(ARM_GaussReplic2D)
	
	//destructor
	virtual ~ARM_GaussReplic2D(){}

	/// what it is for ...
	virtual double Price() const;

	// accessors
	virtual ARM_GP_Matrix GetGLParams2() const {return itsGLParams2;};
	virtual void SetGLParams2(  const ARM_GP_Matrix& glparams ) { itsGLParams2 = glparams; };

	virtual double GetRho() const {return itsRho;};
	virtual void SetRho(  double rho ) { itsRho =rho; };

	ARM_QuantoType GetQuantoType() const {return itsQuantoType;}

    /// Standard ARM object support
	virtual ARM_Object* Clone() const { return new ARM_GaussReplic2D(*this); };
};

CC_END_NAMESPACE()

#endif
//-----------------------------------------------------------------------------
/*---- End of file ----*/



