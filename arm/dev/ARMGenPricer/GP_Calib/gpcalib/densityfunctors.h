/*!
 *
 * Copyright (c) CDC IXIS CM November 2003 Paris
 * 
 * \file numerical.h
 *  \brief Model fitter that calibrates a whole smile numerically
 * 
 *	\author  A. Schauly
 *	\version 1.0
 *	\date August 2005
 */

#ifndef _INGPCALIB_NUMDENSITYFUNCTORS_H
#define _INGPCALIB_NUMDENSITYFUNCTORS_H

#define UNIMPLEMENTED_DENSITYFUNCTOR_FUNCTION  { ARM_THROW( ERR_INVALID_ARGUMENT, "unimplemented function in densityfunctor"); }

#include "gpbase/gpvector.h"
#include "gpbase/assignop.h"
#include "gpbase/rootobject.h"
#include "gpbase/typedef.h"
#include "gpbase/utilityport.h"


#include "gpclosedforms/smile_shiftedlognormal.h"
#include "gpclosedforms/smile_sabr.h"
#include "gpclosedforms/vanilla_bs.h"
#include "gpclosedforms/numerics_interface.h"

#include "gpclosedforms/spreadoption_bisabr_interface.h"
#include "gpclosedforms/bisabr_interface.h"
#include "gpnumlib/gaussiananalytics.h"
#include "gpmodels/qmodelanalytics.h"
#include "gpclosedforms/normal_heston.h"

#include <math.h>

class ARM_Security;
class ARM_Model;

CC_BEGIN_NAMESPACE( ARM )

/// --------------------------------
///	--- Density Functor Base Class
/// --------------------------------
class ARM_DensityFunctor : public ARM_RootObject
{
private:
	bool itsIsDirect; /// to calculate the direct density/ reverse density default(true)
public: 
	/// constructor/destructor
	ARM_DensityFunctor(bool isDirect = true); 
	ARM_DensityFunctor( const ARM_DensityFunctor& rhs ); 
	ARM_DensityFunctor& operator = ( const ARM_DensityFunctor& rhs );
	virtual ~ARM_DensityFunctor() {}

	/// what it is for ...
	virtual ARM_GP_VectorPtr Quantile( const ARM_GP_VectorPtr& x, double fwd, double maturity );
	void MakeIncreasing(ARM_GP_VectorPtr& quantile);

	virtual double Proba(double strike, double fwd, double maturity) const = 0;
	virtual double Quantile( double x, double fwd, double maturity ) const = 0;	
	virtual double Call_Option(double x, double fwd, double maturity ) const = 0;
	/// Accessors

	inline bool GetIsDirect() const {return itsIsDirect;}
	inline void SetIsDirect(bool direct)  {itsIsDirect = direct;}

    /// Standard ARM object support
	virtual ARM_Object* Clone() const = NULL;
	virtual string toString(const string& indent="",const string& nextIndent="") const;
	virtual string ExportShortName() const { return "LDFUN";}

};

/// --------------------------------
///	--- Shifted LN Density functor
/// --------------------------------
class ARM_ShiftedLNDensityFunctor : public ARM_DensityFunctor
{
private: 
	double itsVol;
	double itsShift;
public: 
	/// constructor/destructor
	ARM_ShiftedLNDensityFunctor(double vol = 0.000001, double shift =0.0, bool isDirect = true);
	ARM_ShiftedLNDensityFunctor( const ARM_ShiftedLNDensityFunctor& rhs );
	ASSIGN_OPERATOR(ARM_ShiftedLNDensityFunctor)
	virtual ~ARM_ShiftedLNDensityFunctor() {}

	/// what it is for ...
	virtual double Quantile( double x, double fwd, double maturity ) const;
	virtual double Proba(double strike, double fwd, double maturity) const {return 0.0;};

	virtual double Call_Option(double x, double fwd, double maturity ) const;

	/// Accessors
	inline double getVol() const				{ return itsVol; }
	inline void setVol( double vol )			{ itsVol = vol; }
	inline double getShift() const				{ return itsShift; }
	inline void setShift( double shift )		{ itsShift = shift; }

    /// Standard ARM object support
	virtual ARM_Object* Clone() const { return new ARM_ShiftedLNDensityFunctor(*this); };
	virtual string toString(const string& indent="",const string& nextIndent="") const;

};

/// --------------------------------
///	--- Q Density functor
/// --------------------------------
class ARM_QDensityFunctor : public ARM_DensityFunctor
{
private: 
	double itsVol;
	double itsQ;
public: 
	/// constructor/destructor
	ARM_QDensityFunctor(double vol = 0.000001, double q =1.0,bool isDirect = true);
	ARM_QDensityFunctor( const ARM_QDensityFunctor& rhs );
	ASSIGN_OPERATOR(ARM_QDensityFunctor)
	virtual ~ARM_QDensityFunctor() {}

	/// what it is for ...
	virtual double Quantile( double x, double fwd, double maturity ) const;
	virtual double Proba(double strike, double fwd, double maturity) const {return 0.0;};
	virtual double Call_Option(double x, double fwd, double maturity ) const;

	/// Accessors
	inline double getVol() const			{ return itsVol; }
	inline void setVol( double vol )		{ itsVol = vol; }
	inline double getQ() const				{ return itsQ; }
	inline void setQ( double q )			{ itsQ = q; }

    /// Standard ARM object support
	virtual ARM_Object* Clone() const { return new ARM_QDensityFunctor(*this); };
	virtual string toString(const string& indent="",const string& nextIndent="") const;

};

/// --------------------------------
///	--- Mixture Density functor
/// --------------------------------
class ARM_MixtureDensityFunctor : public ARM_DensityFunctor
{
private:
	double itsVol1;
	double itsVol2;
	double itsAlpha;
	double itsLambda;

public:
	// Basic constructor
	ARM_MixtureDensityFunctor(bool isDirect = true);
	ARM_MixtureDensityFunctor(
		double vol1 , 
		double vol2 ,
		double alpha,
		double lambda,
		bool isDirect = true);

	ARM_MixtureDensityFunctor(
		double fwd,
		double maturity,
		double volATM,
		double decVol,
		double alpha,
		double lambda,
		bool isDirect = true);

	ARM_MixtureDensityFunctor( const ARM_MixtureDensityFunctor& rhs );
	ASSIGN_OPERATOR(ARM_MixtureDensityFunctor)

	virtual ~ARM_MixtureDensityFunctor() {}

	/// what it is for ...
	virtual double Quantile( double x, double fwd, double maturity ) const ;
	virtual double Call_Option(double x, double fwd, double maturity ) const;
	double Proba(double strike, double fwd, double maturity) const;
	static double Call(double strike, double fwd, double maturity, double vol1, double vol2, double alpha, double lambda);
	
	/// Accessors
	inline double getVol1() const			{ return itsVol1; }
	inline void setVol1( double vol1 )		{ itsVol1 = vol1; }
	inline double getVol2() const			{ return itsVol2; }
	inline void setVol2( double vol2 )		{ itsVol2 = vol2; }
	inline double getAlpha() const			{ return itsAlpha; }
	inline void setAlpha( double alpha )	{ itsAlpha = alpha; }
	inline double getLambda() const			{ return itsLambda; }
	inline void setLambda(double lambda)	{ itsLambda = lambda; }

    /// Standard ARM object support
	virtual ARM_Object* Clone() const { return new ARM_MixtureDensityFunctor(*this); };
	virtual string toString(const string& indent="",const string& nextIndent="") const;
};

/// --------------------------------
///	--- Heston Density functor
/// --------------------------------
class ARM_HestonDensityFunctor : public ARM_DensityFunctor
{
private:
	double itsV0;
	double itsKappa;
	double itsTheta;
	double itsRho;
	double itsNu;
	double itsShift;
	double itsLevel;
	double itsSigma;

public:
	// Basic constructor
	ARM_HestonDensityFunctor(bool isDirect = true)
			:
	ARM_DensityFunctor(isDirect)
	{}
	ARM_HestonDensityFunctor(
		double v0, 
		double kappa,
		double theta,
		double rho,
		double nu,
		double shift,
		double level,
		double sigma = 0.0, 
		bool isDirect = true)
		:
	ARM_DensityFunctor(isDirect),
	itsV0(v0),
	itsKappa(kappa),
	itsTheta(theta),
	itsRho(rho),
	itsNu(nu),
	itsShift(shift),
	itsLevel(level),
	itsSigma(sigma)
	{
	}

	ARM_HestonDensityFunctor( const ARM_HestonDensityFunctor& rhs )
	:ARM_DensityFunctor(rhs),
	itsV0(rhs.itsV0),
	itsKappa(rhs.itsKappa),
	itsTheta(rhs.itsTheta),
	itsRho(rhs.itsRho),
	itsNu(rhs.itsNu),
	itsShift(rhs.itsShift),
	itsLevel(rhs.itsLevel),
	itsSigma(rhs.itsSigma)
	{
	}

	ASSIGN_OPERATOR(ARM_HestonDensityFunctor)

	virtual ~ARM_HestonDensityFunctor() {}

	/// what it is for ...
	virtual double Quantile( double x, double fwd, double maturity ) const ;
	virtual double Proba(double strike, double fwd, double maturity) const;
	virtual double Call_Option(double x, double fwd, double maturity ) const;
		
	static double Call(double strike, double fwd, double maturity,
			double v0,double kappa, double theta, double rho, double nu, double shift, double level, double sigma);
	
	/// Accessors
	inline double getV0() const { return itsV0; }
	inline void setV0( double v0 ) { itsV0 = v0; }
	inline double getKappa() const { return itsKappa; }
	inline void setKappa( double kappa ) { itsKappa = kappa; }
	inline double getTheta() const { return itsTheta; }
	inline void setTheta( double theta ) { itsTheta = theta; }
	inline double getRho() const { return itsRho; }
	inline void setRho( double rho ) { itsRho = rho; }
	inline double getNu() const { return itsNu; }
	inline void setNu(double nu) { itsNu = nu; }
	inline double getShift() const { return itsShift; }
	inline void setShift( double shift ) { itsShift = shift; }
	inline double getLevel() const { return itsLevel; }
	inline void setLevel(double level) { itsLevel = level; }
	inline double getSigma() const { return itsSigma; }
	inline void setSigma(double sigma)  { itsSigma = sigma; }

    /// Standard ARM object support
	virtual ARM_Object* Clone() const { return new ARM_HestonDensityFunctor(*this); };
	virtual string toString(const string& indent="",const string& nextIndent="") const;

};

	
/// ---------------------------
///	--- No Vol Density functor
/// ---------------------------
class ARM_NoVolDensityFunctor : public ARM_DensityFunctor
{
public: 
	/// constructor/destructor
	ARM_NoVolDensityFunctor(bool isDirect = true): ARM_DensityFunctor(isDirect) {}
	ARM_NoVolDensityFunctor( const ARM_NoVolDensityFunctor& rhs )
	: ARM_DensityFunctor(rhs) 
	{}
	ASSIGN_OPERATOR(ARM_NoVolDensityFunctor)
	virtual ~ARM_NoVolDensityFunctor() {}

	/// what it is for ...
	inline double Quantile( double x, double fwd, double maturity ) const { return fwd;}
	virtual double Proba(double strike, double fwd, double maturity) const {return 0.0;};
	inline virtual double Call_Option(double x, double fwd, double maturity ) const { return MAX( fwd - x, 0.); }
	
    /// Standard ARM object support
	virtual ARM_Object* Clone() const { return new ARM_NoVolDensityFunctor(*this); };
	virtual string toString(const string& indent="",const string& nextIndent="") const;

};

/// ------------------------
///	--- SABR Density functor
/// ------------------------
class ARM_SABRDensityFunctor : public ARM_DensityFunctor
{
private:
	double	itsAlpha, /// initial vol
			itsBeta,  /// backbone
			itsRho,   /// correl
			itsNu ;   /// vol of vol
	int		itsSabrType;

	/// Quantile precomputation and storage
	size_t		     itsGridSize;
	ARM_GP_VectorPtr itsX;
	ARM_GP_VectorPtr itsQx;
	
public:
	/// constructor/destructor
	ARM_SABRDensityFunctor(bool isDirect = true);
	ARM_SABRDensityFunctor(double alpha, double beta, double rho, double nu, int sabrType, size_t gridSize = 251);
	ARM_SABRDensityFunctor( const ARM_SABRDensityFunctor& rhs );
	ASSIGN_OPERATOR(ARM_SABRDensityFunctor)
	virtual ~ARM_SABRDensityFunctor() {}

	/// what it is for ...
	/// inline = perf issue
	virtual double Quantile( double x, double fwd, double maturity ) const;
	virtual ARM_GP_VectorPtr Quantile( const ARM_GP_VectorPtr& x, double fwd, double maturity );
	virtual double Proba(double strike, double fwd, double maturity) const;
	inline virtual double Call_Option(double x, double fwd, double maturity ) const 
	{ 
		return SABR_smile::call_option(fwd, x, maturity, itsAlpha, itsBeta, itsRho, itsNu, itsSabrType, 120,0.01,1.5,fwd/2); 
		return SABR_smile::call_option(fwd, x, maturity, itsAlpha, itsBeta, itsRho, itsNu, itsSabrType, 120,0.01,1.5,0.025); 
	}

	ARM_ShiftedLNDensityFunctor* toShiftedLN_2strikes(double fwd, double maturity, double strike1, double strike2, double shiftinf, double shiftsup) const;
	ARM_ShiftedLNDensityFunctor* toShiftedLN_1strike1shift(double fwd, double maturity, double strike1, double shift) const;
	
	/// Standard ARM object support
	virtual ARM_Object* Clone() const { return new ARM_SABRDensityFunctor(*this); };
	virtual string toString(const string& indent="",const string& nextIndent="") const;

	// access
	double	GetAlpha() const	{return itsAlpha;};
	inline void SetAlpha( double alpha ) {itsAlpha = alpha;}
	double	GetBeta() const		{return itsBeta;};
	inline void SetBeta( double beta ) {itsBeta = beta; }
	double	GetRho() const		{return itsRho;};
	inline void SetRho( double rho ) { itsRho = rho; }
	double	GetNu() const		{return itsNu;};
	inline void SetNu( double nu ) { itsNu = nu; }
	int GetSabrType() const		{return itsSabrType;};
	inline void SetSabrType( int sabrType ) { itsSabrType = sabrType; }
};


/// --------------------------------
///	--- BISABR Density functor
/// --------------------------------
class ARM_BiSABRDensityFunctor : public ARM_DensityFunctor
{
private:
	double	itsFwd1,
			itsAlpha1,	/// initial vol
			itsBeta1,	/// backbone
			itsRho1,	/// correl
			itsNu1,		/// vol of vol
			itsFwd2,
			itsAlpha2,
			itsBeta2,
			itsRho2,
			itsNu2,
			itsRhoS1S2,
			itsRhoS1V2,
			itsRhoS2V1,
			itsRhoV1V2;
	int		itsSabrType;

	/// Quantile precomputation and storage
	size_t		     itsGridSize;
	ARM_GP_VectorPtr itsX;
	ARM_GP_VectorPtr itsQx;
	
public:
	/// constructor/destructor
	ARM_BiSABRDensityFunctor(bool isDirect = true)
		:	ARM_DensityFunctor(isDirect), 
			itsFwd1(0), 
			itsAlpha1(0), 
			itsBeta1(0), 
			itsRho1(0),
			itsNu1(0), 
			itsFwd2(0),
			itsAlpha2(0), 
			itsBeta2(0),
			itsRho2(0), 
			itsNu2(0), 
			itsRhoS1S2(0),
			itsRhoS1V2(0), 
			itsRhoS2V1(0), 
			itsRhoV1V2(0), 
			itsSabrType (0),
			itsGridSize(0), 
			itsX(ARM_GP_VectorPtr(NULL)),
			itsQx(ARM_GP_VectorPtr(NULL)){} 
	
	ARM_BiSABRDensityFunctor(double fwd1, double alpha1, double beta1, double rho1, double nu1,
						   double fwd2, double alpha2, double beta2, double rho2, double nu2,
						   double rhoS1S2, double rhoS1V2, double rhoS2V1, double rhoV1V2,
						   int sabrType, size_t gridSize = 251,bool isDirect = true)
		:	ARM_DensityFunctor(isDirect), 
			itsFwd1(fwd1), 
			itsAlpha1(alpha1),
			itsBeta1(beta1), 
			itsRho1(rho1), 
			itsNu1(nu1),
			itsFwd2(fwd2), 
			itsAlpha2(alpha2), 
			itsBeta2(beta2), 
			itsRho2(rho2), 
			itsNu2(nu2),
			itsRhoS1S2(rhoS1S2), 
			itsRhoS1V2(rhoS1V2),
			itsRhoS2V1(rhoS2V1), 
			itsRhoV1V2(rhoV1V2),
			itsSabrType (sabrType),
			itsGridSize(gridSize),
			itsX(ARM_GP_VectorPtr(NULL)), 
			itsQx(ARM_GP_VectorPtr(NULL)){} 
	
	ARM_BiSABRDensityFunctor(double alpha1, double beta1, double rho1, double nu1,
						   double alpha2, double beta2, double rho2, double nu2,
						   double rhoS1S2, double rhoS1V2, double rhoS2V1, double rhoV1V2,
						   int sabrType, size_t gridSize = 251,bool isDirect = true)
		:	ARM_DensityFunctor(isDirect), 
			itsFwd1(0), 
			itsAlpha1(alpha1), 
			itsBeta1(beta1),
			itsRho1(rho1), 
			itsNu1(nu1),
			itsFwd2(0),
			itsAlpha2(alpha2), 
			itsBeta2(beta2), 
			itsRho2(rho2),
			itsNu2(nu2),
			itsRhoS1S2(rhoS1S2), 
			itsRhoS1V2(rhoS1V2), 
			itsRhoS2V1(rhoS2V1), 
			itsRhoV1V2(rhoV1V2),
			itsSabrType (sabrType),
			itsGridSize(gridSize), 
			itsX(ARM_GP_VectorPtr(NULL)),
			itsQx(ARM_GP_VectorPtr(NULL)){} 

	ARM_BiSABRDensityFunctor( const ARM_BiSABRDensityFunctor& rhs )
		:	ARM_DensityFunctor(rhs), 
			itsFwd1(rhs.itsFwd1), 
			itsAlpha1(rhs.itsAlpha1), 
			itsBeta1(rhs.itsBeta1), 
			itsRho1(rhs.itsRho1), 
			itsNu1(rhs.itsNu1),
			itsFwd2(rhs.itsFwd2), 
			itsAlpha2(rhs.itsAlpha2), 
			itsBeta2(rhs.itsBeta2), 
			itsRho2(rhs.itsRho2), 
			itsNu2(rhs.itsNu2),
			itsRhoS1S2(rhs.itsRhoS1S2), 
			itsRhoS1V2(rhs.itsRhoS1V2), 
			itsRhoS2V1(rhs.itsRhoS2V1), 
			itsRhoV1V2(rhs.itsRhoV1V2), 
			itsSabrType(rhs.itsSabrType),
			itsGridSize(rhs.itsGridSize), 
			itsX(ARM_GP_VectorPtr(NULL)), 
			itsQx(ARM_GP_VectorPtr(NULL)){} 
	
	ASSIGN_OPERATOR(ARM_BiSABRDensityFunctor)
	virtual ~ARM_BiSABRDensityFunctor() {}

	/// what it is for ...
	/// inline = perf issue	
	virtual double Call_Option(double x, double fwd, double maturity ) const; 
	virtual double Proba(double strike, double fwd, double maturity) const {return 0.0;};
	virtual double Quantile( double x,double fwd, double maturity ) const;
	virtual ARM_GP_VectorPtr Quantile( const ARM_GP_VectorPtr& x, double fwd, double maturity );

	/// Standard ARM object support
	virtual ARM_Object* Clone() const { return new ARM_BiSABRDensityFunctor(*this); };
	virtual string toString(const string& indent="",const string& nextIndent="") const;

	// access
	void	SetForward1(double forward)		{itsFwd1 = forward;};
	void	SetForward2(double forward)		{itsFwd2 = forward;};

	double	GetAlpha(int iEgale1ou2) const	{return iEgale1ou2 = 1 ? itsAlpha1 : itsAlpha2;};
	double	GetBeta(int iEgale1ou2) const	{return iEgale1ou2 = 1 ? itsBeta1 : itsBeta2;};
	double	GetRho(int iEgale1ou2) const	{return iEgale1ou2 = 1 ? itsRho1 : itsRho2;};
	double	GetNu(int iEgale1ou2) const		{return iEgale1ou2 = 1 ? itsNu1 : itsNu2;};
};


class SLNDensity_Approx : public ARM_GP::UnaryFunc<double,double> 
{
public: 
		SLNDensity_Approx(double fwd,double time,double strike1,double target1,double strike2,double target2) :
		itsFwd(fwd),
		itsTime(time),
		itsStrike1(strike1),
		itsTarget1(target1),
		itsStrike2(strike2),
		itsTarget2(target2)
		{
		};

		inline virtual double operator() (double shift) const
		{
			double price2	= BS(itsFwd + shift, itsStrike2 + shift, itsTime, GetVol(shift));
			return price2-itsTarget2;
		};
		inline double GetVol(double shift) const
		{
			bool success(true);
			double sigma = VanillaImpliedVol_BS(itsFwd + shift, itsStrike1 + shift, itsTime, itsTarget1, 1, NULL, &success);
			return sigma;
		};
private:
	double						itsFwd;
	double						itsTime;
	double						itsStrike1;
	double						itsTarget1;
	double						itsStrike2;
	double						itsTarget2;
};



class ARM_SmileViewer : public ARM_RootObject
{
public:
	enum MoneyType
	{
		GAUSS = 0,
		BLACK
	};
private:
	ARM_GP_Vector	itsMoney;
	ARM_GP_Vector	itsStrike;
	ARM_GP_Vector	itsVol;
	MoneyType		itsType;
	ARM_GP_Vector	itsExtraStrikes;
public:
    ARM_SmileViewer(const ARM_GP_Vector& money,MoneyType type = GAUSS,const ARM_GP_Vector& strikes=ARM_GP_Vector());
    ARM_SmileViewer(const ARM_SmileViewer& rhs);
	ASSIGN_OPERATOR(ARM_SmileViewer)
	virtual ~ARM_SmileViewer();

	void Compute(ARM_Security* pSec, ARM_Model* pMod);
	
	inline const ARM_GP_Vector&	GetMoney()		{ return itsMoney;};
	inline const ARM_GP_Vector&	GetStrike()		{ return itsStrike;};
	inline const ARM_GP_Vector&	GetVol()		{ return itsVol;};
	inline const MoneyType&	GetVolType()		{ return itsType;};

	/// Standard ARM object support
	virtual string ExportShortName() const { return "LSMVW";}
	virtual ARM_Object* Clone() const;
    virtual string toString(const string& indent="", const string& nextIndent="") const;
    virtual void View(char* id = NULL, FILE* ficOut = NULL) const;
};

/// --------------------------------
///	--- Spline Density functor
/// --------------------------------
class ARM_SplineDensityFunctor : public ARM_DensityFunctor
{
private: 
	ARM_GP_Vector				itsVar;
	ARM_GP_Vector				itsMoney;
	ARM_GP_Vector				itsDer2;
	size_t						itsSize;
	ARM_SmileViewer::MoneyType  itsType;
public: 
	/// constructor/destructor
	ARM_SplineDensityFunctor( const ARM_GP_Vector& money, const ARM_GP_Vector& vol, ARM_SmileViewer::MoneyType type);
	ARM_SplineDensityFunctor( const ARM_SplineDensityFunctor& rhs );
	ASSIGN_OPERATOR(ARM_SplineDensityFunctor)
	virtual ~ARM_SplineDensityFunctor() {};

	/// what it is for ...
	virtual double Proba(double strike, double fwd, double maturity) const {return 0.0;};
	virtual double Quantile(double x, double fwd, double mat) const;

	inline virtual double Call_Option(double x, double fwd, double maturity ) const;

	void BuildSpline(double yp1, double yp2);

	static double ProbaGauss(double strike, double fwd, double vol, double dvoldm);
	static double ProbaBlack(double strike, double fwd, double vol, double dvoldm);
	static double Proba(double strike, double fwd, const ARM_GP_Vector& money, const ARM_GP_Vector& var, const ARM_GP_Vector& der2,size_t type);
	static double ComputeVol(double strike, double fwd, const ARM_GP_Vector& money, const ARM_GP_Vector& var, const ARM_GP_Vector& der2,size_t type);
	static double ComputeVol(double eps, const ARM_GP_Vector& money, const ARM_GP_Vector& var, const ARM_GP_Vector& der2);
	static double ComputeDVol(double strike, double fwd, const ARM_GP_Vector& money, const ARM_GP_Vector& var, const ARM_GP_Vector& der2,size_t type);
	static double ComputeDVol(double eps, const ARM_GP_Vector& money, const ARM_GP_Vector& var, const ARM_GP_Vector& der2);
	
	double Quantile(double x, double fwd);
	
    /// Standard ARM object support
	virtual ARM_Object* Clone() const { return new ARM_SplineDensityFunctor(*this); };
	virtual string toString(const string& indent="",const string& nextIndent="") const;

};

class ARM_NormalHestonDensityFunctor : public ARM_DensityFunctor
{
private:
	double	itsFwd;
	double	itsV0;
	double	itsKappa;
	double	itsTheta;
	double	itsVVol;
	double	itsRho;
	double	itsLevel;

	/// Quantile precomputation and storage
	size_t		     itsGridSize;
	ARM_GP_VectorPtr itsX;
	ARM_GP_VectorPtr itsQx;
	
public:
	/// constructor/destructor
	ARM_NormalHestonDensityFunctor()
		:	ARM_DensityFunctor(), 
			itsFwd(0), itsV0(0), itsKappa(0), itsTheta(0), itsVVol(0), itsRho(0), itsLevel(0),
			itsGridSize(0), itsX(ARM_GP_VectorPtr(NULL)), itsQx(ARM_GP_VectorPtr(NULL)){} 
	
	ARM_NormalHestonDensityFunctor(double fwd, double v0, double kappa, double theta, double vvol, 
			double rho, double level, size_t gridSize = 251)
		:	ARM_DensityFunctor(), 
			itsFwd(fwd), itsV0(v0), itsKappa(kappa), itsTheta(theta), itsVVol(vvol),
			itsRho(rho), itsLevel(level),
			itsGridSize(gridSize), itsX(ARM_GP_VectorPtr(NULL)), itsQx(ARM_GP_VectorPtr(NULL)){} 
	
	ARM_NormalHestonDensityFunctor( const ARM_NormalHestonDensityFunctor& rhs )
		:	ARM_DensityFunctor(rhs), 
			itsFwd(rhs.itsFwd), itsV0(rhs.itsV0), itsKappa(rhs.itsKappa), itsTheta(rhs.itsTheta), 
			itsVVol(rhs.itsVVol), itsRho(rhs.itsRho), itsLevel(rhs.itsLevel),
			itsGridSize(rhs.itsGridSize), itsX(ARM_GP_VectorPtr(NULL)), itsQx(ARM_GP_VectorPtr(NULL)){} 
	
	ASSIGN_OPERATOR(ARM_NormalHestonDensityFunctor)
	virtual ~ARM_NormalHestonDensityFunctor() {}

	/// what it is for ...
	/// inline = perf issue
	inline double Quantile( double x,double fwd, double maturity ) const
	{
		ARM_THROW( ERR_INVALID_ARGUMENT, "ARM_NormalHestonDensityFunctor::Quantile not defined");
	}
	virtual double Proba(double strike, double fwd, double maturity) const {return 0.0;};
	
	inline virtual double Call_Option(double x, double fwd, double maturity ) const 
	{
		return NormalHeston(itsRho, itsKappa, itsTheta * SQR(itsLevel), itsVVol * itsLevel, itsV0 * SQR(itsLevel), itsFwd, x, maturity, 1);
	}

	/// Standard ARM object support
	virtual ARM_Object* Clone() const { return new ARM_NormalHestonDensityFunctor(*this); };
	virtual string toString(const string& indent="",const string& nextIndent="") const;

	// access
	double	GetFwd() const		{return itsFwd;};
	double	GetV0() const		{return itsV0;};
	double	GetKappa() const	{return itsKappa;};
	double	GetTheta() const	{return itsTheta;};
	double	GetVVol() const		{return itsVVol;};
	double	GetRho() const		{return itsRho;};
	double	GetLevel() const	{return itsLevel;};
};

/// --------------------------------
///	--- Generic density functor
/// --------------------------------
class ARM_GenDensityFunctor : public ARM_DensityFunctor
{
private:
	ARM_Security*		itsSecurity;
	ARM_Model*			itsModel;
	double				itsDecStrike;
	ARM_GP_Vector		itsStrikes;
	ARM_GP_Vector		itsProbas;
	bool				itsIsSO;
	bool				itsIsBS;
public: 
	/// constructor/destructor
	ARM_GenDensityFunctor( ARM_Security* pSec, ARM_Model* pModel, double decStrike, double minProba, double maxProba, bool isDirect = true); 
	ARM_GenDensityFunctor( const ARM_GenDensityFunctor& rhs ); 
	ASSIGN_OPERATOR(ARM_GenDensityFunctor)
	virtual ~ARM_GenDensityFunctor();

	/// what it is for ...
	virtual double Proba(double strike, double fwd, double maturity) const;
	virtual double Call_Option(double x, double fwd, double maturity ) const;
	virtual double Quantile( double x, double fwd, double maturity ) const;	
	
	double Quantile_Direct( double x, double fwd, double maturity ) const;
	double Quantile_Precomputed( double x, double fwd, double maturity ) const;
	void CheckTypes();
	
	inline double DecStrike() const				{return itsDecStrike;};
	inline ARM_Security* GetSecurity() const	{return itsSecurity;};
	inline ARM_Model* GetModel() const			{return itsModel;};

	inline bool IsSO() const					{return itsIsSO;};
	inline bool IsBS() const					{return itsIsBS;};

    /// Standard ARM object support
	virtual ARM_Object* Clone() const {return new ARM_GenDensityFunctor(*this);};
	virtual string toString(const string& indent="",const string& nextIndent="") const;
	virtual string ExportShortName() const { return "LDFUN";}

};

CC_END_NAMESPACE()

#endif
//-----------------------------------------------------------------------------
/*---- End of file ----*/



