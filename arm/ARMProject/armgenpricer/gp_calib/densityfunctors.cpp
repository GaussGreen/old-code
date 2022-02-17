/*!
 *
 * Copyright (c) CDC IXIS CM November 2003 Paris
 * 
 * \file densityfunctors.cpp
 *  \brief Density Functors for Markiv Functional Calibration
 * 
 *	\author  A. Schauly
 *	\version 1.0
 *	\date August 2005
 */

#include "gpcalib/densityfunctors.h"
#include "gpclosedforms\smile_shiftedlognormal.h"
#include "gpnumlib/gaussiananalytics.h"
#include "gpbase/gpmatrix.h"
#include "gpclosedforms/normal.h"

/// gpnumlib
#include "gpnumlib/numfunction.h"
#include "gpnumlib/dichotomy.h"
#include "gpnumlib/newtonraphson.h"

/// gpclosedforms
#include "gpclosedforms/inverse.h"
#include "gpclosedforms/sabrvanilla.h"
#include "gpclosedforms/vanilla_bs.h"
#include "gpclosedforms/vanilla_normal.h"
#include "gpclosedforms/heston_pricer.h"
#include "gpclosedforms/distribution_interface.h"
#include "gpclosedforms/nonparametric_spline.h"


/// gpinfra
//#include <inst/spreadoption.h>

const double DEFAULT_PRECISION				= 1.e-6;
const double CALL_SPREAD_SHIFT_FX			=	0.0001;


CC_BEGIN_NAMESPACE( ARM )

///$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
/// --- Class ARM_DensityFunctor (base)  ---
///$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$

//--------------------------------------------------------------------------------------------------------------
string ARM_DensityFunctor::toString(const string& indent,const string& nextIndent) const
{
	return "Density Functor";
}

ARM_DensityFunctor::ARM_DensityFunctor(bool isDirect) 
		:
	itsIsDirect(isDirect)
{}

ARM_DensityFunctor::ARM_DensityFunctor( const ARM_DensityFunctor& rhs ) 
		:
	ARM_RootObject(rhs),
		itsIsDirect(rhs.itsIsDirect)

{}

ARM_DensityFunctor& ARM_DensityFunctor:: operator = ( const ARM_DensityFunctor& rhs )
{
	if( this != &rhs ){
	ARM_RootObject::operator=( rhs );
	}
	return *this;
}

void ARM_DensityFunctor::MakeIncreasing(ARM_GP_VectorPtr& quantile)
{
	int i, size = quantile->size(), mid = size/2;

	for (i = mid; i < size-1; ++i)
		if ((*quantile)[i+1] < (*quantile)[i])
			(*quantile)[i+1] = (*quantile)[i];

	for (i = mid; i > 0; --i)
		if ((*quantile)[i-1] > (*quantile)[i])
			(*quantile)[i-1] = (*quantile)[i];
}

ARM_GP_VectorPtr ARM_DensityFunctor::Quantile( const ARM_GP_VectorPtr& x, double fwd, double maturity )
{
	std::vector<double>::iterator iter1, iter2, iterEnd;
	ARM_GP_VectorPtr result( new ARM_GP_Vector( x->size() ) );

	iter1 = result->begin();
	iter2 = x->begin(); 
	iterEnd = x->end();

	for( ; iter2 != iterEnd ; ++iter2, ++iter1 )
		(*iter1) = Quantile( *iter2, fwd, maturity );

	MakeIncreasing(result);
	return result;
}

///$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
/// --- Class ARM_ShiftedLNDensityFunctor ---
///$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$

//--------------------------------------------------------------------------------------------------------------

ARM_ShiftedLNDensityFunctor::ARM_ShiftedLNDensityFunctor(double vol, double shift, bool isDirect)
	: 
ARM_DensityFunctor(isDirect), 
	  itsVol(vol), 
	  itsShift(shift)
{
}

ARM_ShiftedLNDensityFunctor::ARM_ShiftedLNDensityFunctor( const ARM_ShiftedLNDensityFunctor& rhs )
	: 
ARM_DensityFunctor(rhs), 
	itsVol(rhs.itsVol), 
	itsShift(rhs.itsShift)
{
}

double ARM_ShiftedLNDensityFunctor::Quantile( double x, double fwd, double maturity ) const 
{ 
	return ShiftedLogNormal_Smile::inverse_distribution( fwd, x, maturity, itsVol, itsShift); 
}

double ARM_ShiftedLNDensityFunctor::Call_Option(double x, double fwd, double maturity ) const
{ 
	return ShiftedLogNormal_Smile::call_option( fwd, x, maturity, itsVol, itsShift);
}

//--------------------------------------------------------------------------------------------------------------
string ARM_ShiftedLNDensityFunctor::toString(const string& indent,const string& nextIndent) const
{	
	CC_Ostringstream os;
	os << indent << "ARM_ShiftedLNDensityFunctor ";
	os << "vol: " << itsVol  << " ";
	os << "shift: " << itsShift;
	os <<  CC_NS(std,endl);
	os << indent << CC_NS(std,endl);
	return os.str();
}


///$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
/// --- Class ARM_QDensityFunctor ---
///$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$

//--------------------------------------------------------------------------------------------------------------


ARM_QDensityFunctor::ARM_QDensityFunctor(double vol, double q,bool isDirect)
	: 
ARM_DensityFunctor(isDirect), 
	  itsVol(vol), 
	  itsQ(q)
{
}

ARM_QDensityFunctor::ARM_QDensityFunctor( const ARM_QDensityFunctor& rhs )
	: 
ARM_DensityFunctor(rhs), 
	itsVol(rhs.itsVol), 
	itsQ(rhs.itsQ)
{
}

double ARM_QDensityFunctor::Quantile( double x, double fwd, double maturity ) const
{
	double quantile;
	if (fabs(itsQ) > K_NEW_DOUBLE_TOL)
		quantile = fwd*(1.0+(exp(-0.5*itsQ*itsQ*itsVol*itsVol*maturity+itsQ*itsVol*sqrt(maturity)*NormalCDFInverse(x))-1.0)/itsQ);
	else
		quantile = fwd*(1.0+itsVol*sqrt(maturity)*NormalCDFInverse(x));

	return quantile;
}

double ARM_QDensityFunctor::Call_Option(double x, double fwd, double maturity ) const
{ 
	return QModelAnalytics::BSQFunction( fwd, x, itsVol,maturity,itsQ,1.0,1,itsQ);
}

//--------------------------------------------------------------------------------------------------------------
string ARM_QDensityFunctor::toString(const string& indent,const string& nextIndent) const
{	
	CC_Ostringstream os;
	os << indent << "ARM_QDensityFunctor ";
	os << "vol: " << itsVol  << " ";
	os << "q: " << itsQ;
	os <<  CC_NS(std,endl);
	os << indent << CC_NS(std,endl);
	return os.str();
}


///$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
/// --- Class ARM_MixtureDensityFunctor -----
///$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
ARM_MixtureDensityFunctor::ARM_MixtureDensityFunctor( const ARM_MixtureDensityFunctor& rhs )
	:
ARM_DensityFunctor(rhs),
	itsVol1(rhs.itsVol1),
	itsVol2(rhs.itsVol2),
	itsAlpha(rhs.itsAlpha),
	itsLambda(rhs.itsLambda)
{
}

ARM_MixtureDensityFunctor::ARM_MixtureDensityFunctor(bool isDirect)
	:
ARM_DensityFunctor(isDirect),
	itsVol1(0.00001),
	itsVol2(0.00001),
	itsAlpha(0.0),
	itsLambda(1.0)
{
}

ARM_MixtureDensityFunctor::ARM_MixtureDensityFunctor(double vol1 , 
		double vol2,
		double alpha,
		double lambda ,
		bool isDirect )
	:
ARM_DensityFunctor(isDirect),
	itsVol1(vol1),
	itsVol2(vol2),
	itsAlpha(alpha),
	itsLambda(lambda)
{
}

ARM_MixtureDensityFunctor::ARM_MixtureDensityFunctor(
		double fwd,
		double maturity,
		double volATM,
		double decVol,
		double alpha,
		double lambda,
		bool isDirect)
		:
ARM_DensityFunctor(isDirect),
	itsVol1(volATM-decVol),
	itsAlpha(alpha),
	itsLambda(lambda)
{
	double ATMPrice = BlackSholes_Formula(fwd,volATM*sqrt(maturity),1.0,fwd,1);
	class Vol2FunctionToInverse  : public ARM_GP::UnaryFunc<double, double> 
	{
	public:
		double m_Fwd;
		double m_Maturity;
		double m_Vol1;
		double m_Alpha;
		double m_Lambda;

		Vol2FunctionToInverse(double _fwd, double _maturity, double _vol1, double _alpha, double _lambda)
			: m_Fwd(_fwd), m_Maturity(_maturity), m_Vol1(_vol1), m_Alpha(_alpha), m_Lambda(_lambda) {}

		virtual double operator()(double vol2) const
		{		
			return ARM_MixtureDensityFunctor::Call(m_Fwd, m_Fwd, m_Maturity, m_Vol1, vol2, m_Alpha, m_Lambda);
		}
	};

	Vol2FunctionToInverse func(fwd,maturity,itsVol1,itsAlpha,itsLambda);
	//double guess = volATM+decVol;
	double sqrtT = sqrt(maturity);
	double stdDev1 = (volATM - decVol)*sqrtT;
	double stdDevATM = volATM*sqrtT;

	UnaryFuncWithNumDerivative<double> funcWithDeriv(func);
	
	T_DichotomySolver<UnaryFuncWithNumDerivative<double> > dicho_solver(funcWithDeriv, ATMPrice);
	
	dicho_solver.setInitialGuess(2*stdDevATM - stdDev1);
	
	double root = dicho_solver.Solve();
	
	T_SmoothNewtonRaphsonSolver<UnaryFuncWithNumDerivative<double> > nr_solver(funcWithDeriv, ATMPrice);
	
	nr_solver.setInitialGuess(root);
	
	itsVol2 = nr_solver.Solve();
	
}

double ARM_MixtureDensityFunctor::Quantile( double proba, double fwd, double maturity) const
{
	class DistributionToInverse : public DoubleToDoubleFunc 
	{
	private:
		ARM_MixtureDensityFunctor itsDensityFunctor;
	public: 
		double m_Fwd;
		double m_Maturity;

 		DistributionToInverse(const ARM_MixtureDensityFunctor& densityFunctor,
			double _fwd, 
			double _maturity)
			: m_Fwd(_fwd), 
			m_Maturity(_maturity), 
			itsDensityFunctor(densityFunctor)
			
		{}
		
		virtual double operator() (double K0)  const
		{
			if(K0<=m_Fwd*1e-9)
			{
				return 0.;
			}
			return itsDensityFunctor.Proba(K0,m_Fwd,m_Maturity);
		}
	};
	
	DistributionToInverse func(*this,fwd,maturity);

	double guess = fwd;
	return Inverse(func,Inverse::REAL)(proba,guess,fwd/5.,1e-12); //my modif 
}

double ARM_MixtureDensityFunctor::Proba(double strike, 
		double fwd, 
		double maturity) const
{
	double bondprice =1.0;
	double callPut = -1.0;

	/// We suppose that strike is non zero
	double sstrike = GetIsDirect() ? strike : 1.0/strike;

	double proba = 0.0;
	if ((1.0-itsLambda) >= K_NEW_DOUBLE_TOL)
	{
		double dig1=itsLambda*DigitalBlackSholes_Formula(fwd-itsAlpha, itsVol1*sqrt(maturity), bondprice, sstrike, callPut);
		double dig2=(1.0-itsLambda)*DigitalBlackSholes_Formula(fwd+itsAlpha*itsLambda/(1.0-itsLambda), itsVol2*sqrt(maturity), bondprice, sstrike, callPut);
		proba=dig1+dig2;
	}
	else
	{
		proba = DigitalBlackSholes_Formula(fwd + itsAlpha, itsVol1*sqrt(maturity), bondprice, sstrike, callPut);
	}

	if(!GetIsDirect())
	{
		double call  = Call(sstrike, fwd,  maturity, itsVol1, itsVol2, itsAlpha, itsLambda);
		proba = 1/fwd*( call + sstrike*(1-proba) );	
	}

	return proba;
};


double ARM_MixtureDensityFunctor::Call(double strike, 
		double fwd, 
		double maturity, 
		double vol1, 
		double vol2, 
		double alpha, 
		double lambda)
{
	if ((1.0-lambda) >= K_NEW_DOUBLE_TOL)
		return lambda*BlackSholes_Formula(fwd-alpha ,vol1*sqrt(maturity), 1.,strike,1)+(1.0-lambda)*BlackSholes_Formula(fwd+lambda/(1.0-lambda)*alpha,vol2*sqrt(maturity), 1.,strike,1);
	else
		return BlackSholes_Formula(fwd ,vol1*sqrt(maturity), 1.,strike,1);
}

double ARM_MixtureDensityFunctor::Call_Option(double strike, double fwd, double maturity) const
{
	return Call(strike, fwd, maturity, itsVol1, itsVol2, itsAlpha, itsLambda);
};

string ARM_MixtureDensityFunctor::toString(const string& indent,const string& nextIndent) const
{
	CC_Ostringstream os;
    os << "\n\n";
    os << indent << "ARM_MixtureDensityFunctor";
    os << indent << "----------------------------\n";
	os << indent << "vol1 : " << itsVol1 << endl;
	os << indent << "vol2 : " << itsVol2 << endl;
	os << indent << "alpha : " << itsAlpha << endl;
	os << indent << "lambda : " << itsLambda << endl;
	os << "\n\n";

	return os.str();	
}

///$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
/// --- Class ARM_HestonDensityFunctor -----
///$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$


double ARM_HestonDensityFunctor::Quantile( double proba, double fwd, double maturity) const
{
	class DistributionToInverse : public DoubleToDoubleFunc 
	{
	private:
		ARM_HestonDensityFunctor itsDensityFunctor;
	public: 
		double m_Fwd;
		double m_Maturity;

 		DistributionToInverse(const ARM_HestonDensityFunctor& densityFunctor,
			double _fwd, 
			double _maturity)
			: m_Fwd(_fwd), 
			m_Maturity(_maturity), 
			itsDensityFunctor(densityFunctor)
			
		{}
		
		virtual double operator() (double K0)  const
		{
			if(K0<=m_Fwd*1e-9)
			{
				return 0.;
			}
			return itsDensityFunctor.Proba(K0,m_Fwd,m_Maturity);
		}
	};
	
	DistributionToInverse func(*this,fwd,maturity);

	double guess = fwd;
	return Inverse(func,Inverse::REAL)(proba,guess,fwd/5.,1e-12); //my modif 
}

double ARM_HestonDensityFunctor::Proba(double strike,
									   double fwd,
									   double maturity) const
{
	/// We suppose that strike is non zero
	double sstrike = GetIsDirect() ? strike : 1.0/strike;

	double eps = CALL_SPREAD_SHIFT_FX;
	int callPut =-1;//only the put case interest us here
	double signed_eps = callPut*eps;
	ARM_MixteHestonOptionPricer pricerUp(
		maturity,
		fwd, 
		sstrike + signed_eps, 
		callPut, 
		itsSigma,
		itsV0, 
		itsKappa, 
		itsTheta, 
		itsRho, 
		itsNu, 
		itsShift,
		itsLevel);

	ARM_MixteHestonOptionPricer pricerDown(
		maturity,
		fwd, 
		sstrike - signed_eps, 
		callPut, 
		itsSigma,
		itsV0, 
		itsKappa, 
		itsTheta, 
		itsRho, 
		itsNu, 
		itsShift,
		itsLevel);

	double priceUp = pricerUp.price();
	double priceDown = pricerDown.price();

	double proba = (priceDown - priceUp)/(2*eps);

	if(!GetIsDirect())
	{
		double call  = Call(sstrike, fwd, maturity, itsV0, itsKappa, itsTheta, itsRho, itsNu, itsShift,	itsLevel, itsSigma);
		proba = 1/fwd*( call + sstrike*(1-proba) );	
	}
	
	return proba;
};

double ARM_HestonDensityFunctor::Call(double strike, double fwd, double maturity,
										double v0,double kappa, double theta, double rho, double nu, double shift, double level, double sigma)
{
	int callPut = 1;//only the call case interest us here
	ARM_MixteHestonOptionPricer pricer(
			maturity,
			fwd, 
			strike, 
			callPut, 
			v0, 
			kappa, 
			theta, 
			rho, 
			nu, 
			shift,
			level,
			sigma);
	return pricer.price();
}

double ARM_HestonDensityFunctor::Call_Option(double strike, double fwd, double maturity) const
{
	return Call(strike, fwd, maturity, itsV0, itsKappa, itsTheta, itsRho, itsNu, itsShift, itsLevel, itsSigma);
};

string ARM_HestonDensityFunctor::toString(const string& indent,const string& nextIndent) const
{
	CC_Ostringstream os;
    os << "\n\n";
    os << indent << "ARM_HestonDensityFunctor";
    os << indent << "----------------------------\n";
	os << indent << "v0 : " << itsV0 << endl;
	os << indent << "kappa : " << itsKappa << endl;
	os << indent << "theta : " << itsTheta << endl;
	os << indent << "rho : " << itsRho << endl;
	os << indent << "nu : " << itsNu << endl;
	os << indent << "shift : " << itsShift << endl;
	os << indent << "level : " << itsLevel << endl;
	os << indent << "sigma : " << itsSigma << endl;
	os << "\n\n";

	return os.str();	
}

///$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
/// --- Class ARM_SplineDensityFunctor ---
///$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$

ARM_SplineDensityFunctor::ARM_SplineDensityFunctor( const std::vector<double>& money, const std::vector<double>& vol, ARM_SmileViewer::MoneyType type)
: ARM_DensityFunctor(),
	itsSize(vol.size()),
	itsVar(0),
	itsMoney(0),
	itsDer2(0),
	itsType(type)
{
	itsVar.resize(itsSize);
	itsMoney.resize(itsSize);
	itsDer2.resize(itsSize);
	for (size_t k=0;k<itsSize;k++)
	{
		itsVar[k]=vol[k]*vol[k];
		itsMoney[k]=money[k];
	}
	if (itsSize>1)
		BuildSpline(1e30,1e30);
}

ARM_SplineDensityFunctor::ARM_SplineDensityFunctor( const ARM_SplineDensityFunctor& rhs )
: ARM_DensityFunctor(rhs),
	itsSize(rhs.itsSize),
	itsVar(rhs.itsVar), 
	itsMoney(rhs.itsMoney),
	itsDer2(rhs.itsDer2),
	itsType(rhs.itsType)
{
}

void ARM_SplineDensityFunctor::BuildSpline(double yp1, double ypn)
{
	//cubicspline_precompute(&itsMoney,&itsVar,yp1,ypn,itsDer2);
	
	/*
	int i,k;
	double p,qn,sig,un;
	std::vector<double> u(itsSize-1);
	if (yp1 > 0.99e30)
		itsDer2[0]=u[0]=0.0;
	else {
		itsDer2[0] = -0.5;
		u[0]=(3.0/(itsMoney[1]-itsMoney[0]))*((itsVar[1]-itsVar[0])/(itsMoney[1]-itsMoney[0])-yp1);
	}
	for (i=1;i<itsSize-1;i++) {
		sig=(itsMoney[i]-itsMoney[i-1])/(itsMoney[i+1]-itsMoney[i-1]);
		p=sig*itsDer2[i-1]+2.0;
		itsDer2[i]=(sig-1.0)/p;
		u[i]=(itsVar[i+1]-itsVar[i])/(itsMoney[i+1]-itsMoney[i]) - (itsVar[i]-itsVar[i-1])/(itsMoney[i]-itsMoney[i-1]);
		u[i]=(6.0*u[i]/(itsMoney[i+1]-itsMoney[i-1])-sig*u[i-1])/p;
	}
	if (ypn > 0.99e30)
		qn=un=0.0;
	else {
		qn=0.5;
		un=(3.0/(itsMoney[itsSize-1]-itsMoney[itsSize-2]))*(ypn-(itsVar[itsSize-1]-itsVar[itsSize-2])/(itsMoney[itsSize-1]-itsMoney[itsSize-2]));
	}
	itsDer2[itsSize-1]=(un-qn*u[itsSize-2])/(qn*itsDer2[itsSize-2]+1.0);
	for (k=itsSize-2;k>=0;k--)
		itsDer2[k]=itsDer2[k]*itsDer2[k+1]+u[k];
	*/
}

double ARM_SplineDensityFunctor::Quantile( double proba, double fwd, double mat) const
{
	class DistributionToInverse : public DoubleToDoubleFunc 
	{
	public: 
		double m_fwd;
		std::vector<double> m_pMoney;
		std::vector<double> m_pVar;
		std::vector<double> m_pDer2;
		size_t		   m_type;
		DistributionToInverse(double fwd,const std::vector<double>& pMoney, 
			const std::vector<double>& pVar, 
			const std::vector<double>& pDer2,
			size_t type):
		m_fwd(fwd),m_pVar(pVar),m_pMoney(pMoney),m_pDer2(pDer2),m_type(type){};
		
		virtual double operator() (double K0)  const
		{
			return ARM_SplineDensityFunctor::Proba(K0,m_fwd,m_pMoney,m_pVar,m_pDer2,m_type);
		}
	};
	
	DistributionToInverse func(fwd,itsMoney,itsVar,itsDer2,itsType);

	double guess = fwd;
	
	if ( itsType == ARM_SmileViewer::BLACK)
		return Inverse(func,Inverse::ALWAYSPOSITIVE)(proba,guess,fwd/15.0,1e-12); 
	else
		return Inverse(func,Inverse::REAL)(proba,guess,fwd/15.0,1e-12); 
}

double ARM_SplineDensityFunctor::Proba(double strike, double fwd, const std::vector<double>& money, const std::vector<double>& var, const std::vector<double>& der2, size_t type)
{
	double eps = 1e-6;
	double volplus = ComputeVol(strike+eps,fwd,money,var,der2,type);
	double volminus = ComputeVol(strike-eps,fwd,money,var,der2,type);
	if (type==0)
		return 1. + (VanillaOption_N(fwd ,volplus,strike+eps,1.,1)-VanillaOption_N(fwd ,volminus,strike-eps,1.,1))/2./eps;
	else
		return 1. + (BlackSholes_Formula(fwd ,volplus, 1.,strike+eps,1)-BlackSholes_Formula(fwd ,volminus, 1.,strike-eps,1))/2./eps;

	double vol = ComputeVol(strike,fwd,money,var,der2,type);
	double dvol = ComputeDVol(strike,fwd,money,var,der2,type);

	double volatm = ComputeVol(0,money,var,der2);
	double dmdk = (type==0?-1./volatm:-1./fwd/volatm);
	dvol=0.5*dvol/vol*dmdk;

	if(type==0)
		return ProbaGauss(strike,fwd,vol,dvol);
	else
		return ProbaBlack(strike,fwd,vol,dvol);
};

double ARM_SplineDensityFunctor::ComputeVol(double strike, double fwd, const std::vector<double>& money, const std::vector<double>& var, const std::vector<double>& der2, size_t type)
{
	
	double volatm1 = sqrt(cubicspline_inter(money,var,der2,0.));
	double moneyness1 = type == 0 ? (strike-fwd)/volatm1 : log(strike/fwd)/volatm1;
	return sqrt(cubicspline_inter(money,var,der2,moneyness1));
	
	/*
	double volatm=ComputeVol(0,money,var,der2);
	double moneyness = (type==0?(strike-fwd)/volatm:log(strike/fwd)/volatm);
	
	int n = money.size();
	double res;
	if (moneyness>=money[0] && moneyness<=money[n-1])
		res=ComputeVol(moneyness,money,var,der2);
	else {
		if (moneyness<money[0]){
			double leftder = (var[1]-var[0])/(money[1]-money[0])-(money[1]-money[0])/6.*der2[1];
			double leftval = var[0];
			if (leftder>0){
				res=leftval*(1.+tanh(leftder/leftval*(moneyness-money[0])));
			}else{
				res=leftval+leftder*(moneyness-money[0]);
			}
		}
		else{
			double rightder = (var[n-1]-var[n-2])/(money[n-1]-money[n-2])+(money[n-1]-money[n-2])/6.*der2[n-2];
			double rightval = var[n-1];
			if (rightder<0){
				res=rightval*(1+tanh(rightder/rightval*(moneyness-money[n-1])));
			}else{
				res=rightval+rightder*(moneyness-money[n-1]);
			}
		}
		res=sqrt(res);
	}
	return res;
	*/
};

double ARM_SplineDensityFunctor::ComputeDVol(double strike, double fwd, const std::vector<double>& money, const std::vector<double>& var, const std::vector<double>& der2, size_t type)
{
	
	double volatm = cubicspline_inter(money,var,der2,0.);
	double moneyness = type == 0 ? (strike-fwd)/volatm : log(strike/fwd)/volatm;
	return cubicspline_interder(money,var,der2,moneyness);
	
	/*
	double volatm=ComputeVol(0,money,var,der2);
	double moneyness = (type==0?(strike-fwd)/volatm:log(strike/fwd)/volatm);
	
	int n = money.size();
	double res;
	if (moneyness>=money[0] && moneyness<=money[n-1])
		res=ComputeDVol(moneyness,money,var,der2);
	else {
		if (moneyness<money[0]){
			double leftder = (var[1]-var[0])/(money[1]-money[0])-(money[1]-money[0])/6.*der2[1];
			double leftval = var[0];
			if (leftder>0){
				res=cosh(leftder/leftval*(moneyness-money[0]));
				res=leftder/res/res;
			}else{
				res=leftder;
			}
		}
		else{
			double rightder = (var[n-1]-var[n-2])/(money[n-1]-money[n-2])+(money[n-1]-money[n-2])/6.*der2[n-2];
			double rightval = var[n-1];
			if (rightder<0){
				res=cosh(rightder/rightval*(moneyness-money[n-1]));
				res=rightder/res/res;
			}else{
				res=rightder;
			}
		}
	}
	return res;
	*/
	
};

double ARM_SplineDensityFunctor::ProbaGauss(double strike, double fwd, double vol, double dvol)
{
	return 1. - VanillaDigitalOption_N(fwd,vol,strike,1.,1) + dvol * VegaVanillaOption_N(fwd,vol,strike,1.,1);
};

double ARM_SplineDensityFunctor::ProbaBlack(double strike, double fwd, double vol, double dvol)
{
	return 1. - DigitalBlackSholes_Formula(fwd,vol,1.,strike,1.) + dvol * BlackSholes_Derivative_2(fwd,vol,1.,strike,1.);
};

double ARM_SplineDensityFunctor::Call_Option(double strike, double fwd, double maturity) const
{
	double vol = ComputeVol(strike,fwd,itsMoney,itsVar,itsDer2,itsType);
	if (itsType==0)
		return VanillaOption_N(fwd ,vol, strike, 1., 1);
	else
		return BlackSholes_Formula(fwd ,vol, 1.,strike,1);
};

double ARM_SplineDensityFunctor::ComputeVol(double x, const std::vector<double>& money, const std::vector<double>& var, const std::vector<double>& der2)
{
	return sqrt(cubicspline_inter(money,var,der2,x));

	/*
	int klo,khi,k;
	double h,b,a;
	int n=money.size();
	klo=0; 
	khi=n-1;
	while (khi-klo > 1) {
		k=(khi+klo) >> 1;
		if (money[k] > x) khi=k;
		else klo=k;
	}
	h=money[khi]-money[klo];
	if (h == 0.0)
		ARM_THROW( ERR_INVALID_ARGUMENT, "ARM_SplineDensityFunctor::ComputeVol: bad input");
	a=(money[khi]-x)/h;
	b=(x-money[klo])/h;
	double y=a*var[klo]+b*var[khi]+((a*a*a-a)*der2[klo]+(b*b*b-b)*der2[khi])*(h*h)/6.0;
	return sqrt(y);
	*/
	
}

double ARM_SplineDensityFunctor::ComputeDVol(double x, const std::vector<double>& money, const std::vector<double>& var, const std::vector<double>& der2)
{
	return cubicspline_interder(money,var,der2,x);

	/*
	int klo,khi,k;
	double h,b,a;
	int n=money.size();
	klo=0; 
	khi=n-1;
	while (khi-klo > 1) {
		k=(khi+klo) >> 1;
		if (money[k] > x) khi=k;
		else klo=k;
	}
	h=money[khi]-money[klo];
	if (h == 0.0)
		ARM_THROW( ERR_INVALID_ARGUMENT, "ARM_SplineDensityFunctor::ComputeDVol: bad input");
	a=(money[khi]-x)/h;
	b=(x-money[klo])/h;
	double y=1./h*(-var[klo]+var[khi]+(-(3*a*a-1)*der2[klo]+(3*b*b-1)*der2[khi])*(h*h)/6.0);
	return y;
	*/
}

string ARM_SplineDensityFunctor::toString(const string& indent,const string& nextIndent) const
{
	CC_Ostringstream os;
    os << "\n\n";
    os << indent << "ARM_SplineDensityFunctor";
    os << indent << "----------------------------\n";
	os << "\n\n";

	size_t size=itsMoney.size();
	if (size>0)
		for (size_t k=0;k<size;k++)
			os << itsMoney[k] << "\t"<< itsVar[k] << "\t"<< itsDer2[k] << "\n";

	return os.str();	
}

///$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
/// --- Class ARM_NoVolDensityFunctor  -----
///$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$

//--------------------------------------------------------------------------------------------------------------
string ARM_NoVolDensityFunctor::toString(const string& indent,const string& nextIndent) const
{
	return "ARM_NoVolDensityFunctor";
}

///$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
/// --- Class ARM_SABRDensityFunctor -----
///$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$

ARM_SABRDensityFunctor::ARM_SABRDensityFunctor(bool isDirect)
	:	
ARM_DensityFunctor(isDirect),
	itsAlpha(0), 
	itsBeta(0), 
	itsRho(0),
	itsNu(0), 
	itsSabrType (0),
	itsGridSize(0), 
	itsX(ARM_GP_VectorPtr(NULL)), 
	itsQx(ARM_GP_VectorPtr(NULL))
{
} 

ARM_SABRDensityFunctor::ARM_SABRDensityFunctor(double alpha,
	double beta, 
	double rho, 
	double nu, 
	int sabrType,
	size_t gridSize)
	:
ARM_DensityFunctor(), 
	itsAlpha(alpha), 
	itsBeta(beta),
	itsRho(rho), 
	itsNu(nu),
	itsSabrType (sabrType),
	itsGridSize(gridSize),
	itsX(ARM_GP_VectorPtr(NULL)), 
	itsQx(ARM_GP_VectorPtr(NULL))
{
} 

ARM_SABRDensityFunctor::ARM_SABRDensityFunctor( const ARM_SABRDensityFunctor& rhs )
	:	
ARM_DensityFunctor(rhs), 
	itsAlpha(rhs.itsAlpha), 
	itsBeta(rhs.itsBeta),
	itsRho(rhs.itsRho), 
	itsNu(rhs.itsNu),
	itsSabrType(rhs.itsSabrType),
	itsGridSize(rhs.itsGridSize), 
	itsX(ARM_GP_VectorPtr(NULL)), 
	itsQx(ARM_GP_VectorPtr(NULL))
{
} 
 
double ARM_SABRDensityFunctor::Quantile( double x, double fwd, double maturity ) const
{	
	return SABR_smile::inverse_distribution(fwd, 
		x, 
		maturity, 
		itsAlpha, 
		itsBeta, 
		itsRho, 
		itsNu, 
		itsSabrType, 120,fwd/4.,1.5,fwd/2.);
		
}

double ARM_SABRDensityFunctor::Proba(double strike,
									   double fwd,
									   double maturity) const
{
	/// We suppose that strike is non zero
	double sstrike = GetIsDirect() ? strike : 1.0/strike;

	double proba = 0;
	if(strike<=fwd*1e-9)
		proba = 0;
	else
	{
		ArgumentList arg(fwd,0.,maturity,itsAlpha,itsBeta,itsRho,itsNu,K_CALL,itsSabrType,120/*not used*/);
		arg.set_nth(1,strike);
		// ArgumentList a(f0,K0,tex0,alpha0,beta0,rho0,nu0,K_CALL,flag0,nbsteps0);
		double s=ARM_CF_SABR_VanillaOption_Formula::specific_shift(1);
		proba = 1.+ARM_CF_SABR_VanillaOption_Formula::value(1,arg,s);
	}
	if(!GetIsDirect())
	{
		double call  = Call_Option(sstrike, fwd, maturity);
		proba = 1/fwd*( call + sstrike*(1-proba) );	
	}
	return proba;
};

//--------------------------------------------------------------------------------------------------------------
ARM_GP_VectorPtr ARM_SABRDensityFunctor::Quantile( const ARM_GP_VectorPtr& x, double fwd, double maturity )
{
	size_t i, j;
	ARM_GP_VectorPtr result( new ARM_GP_Vector( x->size() ) );

	/// Convention 
	/// --> if N = 0, exact computation of quantiles
	if (itsGridSize == 0)
	{
		std::vector<double>::iterator iter1, iter2, iterEnd;
		iter1 = result->begin();
		iter2 = x->begin(); 
		iterEnd = x->end();

		for( ; iter2 != iterEnd ; ++iter2, ++iter1 )
			(*iter1) = Quantile( *iter2, fwd, maturity );

		MakeIncreasing(result);
		return result;
	}

	if (itsGridSize > 2000)
		ARM_THROW( ERR_INVALID_ARGUMENT, ARM_USERNAME + ": ARM_SABRDensityFunctor::Quantile : too big grid size");

	/// Quantiles have not yet been computed
	/// NB : specific x sampling
	if ( itsX.IsNull() )
	{
		double vol  = pow(fwd, itsBeta - 1.) * itsAlpha;
		double xmin = (*x)[0];
		double xmax = (*x)[x->size()-1];
		double Kmin = fwd * exp( - 0.5 * vol * vol * maturity + vol * sqrt( maturity ) * NormalCDFInverse(xmin) ) ;
		double Kmax = fwd * exp( - 0.5 * vol * vol * maturity + vol * sqrt( maturity ) * NormalCDFInverse(xmax) ) ;
		double dK   = (Kmax - Kmin) / (itsGridSize - 1);
		
		itsX  = ARM_GP_VectorPtr( new ARM_GP_Vector(itsGridSize) );
		itsQx = ARM_GP_VectorPtr( new ARM_GP_Vector(itsGridSize) );

		double K, d;
		double volsqrt = vol * sqrt(maturity);

		for (i=0; i<itsGridSize; i++)
		{
			K = Kmin + i * dK;
			d = log(K/fwd) / volsqrt + 0.5 * volsqrt;
			(*itsX)[i]  = NormalCDF(d);
			(*itsQx)[i] = Quantile( (*itsX)[i], fwd, maturity );
		}

		/// to be sure that x[0] and x[n-1] are in the sample
		(*itsX)[0]	    = xmin;
		(*itsX)[itsGridSize-1] = xmax;
		if ( (*itsX)[0]>(*itsX)[1] || (*itsX)[itsGridSize-1]<(*itsX)[itsGridSize-2] )
			ARM_THROW( ERR_INVALID_ARGUMENT, ARM_USERNAME + ": ARM_SABRDensityFunctor::Quantile : sampling problem");
	}


	/// *** Linear Interpolation ***
	i = 0;
	size_t size = x->size() - 1;
	while (/*(i < size-1) &&*/ ((*x)[i] <= (*itsX)[0]))
	{	
		(*result)[i] = (*itsQx)[0]; 
		i++;
	}
	size_t firstIdx = i;

	i = x->size() - 1;
	while (/*(i > 0) && */((*x)[i] >= (*itsX)[itsGridSize-1]))
	{	
		(*result)[i] = (*itsQx)[itsGridSize-1]; 
		i--;
	}
	size_t lastIdx = i;
	
	////
	/// std linear interp
	j = 0;
	double curX, nextX, prevX;
	for( i = firstIdx; i<=lastIdx; ++i)
	{
		curX = (*x)[i];

		while ( (*itsX)[j] < curX)
			++j;
		
		nextX = (*itsX)[j];
		prevX = (*itsX)[j-1];

		(*result)[i] =   (  (curX - prevX) * (*itsQx)[j]
						  + (nextX - curX) * (*itsQx)[j-1] ) / (nextX - prevX);
		
	}

	MakeIncreasing(result);
	
	return result;
}


//--------------------------------------------------------------------------------------------------------------
string ARM_SABRDensityFunctor::toString(const string& indent,const string& nextIndent) const
{
	CC_Ostringstream os;
	os << indent << "ARM_SABRDensityFunctor ";
	os << "alpha : " << itsAlpha << " "; 
	os << "rho : " << itsRho << " "; 
	os << "Nu : " << itsNu << " "; 
	os << "Beta : " << itsBeta; 
	os <<  CC_NS(std,endl);
	os << indent << CC_NS(std,endl);
	return os.str();
}

double ARM_BiSABRDensityFunctor::Call_Option(double x, double fwd, double maturity ) const 
{
	return Export_BiSABR_SpreadOption(itsFwd1, 
		itsAlpha1, 
		itsBeta1, 
		itsRho1,
		itsNu1,
		itsFwd2, 
		itsAlpha2,
		itsBeta2, 
		itsRho2, 
		itsNu2,
		x, 
		maturity, 
		1, 
		itsRhoS1S2, 
		itsRhoV1V2, 
		itsRhoS1V2, 
		itsRhoS2V1);
}
double ARM_BiSABRDensityFunctor::Quantile( double x,double fwd, double maturity ) const
{	
	return Export_BiSABR_Quantile(itsFwd1,
		itsAlpha1,
		itsBeta1, 
		itsRho1, 
		itsNu1, 
		itsFwd2,
		itsAlpha2, 
		itsBeta2, 
		itsRho2, 
		itsNu2,
		itsRhoS1S2, 
		itsRhoV1V2,
		itsRhoS1V2, 
		itsRhoS2V1, 
		x, 
		maturity,
		ARM_CF_SABR_ImplicitVol_Formula::SABR_IMPLNVOL );
}

ARM_GP_VectorPtr ARM_BiSABRDensityFunctor::Quantile( const ARM_GP_VectorPtr& x, double fwd, double maturity )
{
	size_t i, j;
	ARM_GP_VectorPtr result( new ARM_GP_Vector( x->size() ) );

	/// Convention 
	/// --> if N = 0, exact computation of quantiles
	if (itsGridSize == 0)
	{
		std::vector<double>::iterator iter1, iter2, iterEnd;
		iter1 = result->begin();
		iter2 = x->begin(); 
		iterEnd = x->end();

		for( ; iter2 != iterEnd ; ++iter2, ++iter1 )
			(*iter1) = Quantile( *iter2, fwd, maturity );

		MakeIncreasing(result);
		return result;
	}

	if (itsGridSize > 2000)
		ARM_THROW( ERR_INVALID_ARGUMENT, ARM_USERNAME + ": ARM_BiSABRDensityFunctor::Quantile : too big grid size");

	/// Quantiles have not yet been computed
	/// NB : specific x sampling
	if ( itsX.IsNull() )
	{
		ARM_THROW( ERR_INVALID_ARGUMENT, ARM_USERNAME + ": ARM_BiSABRDensityFunctor::Quantile : need x sampling");

	}


	/// *** Linear Interpolation ***
	i = 0;
	size_t size = x->size() - 1;
	while (/*(i < size-1) &&*/ ((*x)[i] <= (*itsX)[0]))
	{	
		(*result)[i] = (*itsQx)[0]; 
		i++;
	}
	size_t firstIdx = i;

	i = x->size() - 1;
	while (/*(i > 0) && */((*x)[i] >= (*itsX)[itsGridSize-1]))
	{	
		(*result)[i] = (*itsQx)[itsGridSize-1]; 
		i--;
	}
	size_t lastIdx = i;
	
	////
	/// std linear interp
	j = 0;
	double curX, nextX, prevX;
	for( i = firstIdx; i<=lastIdx; ++i)
	{
		curX = (*x)[i];

		while ( (*itsX)[j] < curX)
			++j;
		
		nextX = (*itsX)[j];
		prevX = (*itsX)[j-1];

		(*result)[i] =   (  (curX - prevX) * (*itsQx)[j]
						  + (nextX - curX) * (*itsQx)[j-1] ) / (nextX - prevX);
		
	}

	MakeIncreasing(result);
	
	return result;
}

string ARM_BiSABRDensityFunctor::toString(const string& indent,const string& nextIndent) const
{
	CC_Ostringstream os;
	os << indent << "ARM_SABRDensityFunctor ";
	os << "1st underlying : \n";
	os << "alpha : " << itsAlpha1 << " "; 
	os << "rho1 : " << itsRho1 << " "; 
	os << "Nu : " << itsNu1 << " "; 
	os << "Beta : " << itsBeta1 << "\n\n"; 
	os << "2nd underlying : \n";
	os << "alpha : " << itsAlpha2 << " "; 
	os << "rho1 : " << itsRho2 << " "; 
	os << "Nu : " << itsNu2 << " "; 
	os << "Beta : " << itsBeta2 << "\n\n"; 
	os << "cross correlation : \n";
	os << "rho S1,S2 : " << itsRhoS1S2 << " ";
	os << "rho V1,V2 : " << itsRhoV1V2 << " ";
	os << "rho S1,V2 : " << itsRhoS1V2 << " ";
	os << "rho S2,V1 : " << itsRhoS2V1 << " ";

	os <<  CC_NS(std,endl);
	os << indent << CC_NS(std,endl);
	return os.str();
}

ARM_ShiftedLNDensityFunctor* ARM_SABRDensityFunctor::toShiftedLN_2strikes(double fwd, double maturity, double strike1, double strike2, double shiftinf, double shiftsup) const
{
	double mkt1			= Call_Option(strike1, fwd, maturity);
	double mkt2			= Call_Option(strike2, fwd, maturity);
	
	SLNDensity_Approx func(fwd,maturity,strike1,mkt1,strike2,mkt2);

	UnaryFuncWithNumDerivative<double> funcWithDev(func);

	T_DichotomySolver< UnaryFuncWithNumDerivative<double> > solver(funcWithDev,0,DEFAULT_PRECISION,DEFAULT_PRECISION);
	solver.setInitialGuess(0.,shiftinf,shiftsup);

	double result_shift = solver.Solve();
	double result_vol   = func.GetVol(result_shift);
	return new ARM_ShiftedLNDensityFunctor( result_vol , result_shift );
}
	
ARM_ShiftedLNDensityFunctor* ARM_SABRDensityFunctor::toShiftedLN_1strike1shift(double fwd, double maturity, double strike1, double shift) const
{
	double mkt1			= Call_Option(strike1, fwd, maturity);
	
	SLNDensity_Approx func(fwd,maturity,strike1,mkt1,strike1,mkt1);

	double result_shift = shift;
	double result_vol   = func.GetVol(result_shift);
	return new ARM_ShiftedLNDensityFunctor( result_vol , result_shift );
}

//-----------------------------------------------------------------------------------------
ARM_SmileViewer::ARM_SmileViewer(const std::vector<double>& money,MoneyType type,const std::vector<double>& strikes)
:	itsMoney(money),
	itsStrike(0),
	itsVol(0),
	itsType(type),
	itsExtraStrikes(strikes)
{
}

ARM_SmileViewer::ARM_SmileViewer(const ARM_SmileViewer& rhs)
:	itsMoney(rhs.itsMoney),
	itsStrike(rhs.itsStrike),
	itsVol(rhs.itsVol),
	itsType(rhs.itsType)
{
}

ARM_SmileViewer::~ARM_SmileViewer()
{
}

ARM_Object* ARM_SmileViewer::Clone() const
{
	return new ARM_SmileViewer( *this );
}

string ARM_SmileViewer::toString(const string& indent, const string& nextIndent) const
{
	CC_Ostringstream os;
    os << "\n\n";
    os << indent << "Smile Viewer\n";
    return os.str();
}

void ARM_SmileViewer::View(char* id, FILE* ficOut) const
{
	FILE* fOut;
    char fOutName[200];
	
    if ( ficOut == NULL )
    {
		ARM_GetViewFile(ARM_VIEW_FILE, id, fOutName);
		(void) unlink(fOutName);
		fOut = fopen(fOutName, "w");
    }
    else
		fOut = ficOut;   


	CC_Ostringstream os;
	size_t n=itsMoney.size();
	
	os << "\t" << "Moneyness " << "\t" << "Strike    " << "\t" << "Vol       " << "\n";
	for (size_t k=0;k<n;k++)
	{
		os << "\t" << CC_NS(std,fixed) << itsMoney[k] << "\t" << itsStrike[k] << "\t" << itsVol[k] << "\n";
	}
	os <<"\n";

	fprintf(fOut, "\n INFOS FOR SMILE VIEWER\n\n" );
	fprintf(fOut, "%s", (itsType==GAUSS?"VolType : GAUSS":"VolType : BLACK"));
	fprintf(fOut, "\n\n");
	fprintf(fOut, "%s", os.str().c_str());
	fprintf(fOut, " ======> END OF INFOS SMILE VIEWER <================== \n\n" );
    
    if ( ficOut == NULL )
       fclose(fOut);
}

void ARM_SmileViewer::Compute(ARM_Security* pSec, ARM_Model* pMod)
{
	//ARM_BSModel* BSModel = dynamic_cast< ARM_BSModel* >(pMod);
	//if (BSModel)
	//{
	//	//ARM_ZeroCurve* pCurve = pMod->GetZeroCurve();
	//	double asOfDate	  = 0;//pCurve->GetAsOfDateJul();
	//		
	//	ARM_SpreadOption* spreadOption = dynamic_cast<ARM_SpreadOption*>(pSec);
	//	if (spreadOption)
	//	{
	//		spreadOption->SetModel(BSModel);
	//		double price = spreadOption->ComputePrice();
	//		double fwd1 = (*(spreadOption->GetSpreadLeg()->GetFirstLeg()->GetFwdRates()))[0];
	//		double fwd2 = (*(spreadOption->GetSpreadLeg()->GetSecondLeg()->GetFwdRates()))[0];
	//		double weight1 = spreadOption->GetWeight1();
	//		double weight2 = spreadOption->GetWeight2();

	//		//if (weight1!=1 || weight2!=1)
	//		//	ARM_THROW( ERR_INVALID_ARGUMENT, "ARM_SmileViewer : weights should be 1." );

	//		int callPut = spreadOption->IsCap()-spreadOption->IsFloor();
	//		double pay = spreadOption->GetSpreadLeg()->GetFirstLeg()->GetPaymentDates()->Elt(0);
	//		double zcPay = pCurve->DiscountPrice((pay-asOfDate)/K_YEAR_LEN);
	//		
	//		double theta = spreadOption->GetSpreadLeg()->GetFirstLeg()->GetInterestTerms()->Elt(0);
	//		double notio = spreadOption->GetSpreadLeg()->GetFirstLeg()->GetAmount()->CptReferenceValue(pay);
	//		
	//		double mat = (spreadOption->GetSpreadLeg()->GetFirstLeg()->GetResetDates()->Elt(0)-asOfDate)/K_YEAR_LEN;
	//		double fwd = (weight2*fwd2 - weight1*fwd1)/100.;
	//		double volinit = 0.003;
	//		double vol;
	//		double strike;

	//		strike = fwd;
	//		spreadOption->SetStrike(strike*100);
	//		price = spreadOption->ComputePrice();
	//		price/= theta * notio * zcPay;
	//		vol = VanillaImpliedVol_N (fwd, price, strike, mat, callPut, &volinit) * sqrt(mat);
	//		double volatm = vol;
	//		for (size_t i=0;i<itsExtraStrikes.size();i++)
	//			itsMoney.push_back((itsExtraStrikes[i]-fwd)/volatm);
	//		
	//		itsMoney.sort();
	//		itsStrike.resize(itsMoney.size());
	//		itsVol.resize(itsMoney.size());
	//		for (i=0;i<itsMoney.size();i++)
	//		{
	//			strike = fwd+itsMoney[i]*volatm;
	//			spreadOption->SetStrike(strike*100);
	//			price = spreadOption->ComputePrice();
	//			price/= theta * notio * zcPay;
	//			vol = VanillaImpliedVol_N (fwd, price, strike, mat, callPut, &volinit) * sqrt(mat);
	//			itsStrike[i]=strike;
	//			itsVol[i]=vol;
	//		}
	//		
	//	}
	//	else
	//	{
	//		ARM_THROW( ERR_INVALID_ARGUMENT, "ARM_SmileViewer : only SO supported" );
	//	}
	//}
	//else
	//	ARM_THROW( ERR_INVALID_ARGUMENT, "ARM_SmileViewer::Compute : check model type" );
}
	
string ARM_NormalHestonDensityFunctor::toString(const string& indent,const string& nextIndent) const
{
	CC_Ostringstream os;
	
	os << indent << "ARM_NormalHestonDensityFunctor ";

	os << "fwd : " << itsFwd << "\n";
	os << "v0 : " << itsV0 << "\n";
	os << "kappa : " << itsKappa << "\n";
	os << "theta : " << itsTheta << "\n";
	os << "vvol : " << itsVVol << "\n";
	os << "rho : " << itsRho << "\n";
	os << "level : " << itsLevel << " ";
	os <<  CC_NS(std,endl);
	os << indent << CC_NS(std,endl);
	
	return os.str();
}

ARM_GenDensityFunctor::ARM_GenDensityFunctor( ARM_Security* pSec, ARM_Model* pModel, double decStrike, double minProba, double maxProba, bool isDirect)
: ARM_DensityFunctor(isDirect),itsSecurity(pSec),itsModel(pModel),itsDecStrike(decStrike),itsStrikes(0),itsProbas(0),itsIsBS(false),itsIsSO(false)
{
	CheckTypes();
	if (!GetIsDirect())
	{
		double minStrike = Quantile_Direct(minProba,0.01,0.001);
		double maxStrike = Quantile_Direct(maxProba,0.01,0.001);
		int nbStrikes = (int)((maxStrike-minStrike)/decStrike);
		double dx = (maxStrike-minStrike)/(nbStrikes-1.);

		itsStrikes.resize(nbStrikes);
		itsProbas.resize(nbStrikes);
		int i;
		for (i=0;i<nbStrikes;i++)
			itsStrikes[i] =   minStrike + i*dx;
		double csup,cinf = Call_Option(minStrike-dx,0.,0.);
		for (i=0;i<nbStrikes-1;i++)
		{
			csup = Call_Option(itsStrikes[i+1],0.,0.);
			itsProbas[i]=1.-(cinf-csup)/2./dx;
			cinf = Call_Option(itsStrikes[i],0.,0.);
		}
		csup = Call_Option(itsStrikes[i]+dx,0.,0.);
		itsProbas[i]=1.-(cinf-csup)/2./dx;
	}
}

ARM_GenDensityFunctor::~ARM_GenDensityFunctor()
{
	delete itsSecurity;
	delete itsModel;
}

ARM_GenDensityFunctor::ARM_GenDensityFunctor(const ARM_GenDensityFunctor& rhs)
:	ARM_DensityFunctor(rhs),
	//itsSecurity((ARM_Security*)rhs.itsSecurity->Clone()),
	//itsModel((ARM_Model*)rhs.itsSecurity->Clone()),
	itsStrikes(rhs.itsStrikes),
	itsProbas(rhs.itsProbas)
{
}

void ARM_GenDensityFunctor::CheckTypes()
{
	// FIXME
	/*ARM_BSModel* BSModel = dynamic_cast< ARM_BSModel* >(GetModel());
	if (BSModel)
	{
		itsIsBS = true;
		ARM_SpreadOption* spreadOption = dynamic_cast<ARM_SpreadOption*>(GetSecurity());
		if (spreadOption)
			itsIsSO = true;
		else
		{
			ARM_THROW( ERR_INVALID_ARGUMENT, "ARM_GenDensityFunctor::CheckTypes : only SO right now" );
		}
	}
	else
		ARM_THROW( ERR_INVALID_ARGUMENT, "ARM_GenDensityFunctor::CheckTypes : check model type" );*/
}
	
string ARM_GenDensityFunctor::toString(const string& indent,const string& nextIndent) const
{
	CC_Ostringstream os;
	
	os << indent << "ARM_GenDensityFunctor" << endl;
	if (IsSO())
		os << indent << "for a spread option" << endl;
	if (IsBS())
		os << indent << "computed with a BSModel" << endl << endl;
	if (GetIsDirect())
		os << indent << "with direct computation of quantiles (quite long...)" << endl;
	else
	{
		os << indent << "quantiles have been precomputed on a grid" << endl;
		os << indent << "size : " << itsStrikes.size() << endl;
		os << indent << "min strike : " << itsStrikes[0] << endl;
		os << indent << "max strike : " << itsStrikes[itsStrikes.size()-1] << endl;
		os << indent << "min proba : " << itsProbas[0] << endl;
		os << indent << "max proba : " << itsProbas[itsProbas.size()-1] << endl;
	}

	return os.str();
}

double ARM_GenDensityFunctor::Call_Option(double x, double fwd, double maturity) const
{
	return 0;
	// à faire proprement un jour
	/*if (IsBS())
	{
		ARM_BSModel* BSModel = static_cast< ARM_BSModel* >(GetModel());
		ARM_ZeroCurve* pCurve = BSModel->GetZeroCurve();
		double asOfDate	  = pCurve->GetAsOfDateJul();
			
		if (IsSO())
		{
			ARM_SpreadOption* spreadOption = static_cast<ARM_SpreadOption*>(itsSecurity);
			spreadOption->SetCapFloorType(K_CAP);
			spreadOption->SetStrike(x*100.);
			spreadOption->SetModel(BSModel);
			double price = spreadOption->ComputePrice();
			double pay = spreadOption->GetSpreadLeg()->GetFirstLeg()->GetPaymentDates()->Elt(0);
			double zcPay = pCurve->DiscountPrice((pay-asOfDate)/K_YEAR_LEN);
			double theta = spreadOption->GetSpreadLeg()->GetFirstLeg()->GetInterestTerms()->Elt(0);
			double notio = spreadOption->GetSpreadLeg()->GetFirstLeg()->GetAmount()->CptReferenceValue(pay);
			price/= (theta * notio * zcPay);
			return price;
		}
		else
		{
			ARM_THROW( ERR_INVALID_ARGUMENT, "ARM_GenDensityFunctor : only SO supported" );
		}
	}
	else
		ARM_THROW( ERR_INVALID_ARGUMENT, "ARM_GenDensityFunctor::Compute : check model type" );*/
}

double ARM_GenDensityFunctor::Proba(double strike, double fwd, double maturity) const
{
	// à faire proprement un jour
	/*if (IsBS())
	{
		ARM_BSModel* BSModel = static_cast< ARM_BSModel* >(GetModel());
		ARM_ZeroCurve* pCurve = BSModel->GetZeroCurve();
		double asOfDate	  = pCurve->GetAsOfDateJul();
			
		if (IsSO())
		{
			ARM_SpreadOption* spreadOption = static_cast<ARM_SpreadOption*>(GetSecurity());
			spreadOption->SetCapFloorType(K_FLOOR);
			spreadOption->SetModel(BSModel);
			double pay = spreadOption->GetSpreadLeg()->GetFirstLeg()->GetPaymentDates()->Elt(0);
			double zcPay = pCurve->DiscountPrice((pay-asOfDate)/K_YEAR_LEN);
			double theta = spreadOption->GetSpreadLeg()->GetFirstLeg()->GetInterestTerms()->Elt(0);
			double notio = spreadOption->GetSpreadLeg()->GetFirstLeg()->GetAmount()->CptReferenceValue(pay);
			
			spreadOption->SetStrike((strike+DecStrike())*100);
			double priceUp = spreadOption->ComputePrice();
			spreadOption->SetStrike((strike-DecStrike())*100);
			double priceDown = spreadOption->ComputePrice();
			return (priceUp-priceDown)/2./DecStrike()/(theta*notio*zcPay);
		}
		else
		{
			ARM_THROW( ERR_INVALID_ARGUMENT, "ARM_GenDensityFunctor : only SO supported" );
		}
	}
	else
		ARM_THROW( ERR_INVALID_ARGUMENT, "ARM_GenDensityFunctor::Compute : check model type" );*/
	return 0;
}

double ARM_GenDensityFunctor::Quantile_Direct( double proba, double fwd, double mat) const
{
	class DistributionToInverse : public DoubleToDoubleFunc 
	{
	private:
		const ARM_GenDensityFunctor& _functor;
	public: 
		DistributionToInverse(const ARM_GenDensityFunctor& functor):
		_functor(functor){};
		
		virtual double operator() (double K0)  const
		{
			return _functor.Proba(K0,0.,0.);
		}
	};
	
	DistributionToInverse func(*this);

	double guess = fwd;
	
	return Inverse(func,Inverse::REAL)(proba,guess,fwd/15.0,1e-12); 
}

double ARM_GenDensityFunctor::Quantile_Precomputed( double proba, double fwd, double mat) const
{
	int n = itsProbas.size();
	int i = 0;

	if (proba<itsProbas[0])
		i=1;
	else if (proba>=itsProbas[n-1])
		i=n-1;
	else
		while (proba>=itsProbas[i]) i++;

	return ((itsProbas[i]-proba)*itsStrikes[i-1]+(proba-itsProbas[i-1])*itsStrikes[i])/(itsProbas[i]-itsProbas[i-1]);
}

double ARM_GenDensityFunctor::Quantile( double proba, double fwd, double mat) const
{
	if (GetIsDirect())
		return Quantile_Direct(proba,fwd,mat);
	else
		return Quantile_Precomputed(proba,fwd,mat);
}

CC_END_NAMESPACE()

//-----------------------------------------------------------------------------
/*---- End of file ----*/



