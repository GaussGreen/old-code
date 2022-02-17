/*!
 *
 * Copyright (c) CDC IXIS CM July 2003 Paris
 *
 *	\file implicitcorel.cpp
 *
 *  \brief object to check coherence test
 *		between inflation data
 *	\author  N. Belgrade
 *	\version 1.0
 *	\date January 2004
 */


#include "gpinflation/implicitcorrel.h"

/// gpbase
#include "gpbase/checkarg.h"
#include "gpbase/gpvector.h"
#include "gpbase\utilityport.h"


/// kernel
#include <glob/expt.h>

/// standard libraries
#include <string>
#include <cmath>

CC_BEGIN_NAMESPACE( ARM )



///fonction convertissant une chaine de caractère en type tool

void ARM_MarketDataValidator::Convert(const string& x, tool& y )
{
	if (x == "ZC") {y = Case1;} else if (x == "YtY") {y = Case2;} else if(x == "Cor") {y = Case3;};
}


//fonction permettant l'interpolation linéaire d'une fonction

void ARM_MarketDataValidator::Interpol(const double& a, const double& b, const double& c, const double& f_a, double& f_b, const double& f_c)
{
	if (c!=a) f_b=f_a+(b-a)/(c-a)*(f_c-f_a);
}


void ARM_MarketDataValidator::HmgYtYCor_to_ZC( 
		ARM_GP_Vector* pYtYVol, 
		ARM_GP_Vector* pCor, 
		ARM_GP_Vector* pMaturity,
		const int length,
		ARM_GP_Vector& ZCVol )
{
	throw Exception(__LINE__, __FILE__, ERR_INVALID_DATA, "unimplemented HmgYtYCor_to_ZC!" );
}



/// this function computes the corresponding correlation for a set of given volatilities in Black Scholes model

void ARM_MarketDataValidator::Vol_to_Cor(
	const ARM_GP_Vector& ZCVol, 
	const ARM_GP_Vector& YtYVol, 
	const ARM_GP_Vector& Maturity,
	ARM_GP_Vector& Cor )
{
	CC_NS(ARM_Check,CheckSameArgSize)< ARM_GP_Vector, ARM_GP_Vector >( ZCVol, YtYVol,		"Zero coupon vol", "Year to Year Vol" );
	CC_NS(ARM_Check,CheckSameArgSize)< ARM_GP_Vector, ARM_GP_Vector >( ZCVol, Maturity,	"Zero coupon vol", "Maturity" );
	Cor.resize( ZCVol.size()-1 );


	for (int i=2;i<Cor.size()+2;i++)
	{
		Cor.Elt(i-2)=(Maturity.Elt(i-2)*pow(ZCVol.Elt(i-2),2)+Maturity.Elt(i-1)*pow(ZCVol.Elt(i-1),2)-(Maturity.Elt(i-1)-Maturity.Elt(i-2))*pow(YtYVol.Elt(i-1),2))/(2*Maturity.Elt(i-2)*ZCVol.Elt(i-1)*ZCVol.Elt(i-2));
	}
	
}


/// this function computes zero coupon volatilities from the year to year volatilities and correlations in Black Scholes model

void ARM_MarketDataValidator::YtYCor_to_ZC( 
	const ARM_GP_Vector& YtYVol, 
	const ARM_GP_Vector& Cor, 
	const ARM_GP_Vector& Maturity,
	ARM_GP_Vector& ZCVol )
{
	CC_NS(ARM_Check,CheckSameArgSize)< ARM_GP_Vector, ARM_GP_Vector >( YtYVol, Maturity,	"Year to Year Vol", "Maturity" );
	CC_NS(ARM_Check,CheckArgSize)< ARM_GP_Vector >( Cor, "Correlation",  YtYVol.size()-1 );

	ZCVol.resize( YtYVol.size() );
	double delta=0.0;
	
	ZCVol.Elt(0)=YtYVol.Elt(0);
	
	for (int i=1;i<ZCVol.size();i++)
	{
		/// solve the second order equation
		delta = Maturity.Elt(i-1)*pow(ZCVol.Elt(i-1),2) * (pow(Cor.Elt(i-1),2)*Maturity.Elt(i-1)-Maturity.Elt(i)) + Maturity.Elt(i)*(Maturity.Elt(i)-Maturity.Elt(i-1)) * pow(YtYVol.Elt(i),2);
		
		if (delta >=0)
			ZCVol.Elt(i)=(Cor.Elt(i-1)*Maturity.Elt(i-1)*ZCVol.Elt(i-1)+sqrt(delta))/(Maturity.Elt(i));
						 
		else ZCVol.Elt(i)=0.0;
		delta = 0.0;
	}
}


/// this function computes year to year volatilities from zero coupon volatilities and correlations in Black Scholes model

void ARM_MarketDataValidator::ZCCor_to_YtY( 
	const ARM_GP_Vector& ZCVol, 
	const ARM_GP_Vector& Cor, 
	const ARM_GP_Vector& Maturity,
	ARM_GP_Vector& YtYVol )
{
	CC_NS(ARM_Check,CheckSameArgSize)< ARM_GP_Vector, ARM_GP_Vector >( ZCVol, Maturity,	"Zero coupon Vol", "Maturity" );
	CC_NS(ARM_Check,CheckArgSize)< ARM_GP_Vector >( Cor, "Correlation",  ZCVol.size()-1 );
	YtYVol.resize( ZCVol.size() );
	
	YtYVol.Elt(0)=ZCVol.Elt(0);
	
	for (int i=1;i<YtYVol.size();i++)
	{
		YtYVol.Elt(i)=sqrt((Maturity.Elt(i)*pow(ZCVol.Elt(i),2)+Maturity.Elt(i-1)*pow(ZCVol.Elt(i-1),2)-2*Cor.Elt(i-1)*Maturity.Elt(i-1)*ZCVol.Elt(i-1)*ZCVol.Elt(i))/(Maturity.Elt(i)-Maturity.Elt(i-1)));

	}
}


/// this function computes le covariances between CPI fwds from zero coupon volatlities

void ARM_MarketDataValidator::HmgCov( ARM_GP_Vector* pZCVol, 
		ARM_GP_Vector* pMaturity,
		const int length,
		ARM_GP_Vector& Cov )
{
	CC_NS(ARM_Check,CheckSameArgSize)< ARM_GP_Vector, ARM_GP_Vector >( *pZCVol, *pMaturity,	"Zero coupon Vol",	"Maturity" );
	ARM_GP_Vector tempo;
	ARM_GP_Vector f;
	ARM_GP_Vector g;

	double step = 1.0/length;
	
	tempo.resize( pZCVol->size()*length);
	f.resize	( pZCVol->size()*length);
	g.resize	( f.size()-length);
	Cov.resize	( pZCVol->size()-1);
	
// interpolation of zero coupon variances of fwd CPI

	int i=0, j=0;
	double interm = 0.0;
	double spy1 = 0.0;
	double spy2 = (pMaturity->Elt(0))*pow(pZCVol->Elt(0),2);

	double spy3 = 0.0;


// from 0 to the first zero coupon volatlity data
	
	for ( i=0; i<length; i++)
	{
		interm = (i+1)*step;
		ARM_MarketDataValidator::Interpol(0, interm, pMaturity->Elt(0), spy1, tempo.Elt(i), spy2);
		spy3 = tempo.Elt(i);
	}

	interm=0.0;

	int k=length;


// for the rest of zero coupon volatlity datas

	for (i=1; i<pZCVol->size(); i++)
	{
		for (int j=0; j<length; j++)
		{

			interm = pMaturity->Elt(i-1) + (j+1)*step;

			spy1 = (pMaturity->Elt(i-1))*pow(pZCVol->Elt(i-1),2);

			spy2 = (pMaturity->Elt(i))*pow(pZCVol->Elt(i),2);
			
			ARM_MarketDataValidator::Interpol(pMaturity->Elt(i-1), interm, pMaturity->Elt(i), spy1, tempo.Elt(k), spy2);

			spy3 = tempo.Elt(k);
			//spy3 = tempo.Elt(j+i*length);
			
			k+=1;

		}
	}

	interm = 0.0;

// approximation of the volatility function by deriving the zero coupon variances interpoled


	f.Elt(0)= sqrt(tempo.Elt(0));

	for (i=1; i<f.size(); i++)
	{
		spy1 = tempo.Elt(i);
		spy2 = tempo.Elt(i-1);

		f.Elt(i)= sqrt((tempo.Elt(i)-tempo.Elt(i-1))/step);

		spy3 = f.Elt(i);

	}
	
	interm=0.0;


// approximation of the whole covariance

	//g.Elt(0)= 1.0;

	for (i=0; i<g.size(); i++)
	{

		interm += f.Elt(i)*f.Elt(i+length)*step;

		g.Elt(i)= interm;
		
	}
	
	interm=0.0;


// deduction of covariances 

	j=0;

	for (i=length-1; i<g.size(); i+=length)
	{
			
		Cov.Elt(j)= g.Elt(i);
		//spy3= f.Elt(i);
		j += 1;
		
	}

}

/// this function computes the corresponding correlation for a set of given volatilities in homogeneous case

void ARM_MarketDataValidator::HmgVol_to_Cor(
		ARM_GP_Vector* pZCVol, 
		ARM_GP_Vector* pYtYVol, 
		ARM_GP_Vector* pMaturity,
		const int length,
		ARM_GP_Vector& Cor )
{
	CC_NS(ARM_Check,CheckSameArgSize)< ARM_GP_Vector, ARM_GP_Vector >( *pZCVol, *pMaturity,	"Zero coupon Vol",	"Maturity" );
	CC_NS(ARM_Check,CheckSameArgSize)< ARM_GP_Vector, ARM_GP_Vector >( *pYtYVol,*pMaturity,	"Year to Year Vol",	"Maturity" );
	ARM_GP_Vector Cov;
	Cov.resize( pZCVol->size()-1 );
	Cor.resize( pZCVol->size()-1 );

	ARM_MarketDataValidator::HmgCov( pZCVol, pMaturity, length, Cov );

	double spy = 0.0;

	for (int i=0; i<Cor.size(); i++)
	{
		Cor.Elt(i)=( (pMaturity->Elt(i+1)*pow(pZCVol->Elt(i+1),2) + pMaturity->Elt(i)*pow(pZCVol->Elt(i),2) - (pMaturity->Elt(i+1)-pMaturity->Elt(i))*pow(pYtYVol->Elt(i+1),2)) )/(2*Cov.Elt(i));
		spy = Cor.Elt(i);
	}

}



/// this function computes year to year volatilities from zero coupon volatilities and correlations in homogeneous case

void ARM_MarketDataValidator::HmgZCCor_to_YtY( 
	ARM_GP_Vector* pZCVol, 
	ARM_GP_Vector* pCor, 
	ARM_GP_Vector* pMaturity,
	const int length,
	ARM_GP_Vector& YtYVol )
{
	CC_NS(ARM_Check,CheckSameArgSize)< ARM_GP_Vector, ARM_GP_Vector >( *pZCVol, *pMaturity,	"Zero coupon Vol", "Maturity" );
	CC_NS(ARM_Check,CheckArgSize)< ARM_GP_Vector >( *pCor, "Correlation",  pZCVol->size()-1 );
	ARM_GP_Vector Cov;
	Cov.resize( pZCVol->size()-1 );
	YtYVol.resize( pZCVol->size() );

	ARM_MarketDataValidator::HmgCov( pZCVol, pMaturity, length, Cov );
	
	YtYVol.Elt(0)=pZCVol->Elt(0);

	double spy=0.0;
	
	for (int i=1;i<YtYVol.size();i++)
	{
		spy = (pMaturity->Elt(i)*pow(pZCVol->Elt(i),2)+pMaturity->Elt(i-1)*pow(pZCVol->Elt(i-1),2)-2*pCor->Elt(i-1)*Cov.Elt(i-1))/(pMaturity->Elt(i)-pMaturity->Elt(i-1));
		
		if (spy>=0) {YtYVol.Elt(i)=sqrt(spy);} else YtYVol.Elt(i) = 0.0;
		

	}
}







/// computes the confidence intervals for the zero coupon, year to year volatilities and correlations

void ARM_MarketDataValidator::Bounds( 
	ARM_GP_Vector* pZCVol,
	ARM_GP_Vector* pYtYVol,
	ARM_GP_Vector* pCor,
	ARM_GP_Vector* pMaturity,
	ARM_GP_Vector& UBound,
	ARM_GP_Vector& LBound,
	const string& choice,
	string& TBound )
{
	CC_NS(ARM_Check,CheckSameArgSize)< ARM_GP_Vector, ARM_GP_Vector >( *pZCVol, *pMaturity,	"Zero coupon Vol",	"Maturity" );
	CC_NS(ARM_Check,CheckSameArgSize)< ARM_GP_Vector, ARM_GP_Vector >( *pYtYVol,*pMaturity,	"Year to Year Vol",	"Maturity" );
	CC_NS(ARM_Check,CheckArgSize)< ARM_GP_Vector >( *pCor, "Correlation",  pYtYVol->size()-1 );

	UBound.resize( pYtYVol->size()-1);
	LBound.resize( pYtYVol->size()-1);

	int i=0, j=0;
	double spy1 = 0.0, spy2 = 0.0;
	
	tool intern_choice = Case1;
	
	Convert( choice, intern_choice);

	switch( intern_choice )
	{
	case Case1:
		for (i=1;i< pYtYVol->size();i++)
		{	
			
			j=(int) ( pMaturity->Elt(i) - pMaturity->Elt(i-1))-1;
			
			if ( pCor->Elt(i-1) != 0.5){
				
				if ( pCor->Elt(i-1)<0.5)
				{
					TBound="internal";
				}
				else
				{
					TBound="external";
				};
				
				spy1 = (pMaturity->Elt(i)-pMaturity->Elt(i-1)) *	pow(pYtYVol->Elt(i),2) -(1.0-2.0*pCor->Elt(i-1)) * pMaturity->Elt(i-1)
					*pow(pZCVol->Elt(i-1),2);
				
				spy2 = ((pMaturity->Elt(i)-pMaturity->Elt(i-1)) * pow(pYtYVol->Elt(i),2) - pMaturity->Elt(i-1) * pow(pZCVol->Elt(i-1),2)-2*pCor->Elt(i-1)*pMaturity->Elt(j)*pow(pZCVol->Elt(j),2)) / (1.0-2*pCor->Elt(i-1));
				
				
			} else {
				spy1=0.0;
				spy2=0.0;
			};

			if (spy1<0) {spy1=0;};
			if (spy2<0) {spy2=0;}
			
			spy1=sqrt(spy1/pMaturity->Elt(i));
			
			spy2=sqrt(spy2/pMaturity->Elt(i));
			
			LBound.Elt(i-1)=CC_Min(spy1,spy2);
			UBound.Elt(i-1)=CC_Max(spy1,spy2);
			
			spy1=0.0;
			spy2=0.0;
			
		}
		
		break;
		
	case Case2:
		
		TBound = "internal";
		
		for (i=1;i< pYtYVol->size();i++)
		{
			j=(int) (pMaturity->Elt(i)-pMaturity->Elt(i-1))-1;
			
			spy1 = (1.0-2.0*pCor->Elt(i-1))*pMaturity->Elt(i)*pow(pZCVol->Elt(i),2)+pMaturity->Elt(i-1)*pow(pZCVol->Elt(i-1),2)+2.0*pCor->Elt(i-1)*pMaturity->Elt(j)*pow(pZCVol->Elt(j),2);
			
			spy2 = pMaturity->Elt(i)*pow(pZCVol->Elt(i),2)+(1.0-2.0*pCor->Elt(i-1))*pMaturity->Elt(i-1)*pow(pZCVol->Elt(i-1),2);
			
			if (spy1<0) {spy1=0;};
			if (spy2<0) {spy2=0;};
			
			spy1=sqrt(spy1/pMaturity->Elt(j));
			
			spy2=sqrt(spy2/pMaturity->Elt(j));
			
			LBound.Elt(i-1)=CC_Min(spy1,spy2);
			UBound.Elt(i-1)=CC_Max(spy1,spy2);
			
			spy1=0.0;
			spy2=0.0;
		}
		
		break;
		
	case Case3:
		
		TBound = "internal";
		
		for (i=1;i< pYtYVol->size();i++)
		{
			j=(int) (pMaturity->Elt(i)-pMaturity->Elt(i-1))-1;
			
			spy1 = 0.5 * (pMaturity->Elt(i)*pow(pZCVol->Elt(i),2)+ pMaturity->Elt(i-1)*pow(pZCVol->Elt(i-1),2)-(pMaturity->Elt(i)-pMaturity->Elt(i-1))*pow(pYtYVol->Elt(i),2))/(pMaturity->Elt(i)*pow(pZCVol->Elt(i),2)-pMaturity->Elt(j)*pow(pZCVol->Elt(j),2));
			
			spy2 = 0.5 * (1.0+(pMaturity->Elt(i)*pow(pZCVol->Elt(i),2)-(pMaturity->Elt(i)-pMaturity->Elt(i-1))*pow(pYtYVol->Elt(i),2))/(pMaturity->Elt(i-1)*pow(pZCVol->Elt(i-1),2)));
			
			LBound.Elt(i-1)=CC_Min(spy1,spy2);
			UBound.Elt(i-1)=CC_Max(spy1,spy2);
			
			spy1=0.0;
			spy2=0.0;
		}
		
		break;
		
	default:
		throw Exception(__LINE__, __FILE__, ERR_INVALID_DATA,"Invalid choice");
	}
} 




void ARM_Convert_Vol::VolYoY_to_VolSwp(
	ARM_GP_Vector* pDFactor, 
	ARM_GP_Vector* pFwdCPI, 
	ARM_GP_Vector* pVol_DF,
	ARM_GP_Vector* pVol_YoY, 
	ARM_GP_Vector* pAvgCor, 
	ARM_GP_Vector* pDates, 
	ARM_GP_Vector* pTenors, 
	const double SwpRate, 
	double &Vol_Swp)

{
	CC_NS(ARM_Check,CheckSameArgSize)< ARM_GP_Vector, ARM_GP_Vector >( *pDFactor, *pFwdCPI,	"DFactor",	"FwdCPI" );
	CC_NS(ARM_Check,CheckSameArgSize)< ARM_GP_Vector, ARM_GP_Vector >( *pDFactor, *pVol_DF,	"DFactor",	"Vol_DF" );
	CC_NS(ARM_Check,CheckSameArgSize)< ARM_GP_Vector, ARM_GP_Vector >( *pDFactor, *pVol_YoY,	"DFactor",	"Vol_YoY");
	CC_NS(ARM_Check,CheckSameArgSize)< ARM_GP_Vector, ARM_GP_Vector >( *pDFactor, *pDates,	"DFactor",	"Dates"  );
	CC_NS(ARM_Check,CheckSameArgSize)< ARM_GP_Vector, ARM_GP_Vector >( *pDFactor, *pTenors,	"DFactor",	"Tenors" );
	CC_NS(ARM_Check,CheckSameArgSize)< ARM_GP_Vector, ARM_GP_Vector >( *pDFactor, *pAvgCor,	"DFactor",  "AvgCor" );
	
	int i = 0, Maturity = 0;
	double FixedLeg = 0., sum = 0.;
	Maturity = (int) pDates->Elt(pDates->size()-1);

	for (i=0; i<Maturity; i++)
	{
		FixedLeg += pDFactor->Elt(i);
	}

	for (i=0; i<Maturity; i++)
	{
		sum += pow( pDFactor->Elt(i) * pFwdCPI->Elt(i) / SwpRate, 2 ) * pTenors->Elt(i) * pow( pVol_YoY->Elt(i), 2 );
		sum += pow( pFwdCPI->Elt(i) / SwpRate - pTenors->Elt(i)	, 2 ) * pow ( pDFactor->Elt(i), 2 ) * pDates->Elt(i) * pow( pVol_DF->Elt(i), 2);
		sum += 2. * ( pFwdCPI->Elt(i) / SwpRate - pTenors->Elt(i) ) * pow( pDFactor->Elt(i), 2 ) * pFwdCPI->Elt(i) / SwpRate 
			* pAvgCor->Elt(i) * pVol_YoY->Elt(i) * pVol_DF->Elt(i);
	}

	Vol_Swp = pow( sum / pDates->size(), 0.5 ) / FixedLeg;
}


CC_END_NAMESPACE()

/*---------------------------------------------------------------*/
/*---- End Of File ----*/

