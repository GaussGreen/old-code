/*!
 *
 * Copyright (c) CDC IXIS CM November 2003 Paris
 * 
 * \file forexvanilla.h
 *  brief to factorize all vanilla/semi-vanilla fx payoffs
 * 
 *	\author  K. Belkheir
 *	\version 1.0
 *	\date January 2007
 */

//#include "gpbase/removeidentifiedwarning.h"
#include "gpcalculators/forexvanilla.h"
#include "gpbase/numericconstant.h"
#include "gpbase/utilityport.h"  /// for CC_Max


const double DEFAULT_PRECISION				= 1.e-6;

///$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
/// --- Class ARM_FXCall ---
///$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$

//--------------------------------------------------------------------------------------------------------------
CC_BEGIN_NAMESPACE( ARM )

/// standard constructor
ARM_FXCall::ARM_FXCall(	double strike,
		int callput,
		ARM_FXVanilla2DType  fxVanilla2DType,		
		bool isInvG1,
		bool isInvG2)
	:	
ARM_FXVanilla2D (strike,callput, fxVanilla2DType,isInvG1, isInvG2)
{
}
	
	//copy constructor
ARM_FXCall::ARM_FXCall( const ARM_FXCall& rhs )
	:
ARM_FXVanilla2D(rhs)
{
} 

double ARM_FXCall::Payoff( double x ) const
{
	double value;
	switch (itsVanilla2dType)
	{
	case Fx1_Fx2:
	case Fx1_InvFx2:
		{
			value = itsCallPut*(x-itsStrike);
			break;
		}
	case InvFx1_InvFx2:	
	case InvFx1_Fx2:
		{
			value = itsCallPut*(1/x-itsStrike);
			break;
		}
	}
	
	return CC_Max(value,0.0);
}

//--------------------------------------------------------------------------------------------------------------
///$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
/// --- Class ARM_FXSpreadFx1Fx2 ---
///$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$

//--------------------------------------------------------------------------------------------------------------

/// standard constructor
ARM_FXSpread::ARM_FXSpread(double strike,
		int callput,
		ARM_FXVanilla2DType  fxVanilla2Dype, 
		double alpha,
		double beta,
		bool isInvG1,
		bool isInvG2)
	:
ARM_FXVanilla2D (strike,callput, fxVanilla2Dype ,isInvG1, isInvG2),
	itsAlpha(alpha),
	itsBeta(beta)
{
}

ARM_FXSpread::ARM_FXSpread( const ARM_FXSpread& rhs )
	: 
ARM_FXVanilla2D(rhs), 
	itsAlpha(rhs.itsAlpha),
	itsBeta(rhs.itsBeta)
{
}

double ARM_FXSpread::Payoff( double x1, double x2) const
{
	double value;
	switch (itsVanilla2dType)
	{
	case Fx1_Fx2:
		{
			value = itsCallPut*(itsAlpha*x1 - itsBeta*x2 -itsStrike);
			break;
		}
	case InvFx1_InvFx2:
		{
			value = itsCallPut*(itsAlpha/x1 - itsBeta/x2 -itsStrike);
			break;
		}
	case Fx1_InvFx2:
		{
			value = itsCallPut*(itsAlpha*x1 - itsBeta/x2 -itsStrike);
			break;
		}
	case InvFx1_Fx2:
		{
			value = itsCallPut*(itsAlpha/x1 - itsBeta*x2 -itsStrike);
			break;
		}
	}
	
	return CC_Max(value,0.0);
}
///$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
/// --- Class ARM_FXBasket ---
///$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$

//--------------------------------------------------------------------------------------------------------------

double ARM_FXBasket::Payoff( double x1, double x2) const
{
	double opt1 = 0.0;
	double opt2 = 0.0;
	double und1 = itsCallPut*itsAlpha*(x1 - itsStrike);
	double und2 = 0.0;//sCoefficients[3]*itsBeta*(x2 - itsCoefficients[2]);
	if( und1 > 0)
		opt1=und1;
	if( und2 > 0)
		opt2=und2;
	double opt = und1;
	double und = und2-und1;
	//if( itsCoefficients[4]*und > 0)
	//	opt=und1 + und;
	return opt;
}

///$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
/// --- Class ARM_FXDigital ---
///$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$

//--------------------------------------------------------------------------------------------------------------

/// standard constructor
ARM_FXDigital::ARM_FXDigital(double strike,
		int callput,
		ARM_FXVanilla2DType  fxVanilla2DType,		
		bool isInvG1,
		bool isInvG2,
		ARM_DigitType  digitType,
		double epsilon)
	:	
ARM_FXVanilla2D (strike,callput,fxVanilla2DType,isInvG1,isInvG2),
itsDigitType(digitType),
itsEpsilon(epsilon)
{
}
	
	//copy constructor
ARM_FXDigital::ARM_FXDigital( const ARM_FXDigital& rhs )
	:
ARM_FXVanilla2D(rhs),
itsDigitType(rhs.itsDigitType),
itsEpsilon(rhs.itsEpsilon)
{
} 

double ARM_FXDigital::Payoff( double x ) const
{
	double value_down, value_up;
	double signed_eps = itsCallPut*itsEpsilon;
	double signed_norm = itsCallPut/itsEpsilon;
	double strike_down, strike_up;

	switch (itsDigitType)
	{
	case ARM_FXDigitType::analytic:
		{
			return 1;
			break;
		}
	case ARM_FXDigitType::backward:
		{
			strike_down = itsStrike - signed_eps;
			strike_up = itsStrike;
			break;
		}
	case ARM_FXDigitType::forward:
		{
			strike_down = itsStrike;
			strike_up = itsStrike + signed_eps;
			break;
		}
	case ARM_FXDigitType::centred:
		{
			signed_norm = itsCallPut/(2*itsEpsilon);
			strike_down = itsStrike - signed_eps;
			strike_up = itsStrike + signed_eps;
			break;
		}
	}

	switch (itsVanilla2dType)
	{
	case Fx1_Fx2:
	case Fx1_InvFx2:
		{
			value_down = itsCallPut*(x - strike_down );
			value_up   = itsCallPut*(x - strike_up );
			break;
		}
	case InvFx1_InvFx2:	
	case InvFx1_Fx2:
		{
			value_down = itsCallPut*(1/x-strike_down);
			value_up = itsCallPut*(1/x-strike_up);
			break;
		}
	}
	return signed_norm*( CC_Max(value_down,0.0) - CC_Max(value_up,0.0) );
}

///$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
/// --- Class ARM_FXQuotient ---
///$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$

/// standard constructor
ARM_FXQuotient::ARM_FXQuotient(	double strike,
		int callPut,
		ARM_FXVanilla2DType  fxVanilla2DType,
		double alpha,
		double beta,
		double strike2,
		bool isInvG1,
		bool isInvG2)
	:	
ARM_FXVanilla2D (strike,callPut, fxVanilla2DType,isInvG1, isInvG2),
itsAlpha(alpha),
itsBeta(beta),
itsStrike2(strike2)
{
}
	
//copy constructor
ARM_FXQuotient::ARM_FXQuotient( const ARM_FXQuotient& rhs )
	:
ARM_FXVanilla2D(rhs),
itsAlpha(rhs.itsAlpha),
itsBeta(rhs.itsBeta),
itsStrike2(rhs.itsStrike2)
{
} 

double ARM_FXQuotient::Payoff( double x1, double x2) const
{
	double value;
	switch (itsVanilla2dType)
	{
	case Fx1_Fx2:
		{
			value = itsCallPut*( itsStrike2*(x1+itsAlpha)/(x2+itsBeta) - itsStrike );
			break;
		}
	case InvFx1_InvFx2:
		{
			value = itsCallPut*(itsStrike2*( 1.0/x1+itsAlpha)/(1.0/x2+itsBeta) - itsStrike );
			break;
		}
	case Fx1_InvFx2:
		{
			value = itsCallPut*( itsStrike2*(x1+itsAlpha)/(1.0/x2+itsBeta) - itsStrike );
			break;
		}
	case InvFx1_Fx2:
		{
			value = itsCallPut*( itsStrike2*(1.0/x1+itsAlpha)/(x2+itsBeta) - itsStrike );
			break;
		}
	}
	return CC_Max(value,0.0);
}

///$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
/// --- Class ARM_FXBasketInvFx1Fx2 ---
///$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$

//--------------------------------------------------------------------------------------------------------------
/*
double ARM_FXBasketInvFx1Fx2::Payoff( double x1, double x2) const
{
	double opt1 = 0.0;
	double opt2 = 0.0;
	double und1 = itsCallPut*itsAlpha*(1/x1 - itsStrike);
	double und2 = 0.0;//itsCoefficients[3]*itsBeta*(x2 - itsCoefficients[2]);
	if( und1 > 0)
		opt1=und1;
	if( und2 > 0)
		opt2=und2;
	double opt = und1;
	double und = und2-und1;
//	if( itsCoefficients[4]*und > 0)
		//opt=und1 + und;
	return opt;
}

///$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
/// --- Class ARM_FXBasketFx1InvFx2 ---
///$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$

//--------------------------------------------------------------------------------------------------------------

double ARM_FXBasketFx1InvFx2::Payoff( double x1, double x2) const
{
	double opt1 = 0.0;
	double opt2 = 0.0;
	double und1 = itsCallPut*itsAlpha*(x1 - itsStrike);
	double und2 = 0.0;//itsCoefficients[3]*itsBeta*(1/x2 - itsCoefficients[2]);
	if( und1 > 0)
		opt1=und1;
	if( und2 > 0)
		opt2=und2;
	double opt = und1;
	double und = und2-und1;
	//if( itsCoefficients[4]*und > 0)
	//	opt=und1 + und;
	return opt;
}

///$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
/// --- Class ARM_FXBasketInvFx1InvFx2 ---
///$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$

//--------------------------------------------------------------------------------------------------------------

double ARM_FXBasketInvFx1InvFx2::Payoff( double x1, double x2) const
{
	double opt1 = 0.0;
	double opt2 = 0.0;
	double und1 = itsCallPut*itsAlpha*(1/x1 - itsStrike);
	double und2 = 0.0;//itsCoefficients[3]*itsBeta*(1/x2 - itsCoefficients[2]);
	if( und1 > 0)
		opt1=und1;
	if( und2 > 0)
		opt2=und2;
	double opt = und1;
	double und = und2-und1;
	//if( itsCoefficients[4]*und > 0)
	//	opt=und1 + und;
	return opt;
}
*/

///$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
/// --- Class ARM_GaussReplic1D ---
///$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$

//--------------------------------------------------------------------------------------------------------------

/// standard constructor
ARM_GaussReplic1D::ARM_GaussReplic1D(ARM_FXVanilla* vanilla,
					const ARM_GP_Matrix& glParams,
					ARM_1DQuantoType type)
	:
	ARM_GaussReplic( vanilla,glParams ),
		itsQuanto1DType(type)
{
}

ARM_GaussReplic1D::ARM_GaussReplic1D( const ARM_GaussReplic1D& rhs )
	:
ARM_GaussReplic(rhs),
	itsQuanto1DType(rhs.itsQuanto1DType)
		
{
}

double ARM_GaussReplic1D::Price() const
{
	double sum=0.0;
	double dSimple_Int,pay,fwd,xj;
	size_t size = itsGLParams.cols();
	for( int j=0; j<size; ++j )
	{
		xj = itsGLParams.Elt(0,j);
		fwd = itsGLParams.Elt(2,j);
		pay = itsFXVanilla->Payoff( fwd );
		double fwdquanto = itsQuanto1DType == ARM_GaussReplic1D::InvQuanto ? fwd : 1.0;
		dSimple_Int = fwdquanto*pay*exp(-0.5*xj*xj);
		sum += itsGLParams.Elt(1,j)*dSimple_Int;
	}
	return ARM_NumericConstants::ARM_INVSQRT2PI*sum;
}

///$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
/// --- Class ARM_GaussReplic2D ---
///$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$

//--------------------------------------------------------------------------------------------------------------

/// standard constructor
ARM_GaussReplic2D::ARM_GaussReplic2D(ARM_FXVanilla* vanilla,
				double rho,
				ARM_QuantoType  quantoType,
				const ARM_GP_Matrix& glmatrix, 
				const ARM_GP_Matrix& glmatrix2)
	:
ARM_GaussReplic(vanilla, glmatrix),
	 itsGLParams2(glmatrix2),
	 itsRho(rho),
	 itsQuantoType(quantoType)
{
}
	//copy constructor
ARM_GaussReplic2D::ARM_GaussReplic2D ( const ARM_GaussReplic2D& rhs )
	:
ARM_GaussReplic(rhs),
	 itsGLParams2(rhs.itsGLParams2),
	 itsRho(rhs.itsRho),
	 itsQuantoType(rhs.itsQuantoType)
{
}

double ARM_GaussReplic2D::Price() const
{
	double xi, xj, fwd1, fwd2;
	double Double_Int = 0;
	double Simple_Int = 0;
	double dSimple_Int, dDouble_Int;
	double pay=0;
	int i,j;
	size_t size = itsGLParams.cols();
	size_t size2 = itsGLParams2.cols();
	for(i=0 ; i<size; ++i )
	{
		xi = itsGLParams.Elt(0,i);
		fwd1 = itsGLParams.Elt(2,i);
		//Calculation of the second integral
		Simple_Int=0.0;
		for( j=0; j<size2; ++j )
		{
			fwd2 = itsGLParams2.Elt(2,j);
			xj = itsGLParams2.Elt(0,j);
			pay = itsFXVanilla->Payoff(fwd1, fwd2);
			double fwd = itsQuantoType == Quanto1 ? fwd1 : (itsQuantoType == Quanto2 ? fwd2 : 1.0);
			dSimple_Int = fwd*pay*exp(-1/(2*(1-itsRho*itsRho))*(xi*xi+xj*xj-2*itsRho*xi*xj));
			Simple_Int += itsGLParams2.Elt(1,j)*dSimple_Int;
		}
		dDouble_Int=Simple_Int;
		Double_Int += itsGLParams.Elt(1,i)*dDouble_Int;
	}
	return 1/(ARM_NumericConstants::ARM_2_PI*sqrt(1-itsRho*itsRho))*Double_Int;
}

CC_END_NAMESPACE()

//-----------------------------------------------------------------------------
/*---- End of file ----*/



