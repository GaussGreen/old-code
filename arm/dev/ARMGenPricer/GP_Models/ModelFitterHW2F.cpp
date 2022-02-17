/*!
 *
 * Copyright (c) CDC IXIS CM November 2003 Paris
 *
 *	\file ModelFitterHW2F.cpp
 *
 *	\author  J. Messines
 *	\version 1.0
 *	\date June 2006
 */


/// gpmodels
#include "gpmodels/ModelFitterHW2F.h"
#include "gpmodels/HW2F.h"
#include "gpmodels/ModelParamsHW2F.h"

#include "gpcalib/vanillaspreadoption.h"
#include "gpcalib/vanillaswaption.h"
#include "gpcalib/kerneltogp.h"

/// gpclosedforms
#include "gpclosedforms/vanilla_normal.h"

/// gpinfra
#include "gpinfra/curvemodelparam.h"
#include "gpinfra/irrate.h"

#include "gpbase/curveconvert.h"



CC_BEGIN_NAMESPACE( ARM )


////////////////////////////////////////////////////
///	Class  : ARM_ModelFitterHW2F
///	Routine: Constructor (default)
///	Returns: 
///	Action : 
////////////////////////////////////////////////////
ARM_ModelFitterHW2F::ARM_ModelFitterHW2F(ARM_PricingModel* model, const ARM_StdPortfolioPtr portfolio, 
		 const ARM_ModelParamVector& calibParam, size_t nbEquations,ARM_ModelFitterPtr& linkedModelFitter)
: ARM_ModelFitter(model,portfolio,calibParam,linkedModelFitter),
	itsNbEquations(nbEquations%10),itsWithOptim(nbEquations>10)
{
}



////////////////////////////////////////////////////
///	Class  : ARM_ModelFitterHW2F
///	Routine: Constructor
///	Returns: 
///	Action : copy-const
////////////////////////////////////////////////////
ARM_ModelFitterHW2F::ARM_ModelFitterHW2F( const ARM_ModelFitterHW2F& rhs )
:	ARM_ModelFitter(rhs),
	itsNbEquations(rhs.itsNbEquations),
	itsWithOptim(rhs.itsWithOptim)
{
}


////////////////////////////////////////////////////
///	Class  : ARM_ModelFitterHW2F
///	Routine: Destructor
///	Returns: 
///	Action : 
////////////////////////////////////////////////////
ARM_ModelFitterHW2F::~ARM_ModelFitterHW2F()
{
}



////////////////////////////////////////////////////
///	Class  : ARM_ModelFitterHW2F
///	Routine: operator =
///	Returns: this
///	Action : Affectation of a rhs object
////////////////////////////////////////////////////
ARM_ModelFitterHW2F& ARM_ModelFitterHW2F::operator=(const ARM_ModelFitterHW2F& rhs)
{
	if (&rhs != this)
	{ 
		this->~ARM_ModelFitterHW2F();
		new (this) ARM_ModelFitterHW2F (rhs);
	}
	return *this;
}


////////////////////////////////////////////////////
///	Class  : ARM_ModelFitterHW2F
///	Routine: Calibrate
///	Returns: void
///	Action : calibrates the model
////////////////////////////////////////////////////
void ARM_ModelFitterHW2F::Calibrate()
{
	ARM_PricingModel* model = GetPricingModel();

	ARM_HullWhite2F* hw2Model = dynamic_cast< ARM_HullWhite2F* > (model);
	if (!hw2Model)
		throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT, " ARM_ModelFitterHW2F : only for HW2F" );

	ARM_ModelParamsHW2FExt* hw2Params     = dynamic_cast<ARM_ModelParamsHW2FExt* > (hw2Model->GetModelParams());
	if (!hw2Params)
		throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT, " ARM_ModelFitterHW2F : needs ARM_ModelParamsHW2FExt" );

	if(hw2Model->GetNbStored()==GetNbEquations())
		hw2Model->FreeFastCalib();

	StoreDatas();

	if(hw2Model->GetNbStored()==GetNbEquations())
	{
		CheckCoherence();
		SolveEquations();
	}
}

void ARM_ModelFitterHW2F::StoreDatas()
{
	ARM_PricingModel* model = GetPricingModel();

	ARM_HullWhite2F* hw2Model = dynamic_cast< ARM_HullWhite2F* > (model);
	if (!hw2Model)
		throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT, " ARM_ModelFitterHW2F : only for HW2F" );

	double asOfDate = model->GetAsOfDate().GetJulian();
	ARM_ZeroCurvePtr ZcCurve = model->GetZeroCurve();
    
	ARM_ModelParamsHW2FExt* hw2Params     = dynamic_cast<ARM_ModelParamsHW2FExt* > (hw2Model->GetModelParams());
	if (!hw2Params)
		throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT, " ARM_ModelFitterHW2F : needs ARM_ModelParamsHW2FExt" );

	ARM_CurveModelParam* mrsParam		= dynamic_cast<ARM_CurveModelParam*>(&hw2Params->GetModelParam(ARM_ModelParamType::MeanReversion));
	ARM_CurveModelParam* mrsSpreadParam = dynamic_cast<ARM_CurveModelParam*>(&hw2Params->GetModelParam(ARM_ModelParamType::MeanReversionSpread));
	double mrs1 = mrsParam->GetCurve()->Interpolate(0);
	double mrs2 = mrs1 + mrsSpreadParam->GetCurve()->Interpolate(0);

	ARM_CurveModelParam* volParam			= dynamic_cast<ARM_CurveModelParam*>(&hw2Params->GetModelParam(ARM_ModelParamType::Volatility));
	ARM_CurveModelParam* volratioParam		= dynamic_cast<ARM_CurveModelParam*>(&hw2Params->GetModelParam(ARM_ModelParamType::VolatilityRatio));
	ARM_CurveModelParam* correlParam		= dynamic_cast<ARM_CurveModelParam*>(&hw2Params->GetModelParam(ARM_ModelParamType::Correlation));

	size_t nopt = GetPortfolio()->size();
	ARM_GP_Vector breakPointTimes(nopt+1,0.);
	
	for (size_t index=0;index<nopt;index++)
	{
		csecurity* sec;
		ARM_VanillaArg* argVanilla			= ARM_ConverterFromKernel::ConvertSecuritytoArgObject(GetPortfolio()->GetAsset(index), asOfDate);
		ARM_VanillaSpreadOptionArg*  argSO	= dynamic_cast< ARM_VanillaSpreadOptionArg* > (argVanilla);
		if (argSO)
		{
			ARM_Currency* ccy = hw2Model->GetCurrency(hw2Model->GetModelName());
			sec = new cspreadoption(model,(*GetPortfolio()->GetMktPrices())[index],argSO,mrs1,mrs2,ccy);
		}
		else
		{
			ARM_VanillaSwaptionArg*  argSW	= dynamic_cast< ARM_VanillaSwaptionArg* > (argVanilla);
			if (argSW)
			{
				sec = new cswaption(model,(*GetPortfolio()->GetMktPrices())[index],argSW,mrs1,mrs2);
			}
			else
				throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT, " ARM_ModelFitterHW2F : only swaptions and spreadoptions" );
		}
		hw2Model->StoreSec(sec);
		breakPointTimes[index+1] = sec->GetResetTime();
		delete argVanilla;
	}
	
	ARM_ModelParamType::ParamNb typeParam;
	if( ARM_CurveModelParam* curveCalibParam = dynamic_cast<ARM_CurveModelParam*>(GetCalibParam()) )
	{
		typeParam = curveCalibParam->GetType();
	}
	else
		throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT, " ARM_ModelFitterHW2F::StoreDatas: wrong calibparam type" );

	volParam->UpdateValues(&breakPointTimes);
	volratioParam->UpdateValues(&breakPointTimes);
	correlParam->UpdateValues(&breakPointTimes);

	hw2Model->StoreParamType(typeParam);
	hw2Model->StoreSize(nopt);
}

void ARM_ModelFitterHW2F::CheckCoherence()
{
	ARM_PricingModel* model = GetPricingModel();

	ARM_HullWhite2F* hw2Model = dynamic_cast< ARM_HullWhite2F* > (model);
	if (!hw2Model)
		throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT, " ARM_ModelFitterHW2F : only for HW2F" );

	size_t nopt = GetPortfolio()->size();
	for (size_t n=0;n<GetNbEquations();n++)
	{
		for (size_t i=0;i<n;i++)
			if (hw2Model->GetParamType(n)==hw2Model->GetParamType(i))
				throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT, "ModelFitterHW2F : redundant model param" );
		if (hw2Model->GetSize(n)==nopt)
		{
			for (size_t index=0;index<nopt;index++)
				if (fabs(hw2Model->GetSec(n,index)->GetResetTime()-hw2Model->GetSec(0,index)->GetResetTime())>7)
					throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT, "ModelFitterHW2F : differences in date structure" );
		}
		else
			throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT, "ModelFitterHW2F : incoherent portfolio sizes" );
	}

}

void ARM_ModelFitterHW2F::SolveEquations()
{
	ARM_PricingModel* model = GetPricingModel();

	ARM_HullWhite2F* hw2Model = dynamic_cast< ARM_HullWhite2F* > (model);
	if (!hw2Model)
		throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT, " ARM_ModelFitterHW2F : only for HW2F" );

	ARM_ModelParamsHW2FExt* hw2Params     = dynamic_cast<ARM_ModelParamsHW2FExt* > (hw2Model->GetModelParams());
	if (!hw2Params)
		throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT, " ARM_ModelFitterHW2F : needs ARM_ModelParamsHW2FExt" );

	ARM_CurveModelParam* mrsParam		= dynamic_cast<ARM_CurveModelParam*>(&hw2Params->GetModelParam(ARM_ModelParamType::MeanReversion));
	ARM_CurveModelParam* mrsSpreadParam = dynamic_cast<ARM_CurveModelParam*>(&hw2Params->GetModelParam(ARM_ModelParamType::MeanReversionSpread));
	double mrs1 = mrsParam->GetCurve()->Interpolate(0);
	double mrs2 = mrs1 + mrsSpreadParam->GetCurve()->Interpolate(0);

	ARM_CurveModelParam* volParam			= dynamic_cast<ARM_CurveModelParam*>(&hw2Params->GetModelParam(ARM_ModelParamType::Volatility));
	ARM_CurveModelParam* volratioParam		= dynamic_cast<ARM_CurveModelParam*>(&hw2Params->GetModelParam(ARM_ModelParamType::VolatilityRatio));
	ARM_CurveModelParam* correlParam		= dynamic_cast<ARM_CurveModelParam*>(&hw2Params->GetModelParam(ARM_ModelParamType::Correlation));

	size_t nopt = GetPortfolio()->size();
	const ARM_GP_Vector& breakPointTimes = volParam->GetCurve()->GetAbscisses();

	ARM_GP_Vector scale(3);
	scale[0] = 2.0*mrs1;
	scale[1] = 2.0*mrs2;
	scale[2] = mrs1 + mrs2;

	size_t nbEq = GetNbEquations();
	
	for (size_t index=0;index<nopt;index++)
	{
		for (size_t niter=0;niter<5;niter++)
		{
			ARM_GP_Vector* localVarianceToPrec = hw2Params->Variance(0,breakPointTimes[index],scale);
			ARM_GP_Vector* localVarianceToNext = hw2Params->Variance(0,breakPointTimes[index+1],scale);
			ARM_GP_Vector localVariance(3);
			localVariance[0] = xi(2*mrs1,breakPointTimes[index]/K_YEAR_LEN,breakPointTimes[index+1]/K_YEAR_LEN);
			localVariance[1] = xi(2*mrs2,breakPointTimes[index]/K_YEAR_LEN,breakPointTimes[index+1]/K_YEAR_LEN);
			localVariance[2] = xi(mrs1+mrs2,breakPointTimes[index]/K_YEAR_LEN,breakPointTimes[index+1]/K_YEAR_LEN);

			ARM_GP_Vector a(nbEq);
			ARM_GP_Vector b(nbEq);
			ARM_GP_Vector c(nbEq);

			ARM_GP_Vector d(nbEq);
			bool ok=true;
			for (size_t i=0;i<nbEq;i++)
			{
				hw2Model->GetSec(i,index)->BuildCoeffs(localVariance,*localVarianceToPrec,*localVarianceToNext);
				a[i] = hw2Model->GetSec(i,index)->GetA();
				b[i] = hw2Model->GetSec(i,index)->GetB();
				c[i] = hw2Model->GetSec(i,index)->GetC();
				d[i] = hw2Model->GetSec(i,index)->GetD();
			}
			delete localVarianceToNext;
			delete localVarianceToPrec;

			ARM_GP_Vector sv(3);
			sv[0] = volParam->GetValueAtPoint(index);
			sv[1] = volratioParam->GetValueAtPoint(index);
			sv[2] = correlParam->GetValueAtPoint(index);

			ccalib* calib;

			if (nbEq==3)
				calib = new ccalib3eq(a,b,c,d,sv);
			else if (nbEq==2)
			{
				if (hw2Model->IsMissingCalibParam(ARM_ModelParamType::Volatility))
					calib = new ccalib2eqVol(a,b,c,d,sv,0,1,volParam->GetValueAtPoint(index));
				if (hw2Model->IsMissingCalibParam(ARM_ModelParamType::VolatilityRatio))
					calib = new ccalib2eqRatio(a,b,c,d,sv,0,1,volratioParam->GetValueAtPoint(index));
				if (hw2Model->IsMissingCalibParam(ARM_ModelParamType::Correlation))
					calib = new ccalib2eqRho(a,b,c,d,sv,0,1,correlParam->GetValueAtPoint(index));
			}
			else if (nbEq==1)
			{
				if (hw2Model->IsPresentCalibParam(ARM_ModelParamType::Volatility))
					calib = new ccalib1eqVol(a,b,c,d,sv,0);
				if (hw2Model->IsPresentCalibParam(ARM_ModelParamType::VolatilityRatio))
					calib = new ccalib1eqRatio(a,b,c,d,sv,0);
				if (hw2Model->IsPresentCalibParam(ARM_ModelParamType::Correlation))
					calib = new ccalib1eqRho(a,b,c,d,sv,0);
			}
			else
				throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT, " ARM_ModelFitterHW2F::SolveEquations: wrong nb of eq" );

			calib->setOptim(itsWithOptim);
			calib->solve();
			for (size_t k=index+1;k<nopt+1;k++)
			{
				volParam->SetValueAtPoint(k,calib->getVol());
				volratioParam->SetValueAtPoint(k,calib->getRatio());
				correlParam->SetValueAtPoint(k,calib->getCorrel());
			}
			delete calib;

		}
	}
}

////////////////////////////////////////////////////
///	Class  : ARM_ModelFitterHW2F
///	Routine: toString
///	Returns: string
///	Action : 
////////////////////////////////////////////////////
string ARM_ModelFitterHW2F::toString(const string& indent, const string& nextIndent) const
{
	CC_Ostringstream os;
    os << indent <<  "\n";
    os << indent << "ARM_ModelFitterHW2F\n";
    os << indent << "---------------------------\n";
	os << indent << CC_NS(std,endl);
	return os.str();
}

////////////////////////////////////////////////////
///	Class  : ARM_ModelFitterHW2F
///	Routine: ModelFitterName
///	Returns: string
///	Action : 
////////////////////////////////////////////////////
string ARM_ModelFitterHW2F::ModelFitterName() const
{
	return "ARM_ModelFitterHW2F";
}


void ccalib::keepInit()
{
	m_vol=m_sv[0];
	m_ratio=m_sv[1];
	m_correl=m_sv[2];
}
	
void ccalib::switchTo1eqVol()
{
	ccalib1eqVol calib(m_a,m_b,m_c,m_d,m_sv,0);
	calib.solve();
	m_vol = calib.getVol();
	m_ratio = calib.getRatio();
	m_correl = calib.getCorrel();
}

void ccalib::switchTo1eqRatio()
{
	ccalib1eqRatio calib(m_a,m_b,m_c,m_d,m_sv,0);
	calib.solve();
	m_vol = calib.getVol();
	m_ratio = calib.getRatio();
	m_correl = calib.getCorrel();
}

void ccalib::switchTo2eqCorrel()
{
	ccalib2eqRho calib(m_a,m_b,m_c,m_d,m_sv,0,1,m_sv[2]);
	calib.solve();
	m_vol = calib.getVol();
	m_ratio = calib.getRatio();
	m_correl = calib.getCorrel();
}

void ccalib::solveOptim()
{
	OptimCalib func(this);

	double info[LM_INFO_SZ];
	double opts[LM_OPTS_SZ];
	opts[0] = LM_INIT_MU;
	opts[1] = opts[2] = 1E-6;
	opts[3] = 1E-8;	
	opts[4] = 1E-4;

	int dim = m_a.size();
	
	ARM_GP_Vector hx(dim, 0.);
	ARM_GP_Vector fguess(dim,0.);

	initOptim(fguess);

	int status = LEVMARQMinization_WithNumDerivatives(func, fguess, hx, info, 250, opts);
}


void ccalib3eq::solve()
{
	double d = det();
	bool failure = false;

	if(fabs(det())>K_NEW_DOUBLE_TOL)
	{
		double XX = det1()/d;
		double YY = det2()/d;
		double RXY = det3()/d;
		if (XX<=0. || YY<=0.)
			failure = true;
		else
		{
			double rho = RXY/sqrt(XX*YY);
			if (rho<-1. || rho> 1.)
				failure = true;
			else
			{
				m_vol = sqrt(XX);
				m_ratio = sqrt(YY)/m_vol;
				m_correl = RXY/m_vol/m_vol/m_ratio;
			}
		}
	}
	else
		failure = true;

	if (failure)
	{
		keepInit();
		if (m_switchToOptim)
			solveOptim();
		else
			switchTo2eqCorrel();
	}
}

void ccalib3eq::initOptim(ARM_GP_Vector& fguess)
{
	int k=0;
	fguess[k++]=0.005;
	fguess[k++]=1.5;
	fguess[k++]=0.;
}

void ccalib3eq::setParam(double p[])
{
	setVol(p[0]);
	setRatio(p[1]);
	setCorrel(p[2]);
}

void ccalib2eqVol::solve()
{
	bool failure = false;
	double d = det();
	if(fabs(det())>K_NEW_DOUBLE_TOL)
	{
		double RR = det2()/d;
		double CR = det1()/d;
		if (RR<=0.)
		{
			failure = true;
		}
		else
		{
			double rho = CR/sqrt(RR);
			if (rho<-1. || rho> 1.)
				failure = true;
			else
			{
				m_vol = m_sv[0];
				m_ratio = sqrt(RR);
				m_correl = CR/m_ratio;
			}
		}
	}
	else
		failure = true;

	if (failure)
	{
		keepInit();
		if (m_switchToOptim)
			solveOptim();
		else
			switchTo1eqRatio();
	}
}

void ccalib2eqVol::initOptim(ARM_GP_Vector& fguess)
{
	int k=0;
	fguess[k++]=1.5;
	fguess[k++]=0.;
}

void ccalib2eqVol::setParam(double p[])
{
	setRatio(p[0]);
	setCorrel(p[1]);
}


void ccalib2eqRho::solve()
{
	bool failure = false;
	double D = (m_a[m_i]*m_b[m_j]-m_a[m_j]*m_b[m_i]);
	if (fabs(D)>K_NEW_DOUBLE_TOL)
	{
		double A = m_rho*m_rho*(m_b[m_i]*m_c[m_j]-m_b[m_j]*m_c[m_i])*(m_c[m_i]*m_a[m_j]-m_c[m_j]*m_a[m_i])-D*D;
		double B = m_rho*m_rho*((m_b[m_j]*m_d[m_i]-m_b[m_i]*m_d[m_j])*(m_c[m_i]*m_a[m_j]-m_c[m_j]*m_a[m_i])+(m_a[m_i]*m_d[m_j]-m_a[m_j]*m_d[m_i])*(m_b[m_i]*m_c[m_j]-m_b[m_j]*m_c[m_i]));
		double C = m_rho*m_rho*(m_b[m_j]*m_d[m_i]-m_b[m_i]*m_d[m_j])*(m_a[m_i]*m_d[m_j]-m_a[m_j]*m_d[m_i]);

		double delta = B*B-4*A*C;
		double eps = (m_rho>0.?1.:-1.)*(A>0.?1.:-1.);
		
		if (delta>=0.)
		{
			double x1 = (-B+eps*sqrt(delta))/(2*A);
			double x2 = (-B-eps*sqrt(delta))/(2*A);

			double RXY;
			if (x1*m_rho>=0)
			{
				RXY = x1;
				double XX = (m_d[m_i]-m_c[m_i]*RXY)*m_b[m_j]-(m_d[m_j]-m_c[m_j]*RXY)*m_b[m_i];
				double YY = m_a[m_i]*(m_d[m_j]-m_c[m_j]*RXY)-m_a[m_j]*(m_d[m_i]-m_c[m_i]*RXY);
				XX/=D;
				YY/=D;
				if (XX<K_NEW_DOUBLE_TOL || YY<K_NEW_DOUBLE_TOL)
				{
					if (x2*m_rho>=0)
					{
						RXY = x2;
						double XX = (m_d[m_i]-m_c[m_i]*RXY)*m_b[m_j]-(m_d[m_j]-m_c[m_j]*RXY)*m_b[m_i];
						double YY = m_a[m_i]*(m_d[m_j]-m_c[m_j]*RXY)-m_a[m_j]*(m_d[m_i]-m_c[m_i]*RXY);
						XX/=D;
						YY/=D;
						if (XX<K_NEW_DOUBLE_TOL || YY<K_NEW_DOUBLE_TOL)
							failure = true;
						else
						{
							m_vol = sqrt(XX);
							m_ratio = sqrt(YY)/m_vol;
							m_correl = RXY/m_vol/m_vol/m_ratio;
						}
					}
					else
						failure = true;
				}
				else
				{
					m_vol = sqrt(XX);
					m_ratio = sqrt(YY)/m_vol;
					m_correl = RXY/m_vol/m_vol/m_ratio;
				}
			}
			else if (x2*m_rho>=0)
			{
				RXY = x2;
				double XX = (m_d[m_i]-m_c[m_i]*RXY)*m_b[m_j]-(m_d[m_j]-m_c[m_j]*RXY)*m_b[m_i];
				double YY = m_a[m_i]*(m_d[m_j]-m_c[m_j]*RXY)-m_a[m_j]*(m_d[m_i]-m_c[m_i]*RXY);
				XX/=D;
				YY/=D;
				if (XX<K_NEW_DOUBLE_TOL || YY<K_NEW_DOUBLE_TOL)
					failure = true;
				else
				{
					m_vol = sqrt(XX);
					m_ratio = sqrt(YY)/m_vol;
					m_correl = RXY/m_vol/m_vol/m_ratio;
				}
			}
			else
				failure = true;
		}
		else
			failure = true;
	}
	else
		failure = true;

	if (failure)
	{
		keepInit();
		if (m_switchToOptim)
			solveOptim();
		else
			switchTo1eqVol();
	}
}

void ccalib2eqRho::initOptim(ARM_GP_Vector& fguess)
{
	int k=0;
	fguess[k++]=0.005;
	fguess[k++]=1.5;
}

void ccalib2eqRho::setParam(double p[])
{
	setVol(p[0]);
	setRatio(p[1]);
}

void ccalib2eqRatio::solve()
{
	bool failure = false;
	double d = det();
	if(fabs(det())>K_NEW_DOUBLE_TOL)
	{
		double XX = det1()/d;
		double RXX = det2()/d;
		if (XX<=0.)
			failure = true;
		else
		{
			double rho = RXX/XX;
			if (rho<-1. || rho> 1.)
				failure = true;
			else
			{
				m_vol = sqrt(XX);
				m_ratio = m_sv[1];
				m_correl = RXX/XX;
			}
		}
	}
	else
		failure = true;

	if (failure)
	{
		keepInit();
		if (m_switchToOptim)
			solveOptim();
		else
			switchTo1eqVol();
	}
}

void ccalib2eqRatio::initOptim(ARM_GP_Vector& fguess)
{
	int k=0;
	fguess[k++]=0.005;
	fguess[k++]=0;
}

void ccalib2eqRatio::setParam(double p[])
{
	setVol(p[0]);
	setCorrel(p[1]);
}

void ccalib1eqVol::solve()
{
	bool failure = false;
	double coeff = m_a[m_i] + m_b[m_i]*m_sv[1]*m_sv[1] + m_c[m_i]*m_sv[2]*m_sv[1];
	if (coeff!=0.)
	{
		double XX = m_d[m_i]/coeff;
		if (XX<0.)
			failure = true;
		else
		{
			m_vol = sqrt(XX);
			m_ratio = m_sv[1];
			m_correl = m_sv[2];
		}
	}
	else
		failure = true;

	if (failure)
	{
		keepInit();
		if (m_switchToOptim)
			solveOptim();
	}
	
}

void ccalib1eqVol::initOptim(ARM_GP_Vector& fguess)
{
	int k=0;
	fguess[k++]=0.005;
}

void ccalib1eqVol::setParam(double p[])
{
	setVol(p[0]);
}


void ccalib1eqRho::solve()
{
	double XX = m_sv[0]*m_sv[0];
	double YY = m_sv[0]*m_sv[1]*m_sv[0]*m_sv[1];
	if (fabs(m_c[m_i])>K_NEW_DOUBLE_TOL)
	{
		double RXY = (m_d[m_i]-m_a[m_i]*XX-m_b[m_i]*YY)/m_c[m_i];
		m_vol = m_sv[0];
		m_ratio = m_sv[1];
		double rho = RXY/m_sv[0]/m_sv[0]/m_sv[1];
		m_correl = CC_Max<double>(CC_Min<double>(rho,1.),-1.);
	}
}

void ccalib1eqRho::initOptim(ARM_GP_Vector& fguess)
{
	int k=0;
	fguess[k++]=0.;
}

void ccalib1eqRho::setParam(double p[])
{
	setCorrel(p[0]);
}


void ccalib1eqRatio::solve()

{
	bool failure = false;
	double A = m_b[m_i]*m_sv[0]*m_sv[0];
	double B = m_c[m_i]*m_sv[2]*m_sv[0]*m_sv[0];
	double C = m_a[m_i]*m_sv[0]*m_sv[0]-m_d[m_i];
	if (A!=0)
	{
		double delta = B*B-4*A*C;
		if (delta<0)
			failure = true;
		else
		{
			double r1 = (-B-sqrt(delta))/2./A;
			double r2 = (-B+sqrt(delta))/2./A;
			m_vol = m_sv[0];
			m_ratio = CC_Min<double>(CC_Max<double>(r1,0.),CC_Max<double>(r2,0.));
			m_correl = m_sv[2];
		}
	}
	else if(B!=0)
	{
		m_vol = m_sv[0];
		m_ratio = C/B;
		m_correl = m_sv[2];
	}
	else 
		keepInit();

	if (failure)
	{
		keepInit();
		if (m_switchToOptim)
			solveOptim();
	}
}

void ccalib1eqRatio::initOptim(ARM_GP_Vector& fguess)
{
	int k=0;
	fguess[k++]=1.5;
}

void ccalib1eqRatio::setParam(double p[])
{
	setRatio(p[0]);
}

CC_END_NAMESPACE()

/*----------------------------------------------------------------------------*/
/*---- End of file ----*/
