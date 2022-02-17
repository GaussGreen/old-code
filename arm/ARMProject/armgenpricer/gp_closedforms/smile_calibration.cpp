
#include "firsttoinc.h"
#include "gpbase/port.h"
#include "gpbase/gpvector.h"


#include <cmath>
#include <complex>
#include "gpbase/gpmatrix.h"
#include "gpbase/gpvector.h"
#include "gpbase/cloneutilityfunc.h"

#include "gpbase/removenagwarning.h"

#include "gpnumlib/gaussiananalytics.h"
#include "gpbase/numericconstant.h"

#include "gpclosedforms/sabrvanilla.h"
#include "gpclosedforms/sabrimpliedvol.h"
#include "gpclosedforms/smile_sabr.h"

#include "gpclosedforms/sabr_calibration.h"
#include "gpclosedforms/extendedsabrformula.h"
#include "gpclosedforms/optimization1.h"
#include "gpclosedforms/vanille_bs_interface.h"

#include "expt.h"   // for the exceptions

#include "gpclosedforms/smile_calibration.h"

using namespace std;

CC_BEGIN_NAMESPACE(ARM)

ARM_SmileCalibration::~ARM_SmileCalibration()
{
	DeletePointorVector<ARM_GP_Vector>(m_mktvols);
	DeletePointorVector<ARM_GP_Vector>(m_strikes);
}

void ARM_SmileCalibration::Init(double reset, double forward, const std::vector<double>& mktvols, const std::vector<double>& strikes,
								bool calibConstraint, double constraintK, double constraintVol, 
								ARM_SmileCalibration_Params *  params)
{
	m_resets.resize(1);
	m_resets[0] = reset;
	m_forwards.resize(1);
	m_forwards[0] = forward;
	m_constraintStrikes.resize(1);
	m_constraintStrikes[0] = constraintK;
	m_constraintVols.resize(1);
	m_constraintVols[0] = constraintVol;

	m_weights.resize(1,1.);

	m_mktvols.resize(1);
	m_mktvols[0] = new ARM_GP_T_Vector<double>(mktvols.size());
	for(int k = 0; k < mktvols.size(); k++) m_mktvols[0]->Elt(k)= mktvols[k];

	m_strikes.resize(1);
	m_strikes[0] = new ARM_GP_Vector(strikes.size());
	for(int k = 0; k < strikes.size(); k++) (*(m_strikes[0]))[k] = strikes[k];

	m_calibConstraint = calibConstraint;

	m_idx = 0;

	InitOneReset(params);
}

void ARM_SmileCalibration::Init(const std::vector<double>& resets, const std::vector<double>& forwards,
								const ARM_VectorVector& mktvols, const ARM_VectorVector& strikes,
								bool calibConstraint, const std::vector<double>& constraintStrikes, const std::vector<double>& constraintVols,
								const std::vector<double>& weights,
								ARM_SmileCalibration_Params * params)
{
	m_resets			= resets;
	m_forwards			= forwards;
	
	DuplicateCloneablePointorVectorInPlace<ARM_GP_Vector>(mktvols, m_mktvols);
	DuplicateCloneablePointorVectorInPlace<ARM_GP_Vector>(strikes, m_strikes);

	m_calibConstraint	= calibConstraint;
	m_constraintStrikes	= constraintStrikes;
	m_constraintVols	= constraintVols;
	m_weights			= weights;

	if(m_weights.size() == 0) m_weights.resize(m_resets.size(), 1.);

	m_idx = 0;

	InitMultipleResets(params);
}


void ARM_SmileCalibration::CalibConstraint()
{
	CalibConstraintFunc func(this);

	double best;

	SetConstrainedParam(brentSolve(func, m_constraintVols[m_idx], m_constrainedParamLBound, m_constrainedParamUBound, 1e-8, 50, 0, &best));
}

void ARM_SmileCalibration::SetIdx(int k)
{
	m_idx = k;
}

bool ARM_SmileCalibration::Calibrate()
{
	CalibFunc func(this);

	double info[LM_INFO_SZ];
	double opts[LM_OPTS_SZ];
	opts[0] = LM_INIT_MU;
	opts[1] = opts[2] = 1E-6;
	opts[3] = 1E-8;	
	opts[4] = 1E-4;

	int k, size = 0;
	for(k = 0; k < m_mktvols.size(); k++) size += m_mktvols[k]->size();

	std::vector<double> fguess, hx(size, 0.);

	InitfGuess(fguess);

	int status = LEVMARQMinization_WithNumDerivatives(func, fguess, hx, info, 250, opts);

	m_calibsuccess = info[1] > info[0] / 10. ? true : false;

	m_quaderror = info[1];

	return m_calibsuccess;
}

void ARM_SmileCalibration_SABR::InitOneReset(ARM_SmileCalibration_Params * params)
{
	m_params = dynamic_cast<ARM_SmileCalibration_Params_SABR*>(params);

	m_constrainedParamLBound	= 1e-8;
	m_constrainedParamUBound	= 2.;

	m_params->m_rhos.resize(1, m_params->m_rho);
	m_params->m_shifts.resize(1,m_params->m_shift);
	m_params->m_levels.resize(1,1);
}

void ARM_SmileCalibration_SABR::InitMultipleResets(ARM_SmileCalibration_Params * params)
{
	m_params = dynamic_cast<ARM_SmileCalibration_Params_SABR*>(params);

	m_params->m_levels.resize(m_resets.size(), 1.);
	m_params->m_shifts.resize(m_resets.size(), m_params->m_shift);
	m_params->m_rhos.resize(m_params->m_calibeachrho ? m_resets.size() : 1, m_params->m_rho);

	m_constrainedParamLBound	= 0.001;
	m_constrainedParamUBound	= 5.;
}

void ARM_SmileCalibration_SABR::SetParam(double param[])
{
	int k = 0;

	if(m_params->m_calibalpha)
	{
		if(m_params->m_levels.size() > 0 || (m_params->m_levels.size() == 0 && m_calibConstraint == false)) 
		{
			m_params->m_alpha = param[k] < 1e-8 ? 1e-8 : param[k];
			k++;
		}
	}

	if(m_params->m_calibbeta)
	{
		m_params->m_beta = param[k] < 0.1 ? 0.1 : param[k] > 1. ? 1. : param[k];
		k++;
	}

	if(m_params->m_calibrho && m_params->m_calibeachrho == false)
	{
		m_params->m_rhos[0] = param[k] < -0.999 ? -0.999 : param[k] > 0.999 ? 0.999 : param[k];
		k++;
	}

	if(m_params->m_calibnu) 
	{
		m_params->m_nu = param[k] < 1e-8 ? 1e-8 : param[k];
		k++;
	}

	if(m_params->m_levels.size() == 0)
	{
		m_locidx = (m_params->m_calibalpha && m_calibConstraint == false) + m_params->m_calibrho + m_params->m_calibnu;
	}
	else
	{
		m_locidx = m_params->m_calibalpha + m_params->m_calibnu + (m_params->m_calibrho && m_params->m_calibeachrho == false);
	}
}

void ARM_SmileCalibration_SABR::SetConstrainedParam(double x)
{
	if(m_params->m_levels.size() > 1 || (m_params->m_levels.size() == 1 && m_params->m_calibalpha == false)) 
	{
		for(int k = m_idx; k < m_params->m_levels.size(); k++) m_params->m_levels[k] = x < m_constrainedParamLBound ? m_constrainedParamLBound : x > m_constrainedParamUBound ? m_constrainedParamUBound : x;
	}
	else
	{
		if(m_params->m_calibalpha) m_params->m_alpha = x < m_constrainedParamLBound ? m_constrainedParamLBound : x > m_constrainedParamUBound ? m_constrainedParamUBound : x;
	}
}

void ARM_SmileCalibration_SABR::InitfGuess(std::vector<double>& fguess)
{
	if(m_params->m_levels.size() > 0)
	{
		fguess.resize(m_params->m_calibalpha 
					+ m_params->m_calibbeta 
					+ (m_params->m_calibrho && m_params->m_calibeachrho == false)
					+ (m_params->m_calibrho && m_params->m_calibeachrho) * m_params->m_rhos.size()
					+ m_params->m_calibnu 
					+ m_params->m_calibshift * m_params->m_shifts.size());

		int k = 0;

		if(m_params->m_calibalpha) fguess[k++] = 0.1;
		if(m_params->m_calibbeta) fguess[k++] = 0.5;
		if(m_params->m_calibrho && m_params->m_calibeachrho == false) fguess[k++] = -0.1;
		if(m_params->m_calibnu) fguess[k++] = 0.2;

		if(m_params->m_calibeachrho && m_params->m_calibshift) 
		{
			for(int i = 0; i < m_params->m_rhos.size(); i++) 
			{
				fguess[k++] = -0.1;
				fguess[k++] = 0.8;
			}
		}
		else if(m_params->m_calibeachrho)
		{
			for(int i = 0; i < m_params->m_rhos.size(); i++) 
			{
				fguess[k++] = -0.1;
			}
		}
		else if(m_params->m_calibshift)
		{
			for(int i = 0; i < m_params->m_shifts.size(); i++) 
			{
				fguess[k++] = 0.8;
			}
		}
		if(m_params->m_calibshift == false)
		{
			for(int i = 0; i < m_params->m_shifts.size(); i++) m_params->m_shifts[i] = m_params->m_shift;
		}
	}
	else
	{
		fguess.resize((m_params->m_calibalpha && m_calibConstraint == false) + m_params->m_calibbeta + m_params->m_calibrho + m_params->m_calibnu + m_params->m_calibshift);

		int k = 0;

		if(m_params->m_calibalpha && m_calibConstraint == false) fguess[k++] = 0.1;
		if(m_params->m_calibbeta) fguess[k++] = 0.5;
		if(m_params->m_calibrho) fguess[k++] = -0.1;
		if(m_params->m_calibnu) fguess[k++] = 0.2;
		if(m_params->m_calibshift) fguess[k++] = 0.8;
	}

}

void ARM_SmileCalibration_SABR::SetLocalParam(double param[])
{
	if(m_params->m_calibeachrho && m_params->m_calibshift)
	{
		double rho = param[m_locidx] < -0.999 ? -0.999 : param[m_locidx] > 0.999 ? 0.999 : param[m_locidx];
		m_locidx ++;
		double sh = fabs(param[m_locidx]) < K_DOUBLE_TOL ? 1e-3 : param[m_locidx] < -5. ? -5. : param[m_locidx] > 1. ? 1. : param[m_locidx];
		m_locidx ++;

		for(int k = m_idx; k < m_params->m_rhos.size(); k++)
		{
			m_params->m_rhos[k] = rho;
			m_params->m_shifts[k] = sh;
		}
	}
	else if(m_params->m_calibeachrho)
	{
		double rho = param[m_locidx] < -0.999 ? -0.999 : param[m_locidx] > 0.999 ? 0.999 : param[m_locidx];
		m_locidx ++;

		for(int k = m_idx; k < m_params->m_rhos.size(); k++)
		{
			m_params->m_rhos[k] = rho;
		}
	}
	else if(m_params->m_calibshift)
	{
		double sh = fabs(param[m_locidx]) < K_DOUBLE_TOL ? 1e-3 : param[m_locidx] < -5. ? -5. : param[m_locidx] > 1. ? 1. : param[m_locidx];
		m_locidx ++;

		for(int k = m_idx; k < m_params->m_shifts.size(); k++)
		{
			m_params->m_shifts[k] = sh;
		}
	}
}

void ARM_SmileCalibration_Heston::InitOneReset(ARM_SmileCalibration_Params * params)
{
	m_params = dynamic_cast<ARM_SmileCalibration_Params_Heston*>(params);

	m_params->m_shifts.resize(1, m_params->m_shift);
	m_params->m_levels.resize(1, m_params->m_level);
	m_params->m_rhos.resize(1, m_params->m_rho);

	if(m_params->m_sigmas.size() == 0)
	{
		ARM_THROW(ERR_INVALID_ARGUMENT, " sigma size should be at least one when calibrating mixte heston");
	}

	m_params->m_calibeachrho	= false;
	m_params->m_bootstraplevel	= false;

	m_constrainedParamLBound	= 0.001;
	m_constrainedParamUBound	= 5.;

	m_thetaLBound	= 1e-5;
	m_nuLBound		= 1e-5;
}

void ARM_SmileCalibration_Heston::InitMultipleResets(ARM_SmileCalibration_Params * params)
{
	m_params = dynamic_cast<ARM_SmileCalibration_Params_Heston*>(params);

	m_params->m_shifts.resize(m_resets.size(), m_params->m_shift);
	m_params->m_levels.resize(m_resets.size(), 1.);
	m_params->m_rhos.resize(m_params->m_calibeachrho ? m_resets.size() : 1, m_params->m_rho);

	if(m_params->m_sigmas.size() == 0)
	{
		ARM_THROW(ERR_INVALID_ARGUMENT, " sigma size should be at least one when calibrating mixte heston");
	}
	if(m_params->m_sigmas.size() == 1)
	{
		m_params->m_sigmas.resize(m_resets.size(), m_params->m_sigmas[0]);
	}
	else if(m_params->m_sigmas.size() < m_resets.size())
	{
		ARM_THROW(ERR_INVALID_ARGUMENT, " sigma size should be at least equal to reset size when calibrating mixte heston");
	}

	m_constrainedParamLBound	= 0.001;
	m_constrainedParamUBound	= 5.;

	m_thetaLBound	= 1e-5;
	m_nuLBound		= 1e-5;
}

void ARM_SmileCalibration_Heston::SetParam(double param[])
{
	int k = 0;

	if(m_params->m_calibkappa) 
	{
		m_params->m_kappa = param[k] < 1e-5 ? 1e-5 : param[k];
		k++;
	}

	if(m_params->m_calibtheta)
	{
		m_params->m_theta = param[k] < m_thetaLBound ? m_thetaLBound : param[k];
		k++;
	}

	if(m_params->m_calibrho && m_params->m_calibeachrho == false)
	{
		m_params->m_rhos[0] = param[k] < -0.999 ? -0.999 : param[k] > 0.999 ? 0.999 : param[k];
		k++;
	}

	if(m_params->m_calibnu)
	{
		m_params->m_nu = param[k] < m_nuLBound ? m_nuLBound : param[k];
		k++;
	}

	m_locidx = m_params->m_calibkappa + m_params->m_calibtheta + (m_params->m_calibrho && m_params->m_calibeachrho == false) + m_params->m_calibnu;
}

void ARM_SmileCalibration_Heston::SetConstrainedParam(double x)
{
	for(int k = m_idx; k < m_params->m_levels.size(); k++) m_params->m_levels[k] = x < m_constrainedParamLBound ? m_constrainedParamLBound : x > m_constrainedParamUBound ? m_constrainedParamUBound : x; 
}

void ARM_SmileCalibration_Heston::SetLocalParam(double param[])
{
	if(m_params->m_calibeachrho && m_params->m_calibshift)
	{
		double rho = param[m_locidx] < -0.999 ? -0.999 : param[m_locidx] > 0.999 ? 0.999 : param[m_locidx];
		m_locidx ++;
		double sh = fabs(param[m_locidx]) < K_DOUBLE_TOL ? 1e-3 : param[m_locidx] < -5. ? -5. : param[m_locidx] > 1. ? 1. : param[m_locidx];
		m_locidx ++;

		for(int k = m_idx; k < m_params->m_rhos.size(); k++)
		{
			m_params->m_rhos[k] = rho;
			m_params->m_shifts[k] = sh;
		}
	}
	else if(m_params->m_calibeachrho)
	{
		double rho = param[m_locidx] < -0.999 ? -0.999 : param[m_locidx] > 0.999 ? 0.999 : param[m_locidx];
		m_locidx ++;

		for(int k = m_idx; k < m_params->m_rhos.size(); k++)
		{
			m_params->m_rhos[k] = rho;
		}
	}
	else if(m_params->m_calibshift)
	{
		double sh = fabs(param[m_locidx]) < K_DOUBLE_TOL ? 1e-3 : param[m_locidx] < -5. ? -5. : param[m_locidx] > 1. ? 1. : param[m_locidx];
		m_locidx ++;

		for(int k = m_idx; k < m_params->m_shifts.size(); k++)
		{
			m_params->m_shifts[k] = sh;
		}
	}
}

void ARM_SmileCalibration_Heston::InitfGuess(std::vector<double>& fguess)
{
	fguess.resize(m_params->m_calibkappa 
				+ m_params->m_calibtheta 
				+ (m_params->m_calibrho && m_params->m_calibeachrho == false)
				+ (m_params->m_calibrho && m_params->m_calibeachrho) * m_params->m_rhos.size()
				+ m_params->m_calibnu 
				+ m_params->m_calibshift * m_params->m_shifts.size());

	int k = 0;

	if(m_params->m_calibkappa) fguess[k++] = 0.1;
	if(m_params->m_calibtheta) fguess[k++] = m_params->m_v0;
	if(m_params->m_calibrho && m_params->m_calibeachrho == false) fguess[k++] = -0.1;
	if(m_params->m_calibnu) fguess[k++] = 0.2;
	if(m_params->m_calibeachrho && m_params->m_calibshift) 
	{
		for(int i = 0; i < m_params->m_rhos.size(); i++) 
		{
			fguess[k++] = -0.1;
			fguess[k++] = 0.8;
		}
	}
	else if(m_params->m_calibeachrho)
	{
		for(int i = 0; i < m_params->m_rhos.size(); i++) 
		{
			fguess[k++] = -0.1;
		}
	}
	else if(m_params->m_calibshift)
	{
		for(int i = 0; i < m_params->m_shifts.size(); i++) 
		{
			fguess[k++] = 0.8;
		}
	}
	if(m_params->m_calibshift == false)
	{
		for(int i = 0; i < m_params->m_shifts.size(); i++) m_params->m_shifts[i] = m_params->m_shift;
	}
}

void ARM_SmileCalibration_NormalHeston::InitOneReset(ARM_SmileCalibration_Params * params)
{
	m_params = dynamic_cast<ARM_SmileCalibration_Params_Heston*>(params);

	m_constrainedParamLBound	= 1e-6;
	m_constrainedParamUBound	= 5.;

	AdjustVolToBP();

	m_thetaLBound	= 1e-8;
	m_nuLBound		= 1e-8;

	m_params->m_levels.resize(1,m_params->m_level);
	m_params->m_rhos.resize(1, m_params->m_rho);

	m_params->m_calibeachrho	= false;
	m_params->m_calibshift		= false;
	m_params->m_bootstraplevel	= false;
}

void ARM_SmileCalibration_NormalHeston::InitMultipleResets(ARM_SmileCalibration_Params * params)
{
	m_params = dynamic_cast<ARM_SmileCalibration_Params_Heston*>(params);

	m_constrainedParamLBound	= 1e-6;
	m_constrainedParamUBound	= 5.;

	AdjustVolToBP();

	m_thetaLBound	= 1e-8;
	m_nuLBound		= 1e-8;

	m_params->m_levels.resize(m_resets.size(), m_params->m_level);
	m_params->m_rhos.resize(m_params->m_calibeachrho ? m_resets.size() : 1, m_params->m_rho);

	m_params->m_calibshift		= false;
}

void ARM_SmileCalibration_NormalHeston::InitfGuess(std::vector<double>& fguess)
{
	fguess.resize(m_params->m_calibkappa 
				+ m_params->m_calibtheta 
				+ (m_params->m_calibrho && m_params->m_calibeachrho == false)
				+ (m_params->m_calibrho && m_params->m_calibeachrho) * m_params->m_rhos.size()
				+ m_params->m_calibnu 
				+ m_params->m_calibshift * m_params->m_shifts.size());

	int k = 0;

	if(m_params->m_calibkappa) fguess[k++] = 0.1;
	if(m_params->m_calibtheta) fguess[k++] = m_params->m_v0;
	if(m_params->m_calibrho && m_params->m_calibeachrho == false) fguess[k++] = 0.1;
	if(m_params->m_calibnu) fguess[k++] = 0.01;
	if(m_params->m_calibeachrho) for(int i = 0; i < m_params->m_rhos.size(); i++) fguess[k+i] = 0.1;
}

void ARM_SmileCalibration_BiSABR::InitOneReset(ARM_SmileCalibration_Params * params)
{
	m_params = dynamic_cast<ARM_SmileCalibration_Params_BiSABR*>(params);

	AdjustVolToBP();

	m_constrainedParamLBound	= -0.999;
	m_constrainedParamUBound	= 0.999;
}

void ARM_SmileCalibration_BiSABR::InitMultipleResets(ARM_SmileCalibration_Params * params)
{
	m_params = dynamic_cast<ARM_SmileCalibration_Params_BiSABR*>(params);

	AdjustVolToBP();

	m_constrainedParamLBound	= -0.999;
	m_constrainedParamUBound	= 0.999;	
}

void ARM_SmileCalibration_BiSABR::SetParam(double param[])
{
	int k = 0;

	if(m_params->m_calibrhoSS && m_calibConstraint == false) 
	{
		m_params->m_rhoS1S2 = param[k] < -0.999 ? -0.999 : param[k] > 0.999 ? 0.999 : param[k];
		k++;
	}

	if(m_params->m_calibrhoVV)
	{
		m_params->m_rhoV1V2 = param[k] < -0.999 ? -0.999 : param[k] > 0.999 ? 0.999 : param[k];
		k++;
	}

	if(m_params->m_calibrhoSV)
	{
		m_params->m_rhoS1V2 = param[k] < -0.999 ? -0.999 : param[k] > 0.999 ? 0.999 : param[k];
		k++;
		m_params->m_rhoS2V1 = param[k] < -0.999 ? -0.999 : param[k] > 0.999 ? 0.999 : param[k];
		k++;
	}
}

void ARM_SmileCalibration_BiSABR::SetConstrainedParam(double x)
{
	if(m_params->m_calibrhoSS == false) return;

	m_params->m_rhoS1S2 = x < m_constrainedParamLBound ? m_constrainedParamLBound : x > m_constrainedParamUBound ? m_constrainedParamUBound : x; 
}

void ARM_SmileCalibration_BiSABR::InitfGuess(std::vector<double>& fguess)
{
	fguess.resize((m_params->m_calibrhoSS && m_calibConstraint == false) + m_params->m_calibrhoVV + 2*m_params->m_calibrhoSV);

	int k = 0;

	if(m_params->m_calibrhoSS && m_calibConstraint == false) fguess[k++] = 0.8;

	if(m_params->m_calibrhoVV) fguess[k++] = 0.5;

	if(m_params->m_calibrhoSV) 
	{
		fguess[k++] = -0.5;
		fguess[k++] = -0.3;
	}
}

void ARM_SmileCalibration_Heston2b::InitOneReset(ARM_SmileCalibration_Params * params)
{
	m_params = dynamic_cast<ARM_SmileCalibration_Params_Heston2b*>(params);

	m_params->m_shifts.resize(1, m_params->m_shift);
	m_params->m_levels.resize(1, m_params->m_level);
	m_params->m_rhos1.resize(1, m_params->m_rho1);
	m_params->m_rhos2.resize(1, m_params->m_rho2);

	m_params->m_calibeachrho1	= false;
	m_params->m_calibeachrho2	= false;
	m_params->m_bootstraplevel	= false;

	m_constrainedParamLBound	= 0.001;
	m_constrainedParamUBound	= 5.;

	m_thetaLBound	= 1e-5;
	m_nuLBound		= 1e-5;
}

void ARM_SmileCalibration_Heston2b::InitMultipleResets(ARM_SmileCalibration_Params * params)
{
	m_params = dynamic_cast<ARM_SmileCalibration_Params_Heston2b*>(params);

	m_params->m_shifts.resize(m_resets.size(), m_params->m_shift);
	m_params->m_levels.resize(m_resets.size(), 1.);
	m_params->m_rhos1.resize(m_params->m_calibeachrho1 ? m_resets.size() : 1, m_params->m_rho1);
	m_params->m_rhos2.resize(m_params->m_calibeachrho2 ? m_resets.size() : 1, m_params->m_rho2);

	m_constrainedParamLBound	= 0.001;
	m_constrainedParamUBound	= 5.;

	m_thetaLBound	= 1e-5;
	m_nuLBound		= 1e-5;
}

void ARM_SmileCalibration_Heston2b::SetParam(double param[])
{
	int k = 0;

	if(m_params->m_calibkappa1) 
	{
		m_params->m_kappa1 = param[k] < 1e-5 ? 1e-5 : param[k];
		k++;
	}

	if(m_params->m_calibtheta1)
	{
		m_params->m_theta1 = param[k] < m_thetaLBound ? m_thetaLBound : param[k];
		k++;
	}

	if(m_params->m_calibrho1 && m_params->m_calibeachrho1 == false)
	{
		m_params->m_rhos1[0] = param[k] < -0.999 ? -0.999 : param[k] > 0.999 ? 0.999 : param[k];
		k++;
	}

	if(m_params->m_calibnu1)
	{
		m_params->m_nu1 = param[k] < m_nuLBound ? m_nuLBound : param[k];
		k++;
	}

	if(m_params->m_calibkappa2) 
	{
		m_params->m_kappa2 = param[k] < 1e-5 ? 1e-5 : param[k];
		k++;
	}

	if(m_params->m_calibtheta2)
	{
		m_params->m_theta2 = param[k] < m_thetaLBound ? m_thetaLBound : param[k];
		k++;
	}

	if(m_params->m_calibrho2 && m_params->m_calibeachrho2 == false)
	{
		m_params->m_rhos2[0] = param[k] < -0.999 ? -0.999 : param[k] > 0.999 ? 0.999 : param[k];
		k++;
	}

	if(m_params->m_calibnu2)
	{
		m_params->m_nu2 = param[k] < m_nuLBound ? m_nuLBound : param[k];
		k++;
	}

	m_locidx = m_params->m_calibkappa1 + m_params->m_calibtheta1 + (m_params->m_calibrho1 && m_params->m_calibeachrho1 == false) + m_params->m_calibnu1
			+ m_params->m_calibkappa2 + m_params->m_calibtheta2 + (m_params->m_calibrho2 && m_params->m_calibeachrho2 == false) + m_params->m_calibnu2;
}

void ARM_SmileCalibration_Heston2b::SetConstrainedParam(double x)
{
	for(int k = m_idx; k < m_params->m_levels.size(); k++) m_params->m_levels[k] = x < m_constrainedParamLBound ? m_constrainedParamLBound : x > m_constrainedParamUBound ? m_constrainedParamUBound : x; 
}

void ARM_SmileCalibration_Heston2b::SetLocalParam(double param[])
{
	if(m_params->m_calibeachrho1 && m_params->m_calibeachrho2)
	{
		double rho1 = param[m_locidx] < -0.999 ? -0.999 : param[m_locidx] > 0.999 ? 0.999 : param[m_locidx];
		m_locidx ++;
		double rho2 = param[m_locidx] < -0.999 ? -0.999 : param[m_locidx] > 0.999 ? 0.999 : param[m_locidx];
		m_locidx ++;
		double sh;
		if(m_params->m_calibshift)
		{	
			sh = fabs(param[m_locidx]) < K_DOUBLE_TOL ? 1e-3 : param[m_locidx] < -5. ? -5. : param[m_locidx] > 1. ? 1. : param[m_locidx];
			m_locidx ++;
		}

		for(int k = m_idx; k < m_params->m_rhos1.size(); k++)
		{
			m_params->m_rhos1[k] = rho1;
			m_params->m_rhos2[k] = rho2;
			if(m_params->m_calibshift) m_params->m_shifts[k] = sh;
		}
	}
	else if(m_params->m_calibeachrho1)
	{
		double rho = param[m_locidx] < -0.999 ? -0.999 : param[m_locidx] > 0.999 ? 0.999 : param[m_locidx];
		m_locidx ++;

		double sh;
		if(m_params->m_calibshift)
		{	
			sh = fabs(param[m_locidx]) < K_DOUBLE_TOL ? 1e-3 : param[m_locidx] < -5. ? -5. : param[m_locidx] > 1. ? 1. : param[m_locidx];
			m_locidx ++;
		}

		for(int k = m_idx; k < m_params->m_rhos1.size(); k++)
		{
			m_params->m_rhos1[k] = rho;
			if(m_params->m_calibshift) m_params->m_shifts[k] = sh;
		}
	}
	else if(m_params->m_calibeachrho2)
	{
		double rho = param[m_locidx] < -0.999 ? -0.999 : param[m_locidx] > 0.999 ? 0.999 : param[m_locidx];
		m_locidx ++;

		double sh;
		if(m_params->m_calibshift)
		{	
			sh = fabs(param[m_locidx]) < K_DOUBLE_TOL ? 1e-3 : param[m_locidx] < -5. ? -5. : param[m_locidx] > 1. ? 1. : param[m_locidx];
			m_locidx ++;
		}

		for(int k = m_idx; k < m_params->m_rhos2.size(); k++)
		{
			m_params->m_rhos2[k] = rho;
			if(m_params->m_calibshift) m_params->m_shifts[k] = sh;
		}
	}
	else if(m_params->m_calibshift)
	{
		double sh = fabs(param[m_locidx]) < K_DOUBLE_TOL ? 1e-3 : param[m_locidx] < -5. ? -5. : param[m_locidx] > 1. ? 1. : param[m_locidx];
		m_locidx ++;

		for(int k = m_idx; k < m_params->m_shifts.size(); k++)
		{
			m_params->m_shifts[k] = sh;
		}
	}
}

void ARM_SmileCalibration_Heston2b::InitfGuess(std::vector<double>& fguess)
{
	fguess.resize(m_params->m_calibkappa1 
				+ m_params->m_calibtheta1 
				+ (m_params->m_calibrho1 && m_params->m_calibeachrho1 == false)
				+ (m_params->m_calibrho1 && m_params->m_calibeachrho1) * m_params->m_rhos1.size()
				+ m_params->m_calibnu1
				+ m_params->m_calibkappa2 
				+ m_params->m_calibtheta2 
				+ (m_params->m_calibrho2 && m_params->m_calibeachrho2 == false)
				+ (m_params->m_calibrho2 && m_params->m_calibeachrho2) * m_params->m_rhos2.size()
				+ m_params->m_calibnu2
				+ m_params->m_calibshift * m_params->m_shifts.size());

	int k = 0;

	if(m_params->m_calibkappa1) fguess[k++] = 0.1;
	if(m_params->m_calibtheta1) fguess[k++] = m_params->m_v01;
	if(m_params->m_calibrho1 && m_params->m_calibeachrho1 == false) fguess[k++] = -0.1;
	if(m_params->m_calibnu1) fguess[k++] = 0.2;
	if(m_params->m_calibkappa2) fguess[k++] = 0.1;
	if(m_params->m_calibtheta2) fguess[k++] = m_params->m_v02;
	if(m_params->m_calibrho2 && m_params->m_calibeachrho2 == false) fguess[k++] = -0.3;
	if(m_params->m_calibnu2) fguess[k++] = 0.3;

	if(m_params->m_calibeachrho1 && m_params->m_calibeachrho2) 
	{
		for(int i = 0; i < m_params->m_rhos1.size(); i++) 
		{
			fguess[k++] = -0.1;
			fguess[k++] = -0.3;
			if(m_params->m_calibshift) fguess[k++] = 0.8;
		}
	}
	else if(m_params->m_calibeachrho1)
	{
		for(int i = 0; i < m_params->m_rhos1.size(); i++) 
		{
			fguess[k++] = -0.1;
			if(m_params->m_calibshift) fguess[k++] = 0.8;
		}
	}
	else if(m_params->m_calibeachrho2)
	{
		for(int i = 0; i < m_params->m_rhos2.size(); i++) 
		{
			fguess[k++] = -0.3;
			if(m_params->m_calibshift) fguess[k++] = 0.8;
		}
	}
	else if(m_params->m_calibshift)
	{
		for(int i = 0; i < m_params->m_shifts.size(); i++) 
		{
			fguess[k++] = 0.8;
		}
	}
	if(m_params->m_calibshift == false)
	{
		for(int i = 0; i < m_params->m_shifts.size(); i++) m_params->m_shifts[i] = m_params->m_shift;
	}
}

void ARM_SmileCalibration_SABR2beta::InitOneReset(ARM_SmileCalibration_Params * params)
{
	m_params = dynamic_cast<ARM_SmileCalibration_Params_SABR2beta*>(params);

	m_constrainedParamLBound	= 1e-8;
	m_constrainedParamUBound	= 2.;

}

void ARM_SmileCalibration_SABR2beta::InitMultipleResets(ARM_SmileCalibration_Params * params)
{
	ARM_THROW( ERR_INVALID_ARGUMENT, ARM_USERNAME + " : _SABR2beta::InitMultipleResets not available");
}

void ARM_SmileCalibration_SABR2beta::SetParam(double param[])
{
	int k = 0;

	if(m_params->m_calibalpha)
	{
		if(m_calibConstraint == false) 
		{
			m_params->m_alpha = param[k] < 1e-8 ? 1e-8 : param[k];
			k++;
		}
	}

	if(m_params->m_calibbeta1)
	{
		m_params->m_beta1 = param[k] < 0.05 ? 0.05 : param[k] > 0.5 ? 0.5 : param[k];
		k++;
	}

	if(m_params->m_calibbeta2)
	{
		m_params->m_beta2 = param[k] < 0.5 ? 0.5 : param[k] > 1. ? 1. : param[k];
		k++;
	}

	if(m_params->m_calibrho) 
	{
		m_params->m_rho = param[k] < -0.999 ? -0.999 : param[k] > 0.999 ? 0.999 : param[k];
		k++;
	}

	if(m_params->m_calibnu) 
	{
		m_params->m_nu = param[k] < 1e-8 ? 1e-8 : param[k];
		k++;
	}

	if(m_params->m_calibzero) 
	{
		m_params->m_zero = param[k];
		k++;
	}

	if(m_params->m_caliblambda) 
	{
		m_params->m_lambda = param[k];
		k++;
	}
}

void ARM_SmileCalibration_SABR2beta::SetConstrainedParam(double x)
{
	if(m_params->m_calibalpha) m_params->m_alpha = x < m_constrainedParamLBound ? m_constrainedParamLBound : x > m_constrainedParamUBound ? m_constrainedParamUBound : x;
}

void ARM_SmileCalibration_SABR2beta::InitfGuess(std::vector<double>& fguess)
{
	fguess.resize(	m_params->m_calibalpha + 
					m_params->m_calibbeta1 + 
					m_params->m_calibbeta2 + 
					m_params->m_calibrho + 
					m_params->m_calibnu + 
					m_params->m_calibzero +
					m_params->m_caliblambda);

	int k = 0;

	if(m_params->m_calibalpha && m_calibConstraint == false) fguess[k++] = 0.1;
	if(m_params->m_calibbeta1) fguess[k++] = 0.25;
	if(m_params->m_calibbeta2) fguess[k++] = 0.55;
	if(m_params->m_calibrho) fguess[k++] = -0.1;
	if(m_params->m_calibnu) fguess[k++] = 0.2;
	if(m_params->m_calibzero) fguess[k++] = -0.25/100.;
	if(m_params->m_caliblambda) fguess[k++] = 100.;
}

void ARM_SmileCalibration_SABR2beta::SetLocalParam(double param[])
{
}

void ARM_SmileCalibration_Merton::InitOneReset(ARM_SmileCalibration_Params * params)
{
	m_params = dynamic_cast<ARM_SmileCalibration_Params_Merton*>(params);

	m_constrainedParamLBound	= 1e-8;
	m_constrainedParamUBound	= 2.;
}

void ARM_SmileCalibration_Merton::InitMultipleResets(ARM_SmileCalibration_Params * params)
{
	m_params = dynamic_cast<ARM_SmileCalibration_Params_Merton*>(params);

	m_constrainedParamLBound	= 1e-8;
	m_constrainedParamUBound	= 2.;
}

void ARM_SmileCalibration_Merton::SetParam(double param[])
{
	int k = 0;

	if(m_params->m_calibsigma && m_calibConstraint == false) 
	{
		m_params->m_sigma = param[k] < 1e-8 ? 1e-8 : param[k] > 2. ? 2. : param[k];
		k++;
	}

	if(m_params->m_caliblambda1)
	{
		m_params->m_lambda1 = param[k] < 0. ? 0. : param[k];
		k++;
	}

	if(m_params->m_calibU1)
	{
		m_params->m_U1 = param[k] < 0. ? 0. : param[k] > 2. ? 2. : param[k];
		k++;
	}

	if(m_params->m_caliblambda2)
	{
		m_params->m_lambda2 = param[k] < 0. ? 0. : param[k];
		k++;
	}

	if(m_params->m_calibU2)
	{
		m_params->m_U2 = param[k] < 0. ? 0. : param[k] > 2. ? 2. : param[k];
		k++;
	}
}

void ARM_SmileCalibration_Merton::SetConstrainedParam(double x)
{
	if(m_params->m_calibsigma)
	{
		m_params->m_sigma = x < m_constrainedParamLBound ? m_constrainedParamLBound : x > m_constrainedParamUBound ? m_constrainedParamUBound : x; 
	}
}

void ARM_SmileCalibration_Merton::SetLocalParam(double param[])
{
}

void ARM_SmileCalibration_Merton::InitfGuess(std::vector<double>& fguess)
{
	fguess.resize((m_params->m_calibsigma && m_calibConstraint == false) 
				+ m_params->m_caliblambda1
				+ m_params->m_calibU1
				+ m_params->m_caliblambda2
				+ m_params->m_calibU2);

	int k = 0;

	if(m_params->m_calibsigma && m_calibConstraint == false) fguess[k++] = 0.1;
	if(m_params->m_caliblambda1) fguess[k++] = 0.75;
	if(m_params->m_calibU1) fguess[k++] = 0.1;
	if(m_params->m_caliblambda2) fguess[k++] = 0.25;
	if(m_params->m_calibU2) fguess[k++] = 0.1;
}

double ARM_SmileCalibration_Merton::ImpliedVol(double strike, int sens)
{
	double price = MertonOption(m_forwards[m_idx], strike, sens, m_resets[m_idx],
							m_params->m_sigma, m_params->m_lambda1, m_params->m_U1, m_params->m_lambda2, m_params->m_U2);

	return LocalBSImpliedVol(m_forwards[m_idx], strike, m_resets[m_idx], sens, price);
}

void ARM_SmileCalibration_Spread2Heston::Init(double reset,
											  double fwd1, const std::vector<double>& mktvols1, const std::vector<double>& strikes1,
											  double constrvol1, double constrK1,
											  double fwd2, const std::vector<double>& mktvols2, const std::vector<double>& strikes2,
											  double constrvol2, double constrK2,
											  const std::vector<double>& mktvolspread, const std::vector<double>& strikespread, 
											  double constrvolspread, double constrKspread,
											  ARM_SmileCalibration_Params_Spread2Heston * params)
{
	m_reset			= reset;
	m_fwd1			= fwd1;
	m_mktvols1		= mktvols1;
	m_strikes1		= strikes1;
	m_constrvol1	= constrvol1;
	m_constrK1		= constrK1;
	m_fwd2			= fwd2;
	m_mktvols2		= mktvols2;
	m_strikes2		= strikes2;
	m_constrvol2	= constrvol2;
	m_constrK2		= constrK2;
	m_mktvolspread	= mktvolspread;
	std::transform(m_mktvolspread.begin(), m_mktvolspread.end(), m_mktvolspread.begin(),
               std::bind1st(std::multiplies<double>(),1/VolToBP));
	m_strikespread	= strikespread;
	m_constrVolSpread	= constrvolspread / VolToBP;
	m_constrKSpread		= constrKspread;
	m_params		= params;
}

bool ARM_SmileCalibration_Spread2Heston::Calibrate()
{
	CalibSpreadFunc func(this);

	double info[LM_INFO_SZ];
	double opts[LM_OPTS_SZ];
	opts[0] = LM_INIT_MU;
	opts[1] = opts[2] = 1E-6;
	opts[3] = 1E-8;	
	opts[4] = 1E-4;

	int size = m_mktvolspread.size();

	std::vector<double> fguess(2), hx(size, 0.);

	fguess[0] = 0.4;
	fguess[1] = 0.55;

	int status = LEVMARQMinization_WithNumDerivatives(func, fguess, hx, info, 250, opts);

	m_calibsuccess = info[1] < info[0] / 10. ? true : false;

	m_quaderror = info[1];

	if(m_calibsuccess == false || fabs(m_params->m_correl - 0.99) < K_DOUBLE_TOL)
	{
		CalibSpreadConstraint();
	}

	return m_calibsuccess;

}

void ARM_SmileCalibration_Spread2Heston::CalibrateCorrel(double reset, double fwd1, double fwd2, 
														 double mktvolspread, double strikespread,
														 ARM_SmileCalibration_Params_Spread2Heston * params)
{
	m_reset			= reset;
	m_fwd1			= fwd1;
	m_fwd2			= fwd2;
	m_constrVolSpread	= mktvolspread / VolToBP;
	m_constrKSpread		= strikespread;
	m_params		= params;
	
	CalibSpreadConstraint();
}

void ARM_SmileCalibration_Spread2Heston::CalibFwd()
{
	ARM_SmileCalibration_Heston calibtool;

	m_params->m_params1.build(m_params->m_params1.v0(), m_params->m_params1.kappa(), m_params->m_params1.theta(),
		0., m_params->m_params1.nu(), m_params->m_shift1, 0.5, false, m_params->m_params1.calibtheta(), true, true, false);

	calibtool.Init(m_reset, m_fwd1, m_mktvols1, m_strikes1, true, m_constrK1, m_constrvol1, &(m_params->m_params1));
	calibtool.Calibrate();

	m_params->m_params2.build(m_params->m_params1.v0(), m_params->m_params1.kappa(), m_params->m_params1.theta(),
		0., m_params->m_params1.nu(), m_params->m_shift2, 0.5, false, false, true, false, false);

	calibtool.Init(m_reset, m_fwd2, m_mktvols2, m_strikes2, true, m_constrK2, m_constrvol2, &(m_params->m_params2));
	calibtool.Calibrate();
}

void ARM_SmileCalibration_Spread2Heston::SetParamsForSpreadCalib(double p[])
{
	m_params->m_shift1 = p[0] < 0.001 ? 0.001 : p[0] > 1. ? 1. : p[0];
	m_params->m_shift2 = p[1] < 0.001 ? 0.001 : p[1] > 1. ? 1. : p[1];
}

void ARM_SmileCalibration_Spread2Heston::CalibSpreadConstraint()
{
	CalibSpreadConstrFunc func(this);

	double best;

	m_params->m_correl = brentSolve(func, m_constrVolSpread, -0.99, 0.99, 1e-8, 50, 0, &best);
}

ARM_Spread2Heston_CAPCalibration_Params::~ARM_Spread2Heston_CAPCalibration_Params()
{
	DeletePointorVector<ARM_GP_Vector>(m_mktvols1);
	DeletePointorVector<ARM_GP_Vector>(m_strikes1);

	DeletePointorVector<ARM_GP_Vector>(m_mktvols2);
	DeletePointorVector<ARM_GP_Vector>(m_strikes2);
}

ARM_Spread2Heston_CAPCalibration_Params::ARM_Spread2Heston_CAPCalibration_Params(const std::vector<double>& fwdReset, 
	const std::vector<double>& fwd1, const ARM_GP_Matrix& mktvols1, const ARM_GP_Matrix& strikes1, 
	const std::vector<double>& fwd2, const ARM_GP_Matrix& mktvols2, const ARM_GP_Matrix& strikes2, 
	double v0, double kappa, double theta, bool calibtheta, bool caliballShift, bool constrainedCorrel) : 
	m_fwdReset(fwdReset), m_fwd1(fwd1), m_fwd2(fwd2),
	m_shift1(fwdReset.size()), m_shift2(fwdReset.size()), m_calibAllShift(caliballShift), m_isCalibrated(fwdReset.size(), 0),
	m_constrainedCorrel(constrainedCorrel)
{
	m_params1 = new ARM_SmileCalibration_Params_Heston [fwdReset.size()];
	m_params2 = new ARM_SmileCalibration_Params_Heston [fwdReset.size()];

	int k, i;

	for(k = 0; k < (int)fwdReset.size(); k++)
	{
		m_params1[k].build(v0, kappa, theta, 0., 0., 1., 1., false, calibtheta, true, true, false);
		m_params2[k].build(v0, kappa, theta, 0., 0., 1., 1., false, calibtheta, true, true, false);
	}

	m_mktvols1.resize(fwdReset.size());
	m_strikes1.resize(fwdReset.size());
	m_mktvols2.resize(fwdReset.size());
	m_strikes2.resize(fwdReset.size());
	
	m_atm1.resize(fwdReset.size());
	m_atm2.resize(fwdReset.size());

	for(k = 0; k < (int)fwdReset.size(); k++)
	{
		m_mktvols1[k] = new ARM_GP_Vector(mktvols1.cols());
		m_strikes1[k] = new ARM_GP_Vector(mktvols1.cols());
		for(i = 0; i < (int)mktvols1.cols(); i++)
		{
			(*m_mktvols1[k])[i] = mktvols1(k,i);
			(*m_strikes1[k])[i] = strikes1(k,i);

			if(fabs(strikes1(k,i) - fwd1[k]) < K_DOUBLE_TOL) m_atm1[k] = mktvols1(k,i);
		}

		m_mktvols2[k] = new ARM_GP_Vector(mktvols2.cols());
		m_strikes2[k] = new ARM_GP_Vector(mktvols2.cols());
		for(i = 0; i < (int)mktvols2.cols(); i++)
		{
			(*m_mktvols2[k])[i] = mktvols2(k,i);
			(*m_strikes2[k])[i] = strikes2(k,i);

			if(fabs(strikes2(k,i) - fwd2[k]) < K_DOUBLE_TOL) m_atm2[k] = mktvols2(k,i);
		}
	}
}

int ARM_Spread2Heston_CAPCalibration_Params::nbFwdToCalibrate(double firstCapReset, double lastCapReset)
{
	int k = 0, size = (int)m_fwdReset.size();

	while(m_fwdReset[k] < firstCapReset && k < size) k++;
	if(k > 0) k--;

	int i, nbFwds = 0;

	m_currFwdCalib.resize(0);

	for(i = k; i < size; i++)
	{
		if(m_fwdReset[i] < lastCapReset && m_isCalibrated[i] == 0)
		{
			m_currFwdCalib.push_back(i);
			nbFwds ++;
		}
		if(m_fwdReset[i] >= lastCapReset && m_isCalibrated[i] == 0)
		{
			m_currFwdCalib.push_back(i);
			nbFwds ++;
			break;
		}
	}

	return nbFwds;
}

void ARM_Spread2Heston_CAPCalibration_Params::Calibrate(double firstCapReset, double lastCapReset)
{
	ARM_SmileCalibration_Heston calibtool;

	int k, i, size = (int)m_currFwdCalib.size();

	for(k = 0; k < size; k++)
	{
		i = m_currFwdCalib[k];

		m_params1[i].build(m_params1[i].v0(), m_params1[i].kappa(), m_params1[i].theta(),
			0., m_params1[i].nu(), m_shift1[i], 0.5, false, m_params1[i].calibtheta(), true, true, false);

		calibtool.Init(m_fwdReset[i], m_fwd1[i], (*m_mktvols1[i]).GetValues(), (*m_strikes1[i]).GetValues(), true, m_fwd1[i], m_atm1[i], &m_params1[i]);
		calibtool.Calibrate();

		m_params2[i].build(m_params1[i].v0(), m_params1[i].kappa(), m_params1[i].theta(),
			0., m_params1[i].nu(), m_shift2[i], 0.5, false, false, true, false, false);

		calibtool.Init(m_fwdReset[i], m_fwd2[i], (*m_mktvols2[i]).GetValues(), (*m_strikes2[i]).GetValues(), true, m_fwd2[i], m_atm2[i], &m_params2[i]);
		calibtool.Calibrate();
	}
}

void ARM_Spread2Heston_CAPCalibration_Params::SetIsCalibrated(double firstCapReset, double lastCapReset)
{
	int k, size = (int)m_currFwdCalib.size();
	for(k = 0; k < size; k++) m_isCalibrated[m_currFwdCalib[k]] = 1;
}

void ARM_Spread2Heston_CAPCalibration_Params::getData(double reset, double& theta_, double& nu_, double& rho1_, double& rho2_,
													  double& shift1_, double& shift2_, double& level1_, double& level2_)
{
	int kinf, ksup;
	getInterpolatePoints(reset, kinf, ksup);
	double rinf = m_fwdReset[kinf];
	double rsup = m_fwdReset[ksup];
	double fac = kinf == ksup ? 0. : (reset - rinf) / (rsup - rinf);

	theta_	= kinf == ksup ? theta(kinf) : theta(kinf) + (theta(ksup) - theta(kinf)) * fac;
	nu_		= kinf == ksup ? nu(kinf) : nu(kinf) + (nu(ksup) - nu(kinf)) * fac;
	rho1_	= kinf == ksup ? rho1(kinf): rho1(kinf) + (rho1(ksup) - rho1(kinf)) * fac;
	shift1_	= kinf == ksup ? shift1(kinf): shift1(kinf) + (shift1(ksup) - shift1(kinf)) * fac;
	level1_	= kinf == ksup ? level1(kinf): level1(kinf) + (level1(ksup) - level1(kinf)) * fac;
	rho2_	= kinf == ksup ? rho2(kinf): rho2(kinf) + (rho2(ksup) - rho2(kinf)) * fac;
	shift2_	= kinf == ksup ? shift2(kinf): shift2(kinf) + (shift2(ksup) - shift2(kinf)) * fac;
	level2_	= kinf == ksup ? level2(kinf): level2(kinf) + (level2(ksup) - level2(kinf)) * fac;
}

void ARM_Spread2Heston_CAPCalibration_Params::getInterpolatePoints(double reset, int& kinf, int& ksup)
{
	int k = 0, size = m_fwdReset.size();
	
	while(m_fwdReset[k] < reset)
	{
		k++;
		if(k >= size) break;
	}

	kinf = k > 0 ? k-1 : k;
	
	k = size-1;
	while(m_fwdReset[k] > reset)
	{
		k--;
		if(k < 0) break;
	}

	ksup = k < size-1 ? k+1 : k;
}

void ARM_SmileCalibration_Spread2HestonCAP::Init(const std::vector<double>& resets, const std::vector<double>& dfs,
												 const std::vector<double>& fwd1, const std::vector<double>& fwd2,
												 const std::vector<double>& capprices, const std::vector<double>& capstrikes, 
												 double constrCapPrice, double constrCapStrike, 
												 ARM_Spread2Heston_CAPCalibration_Params * params)
{
	m_resets			= resets;
	m_dfs				= dfs;
	m_fwd1				= fwd1;
	m_fwd2				= fwd2;

	m_capprices			= capprices;
	m_capstrikes		= capstrikes;
	m_constrCapPrice	= constrCapPrice;
	m_constrCapStrike	= constrCapStrike;

	m_params			= params;

}

void ARM_SmileCalibration_Spread2HestonCAP::Calibrate()
{
	// initialisation
	CalibSpreadFunc func(this);

	double info[LM_INFO_SZ];
	double opts[LM_OPTS_SZ];
	opts[0] = LM_INIT_MU;
	opts[1] = opts[2] = 1E-4;
	opts[3] = 1E-6;	
	opts[4] = 1E-4;

	int k, size = m_capprices.size();

	std::vector<double> fguess, hx(size, 0.);

	int nbfwdcalib = m_params->nbFwdToCalibrate(m_resets[0], m_resets[m_resets.size()-1]);

	if(m_params->m_calibAllShift)
	{
		fguess.resize(nbfwdcalib*2 + (m_params->m_constrainedCorrel == false));
		k = 0;
		if(m_params->m_constrainedCorrel == false) 
			fguess[k++] = 0.9;

		for(; k < nbfwdcalib; k++)
		{
			fguess[2*k] = 0.4;
			fguess[2*k+1] = 0.55;
		}
	}
	else
	{
		fguess.resize(2 + (m_params->m_constrainedCorrel == false));

		k = 0;
		if(m_params->m_constrainedCorrel == false) 
			fguess[k++] = 0.9;

		fguess[k++] = 0.4;
		fguess[k++] = 0.55;
	}

	int status = LEVMARQMinization_WithNumDerivatives(func, fguess, hx, info, 250, opts);

	m_params->SetIsCalibrated(m_resets[0], m_resets[m_resets.size()-1]);
}

void ARM_SmileCalibration_Spread2HestonCAP::SetParamsForSpreadCalib(double p[])
{
	int k = 0;
	if(m_params->m_constrainedCorrel == false)
	{
		m_params->m_correl = p[k] < 0 ? 0. : p[k] > 0.99 ? 0.99 : p[k];
		k++;
	}
	if(m_params->m_calibAllShift)
	{
		for(; k < (int)m_params->m_currFwdCalib.size(); k++)
		{
			m_params->m_shift1[m_params->m_currFwdCalib[k]] = p[2*k] < 0.001 ? 0.001 : p[2*k] > 1. ? 1. : p[2*k];
			m_params->m_shift2[m_params->m_currFwdCalib[k]] = p[2*k+1] < 0.001 ? 0.001 : p[2*k+1] > 1. ? 1. : p[2*k+1];
		}
	}
	else
	{
		for(int i = 0; i < (int)m_params->m_currFwdCalib.size(); i++)
		{
			m_params->m_shift1[m_params->m_currFwdCalib[i]] = p[k] < 0.001 ? 0.001 : p[k] > 1. ? 1. : p[k];
			m_params->m_shift2[m_params->m_currFwdCalib[i]] = p[k+1] < 0.001 ? 0.001 : p[k+1] > 1. ? 1. : p[k+1];
		}
	}
}

void ARM_SmileCalibration_Spread2HestonCAP::CalibFwd()
{
	m_params->Calibrate(m_resets[0], m_resets[m_resets.size()-1]);
}

void ARM_SmileCalibration_Spread2HestonCAP::CalibSpreadConstraint()
{
	if(m_params->m_constrainedCorrel == false) return;

	CalibSpreadConstrFunc func(this);

	double best;

	m_params->m_correl = brentSolve(func, m_constrCapPrice, -0.99, 0.99, 1e-8, 50, 0, &best);
}

double ARM_SmileCalibration_Spread2HestonCAP::cap(double strike)
{
	int k, size = m_resets.size();

	double price = 0., caplet;
	double v0, kappa, theta, nu, shift1, shift2, level1, level2, rho1, rho2, correl;

	v0 = m_params->v0();
	kappa = m_params->kappa();
	correl = m_params->correl();

	for(k = 0; k < size; k++)
	{
		m_params->getData(m_resets[k], theta, nu, rho1, rho2, shift1, shift2, level1, level2);
		
		caplet = Spread2HestonVanilla(m_resets[k], m_fwd1[k], m_fwd2[k], strike, -1, 
					v0, kappa, theta, nu, rho1, rho2, shift1, shift2, level1, level2, correl);

		price += caplet * m_dfs[k];
	}

	return price;
}

void ARM_Spread2Heston_TOTEMCalibration(const std::vector<double>& TOTEMMaturities, const std::vector<double>& TOTEMStrikes,
										const std::vector<double>& TOTEMPrices,
										const std::vector<double>& FullSchedReset,
										const std::vector<double>& FullSchedAnnuity,
										const std::vector<double>& FullSchedFwd1, const std::vector<double>& FullSchedFwd2,
										const std::vector<double>& FwdCalibReset,
										const std::vector<double>& LongFwds, const ARM_GP_Matrix& LongVols, const ARM_GP_Matrix& LongK,
										const std::vector<double>& ShortFwds, const ARM_GP_Matrix& ShortVols, const ARM_GP_Matrix& ShortK,
										double v0, double kappa, double theta, bool constrCorrel,
										std::vector<double>& Mat, std::vector<double>& Theta, std::vector<double>& Nu, 
										std::vector<double>& LongRho, std::vector<double>& LongShift, std::vector<double>& LongLevel,
										std::vector<double>& ShortRho, std::vector<double>& ShortShift, std::vector<double>& ShortLevel,
										std::vector<double>& Correl)
{
	// check des size
	if(TOTEMMaturities.size() != TOTEMStrikes.size() || TOTEMMaturities.size() != TOTEMPrices.size() || TOTEMStrikes.size() != TOTEMPrices.size())
	{
		ARM_THROW(ERR_INVALID_ARGUMENT, ARM_USERNAME + "TOTEM Datas should have the same size");
	}

	if(FullSchedReset.size() != FullSchedAnnuity.size() || FullSchedReset.size() != FullSchedFwd1.size() || FullSchedAnnuity.size() != FullSchedFwd1.size())
	{
		ARM_THROW(ERR_INVALID_ARGUMENT, ARM_USERNAME + "Full Schedule Datas should have the same size");
	}

	if(FullSchedFwd1.size() != FullSchedFwd2.size())
	{
		ARM_THROW(ERR_INVALID_ARGUMENT, ARM_USERNAME + "Full Schedule Datas should have the same size");
	}

	if(LongFwds.size() != FwdCalibReset.size())
	{
		ARM_THROW(ERR_INVALID_ARGUMENT, ARM_USERNAME + "Long Fwds should have same size as Full Schedule");
	}
	if(LongVols.rows() != FwdCalibReset.size())
	{
		ARM_THROW(ERR_INVALID_ARGUMENT, ARM_USERNAME + "Long Vols Rows should have same size as Full Schedule");
	}
	if(LongK.rows() != FwdCalibReset.size())
	{
		ARM_THROW(ERR_INVALID_ARGUMENT, ARM_USERNAME + "Long Strikes Rows should have same size as Full Schedule");
	}
	if(LongVols.cols() != LongK.cols())
	{
		ARM_THROW(ERR_INVALID_ARGUMENT, ARM_USERNAME + "Long Strikes and Long Vols should have same number of cols");
	}

	if(ShortFwds.size() != FwdCalibReset.size())
	{
		ARM_THROW(ERR_INVALID_ARGUMENT, ARM_USERNAME + "Short Fwds should have same size as Full Schedule");
	}
	if(ShortVols.rows() != FwdCalibReset.size())
	{
		ARM_THROW(ERR_INVALID_ARGUMENT, ARM_USERNAME + "Short Vols Rows should have same size as Full Schedule");
	}
	if(ShortK.rows() != FwdCalibReset.size())
	{
		ARM_THROW(ERR_INVALID_ARGUMENT, ARM_USERNAME + "Short Strikes Rows should have same size as Full Schedule");
	}
	if(ShortVols.cols() != ShortK.cols())
	{
		ARM_THROW(ERR_INVALID_ARGUMENT, ARM_USERNAME + "Short Strikes and Short Vols should have same number of cols");
	}

	
	int nbmat = 1;
	int k, size = TOTEMMaturities.size();

	double currTOTEMMat = TOTEMMaturities[0];
	Mat.resize(1,currTOTEMMat/12.);

	for(k = 1; k < size; k++)
	{
		if(fabs(TOTEMMaturities[k] - currTOTEMMat) < K_DOUBLE_TOL) continue;

		currTOTEMMat = TOTEMMaturities[k];
		Mat.push_back(currTOTEMMat/12.);
		nbmat ++;
	}

	Correl.resize(nbmat);

	std::vector<double> currCalibPrices, currTOTEMStrikes, currTOTEMPrices, prevTOTEMStrikes(0), prevTOTEMPrices(0);
	std::vector<double> currReset, currFwd1, currFwd2, currAnnuity;
	double constrCapPrice, constrCapStrike = 0.;
	
	int schedSize = FullSchedReset.size(), fwdSize = FwdCalibReset.size();
	int i = 0, j = 0, idxMat = 0;
	int currTOTEMSize = 0;
	int currCAPSize = 0;

	ARM_Spread2Heston_CAPCalibration_Params params(FwdCalibReset, LongFwds, LongVols, LongK, 
		ShortFwds, ShortVols, ShortK, v0, kappa, theta, fabs(theta) < K_DOUBLE_TOL ? true : false, constrCorrel, true);

	ARM_SmileCalibration_Spread2HestonCAP calibtool;

	for(;;)
	{
		
		// récupération données totem pour une maturité donnée
		currTOTEMSize = 0;

		while(fabs(TOTEMMaturities[i + currTOTEMSize]/12 - Mat[idxMat]) < K_DOUBLE_TOL)
		{
			currTOTEMSize ++;
			if(i + currTOTEMSize == size) break;
		}

		currTOTEMStrikes.resize(currTOTEMSize);
		currTOTEMPrices.resize(currTOTEMSize);
		currCalibPrices.resize(currTOTEMSize);

		for(k = 0; k < currTOTEMSize; k++)
		{
			currTOTEMStrikes[k] = TOTEMStrikes[k+i];
			currTOTEMPrices[k] = TOTEMPrices[k+i];
			currCalibPrices[k] = TOTEMPrices[k+i];

			for(int kk = 0; kk < (int)prevTOTEMStrikes.size(); kk++)
			{
				if(fabs(prevTOTEMStrikes[kk] - currTOTEMStrikes[k]) < K_DOUBLE_TOL)
				{
					currCalibPrices[k] -= prevTOTEMPrices[kk];
					break;
				}
			}

			if(fabs(currTOTEMStrikes[k] - constrCapStrike) < K_DOUBLE_TOL)
				constrCapPrice = currCalibPrices[k];
		}

		// récupération des données pour le cap (on suppose fixing in advance)
		currCAPSize = 0;

		while(FullSchedReset[j + currCAPSize] < Mat[idxMat]-0.1)
		{
			currCAPSize ++;
			if(j + currCAPSize == schedSize) break;
		}

		currReset.resize(currCAPSize);
		currFwd1.resize(currCAPSize);
		currFwd2.resize(currCAPSize);
		currAnnuity.resize(currCAPSize);

		for(k = 0; k < currCAPSize; k++)
		{
			currReset[k]	= FullSchedReset[k+j];
			currFwd1[k]		= FullSchedFwd1[k+j];
			currFwd2[k]		= FullSchedFwd2[k+j];
			currAnnuity[k]	= FullSchedAnnuity[k+j];
		}
		
		// calibration
		calibtool.Init(currReset, currAnnuity, currFwd1, currFwd2, currCalibPrices, currTOTEMStrikes, constrCapPrice, constrCapStrike, &params);

		calibtool.Calibrate();

		// test
		double pr = 0.;
		for(k = 0; k < currTOTEMSize; k++)
		{
			pr = calibtool.cap(currTOTEMStrikes[k]);
		}

		Correl[idxMat] = params.correl();

		i += currTOTEMSize;
		j += currCAPSize;

		prevTOTEMStrikes = currTOTEMStrikes;
		prevTOTEMPrices = currTOTEMPrices;

		idxMat ++;
		if(i == size || j == schedSize) break;
	}

	Theta.resize(FwdCalibReset.size());
	Nu.resize(FwdCalibReset.size());
	LongShift.resize(FwdCalibReset.size());
	LongLevel.resize(FwdCalibReset.size());
	LongRho.resize(FwdCalibReset.size());
	ShortShift.resize(FwdCalibReset.size());
	ShortLevel.resize(FwdCalibReset.size());
	ShortRho.resize(FwdCalibReset.size());

	for(k = 0; k < (int)FwdCalibReset.size(); k++)
	{
		Theta[k] = params.theta(k);
		Nu[k] = params.nu(k);
		LongShift[k] = params.shift1(k);
		LongRho[k] = params.rho1(k);
		LongLevel[k] = params.level1(k);
		ShortShift[k] = params.shift2(k);
		ShortRho[k] = params.rho2(k);
		ShortLevel[k] = params.level2(k);
	}
}


CC_END_NAMESPACE()
