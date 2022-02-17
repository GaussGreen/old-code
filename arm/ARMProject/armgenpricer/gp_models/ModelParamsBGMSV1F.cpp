/*!
 *
 * Copyright (c) CDC IXIS CM July 2003 Paris
*/

/// gpmodels
#include "gpmodels/ModelParamsBGMSV1F.h"

/// gpbase

/// gpinfra



CC_BEGIN_NAMESPACE( ARM )

ARM_ModelParamsBGMSV1F::ARM_ModelParamsBGMSV1F(const ARM_ModelParamVector& params, const std::vector<double>& stddev, 
											   double recorrel, int nbfactors, double minratio, 
											   bool localrhocalib, bool globalCalib) : 
ARM_ModelParams(params), 
itsRecorrel(recorrel), 
itsMinRatio(minratio), 
itsNbFactors(nbfactors), 
itsLocalRhoCalib(localrhocalib), 
itsStdDev(stddev),
itsGlobalCalib(globalCalib), 
itsCalibrated(false), 
itsVVar(0)
{
	if(itsStdDev.size() == 0)
	{
		itsStdDev.resize(11);
		itsStdDev[0] = -2.;
		for(int k = 1; k < 11; k++) itsStdDev[k] = itsStdDev[k-1]+0.4;
	}
}

ARM_ModelParamsBGMSV1F::ARM_ModelParamsBGMSV1F(const ARM_ModelParamsBGMSV1F& rhs) : 
ARM_ModelParams(rhs), 
itsRecorrel(rhs.itsRecorrel), 
itsMinRatio(rhs.itsMinRatio), 
itsNbFactors(rhs.itsNbFactors), 
itsLocalRhoCalib(rhs.itsLocalRhoCalib), 
itsStdDev(rhs.itsStdDev),
itsGlobalCalib(rhs.itsGlobalCalib), 
itsCalibrated(rhs.itsCalibrated), 
itsVVar(rhs.itsVVar)
{
}

ARM_ModelParamsBGMSV1F::~ARM_ModelParamsBGMSV1F()
{
}

void ARM_ModelParamsBGMSV1F::SetParams(const std::vector<double>& resetTimes, double v0, double kappa, double theta,
									   const std::vector<double>& shifts, const std::vector<double>& levels, 
									   const std::vector<double>& vvols, const std::vector<double>& rhos)
{
	( (ARM_CurveModelParam&) GetModelParam(ARM_ModelParamType::VolMeanReversion) ).GetCurve()->SetOrdinates(std::vector<double>(1,kappa));
	( (ARM_CurveModelParam&) GetModelParam(ARM_ModelParamType::InitialVol) ).GetCurve()->SetOrdinates(std::vector<double>(1,v0));
	( (ARM_CurveModelParam&) GetModelParam(ARM_ModelParamType::LongTermVol) ).GetCurve()->SetOrdinates(std::vector<double>(1,theta));

	if(shifts.size() != resetTimes.size() && shifts.size() != 1)
	{
		ARM_THROW(ERR_INVALID_ARGUMENT, "ARM_ModelParamsBGMSV1F::SetParam shift and times should have same size");
	}

	if(DoesModelParamExist(ARM_ModelParamType::Shift))
	{
		( (ARM_CurveModelParam&) GetModelParam(ARM_ModelParamType::Shift) ).GetCurve()->SetAbscisses(resetTimes);

		if(shifts.size() == resetTimes.size())
		{
			( (ARM_CurveModelParam&) GetModelParam(ARM_ModelParamType::Shift) ).GetCurve()->SetOrdinates(shifts);
		}
		else if(shifts.size() == 1)
		{
			( (ARM_CurveModelParam&) GetModelParam(ARM_ModelParamType::Shift) ).GetCurve()->SetOrdinates(std::vector<double>(resetTimes.size(), shifts[0]));
		}
	}

	if(levels.size() != resetTimes.size() && levels.size() != 1)
	{
		ARM_THROW(ERR_INVALID_ARGUMENT, "ARM_ModelParamsBGMSV1F::SetParam levels and times should have same size");
	}

	( (ARM_CurveModelParam&) GetModelParam(ARM_ModelParamType::Alpha) ).GetCurve()->SetAbscisses(resetTimes);

	if(levels.size() == resetTimes.size())
	{
		( (ARM_CurveModelParam&) GetModelParam(ARM_ModelParamType::Alpha) ).GetCurve()->SetOrdinates(levels);
	}
	else if(levels.size() == 1)
	{
		( (ARM_CurveModelParam&) GetModelParam(ARM_ModelParamType::Alpha) ).GetCurve()->SetOrdinates(std::vector<double>(resetTimes.size(),levels[0]));
	}

	if(rhos.size() != resetTimes.size() && rhos.size() != 1)
	{
		ARM_THROW(ERR_INVALID_ARGUMENT, "ARM_ModelParamsBGMSV1F::SetParam rhos and times should have same size");
	}

	( (ARM_CurveModelParam&) GetModelParam(ARM_ModelParamType::QParameter) ).GetCurve()->SetAbscisses(resetTimes);

	if(rhos.size() == resetTimes.size())
	{
		( (ARM_CurveModelParam&) GetModelParam(ARM_ModelParamType::QParameter) ).GetCurve()->SetOrdinates(rhos);
	}
	else if(rhos.size() == 1)
	{
		( (ARM_CurveModelParam&) GetModelParam(ARM_ModelParamType::QParameter) ).GetCurve()->SetOrdinates(std::vector<double>(resetTimes.size(),rhos[0]));
	}


	if(vvols.size() == 1)
	{
		( (ARM_CurveModelParam&) GetModelParam(ARM_ModelParamType::VolOfVol) ).GetCurve()->SetAbscisses(std::vector<double>(1,0.));
		( (ARM_CurveModelParam&) GetModelParam(ARM_ModelParamType::VolOfVol) ).GetCurve()->SetOrdinates(std::vector<double>(1,0.));
	}
	else if(vvols.size() == resetTimes.size())
	{
		( (ARM_CurveModelParam&) GetModelParam(ARM_ModelParamType::VolOfVol) ).GetCurve()->SetAbscisses(resetTimes);

		std::vector<double> bvvols(vvols);

		itsVVar.resize(vvols.size());

		itsVVar[0] = vvols[0]*vvols[0]*resetTimes[0]/365.;
		for(int i = 1; i < vvols.size(); i++)
		{
			itsVVar[i] = vvols[i]*vvols[i]*resetTimes[i]/365.;
			bvvols[i] = itsVVar[i] < itsVVar[i-1] ? 0. : sqrt((itsVVar[i] - itsVVar[i-1])*365./(resetTimes[i]-resetTimes[i-1]));
		}

		( (ARM_CurveModelParam&) GetModelParam(ARM_ModelParamType::VolOfVol) ).GetCurve()->SetOrdinates(bvvols);
	}
	else
	{
		ARM_THROW(ERR_INVALID_ARGUMENT, "ARM_ModelParamsBGMSV1F::SetParam vol of vol and times should have same size");
	}

}

double ARM_ModelParamsBGMSV1F::GetShift(int i) const
{
	return ( (ARM_CurveModelParam&) GetModelParam(ARM_ModelParamType::Shift) ).GetCurve()->GetOrdinates()[i];
}

double ARM_ModelParamsBGMSV1F::GetLevel(int i) const
{
	return ( (ARM_CurveModelParam&) GetModelParam(ARM_ModelParamType::Alpha) ).GetCurve()->GetOrdinates()[i];
}

double ARM_ModelParamsBGMSV1F::GetRho(int i) const
{
	return ( (ARM_CurveModelParam&) GetModelParam(ARM_ModelParamType::QParameter) ).GetCurve()->GetOrdinates()[i];
}

double ARM_ModelParamsBGMSV1F::GetVVol(double time) const
{
	return ( (ARM_CurveModelParam&) GetModelParam(ARM_ModelParamType::VolOfVol) ).GetCurve()->Interpolate(time);
}

double ARM_ModelParamsBGMSV1F::GetVVol(int i) const
{
	if(itsGlobalCalib || ( (ARM_CurveModelParam&) GetModelParam(ARM_ModelParamType::VolOfVol) ).GetCurve()->GetOrdinates().size() == 1)
	{
		return ( (ARM_CurveModelParam&) GetModelParam(ARM_ModelParamType::VolOfVol) ).GetCurve()->GetOrdinates()[0];
	}
	else
	{
		double Ti = ( (ARM_CurveModelParam&) GetModelParam(ARM_ModelParamType::VolOfVol) ).GetCurve()->GetAbscisses()[i]/365.;

		return sqrt(itsVVar[i]/Ti);
	}
}

ARM_Curve * ARM_ModelParamsBGMSV1F::GetVVolt(double evalTime) const
{
	if(evalTime < K_DOUBLE_TOL)
		return ((ARM_CurveModelParam&) GetModelParam(ARM_ModelParamType::VolOfVol)).GetCurve();
	else
	{
		ARM_Curve * initCurve = ((ARM_CurveModelParam&) GetModelParam(ARM_ModelParamType::VolOfVol)).GetCurve();

		int i, size = initCurve->GetAbscisses().size();

		if(size == 1) return initCurve;

		i = 0; 
		while(initCurve->GetAbscisses()[i] < evalTime)
		{
			if(i == size-1) break;
			i++;
		}

		std::vector<double> times(size-i), values(size-i);
		for(int k = 0; k <times.size(); k++)
		{
			times[k] = initCurve->GetAbscisses()[i+k];
			values[k] = initCurve->GetOrdinates()[i+k];
		}

		return new ARM_Curve(times,values);
	}
}

double ARM_ModelParamsBGMSV1F::GetKappa() const
{
	double kappa = ( (ARM_CurveModelParam&) GetModelParam(ARM_ModelParamType::VolMeanReversion) ).GetCurve()->GetOrdinates()[0];
	return kappa < 1e-5 ? 1e-5 : kappa;
}

double ARM_ModelParamsBGMSV1F::GetV0() const
{
	return ( (ARM_CurveModelParam&) GetModelParam(ARM_ModelParamType::InitialVol) ).GetCurve()->GetOrdinates()[0];
}

double ARM_ModelParamsBGMSV1F::GetTheta() const
{
	return ( (ARM_CurveModelParam&) GetModelParam(ARM_ModelParamType::LongTermVol) ).GetCurve()->GetOrdinates()[0];
}

double ARM_ModelParamsBGMSV1F::RateRateCorrel(double t, double ithResetTime, int ithRate, double jthResetTime, int jthRate) const
{
	if(ithRate == jthRate) return 1.;

	double beta		= ((ARM_CurveModelParam&)GetModelParam(ARM_ModelParamType::BetaCorrelation)).GetCurve()->Interpolate(t);
	
	double dt		= (ithResetTime - jthResetTime) / K_YEAR_LEN;
	double delta	= (ithResetTime + jthResetTime - 2. * t) / K_YEAR_LEN;

	return fabs(itsRecorrel) < K_NEW_DOUBLE_TOL || delta < K_NEW_DOUBLE_TOL ? exp(- beta * fabs(dt)) : exp( - beta * fabs(dt) / pow(delta, itsRecorrel));
}

double ARM_ModelParamsBGMSV1F::GetMinRatio() const
{
	return itsMinRatio;
}

void ARM_ModelParamsBGMSV1F::checkFactorCount(int nbResetDates)
{
	itsNbFactors = itsNbFactors <= 0 || itsNbFactors > nbResetDates + 1 ? nbResetDates + 1 : itsNbFactors;
}

void ARM_ModelParamsBGMSV1F::SetFactorCount(int factorsNb)
{
	itsNbFactors = factorsNb;
}

string ARM_ModelParamsBGMSV1F::toString(const string& indent,const string& nextIndent) const
{
	CC_Ostringstream os;

    os << "ARM_ModelParamsBGMSV1F\n";
    os << "----------------------\n\n";
    os << ARM_ModelParams::toString();
	os << "\n\n";
	os << indent << "Effective number of factors : \n";
	os << itsNbFactors << "\n";

	return os.str();
}

void ARM_ModelParamsBGMSV1F::Calibrate(const std::vector<double>& resetTimes, const std::vector<double>& fwdRates,
									   const ARM_VanillaSecDensityPtrVector& CalibSecDensities)
{
	if(itsCalibrated) return;

	ARM_VectorVector calibStrikes, calibVols;
	std::vector<double> atmvols, calibATM, calibTimes, calibFwd, calibWeights;

	GetCalibrationData(resetTimes, fwdRates, CalibSecDensities, calibTimes, calibFwd, calibStrikes, calibVols, calibATM, calibWeights, atmvols);

	if(itsGlobalCalib)
	{
		double v0 = GetV0();
		double kappa = GetKappa();
		double rho = GetRho(0);
		double shift = GetShift(0);
		double theta = GetTheta();

		bool calibtheta = theta < K_DOUBLE_TOL ? true : false;
		bool calibrho = fabs(rho) > 0.999 ? true : false;
		bool calibshift = fabs(shift) < K_DOUBLE_TOL ? true : false;
		
		if(calibrho == false) itsLocalRhoCalib = false;

		ARM_SmileCalibration_Params_Heston params(v0, kappa, theta, rho, 0.1, shift, 1., 
												  false, // calibkappa
												  calibtheta,
												  calibrho, 
												  true,
												  calibshift,
												  itsLocalRhoCalib);

		ARM_SmileCalibration_Heston calibtool;
		calibtool.Init(calibTimes, calibFwd, calibVols, calibStrikes, true, calibFwd, calibATM, std::vector<double>(calibTimes.size(), 1.), 
			&params);
		calibtool.Calibrate();

		// interpolation avec recalibration de la vol monnaie
		finishGlobalCalibration(resetTimes, fwdRates, atmvols, calibWeights, calibTimes, params);	
	}
	else
	{
		double v0 = GetV0();
		double kappa = GetKappa();
		double rho = GetRho(0);
		double shift = GetShift(0);
		double theta = GetV0();

		bool calibrho = fabs(rho) > 0.999 ? true : false;
		bool calibshift = fabs(shift) < K_DOUBLE_TOL ? true : false;

		int k, size = calibTimes.size();

		std::vector<double> shifts(size), levels(size), rhos(size), vvols(size);

		for(k = 0; k < size; k++)
		{
			ARM_SmileCalibration_Params_Heston params(v0, kappa, theta, rho, 0.1, shift, 1.,
													  false,
													  false,
													  calibrho,
													  true,
													  calibshift);

			ARM_SmileCalibration_Heston calibtool;
			calibtool.Init(calibTimes[k], calibFwd[k], *(calibVols[k]), *(calibStrikes[k]), true, calibFwd[k], 
				calibATM[k], &params);
			calibtool.Calibrate();

			shifts[k]	= params.shift();
			levels[k]	= params.level();
			rhos[k]		= params.rho();
			vvols[k]	= params.nu();
		}

		finishLocalCalibration(resetTimes, fwdRates, atmvols, calibWeights, calibTimes, shifts, levels, rhos, vvols);
	}
}

void ARM_ModelParamsBGMSV1F::GetCalibrationData(const std::vector<double>& resetTimes, const std::vector<double>& fwdRates,
												const ARM_VanillaSecDensityPtrVector& CalibSecDensities,
												std::vector<double>& calibTimes, std::vector<double>& calibFwd,
												ARM_VectorVector& calibStrikes, ARM_VectorVector& calibVols,
												std::vector<double>& calibATM, std::vector<double>& calibWeights,
												std::vector<double>& atmvols)
{
	int i, k, size = resetTimes.size(), realCalibSize = 0;

	calibWeights.resize(size,1.);

	double dt, previousReset = 0.;

	for(i = 0; i < size; i++)
	{
		dt = i == 0 ? resetTimes[i] : resetTimes[i] - previousReset;

		if(fabs(CalibSecDensities[i]->getWeight()) < K_DOUBLE_TOL)
		{
			calibWeights[i] = 0.;
		}
		else
		{
			if(dt < 360)
			{
				calibWeights[i] = 0.;
			}
			else
			{
				calibWeights[i] = 1.;
				realCalibSize ++;
				previousReset = resetTimes[i];
			}
		}
	}

	if(realCalibSize == 0)
	{
		realCalibSize = 1;
		calibWeights[size-1] = 1.;
	}

	calibTimes.resize(realCalibSize);
	calibFwd.resize(realCalibSize);
	calibStrikes.resize(realCalibSize);
	calibVols.resize(realCalibSize);
	calibATM.resize(realCalibSize);
	atmvols.resize(size);

	double price;

	for(i = 0, k = 0; i < size; i++)
	{
		price = CalibSecDensities[i]->getDensityFunctor()->Call_Option(fwdRates[i], fwdRates[i], resetTimes[i] / 365.);
		ARM_ImpliedVolBS inverse(fwdRates[i], fwdRates[i], resetTimes[i]/365., 1);
		atmvols[i] = inverse.vol(price);
		
		if(fabs(calibWeights[i]) < K_DOUBLE_TOL) continue;

		calibStrikes[k] = new std::vector<double>(itsStdDev.size());
		calibVols[k] = new std::vector<double>(itsStdDev.size());
		
		calibTimes[k] = resetTimes[i] / 365.;
		calibFwd[k] = fwdRates[i];
		calibATM[k] = atmvols[i];

		for(int n = 0; n < itsStdDev.size(); n++)
		{
			(*calibStrikes[k])[n] = calibFwd[k] * exp(atmvols[i]*itsStdDev[n]*sqrt(calibTimes[k]));
			price = CalibSecDensities[i]->getDensityFunctor()->Call_Option((*calibStrikes[k])[n], fwdRates[i], resetTimes[i] / 365.);
			(*calibVols[k])[n] = inverse.vol(price, (*calibStrikes[k])[n]);
		}

		k++;
	}

}

void ARM_ModelParamsBGMSV1F::finishGlobalCalibration(const std::vector<double>& resetTimes, 
													 const std::vector<double>& fwdRates,
													 const std::vector<double>& atmvols, std::vector<double>& calibWeights, 
													 const std::vector<double>& calibTimes, 
													 const ARM_SmileCalibration_Params_Heston& params)
{
	(*( (ARM_CurveModelParam&) GetModelParam(ARM_ModelParamType::LongTermVol) ).GetCurve())[0] = params.theta();	
	(*( (ARM_CurveModelParam&) GetModelParam(ARM_ModelParamType::VolOfVol) ).GetCurve())[0] = params.nu();	

	int i, size = resetTimes.size();
	std::vector<double> shifts(size), rhos(size), levels(size);

	if(itsLocalRhoCalib == false)
	{
		for(i = 0; i < size; i++) rhos[i] = params.rho();
	}
	else
	{
		ARM_Curve curve(calibTimes, params.rhos(), new ARM_LinInterpCstExtrapolDble);

		for(i = 0; i < size; i++)
		{
			rhos[i] = curve(resetTimes[i]/365.);
		}
	}

	if(params.calibshift() == false)
	{
		for(i = 0; i < size; i++) shifts[i] = params.shift();
	}
	else
	{
		ARM_Curve curve(calibTimes, params.shifts(), new ARM_LinInterpCstExtrapolDble);

		for(i = 0; i < size; i++)
		{
			shifts[i] = curve(resetTimes[i]/365.);
		}
	}

	ARM_Curve curve(calibTimes, params.levels(), new ARM_LinInterpCstExtrapolDble);

	for(i = 0; i < resetTimes.size(); i++)
	{
		levels[i] = curve(resetTimes[i]/365.);

		if(fabs(fabs(rhos[i])-0.999) < K_DOUBLE_TOL && i > 0)
		{
			rhos[i] = rhos[i-1];
			calibWeights[i] = 0.;
		}

		if(fabs(calibWeights[i]) < K_DOUBLE_TOL)
		{
			ARM_SmileCalibration_Params_Heston tmpParams(GetV0(), GetKappa(), GetTheta(), rhos[i], GetVVol(0),
														shifts[i], 1., 
														false, false, false, false, false);
			ARM_SmileCalibration_Heston calibtool;
			calibtool.Init(resetTimes[i]/365., fwdRates[i], std::vector<double>(5,1.), std::vector<double>(5,1.),
				true, fwdRates[i], atmvols[i], &tmpParams);
			
			calibtool.Calibrate();

			levels[i] = tmpParams.levels()[0];
		}
	}

	( (ARM_CurveModelParam&) GetModelParam(ARM_ModelParamType::QParameter) ).GetCurve()->SetAbscisses(resetTimes);
	( (ARM_CurveModelParam&) GetModelParam(ARM_ModelParamType::Shift) ).GetCurve()->SetAbscisses(resetTimes);
	( (ARM_CurveModelParam&) GetModelParam(ARM_ModelParamType::Alpha) ).GetCurve()->SetAbscisses(resetTimes);

	( (ARM_CurveModelParam&) GetModelParam(ARM_ModelParamType::QParameter) ).GetCurve()->SetOrdinates(rhos);
	( (ARM_CurveModelParam&) GetModelParam(ARM_ModelParamType::Shift) ).GetCurve()->SetOrdinates(shifts);
	( (ARM_CurveModelParam&) GetModelParam(ARM_ModelParamType::Alpha) ).GetCurve()->SetOrdinates(levels);
}

void ARM_ModelParamsBGMSV1F::finishLocalCalibration(const std::vector<double>& resetTimes, 
													const std::vector<double>& fwdRates, 
													const std::vector<double>& atmvols, std::vector<double>& calibWeights,
													const std::vector<double>& calibTimes, 
													const std::vector<double>& calibShifts, 
													const std::vector<double>& calibLevels, 
													const std::vector<double>& calibRhos, 
													const std::vector<double>& calibVvols)
{
	(*( (ARM_CurveModelParam&) GetModelParam(ARM_ModelParamType::LongTermVol) ).GetCurve())[0] = GetV0();	

	int i, size = resetTimes.size();
	
	std::vector<double> shifts(size), levels(size), rhos(size), vvols(size);
	itsVVar.resize(size);

	ARM_Curve curve(calibTimes, calibShifts, new ARM_LinInterpCstExtrapolDble);

	for(i = 0; i < size; i++) shifts[i] = curve(resetTimes[i]/365.);

	curve.SetOrdinates(calibRhos);

	for(i = 0; i < size; i++) rhos[i] = curve(resetTimes[i]/365.);

	curve.SetOrdinates(calibVvols);

	for(i = 0; i < size; i++) vvols[i] = curve(resetTimes[i]/365.);


	bool rescale = false;

	if(rescale == false)
	{
		curve.SetOrdinates(calibLevels);

		for(i = 0; i < size; i++)
		{
			levels[i] = curve(resetTimes[i]/365.);

			if(fabs(calibWeights[i]) < K_DOUBLE_TOL)
			{
				ARM_SmileCalibration_Params_Heston tmpParams(GetV0(), GetKappa(), GetTheta(), rhos[i], vvols[i],
															shifts[i], 1., 
															false, false, false, false, false);
				ARM_SmileCalibration_Heston calibtool;
				calibtool.Init(resetTimes[i]/365., fwdRates[i], std::vector<double>(5,1.), std::vector<double>(5,1.),
					true, fwdRates[i], atmvols[i], &tmpParams);
				
				calibtool.Calibrate();

				levels[i] = tmpParams.levels()[0];
			}
		}
	}

	// bootstrap de la vol de vol
	itsVVar[0] = vvols[0]*vvols[0]*resetTimes[0]/365.;

	for(i = 1; i < size; i++)
	{
		itsVVar[i] = vvols[i]*vvols[i]*resetTimes[i]/365.;

		vvols[i] = itsVVar[i] < itsVVar[i-1] ? 0. : sqrt((itsVVar[i]-itsVVar[i-1])/((resetTimes[i]-resetTimes[i-1])/365.));
	}


	if(rescale == true)
	{
		curve.SetAbscisses(resetTimes);
		curve.SetOrdinates(vvols);

		for(i = 0; i < size; i++)
		{
			ARM_HestonOptionPricerVVolt pricer(resetTimes[i]/365., fwdRates[i], fwdRates[i], -1.,
											   GetV0(), GetKappa(), GetTheta(), rhos[i], curve, shifts[i], 1.);
			
			ARM_ImpliedVolBS inverse(fwdRates[i], fwdRates[i], resetTimes[i]/365., -1.);

			CATMLocalConstr func(&pricer, &inverse);

			double best;
			levels[i] = brentSolve(func, atmvols[i], 0.0001, 5., 1e-8, 50, 0, &best);
		}
	}

	( (ARM_CurveModelParam&) GetModelParam(ARM_ModelParamType::QParameter) ).GetCurve()->SetAbscisses(resetTimes);
	( (ARM_CurveModelParam&) GetModelParam(ARM_ModelParamType::Shift) ).GetCurve()->SetAbscisses(resetTimes);
	( (ARM_CurveModelParam&) GetModelParam(ARM_ModelParamType::Alpha) ).GetCurve()->SetAbscisses(resetTimes);
	( (ARM_CurveModelParam&) GetModelParam(ARM_ModelParamType::VolOfVol) ).GetCurve()->SetAbscisses(resetTimes);

	( (ARM_CurveModelParam&) GetModelParam(ARM_ModelParamType::QParameter) ).GetCurve()->SetOrdinates(rhos);
	( (ARM_CurveModelParam&) GetModelParam(ARM_ModelParamType::Shift) ).GetCurve()->SetOrdinates(shifts);
	( (ARM_CurveModelParam&) GetModelParam(ARM_ModelParamType::Alpha) ).GetCurve()->SetOrdinates(levels);
	( (ARM_CurveModelParam&) GetModelParam(ARM_ModelParamType::VolOfVol) ).GetCurve()->SetOrdinates(vvols);
}


CC_END_NAMESPACE()