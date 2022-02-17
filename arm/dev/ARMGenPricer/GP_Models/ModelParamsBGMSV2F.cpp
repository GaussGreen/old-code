#include "gpmodels/ModelParamsBGMSV2F.h"
#include "gpclosedforms/smile_calibration.h"

CC_BEGIN_NAMESPACE( ARM )

ARM_ModelParamsBGMSV2F::ARM_ModelParamsBGMSV2F(const ARM_ModelParamVector& params, 
											   double v01, double kappa1, double rho1,
											   double v02, double kappa2, double rho2, double shift, 
											   const ARM_GP_Vector& stddev,
											   double recorrel, int nbfactor, double minratio, 
											   bool localrho1calib, bool localrho2calib) : 
ARM_ModelParams(params),
itsRecorrel(recorrel),
itsNbFactors(nbfactor),
itsMinRatio(minratio),
itsLocalRho1Calib(localrho1calib),
itsLocalRho2Calib(localrho2calib),
itsStdDev(stddev),
itsCalibrated(false),
itsV01(v01), itsV02(v02), itsKappa1(kappa1), itsKappa2(kappa2), itsTheta1(v01), itsTheta2(v02), 
itsVVol1(0.), itsVVol2(0.),
itsRho1(1,rho1), itsRho2(1,rho2), itsShift(1,shift), itsSigma(0)
{
	if(itsStdDev.size() == 0)
	{
		itsStdDev.resize(11);
		itsStdDev[0] = -2.;
		for(int k = 1; k < 11; k++) itsStdDev[k] = itsStdDev[k-1]+0.4;
	}

	itsCalibRho1 = fabs(rho1) > 0.999 ? true : false;
	itsCalibRho2 = fabs(rho2) > 0.999 ? true : false;
	itsCalibShift = fabs(shift) < K_DOUBLE_TOL ? true : false;
}

ARM_ModelParamsBGMSV2F::ARM_ModelParamsBGMSV2F(const ARM_ModelParamsBGMSV2F& rhs) : 
ARM_ModelParams(rhs),
itsRecorrel(rhs.itsRecorrel),
itsNbFactors(rhs.itsNbFactors),
itsMinRatio(rhs.itsMinRatio),
itsLocalRho1Calib(rhs.itsLocalRho1Calib),
itsLocalRho2Calib(rhs.itsLocalRho2Calib),
itsCalibrated(rhs.itsCalibrated),
itsV01(rhs.itsV01),
itsKappa1(rhs.itsKappa1),
itsTheta1(rhs.itsTheta1),
itsVVol1(rhs.itsVVol1),
itsRho1(rhs.itsRho1),
itsV02(rhs.itsV02),
itsKappa2(rhs.itsKappa2),
itsTheta2(rhs.itsTheta2),
itsVVol2(rhs.itsVVol2),
itsRho2(rhs.itsRho2),
itsShift(rhs.itsShift),
itsSigma(rhs.itsSigma),
itsStdDev(rhs.itsStdDev),
itsCalibRho1(rhs.itsCalibRho1),
itsCalibRho2(rhs.itsCalibRho2),
itsCalibShift(rhs.itsCalibShift)
{
}

ARM_ModelParamsBGMSV2F::~ARM_ModelParamsBGMSV2F()
{
}

string ARM_ModelParamsBGMSV2F::toString(const string& indent, const string& nextIndent) const
{
	CC_Ostringstream os;

    os << "ARM_ModelParamsBGMSV2F\n";
    os << "----------------------\n\n";
    os << ARM_ModelParams::toString();
	os << "\n\n";
	os << indent << "v01 : \t" << itsV01 << " \n";
	os << indent << "Kappa1 : \t" << itsKappa1 << " \n";
	os << indent << "Theta1 : \t" << itsTheta1 << " \n";
	os << indent << "VVol1 : \t" << itsVVol1 << " \n";
	os << indent << "v02 : \t" << itsV02 << " \n";
	os << indent << "Kappa2 : \t" << itsKappa2 << " \n";
	os << indent << "Theta2 : \t" << itsTheta2 << " \n";
	os << indent << "VVol2 : \t" << itsVVol2 << " \n";

	int size = itsSigma.size();

	os << indent << "Rho1 \t" << "Rho2 \t" << "Shift \t" << "Level \n";

	for(int i = 0; i < size; i++)
	{
		os << indent << CC_NS(std,setprecision)(5) << itsRho1[i] << "\t" 
					 << CC_NS(std,setprecision)(5) << itsRho2[i] << "\t" 
					 << CC_NS(std,setprecision)(5) << itsShift[i] << "\t" 
					 << CC_NS(std,setprecision)(5) << itsSigma[i] << "\n";
	}

	return os.str();

}


void ARM_ModelParamsBGMSV2F::checkFactorCount(int nbResetDates)
{
	itsNbFactors = itsNbFactors <= 0 || itsNbFactors > nbResetDates + 1 ? nbResetDates + 1 : itsNbFactors;
}

void ARM_ModelParamsBGMSV2F::SetFactorCount(int factorsNb)
{
	itsNbFactors = factorsNb;
}

double ARM_ModelParamsBGMSV2F::RateRateCorrel(double t, double ithResetTime, int ithRate, double jthResetTime, int jthRate) const
{
	if(ithRate == jthRate) return 1.;

	double beta		= ((ARM_CurveModelParam&)GetModelParam(ARM_ModelParamType::BetaCorrelation)).GetCurve()->Interpolate(t);
	
	double dt		= (ithResetTime - jthResetTime) / K_YEAR_LEN;
	double delta	= (ithResetTime + jthResetTime - 2. * t) / K_YEAR_LEN;

	return fabs(itsRecorrel) < K_NEW_DOUBLE_TOL || delta < K_NEW_DOUBLE_TOL ? exp(- beta * fabs(dt)) : exp( - beta * fabs(dt) / pow(delta, itsRecorrel));
}

void ARM_ModelParamsBGMSV2F::Calibrate(const ARM_GP_Vector& resetTimes, const ARM_GP_Vector& fwdRates, 
									   const ARM_VanillaSecDensityPtrVector& CalibSecDensities)
{
	if(itsCalibrated) return;

	ARM_VectorVector calibStrikes, calibVols;
	ARM_GP_Vector atmvols, calibATM, calibTimes, calibFwd, calibWeights;

	GetCalibrationData(resetTimes, fwdRates, CalibSecDensities, calibTimes, calibFwd, calibStrikes, calibVols, calibATM, calibWeights, atmvols);

	// calibration du smile
	ARM_SmileCalibration_Params_Heston2b params(itsV01, itsKappa1, itsTheta1, itsRho1[0], itsVVol1, 
												itsV02, itsKappa2, itsTheta2, itsRho2[0], itsVVol2, 
												itsShift[0], 1., 
												false, false, itsCalibRho1, true,
												false, false, itsCalibRho2, true,
												itsCalibShift, itsLocalRho1Calib, itsLocalRho2Calib);

	ARM_SmileCalibration_Heston2b calibtool;
	calibtool.Init(calibTimes, calibFwd, calibVols, calibStrikes, true, calibFwd, calibATM, ARM_GP_Vector(calibTimes.size(), 1.), 
		&params);
	calibtool.Calibrate();

	// interpolation avec recalibration de la vol monnaie
	finishCalibration(resetTimes, fwdRates, atmvols, calibWeights, calibTimes, params);	

	itsCalibrated = true;
}

void ARM_ModelParamsBGMSV2F::GetCalibrationData(const ARM_GP_Vector& resetTimes, const ARM_GP_Vector& fwdRates,
												const ARM_VanillaSecDensityPtrVector& CalibSecDensities,
												ARM_GP_Vector& calibTimes, ARM_GP_Vector& calibFwd,
												ARM_VectorVector& calibStrikes, ARM_VectorVector& calibVols,
												ARM_GP_Vector& calibATM, ARM_GP_Vector& calibWeights,
												ARM_GP_Vector& atmvols)
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

		calibStrikes[k] = new ARM_GP_Vector(itsStdDev.size());
		calibVols[k] = new ARM_GP_Vector(itsStdDev.size());
		
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

void ARM_ModelParamsBGMSV2F::finishCalibration(const ARM_GP_Vector& resetTimes, const ARM_GP_Vector& fwdRates,
											   const ARM_GP_Vector& atmvols, ARM_GP_Vector& calibWeights, 
											   const ARM_GP_Vector& calibTimes, const ARM_SmileCalibration_Params_Heston2b& params)
{
	itsSigma.resize(resetTimes.size());

	itsVVol1 = params.nu1();
	itsVVol2 = params.nu2();

	if(itsLocalRho1Calib == false)
	{
		itsRho1.resize(resetTimes.size(), params.rho1());
	}
	else
	{
		ARM_Curve curve(calibTimes, params.rhos1(), new ARM_LinInterpCstExtrapolDble);

		itsRho1.resize(resetTimes.size());
		for(int i = 0; i < resetTimes.size(); i++)
		{
			itsRho1[i] = curve(resetTimes[i]/365.);
		}
	}

	if(itsLocalRho2Calib == false)
	{
		itsRho2.resize(resetTimes.size(), params.rho2());
	}
	else
	{
		ARM_Curve curve(calibTimes, params.rhos2(), new ARM_LinInterpCstExtrapolDble);

		itsRho2.resize(resetTimes.size());
		for(int i = 0; i < resetTimes.size(); i++)
		{
			itsRho2[i] = curve(resetTimes[i]/365.);
		}
	}

	if(itsCalibShift == false)
	{
		itsShift.resize(resetTimes.size(), params.shifts()[0]);
	}
	else
	{
		ARM_Curve curve(calibTimes, params.shifts(), new ARM_LinInterpCstExtrapolDble);

		itsShift.resize(resetTimes.size());
		for(int i = 0; i < resetTimes.size(); i++)
		{
			itsShift[i] = curve(resetTimes[i]/365.);
		}
	}

	ARM_Curve curve(calibTimes, params.levels(), new ARM_LinInterpCstExtrapolDble);

	itsSigma.resize(resetTimes.size());
	for(int i = 0; i < resetTimes.size(); i++)
	{
		itsSigma[i] = curve(resetTimes[i]/365.);

		if(fabs(fabs(itsRho1[i])-0.999) < K_DOUBLE_TOL && i > 0)
		{
			itsRho1[i] = itsRho1[i-1];
			calibWeights[i] = 0.;
		}

		if(fabs(fabs(itsRho2[i])-0.999) < K_DOUBLE_TOL && i > 0)
		{
			itsRho2[i] = itsRho2[i-1];
			calibWeights[i] = 0.;
		}

		if(fabs(calibWeights[i]) < K_DOUBLE_TOL)
		{
			ARM_SmileCalibration_Params_Heston2b tmpParams(itsV01, itsKappa1, itsTheta1, itsRho1[i], itsVVol1, 
														itsV02, itsKappa2, itsTheta2, itsRho2[i], itsVVol2, 
														itsShift[i], 1., 
														false, false, false, false,
														false, false, false, false,
														false);
			ARM_SmileCalibration_Heston2b calibtool;
			calibtool.Init(resetTimes[i]/365., fwdRates[i], ARM_GP_Vector(5,1.), ARM_GP_Vector(5,1.),
				true, fwdRates[i], atmvols[i], &tmpParams);
			
			calibtool.Calibrate();

			itsSigma[i] = tmpParams.levels()[0];
		}
	}


}



CC_END_NAMESPACE()