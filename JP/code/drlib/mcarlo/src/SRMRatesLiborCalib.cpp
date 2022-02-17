//----------------------------------------------------------------------------
//
//   Group       : QR
//
//   Filename    : SRMRatesLiborCalib.cpp
//
//   Description : forward libor model path generation
//
//
//----------------------------------------------------------------------------

#include "edginc/config.hpp"
#include "edginc/QMCRatesDiffuse.hpp"
#include "edginc/SRMRatesLiborDiffuse.hpp"
#include "edginc/SRMRatesUtil.hpp"
#include "edginc/SRMFXDiffuse.hpp"
#include "edginc/SRMUtil.hpp"
#include "edginc/QMCHelperCachingTimeLogic.hpp"
#include "edginc/SwapMaturityVolRequest.hpp"

#include <cassert>
#include <algorithm>
#include "stdlib.h"
#include <limits>

DRLIB_BEGIN_NAMESPACE

const double TEN_BP = 0.001;

// make sure all elementary types are initialized to some unreal values

void SRMRatesLiborDiffuse::analyticCalib(double calibAccuracy)
{
	static const string method = "SRMRatesLiborDiffuse::analyticCalib";

	const long   
		nbRuns = 2,
		maxNbBackSteps = 10,
		nbSteps = 100 ;

	const double 
		volA     = 0.01,        // vol/temperature
		volB     = 0.01,
		volD     = 0.01,    
		deltaA   = 0.001, 
		deltaB   = 0.001,
		deltaD   = 0.001,
		eps      = TEN_BP, 
		largeContract = 0.25,
		smallContract = 0.5,
		corEps = 0.0025;

	double       startDt  = 0.1;


	long factor = -1,
		i = -1,
		run = -1,
		bestRun = -1,
		bestCorRun = -1,
		backStepCounter = -1;

	double 		
		temp = -1.0,
		bestError = -1.,
		bestCorError = -1.,
		dt = -1.,
		sqrtDt = -1.,
		muA = 0.,
		muB = 0.,
		muD = 0.,
		rhoLong = -1.,
		sqrtOneMinusRhoLong = -1.,
		marketBlackVol = 0.,
		modelBlackVol  = 0.,
		totalError = 0.,
		totalErrorNew = 0.,
		totalErrorDiff = 0.,
		errorScale = 0.,
		errorScaleCor = 0.0,
		tempParamChange = 0.,
		gradVolLength2 = 0.,
		gradCorLength2 = 0.,
		gradVolxgradCor = 0.0,
		lambda = 0.,
		shiftLam = 0.,
		gradVolL = 0.,
		gradCorL = 0.,
		scale = 0.,
		scaleA = 0.,
		scaleB = 0.,
		scaleD = 0.,
		noiseA = 0.,
		noiseB = 0.,
		noiseD = 0.,
		corError = 0.0,
		corErrorNew = 0.0,
		corErrorDiff = 0.0,
		gradProduct = 0.0,
		oneOverScale = 0.;

	vector<double> 
		a(m_nbFactors),
		b(m_nbFactors),
		d(m_nbFactors),
		shiftA(m_nbFactors),
		shiftB(m_nbFactors),
		shiftD(m_nbFactors),
		gradVolA(m_nbFactors),
		gradVolB(m_nbFactors),
		gradVolD(m_nbFactors),
		gradCorA(m_nbFactors),
		gradCorB(m_nbFactors),
		gradCorD(m_nbFactors),
		corErrorOptimal(nbRuns),
		corErrorOptimalPrevious(nbRuns),
		totalErrorOptimal(nbRuns),
		totalErrorOptimalPrevious(nbRuns);

	vector< vector<double> > 
		aOptimal(nbRuns),
		bOptimal(nbRuns),
		dOptimal(nbRuns);

	for(i = 0; i < nbRuns; i++)
	{
		aOptimal[i].resize(m_nbFactors);
		bOptimal[i].resize(m_nbFactors);
		dOptimal[i].resize(m_nbFactors);
	}

	// only initialised here
	m_humpParams[0] = Maths::max(m_calibMat, m_calibMatCMS)/12.;
	m_humpParams[1] = 1./30.;

	// set rng and get random numbers for Metropolis
	RngGaussianNRClass rng(123456789);
	vector<double> noiseVector(8 * nbSteps,0.);

	// set maximum total error between model-market swaption vols
	errorScale = eps * eps;
	errorScaleCor = corEps * corEps;

	// determine the constant multiplicative factors
	// automatically calibrates caplet prices
	calibrateCaplets();

	// each run uses a different starting point
	// - the purpose is to make it less likely that the minimisation gets
	// stuck in a local minimum for the error - this is also the reason 
	// for adding noise to the search paths
	for(run = 0; run < nbRuns; run++)
	{
		dt = startDt;
		sqrtDt = sqrt(dt);

		// get new random numbers - not used right now
		rng.getRandVector(&(noiseVector[0]), 8 * nbSteps);	

		// define start point 
		// short/long correlation between -10% and 10%
		rhoLong = 0.2 * run;
		rhoLong = rhoLong / nbRuns - 0.1;

		if (fabs(rhoLong) < 1. - TOO_SMALL)
		{
			sqrtOneMinusRhoLong = sqrt(1. - rhoLong * rhoLong);
		}
		else
		{
			sqrtOneMinusRhoLong = 0.;
		}

		// minimum value ensures some decorrelation
		lambda = m_lambdaInit;

		// initialize parameters 
		for (factor = 0; factor < m_nbFactors; factor++)
		{
			// Note : these settings, and calibrations below,
			//        are for the case of a 2-factor model only

			// choice of a and d vectors ensures that 
			// - a+d and d vectors both have lengths one
			// - dot product between a+d and d vectors equals rho
			a[factor] = (factor == 0) ? 
				m_aInit[factor] - rhoLong : 
			m_aInit[factor] - sqrtOneMinusRhoLong;

			// needs to scale with lambda to ensure that m_b
			// represents largest % deviation (otherwise max 
			// value of b term would scale with lambda)
			b[factor] = (factor == 0) ? lambda * m_bInit[factor] : 0.;

			// this d setting ensures that second factor only, in  
			// case of rhoLong = 0, determines long-end dynamics
			d[factor] = (factor == 0) ? 
				m_dInit[factor] + rhoLong : 
			m_dInit[factor] + sqrtOneMinusRhoLong;

			// these matrices store best results from all runs 
			// - at the end, we then pick the best run
			aOptimal[run][factor] = a[factor];
			bOptimal[run][factor] = b[factor];
			dOptimal[run][factor] = d[factor];
		}

		// initialise
		backStepCounter = 0;
		setParameters(a,b,d,lambda,m_humpParams);
		numericCalibCaplets();
		corError = getCorError();
		corErrorOptimalPrevious[run] = corError;
		totalError = getTotalError();
		totalErrorOptimalPrevious[run] = totalError;

		// this bit of code (i loop) is written for the generic case 
		// (i.e., calibrating all elements of all parameter vectors),
		// but is here used only for the two-factor model, with "mean 
		// reversion corrTwos" m_lambda and m_corrOne as user inputs -
		// hence we calibrate only for a[0], b[0] and b[1]. That's the
		// reason for the commented-out code.

		for (i = 0; i < nbSteps; i++)
		{
			// always comparing with the optimal error
			totalError = totalErrorOptimalPrevious[run];

			//numerically calculate gradient of error function
			gradVolLength2 = 0.;

			for (factor = 0; factor < m_nbFactors; factor++)
			{
				if (fabs(deltaA * a[factor]) > TOO_SMALL)
				{
					shiftA[factor] = deltaA * a[factor];
					scaleA = a[factor];
				}
				else
				{
					shiftA[factor] = 0.01;
					scaleA = shiftA[factor];
				}

				shiftB[factor] = 0.01; // magnitude of lambda * bInit
				scaleB = 1.; // semi-empirical guesstimate of good scaling

				if (fabs(deltaD * d[factor]) > TOO_SMALL)
				{
					shiftD[factor] = deltaD * d[factor];
					scaleD = d[factor];
				}
				else
				{
					shiftD[factor] = 0.01;
					scaleD = shiftD[factor];
				}

				if (factor == 0)
				{
					a[factor] += shiftA[factor];
					setParameters(a,b,d,lambda,m_humpParams);

					totalErrorNew = getTotalError();

					gradVolA[factor] = scaleA * (totalErrorNew - totalError)/shiftA[factor];
					a[factor] -= shiftA[factor];
				}
				else
				{
					gradVolA[factor] = 0.;
				}

				// calculate gradB
				b[factor] += shiftB[factor];
				setParameters(a,b,d,lambda,m_humpParams);				
				totalErrorNew = getTotalError();
				gradVolB[factor] = scaleB * (totalErrorNew - totalError )/shiftB[factor];
				b[factor] -= shiftB[factor];

				// calculate gradD
				/*
				d[factor] += shiftD[factor];
				setParameters(a,b,d,lambda,h);
				totalErrorNew = getTotalError();
				gradVolD[factor] = scaleD * (totalErrorNew - totalError)/shiftD[factor];
				d[factor] -= shiftD[factor];
				*/
				gradVolD[factor] = 0.;

				gradVolLength2 += 
					gradVolA[factor] * gradVolA[factor] +
					gradVolB[factor] * gradVolB[factor] +
					gradVolD[factor] * gradVolD[factor]; 

			} // end factor

			if (gradVolLength2 > TOO_SMALL * errorScale * errorScale)
			{
				scale = totalErrorOptimalPrevious[run] / gradVolLength2;
			}
			else 
			{
				break; // next run
			}

			// update parameters (gradient descent rule)
			muA = 0.;
			muB = 0.;
			muD = 0.;

			for (factor = 0; factor < m_nbFactors; factor++)
			{
				// move down gradient
				muA = scale *  gradVolA[factor] ;
				muB = scale *  gradVolB[factor] ;
				muD = scale *  gradVolD[factor] ;;

				if (factor == 0)
				{
					noiseA = noiseVector[8 * i + factor];
					noiseB = noiseVector[8 * i + 2 + factor];
					noiseD = noiseVector[8 * i + 4 + factor];
				}
				else
				{
					noiseA = rhoLong * noiseVector[8 * i] +
						sqrtOneMinusRhoLong * noiseVector[8 * i + 1];
					noiseB = rhoLong * noiseVector[8 * i + 2] +
						sqrtOneMinusRhoLong * noiseVector[8 * i + 3];
					noiseD = rhoLong * noiseVector[8 * i + 4] +
						sqrtOneMinusRhoLong * noiseVector[8 * i + 5];
				}

				a[factor]      += - muA * dt + muA * volA * sqrtDt * noiseA;
				b[factor]      += - muB * dt + muB * volB * sqrtDt * noiseB;
				d[factor]      += - muD * dt + muD * volD * sqrtDt * noiseD;

			}

			// set new parameter values
			setParameters(a,b,d,lambda,m_humpParams);
			// calculate new total error
			totalError = getTotalError();

			if (fabs(totalError) < errorScale)
			{
				// store parameters which have resulted in best error so far
				for (factor=0; factor<m_nbFactors; factor++)
				{
					aOptimal[run][factor] = a[factor];
					bOptimal[run][factor] = b[factor];
					dOptimal[run][factor] = d[factor];
				}			

				// done - next run!
				break; 
			}

			// check if error decreasing or increasing
			totalErrorDiff = totalError - totalErrorOptimalPrevious[run];

			if ( totalErrorDiff < 0. ) // new best error
			{
				// increase dt 
				dt /= smallContract;
				sqrtDt = sqrt(dt);
				backStepCounter = 0;

				// new best error
				totalErrorOptimalPrevious[run] = totalError; 

				// store parameters which have resulted in best error so far
				for (factor=0; factor<m_nbFactors; factor++)
				{
					aOptimal[run][factor] = a[factor];
					bOptimal[run][factor] = b[factor];
					dOptimal[run][factor] = d[factor];
				}			
			}
			else // worse than best error
			{
				// revert
				for (factor = 0; factor < m_nbFactors; factor++)
				{
					a[factor] = aOptimal[run][factor];
					b[factor] = bOptimal[run][factor];
					d[factor] = dOptimal[run][factor];
				}

				// reduce dt in hope of finding minimum
				dt *= largeContract;
				sqrtDt = sqrt(dt);

				backStepCounter++;

			}  // end backtrack step

			// after being stuck in a plateau for maxCounter steps, we exit 
			if (backStepCounter > maxNbBackSteps)
			{
				break; // next run
			}			

		} // end for i (minimization loop) - vol minimization

		fudgeCalib(m_a, 
			m_b, 
			m_d, 
			m_lambda,
			m_humpParams);

		// minimize correlation error, start from minimum of vol error
		startDt = 0.1;
		dt = startDt;
		sqrtDt = sqrt(dt);

		// get new random numbers 
		rng.getRandVector(&(noiseVector[0]), 8 * nbSteps);

		// initialize parameters : start from vol minimization parameters
		for (factor = 0; factor < m_nbFactors; factor++)
		{
			a[factor] = aOptimal[run][factor];
			b[factor] = bOptimal[run][factor];
			d[factor] = dOptimal[run][factor];
		}

		// initialise
		backStepCounter = 0;
		setParameters(a,b,d,lambda,m_humpParams);

		// "previous" values are for checking that errors below do not 
		// increase, relative to errors above
		totalError = getTotalError();
		totalErrorOptimalPrevious[run] = totalError;

		corError = getCorError();
		corErrorOptimal[run] = corError;

		// needed for proper gradient descent subject to constraint
		// that vol error remains roughly constant
		scaleA = 1.;
		scaleB = 1.;
		scaleD = 1.;

		for (i = 0; i < nbSteps; i++)
		{
			// always comparing with the optimal error
			// used as centre values for calculating gradients only
			totalError = totalErrorOptimal[run];
			corError = corErrorOptimal[run];

			//numerically calculate gradient of error function
			gradVolLength2 = 0.;
			gradCorLength2 = 0.0;
			gradVolxgradCor =  0.0;
			for (factor = 0; factor < m_nbFactors; factor++)
			{
				if (fabs(deltaA * a[factor]) > TOO_SMALL)
				{
					shiftA[factor] = deltaA * a[factor];
				}
				else
				{
					shiftA[factor] = 0.01;
				}

				shiftB[factor] = 0.01; // magnitude of lambda * bInit
				scaleB = 1.; // semi-empirical guesstimate of good scaling

				if (fabs(deltaD * d[factor]) > TOO_SMALL)
				{
					shiftD[factor] = deltaD * d[factor];
				}
				else
				{
					shiftD[factor] = 0.01;
				}

				/*
				a[factor] += shiftA[factor];
				setParameters(a,b,d,lambda,h);

				totalErrorNew = getTotalError();
				corErrorNew = getCorError();

				gradVolA[factor] = scaleA * (totalErrorNew - totalError)/shiftA[factor];
				gradCorA[factor] = scaleA * (corErrorNew - corError)/shiftA[factor];
				a[factor] -= shiftA[factor];
				*/
				gradVolA[factor] = 0.;
				gradCorA[factor] = 0.;

				// calculate gradB
				/*
				b[factor] += shiftB[factor];
				setParameters(a,b,d,lambda,h);

				totalErrorNew = getTotalError();
				corErrorNew = getCorError();

				gradVolB[factor] = scaleB * (totalErrorNew - totalError )/shiftB[factor];
				gradCorB[factor] = scaleB * (corErrorNew - corError)/shiftB[factor];
				b[factor] -= shiftB[factor];
				*/
				gradVolB[factor] = 0.;
				gradCorB[factor] = 0.;

				// calculate gradD
				d[factor] += shiftD[factor];
				setParameters(a,b,d,lambda,m_humpParams);

				totalErrorNew = getTotalError();
				corErrorNew = getCorError();

				gradVolD[factor] = scaleD * (totalErrorNew - totalError)/shiftD[factor];
				gradCorD[factor] = scaleD * (corErrorNew - corError)/shiftD[factor];
				d[factor] -= shiftD[factor];

				/*
				gradVolD[factor] = 0.;
				gradCorD[factor] = 0.;
				*/

				gradVolLength2 += 
					gradVolA[factor] * gradVolA[factor] +
					gradVolB[factor] * gradVolB[factor] +
					gradVolD[factor] * gradVolD[factor];

				gradCorLength2 += 
					gradCorA[factor] * gradCorA[factor] +
					gradCorB[factor] * gradCorB[factor] +
					gradCorD[factor] * gradCorD[factor];

				gradVolxgradCor += 
					gradVolA[factor] * gradCorA[factor] +
					gradVolB[factor] * gradCorB[factor] +
					gradVolD[factor] * gradCorD[factor];

			} // end factor

			// add more checks ....
			if ( gradVolLength2 > TOO_SMALL * errorScale * errorScale  )
			{
				//gradProduct = 0.0;
				gradProduct = gradVolxgradCor;

				gradProduct /= gradVolLength2; 
				oneOverScale = gradCorLength2 
					- gradVolxgradCor * gradVolxgradCor / gradVolLength2;

				if (oneOverScale > TOO_SMALL * errorScale * errorScale ) 
				{
					oneOverScale /= oneOverScale;
				}
			}
			else
			{
				gradProduct = 0.;
				oneOverScale = 1./gradCorLength2;
			}

			scale = corErrorOptimal[run] * oneOverScale ;

			// update parameters (gradient descent rule)
			muA = 0.;
			muB = 0.;
			muD = 0.;

			for (factor = 0; factor < m_nbFactors; factor++)
			{
				if (factor == 0)
				{
					noiseA = noiseVector[8 * i + factor];
					noiseB = noiseVector[8 * i + 2 + factor];
					noiseD = noiseVector[8 * i + 4 + factor];
				}
				else
				{
					noiseA = rhoLong * noiseVector[8 * i] +
						sqrtOneMinusRhoLong * noiseVector[8 * i + 1];
					noiseB = rhoLong * noiseVector[8 * i + 2] +
						sqrtOneMinusRhoLong * noiseVector[8 * i + 3];
					noiseD = rhoLong * noiseVector[8 * i + 4] +
						sqrtOneMinusRhoLong * noiseVector[8 * i + 5];
				}

				// move in direction perpendicular to gradVol
				// costruct direction by projecting gradCor onto
				// the surfice of "minimum" correlation error
				muA = scale * (gradCorA[factor] - gradProduct * gradVolA[factor]);
				muB = scale * (gradCorB[factor] - gradProduct * gradVolB[factor]);
				muD = scale * (gradCorD[factor] - gradProduct * gradVolD[factor]);

				a[factor]      += - muA * dt;// + volA * sqrtDt * noiseA);
				b[factor]      += - muB * dt;// + volB * sqrtDt * noiseB);
				d[factor]      += - muD * dt;// + volD * sqrtDt * noiseD);
			}

			// set new parameter values
			setParameters(a,b,d,lambda,m_humpParams);
			// calculate new total error
			totalError = getTotalError();
			corError = getCorError();

			// check if corr has converged and that there is no too much error
			// in vol, relative to result from above
			corErrorDiff = corError - corErrorOptimal[run];

			// new best error
			if ( (corErrorDiff < 0.) &&
				(totalError < 1.25 * totalErrorOptimalPrevious[run]) ) 
			{
				// increase dt 
				dt /= smallContract;
				//dt = startDt;
				sqrtDt = sqrt(dt);
				backStepCounter = 0;

				// new best error
				corErrorOptimal[run] = corError;
				totalErrorOptimal[run] = totalError;

				// store parameters which have resulted in best error so far
				for (factor=0; factor<m_nbFactors; factor++)
				{
					aOptimal[run][factor] = a[factor];
					bOptimal[run][factor] = b[factor];
					dOptimal[run][factor] = d[factor];
				}			
			}
			else // worse than best error
			{
				// revert
				for (factor = 0; factor < m_nbFactors; factor++)
				{
					a[factor] = aOptimal[run][factor];
					b[factor] = bOptimal[run][factor];
					d[factor] = dOptimal[run][factor];
				}

				// reduce dt in hope of finding minimum
				dt *= largeContract;
				sqrtDt = sqrt(dt);

				backStepCounter++;

			}  // end backtrack step

			// after being stuck in a plateau for maxCounter steps, we exit 
			if (backStepCounter > maxNbBackSteps)
			{
				break; // next run
			}			

		} // end for i (minimization loop)

	} // end for loop over runs

//SET_AND_RETURN :

	// find best run
	bestRun = 0;
	bestError = totalErrorOptimal[0] + corErrorOptimal[0];
	for(run = 1; run < nbRuns; run++)
	{
		if (totalErrorOptimal[run] + corErrorOptimal[run] < bestError)
		{
			bestRun = run;
			bestError = totalErrorOptimal[run] + corErrorOptimal[run];
		}
	}

	for(i = 0; i < m_nbFactors; i++)
	{
		m_a[i] = aOptimal[bestRun][i];
		m_b[i] = bOptimal[bestRun][i];
		m_d[i] = dOptimal[bestRun][i];
	}

	setParameters(m_a, 
		m_b, 
		m_d, 
		m_lambda,
		m_humpParams);

	totalError = getTotalError();
	corError = getCorError();

	numericCalibCaplets();

	totalError = getTotalError();
	corError = getCorError();

	// set correlations today
	setModelCorrelationMatrix(0.,1./365.); // 1 days
	setModelSwaptionGrid();

    printModel();

	//printRateModel();

	//checkCalibErrors(0.005,   // warning if above
	//	0.01,    // failure if above
	//	0.005,
    //	0.01);  

	return;
}


//****************************************************************************//
//  end: analyticCalib                                                        //
//****************************************************************************//

//*****************************************************************//
// returns the total error between market and model swaption vols. //
// for a given set of model's parameters                           //
//*****************************************************************//

/** fudgeCalib */
void SRMRatesLiborDiffuse::fudgeCalib(
							   const vector<double> &a, 
							   const vector<double> &b,
							   const vector<double> &d, 
							   const double &lambda,
							   const vector<double> &h)
{
	static const string method = "SRMRatesLiborDiffuse::fudgeCalib";

	const int nbRuns = 10;

	int i = -1,
		run = -1,
		nbBumps = -1;

	double 
		temp = -1.,
		length = 0.,
		oldError = 0.,
		newError = 0.,
		totalBump = -1.,
		deltaBumpA = -1.,
		newBumpA = -1.,
		firstExpiry = -1,
		lastExpiry = -1;

	vector<double>  
		firstExpiries,
		lastExpiries,
		gradVector,
		scale;

	SRM_VERIFY(a.size() == b.size(),
		"Error : a.size() != b.size() in fudgeCalib", 
		method); 

	SRM_VERIFY(a.size() == d.size(), 
		"Error : a.size() != d.size() in fudgeCalib",
		method); 

	nbBumps = m_calibIndices.size() - 1;
	gradVector.resize(nbBumps);	
	scale.resize(nbBumps);	
	firstExpiries.resize(nbBumps);	
	lastExpiries.resize(nbBumps);	

	// initialise
	deltaBumpA = 0.2;
	for(i = 0; i < nbBumps; i++)
	{
		scale[i] = 1.;
		firstExpiries[i] = today.yearFrac(m_gridSwap[m_calibIndices[i]]->getStartDate());
		lastExpiries[i] = (i < nbBumps - 1) ? 
			today.yearFrac(m_gridSwap[m_calibIndices[i+1]]->getStartDate()) :
		    100.;

		// for CMS spread calibration, this means that 
		// only the first column (lowest tenor) is 
		// modified by bumping the fudge parameter
		if (lastExpiries[i] <= firstExpiries[i]) 
		{
			nbBumps = i + 1;
			lastExpiries[i] = 100.;
			break;
		}
	}

	for(run = 0; run < nbRuns; run++)
	{
		oldError = getTotalError() ;
		length = 0.;
		totalBump = 0.;

		// calculate non-zero elements in gradVector
		for(i = 0; i < nbBumps; i++)
		{
			firstExpiry = firstExpiries[i]; 
			lastExpiry = lastExpiries[i];

			// for CMS spread calibration, this means that 
			// only the first column (lowest tenor) is 
			// modified by bumping the fudge parameter
			if (lastExpiry <= firstExpiry) 
			{
				throw ModelException(method, "Error in BGM calibration (fudgeCalib)");
			}

			bumpFudge(firstExpiry, lastExpiry, deltaBumpA, 0.); 
			setParameters(a,b,d,lambda,h);
			newError = getTotalError();
			// revert to starting point
			bumpFudge(firstExpiry, lastExpiry, - deltaBumpA, 0.); 

			if (oldError > TOO_SMALL)
			{
				temp = (newError - oldError) / oldError;
			}

			newBumpA = - 10. * temp * scale[i];
			bumpFudge(firstExpiry, lastExpiry, newBumpA, 0.);
			setParameters(a,b,d,lambda,h);
			newError = getTotalError();	

			if (1.001 * oldError < newError)
			{
				bumpFudge(firstExpiry, lastExpiry, - newBumpA, 0.); 
				setParameters(a,b,d,lambda,h);
				scale[i] *= 0.5;
			}
			else
			{
				oldError = newError;
				scale[i] *= 1.5;
			}
		}
	}

	return;
}

//****************************************************************************//
//  end: fudgeCalib                                                           //
//****************************************************************************//

/** calibrate caplets */
void SRMRatesLiborDiffuse::calibrateCaplets(void)
{
	static const string method = "SRMRatesLiborDiffuse::calibrateCaplets";

	int row = -1 , libor = -1, firstRow = -1;

	double 
		firstExpiry, lastExpiry,
		lastDate, firstDate, resetDate, previousDate, lastPaymentDate,
		previousExpiry, nextExpiry,
		temp = -1.0, lastCapVol = -1., longVol = -1.,
		volSum = -1.;
	
	int nbExpiries = m_marketExpiries.size();
	int nbTenors = m_marketTenors.size();

	firstExpiry = m_marketExpiries[0];
	lastExpiry = m_marketExpiries[nbExpiries-1];

	// last volatility that is known from 1st column of swaption grid
	lastCapVol = m_marketSwaptionVols[nbExpiries-1][0];

	// note : empirically, a decay to asymptotic libor vol = longVol
	//        over 10 years from last swaption expiry date seems to work well
	//        could be calibrated, in principle
	lastPaymentDate = 100.;

	// generally, this is the smallest caplet vol of interest
	// - it is the approximate vol of libors resetting after lastPaymentDate
	longVol = m_marketSwaptionVols[nbExpiries-1][nbTenors-1];

	row = 0;
	for (libor = 0; libor < m_nbLibors; libor++)
	{
		if ( m_resetDates[libor] <= firstExpiry )
		{
			m_liborVols[libor] = m_marketSwaptionVols[0][0];
		}
		else if ( m_resetDates[libor] >= lastExpiry )
		{
			temp = (m_resetDates[libor] - lastExpiry) /
				   m_marketTenors[nbTenors-1];
			temp = Maths::min(1., temp);

			m_liborVols[libor] = lastCapVol + temp * (longVol - lastCapVol);
		}
		else
		{
			nextExpiry = m_marketExpiries[row];

			while ( 
				     (m_resetDates[libor] > nextExpiry) &&
				     (row < nbExpiries - 1) // should always be satisfied
				  )
			{
				row++;
				nextExpiry = m_marketExpiries[row];
			}

			if (row == 0) 	
			{
				m_liborVols[libor] = m_marketSwaptionVols[row][0];
			}
			else	// row >= 1 && row < nbExpiries - 1
			{
				previousExpiry = m_marketExpiries[row-1];

				temp = (nextExpiry - m_resetDates[libor]) /
				       (nextExpiry - previousExpiry);

				// add check that libor pay date == expiryDate
				// before assigning the m_liborVols
				double vol2Next = m_marketSwaptionVols[row][0],
					   vol2Prev = m_marketSwaptionVols[row-1][0];

				vol2Next *= vol2Next;
				vol2Prev *= vol2Prev;

				m_liborVols[libor] = sqrt(vol2Next + temp * (vol2Prev - vol2Next));
			}

			m_liborVols[libor] = Maths::max(BGM_MIN_VOL, m_liborVols[libor]);
		}


	} // end for libor


	// fix for short-dated libors (say, resets at 1m, 3m, 6m, 12m, 18m, tenor = 12m)
	// libor vols are interpolated linearly
	firstRow = 0;
	while ( (firstRow < nbExpiries - 2) &&
		    (SRMRound(m_marketExpiries[firstRow+1] - m_marketExpiries[firstRow]) < m_liborInterval/12.) ) 
	{
		firstRow++;
	}

	row = firstRow; 
	while ( row >=0 )
	{
		firstDate = m_marketExpiries[row];
		lastDate = m_marketExpiries[row+1]; // ok as row <= nbExpiries - 2

		// find first libor with reset in interval
		libor = m_nbLibors - 1;
		while ( m_resetDates[libor] >= lastDate )
		{
			libor--;
		}

		// we can't modify vols of libors contributing in later intervals,
		// therefore we look only at libors paying in this period, but we
		// still need to estimate contribution from later libors
		double invAnnuity = 1.;
		volSum = 0.;
		resetDate = lastDate;
		while (m_resetDates[libor] > firstDate)
		{
			previousDate = resetDate;
			resetDate = m_resetDates[libor];
			double delta = (previousDate - resetDate) /
			               (lastDate - firstDate);
			invAnnuity *= 1. + m_initialLibor[libor] * delta;
			volSum -= m_initialLibor[libor] * m_liborVols[libor] * delta;
			libor--;
		}

		// this corresponds to assuming the short-end libors
		// perfectly correlated

		previousDate = resetDate;
		resetDate = m_resetDates[libor];

		double delta = (previousDate - resetDate) /
		               (lastDate - firstDate);

		// this is an approximation
		invAnnuity *= 1. + m_initialLibor[libor] * delta;

		volSum += (invAnnuity - 1.) * m_marketSwaptionVols[row][0];

		volSum /= m_initialLibor[libor] * delta;


		if ( volSum > TOO_SMALL )
		{
			m_liborVols[libor] = volSum;
		}
		else
		{
			m_liborVols[libor] = m_marketSwaptionVols[row][0];
		}

		// to avoid too small vols (which may not get modified by the later 
		// numerical calibration)
		m_liborVols[libor] = Maths::max(0.1 * m_marketSwaptionVols[row][0], 
			                            m_liborVols[libor]);

		row--;
	}

	return;
}

//****************************************************************************//
//  end : calibrateCaplets                                                    //
//****************************************************************************//

/** numericCalibCaplets */
// calibrates libor vols for rates with resets < 5.
void SRMRatesLiborDiffuse::numericCalibCaplets()
{
	static const string method = "SRMRatesLiborDiffuse::numericCalibCaplets";

	const int 
		nbSteps = 100,
		maxNbBackSteps = 10;

	const double 
		delta = 0.01,
		eps = 0.0001,
		startDt = 1.,
		largeContract = 0.5,
		smallContract = 0.25;

	int libor = -1,
		row = -1,
		calibRow = -1,
		backStepCounter = 0,
		liborFirst,
		liborLast,
		lastLiborIndex,
		i;

	double 
		lambda = 0.,
		error = -1.0,
		newError = -1.0,
		errorOptimal = -1.0,
		errorDiff = -1.0,
		gradLength   = -1.0,
		errorScale = -1.0,
		bump = -1.,
		dt,
		scale,
		firstDate,
		lastDate,
		expiry;

	vector<double> 
		sigmaBarOptimal,
		gradSigmaBar,
		a(m_nbFactors),
		b(m_nbFactors),
		d(m_nbFactors);

	int nbExpiries = m_marketExpiries.size(),
		nbTenors = m_marketTenors.size();

	lambda = m_lambda;
	// set parameters with current values
	for ( i = 0; i < m_nbFactors; i++)
	{
		a[i] = m_a[i];
		b[i] = m_b[i];
		d[i] = m_d[i];
	}

	calibRow = 0;
    while ( (calibRow < nbExpiries - 1) &&
			(m_marketExpiries[calibRow+1] < 5.) 
		  ) 
	{
		calibRow++;
	}

	lastLiborIndex = 0;
	expiry = m_marketExpiries[calibRow];

	while ( 
		(m_payDates[lastLiborIndex] <= expiry) &&
		(lastLiborIndex < m_nbLibors - 1) )
	{
		lastLiborIndex++;
	}

	sigmaBarOptimal.resize(lastLiborIndex);
	gradSigmaBar.resize(lastLiborIndex);

	// initialize sigmaBar
	for (libor = 0; libor < lastLiborIndex; libor++)
	{
		sigmaBarOptimal[libor] = m_liborVols[libor];
	}

	errorScale = eps * eps;

	for (row = calibRow - 1; row >= 0; row--)
	{
		backStepCounter = 0;
		bump = delta;
		dt = startDt;

		setParameters(a,b,d,lambda,m_humpParams);
		error = getVolError(row,0);
		errorOptimal = error;

		firstDate = m_marketExpiries[row];
		lastDate = m_marketExpiries[row+1];

		// first libor with pay date after firstDate
		liborFirst = 0;
		while ( (m_payDates[liborFirst] <= firstDate) &&
		    	(liborFirst < m_nbLibors - 1) )
		{ 
			liborFirst++;
		}

		// last libor with pay date on or before lastDate
		liborLast = liborFirst;
		while ( (m_payDates[liborLast+1] <= lastDate) &&
			    (liborLast < m_nbLibors -1) )
		{
			liborLast++;
		}

		liborLast = Maths::min(liborLast,lastLiborIndex-1);

		// first check if we need to do anything
		if ( fabs(error) < errorScale )
		{
			for (libor = liborFirst; libor <= liborLast; libor++)
			{
				sigmaBarOptimal[libor] = m_liborVols[libor];
			}

			// done!
			continue;
		}

		// if not good enough, start gradient descent algorithm
		for (i = 0; i < nbSteps; i++ )
		{
			// get current error
			error = errorOptimal;

			// calculate numerically gradient
			gradLength = 0.0;
			for (libor = liborLast; libor >= liborFirst; libor--)
			{
				m_liborVols[libor] += bump;
				setParameters(a,b,d,lambda,m_humpParams);
				newError = getVolError(row,0);
				m_liborVols[libor] -= bump;

				gradSigmaBar[libor] = (newError - error) / bump;
				gradLength += gradSigmaBar[libor] * gradSigmaBar[libor];
			}

			gradLength = sqrt(gradLength);

			// update parameters
			if ( gradLength > TOO_SMALL )
			{
				scale = 0.001 / gradLength; // scale is 0.1%
			}
			else
			{
				break;
			}

			for (libor = liborFirst; libor <= liborLast; libor++)
			{
				m_liborVols[libor] += - scale * gradSigmaBar[libor] * dt;
			    // floor to avoid spurious minima at negative vols
				m_liborVols[libor] = Maths::max(BGM_MIN_VOL, m_liborVols[libor]);
			}

			// get new error
			setParameters(a,b,d,lambda,m_humpParams);
			error = getVolError(row,0);

			if ( fabs(error) < errorScale )
			{
				for (libor = liborFirst; libor <= liborLast; libor++)
				{
					sigmaBarOptimal[libor] = m_liborVols[libor];
				}

				// done!
				break;
			}

			errorDiff = error - errorOptimal;
			if ( errorDiff < 0. )
			{
				// increase dt
				dt /= largeContract;
				bump /= largeContract;
				backStepCounter = 0;

				// new best error
				errorOptimal = error;

				// store parameters
				for (libor = liborFirst; libor <= liborLast; libor++)
				{
					sigmaBarOptimal[libor] = m_liborVols[libor];
				}
			}
			else
			{
				// revert
				for (libor = liborFirst; libor <= liborLast; libor++)
				{
					m_liborVols[libor] = sigmaBarOptimal[libor];
				}

				//setParameters(m_aInit,m_bInit,m_dInit,m_lambdaInit);

				dt *= smallContract;
				bump *= smallContract;
				bump = Maths::max(bump, 0.0001);

				backStepCounter++;
			}

			// after being stuck in a plateau for maxCounter steps, we exit
			if (backStepCounter > maxNbBackSteps)
			{
				break;
			}
		} // end i (over steps)
	}// end row

	for (libor = 0; libor < lastLiborIndex; libor++)
	{
		m_liborVols[libor] = sigmaBarOptimal[libor];
		m_liborVols[libor] = Maths::max(BGM_MIN_VOL, m_liborVols[libor]);
	}

	setParameters(a,b,d,lambda,m_humpParams);

	return;
}
//****************************************************************************//
//  end : numericCalibCaplet                                                  //
//****************************************************************************//

/** bumpFudge */
void SRMRatesLiborDiffuse::bumpFudge(
									 double startDate, 
									 double endDate,
									 double bumpA,
									 double bumpD)
{
	static const string method = "SRMRatesLiborDiffuse::bumpFudge";

	int i = 0;

	// find start
	while (m_resetDates[i] < startDate)
	{
		i++;
	}

	while ((i < m_nbLibors) && 
		(m_resetDates[i] < endDate)
		)
	{
		m_fudgeA[i] += bumpA;
		m_fudgeD[i] += bumpD;
		i++;
	}

	return;
}

//****************************************************************************//
//  end: bumpFudge                                                            //
//****************************************************************************//

/** setCalibDetails */
void SRMRatesLiborDiffuse::setCalibDetails(
	const string& calibStyle,
	double calibMat,
	double calibMatCMS)
{
	static const string method = "SRMRatesLiborDiffuse::setCalibDetails";

	int 
		calibCol  = 0,
		calibCol1 = 0,
		calibCol2 = 0,
		offset1 = 0,
		offset2 = 0,
		swaptionIndex = -1,
		swaptionIndex1 = -1,
		swaptionIndex2 = -1,
		row = -1,
		col = -1,
		expiry = -1;

	m_calibIndices.reserve(100);
	m_marketCalibVol.reserve(100);

	// find five-year expiry row
	for (row = 0; row < (int)m_marketExpiries.size(); row++)
	{
		if (SRMRound(m_marketExpiries[row]) >= 5)
		{
			offset1 = row *  m_marketTenors.size();
			break;
		}
	}

	// find long expiry row
	for (row = 0; row < (int)m_marketExpiries.size(); row++)
	{
		if (SRMRound(m_marketExpiries[row]) >= 20)
		{
			offset2 = row *  m_marketTenors.size();
			break;
		}
	}

	if (calibStyle == SwapMaturityVolRequest::CALIB_CMS)
	{
		// find first tenor above target - if not found, 
		// we use the last column in grid => flat extrapolation
		while (calibCol1 < (int)m_marketTenors.size())
		{
			if (m_marketTenors[calibCol1] >= calibMat)
			{
				break;
			}
			if (calibCol1 == m_marketTenors.size() - 1)
			{
				break;
			}
			calibCol1++;
		}

		while (calibCol2 < (int)m_marketTenors.size())
		{
			if (m_marketTenors[calibCol2] >= calibMatCMS)
			{
				break;
			}
			if (calibCol2 == m_marketTenors.size() - 1)
			{
				break;
			}
			calibCol2++;
		}

		const int 
			firstCol = Maths::min(calibCol1, calibCol2),
			lastCol  = Maths::max(calibCol1, calibCol2);

		swaptionIndex1 = firstCol;
		swaptionIndex2 = lastCol;

		// first column 1
		for (row = 0; row < (int)m_marketExpiries.size(); row++)
		{				
			// only calibrate to swaptions with expiry >= 1Y 
			if (m_marketExpiries[row] >= 1.)
			{
				m_calibIndices.push_back(swaptionIndex1);
				m_marketCalibVol.push_back(m_marketSwaptionVols[row][firstCol]);
			}

			// update swaptionIndex
			swaptionIndex1 += m_marketTenors.size();

		}

		m_corrSwapIndices[0] = offset1 + firstCol;
		m_corrSwapIndices[1] = offset1 + lastCol;

		if (lastCol > firstCol)
		{
			// then column 2
			for (row = 0; row < (int)m_marketExpiries.size(); row++)
			{		
				// only calibrate to swaptions with expiry >= 1Y 
				if (m_marketExpiries[row] >= 1.)
				{
					m_calibIndices.push_back(swaptionIndex2);
					m_marketCalibVol.push_back(m_marketSwaptionVols[row][lastCol]);
				}

				// update swaptionIndex
				swaptionIndex2 += m_marketTenors.size();
			}
		}

		m_corrSwapIndices[2] = offset2 + firstCol;
		m_corrSwapIndices[3] = offset2 + lastCol;

	} //end case CMSspred
	else if (calibStyle == SwapMaturityVolRequest::CALIB_FIX)
	{
		// code for diagonal calibration
		for (row=0; row < (int)m_marketExpiries.size(); row++)
		{		
			// for the given row (i.e expiry) find column with tenor
			// which results expiry+tenor = calibTenor
			calibCol = 0;
			while (calibCol < (int)m_marketTenors.size())
			{
				long temp = m_marketTenors[calibCol] + m_marketExpiries[row];
				if ( (m_marketTenors[calibCol] + m_marketExpiries[row]) >= calibMat )
				{
					break;
				}
				if (calibCol == m_marketTenors.size() - 1)
				{
					break;
				}
				calibCol++;
			}

			// we don't want to use the first column, it does
			// not contribute to the error because it is already calibrated
			if ( calibCol == 0 )
			{
				break;
			}

			// only calibrate to swaptions with expiry >= 1Y 
			if (m_marketExpiries[row] >= 1.)
			{
				swaptionIndex = row * m_marketTenors.size() + calibCol;                   
				m_calibIndices.push_back(swaptionIndex);
				m_marketCalibVol.push_back(m_marketSwaptionVols[row][calibCol]);
			}

		}// end for row

		m_corrSwapIndices[0] = offset1;
		m_corrSwapIndices[1] = offset1 + calibCol;

		m_corrSwapIndices[2] = offset2;
		m_corrSwapIndices[3] = offset2 + calibCol;

	} // end case diagonal
	else
	{
		throw ModelException(method, "Calibration style (CMS/FIX) not specified");
	}  

	return;

}

//****************************************************************************//
//  end: setCalibDetails                                                      //
//****************************************************************************//

double SRMRatesLiborDiffuse::getModelBlackVolatility(const SRMSwapClass &swap)
{
	static const string method = "SRMRatesLiborDiffuse::getModelBlackVolatility";

	int i  = -1,
		j  = -1,
		li = -1,
		lj = -1,
		liStart = -1,
		liEnd = -1,
		ljStart = -1,
		ljEnd = -1,
		factor = -1,
		firstLiborAlive = -1,
		lastLiborInFirstPeriod = -1,
		lastLiborUsed = -1,
		nbFixPay = -1;

	double 
		covLower = -1.,
		startLower = -1.,
		blackVol = 0.,
		weight_i = -1.,
		weight_j = -1.;

	vector<double> swapStep;
	vector<double> weights;
	vector<double> beginDates;
	vector<double> endDates;

	DateTimeArray swapResetDates(swap.getResetDates());
	DateTimeArray swapPayDates(swap.getPayDates()); // fix, that is

	nbFixPay = swapResetDates.size();

	SRM_VERIFY(nbFixPay >= 1,
		"Invalid swap used for calibration", 
		method); 

	double swapStart = today.yearFrac(swap.getStartDate());
	double swapEnd = today.yearFrac(swap.getMatDate());

	weights.resize(nbFixPay);
	swapStep.resize(nbFixPay);

	beginDates.resize(nbFixPay);
	endDates.resize(nbFixPay);

	for(i = 0; i < nbFixPay; i++)
	{
		beginDates[i] = today.yearFrac(swapResetDates[i]);
		endDates[i] = today.yearFrac(swapPayDates[i]);
		swapStep[i] = swapResetDates[i].yearFrac(swapPayDates[i]);
	}

	// find indices of first and last LIBORS of the swap

	// returns index of first libor with reset on or before swapStart
	firstLiborAlive = 0;
	while ( (firstLiborAlive < m_nbLibors - 1) &&
		    (m_resetDates[firstLiborAlive+1] <= swapStart)
		  )
	{
		firstLiborAlive++;
	}

	lastLiborInFirstPeriod = firstLiborAlive;
	while ( (lastLiborInFirstPeriod < m_nbLibors - 1) &&
		    (m_resetDates[lastLiborInFirstPeriod+1] <= swapStart)
		  )
	{
		lastLiborInFirstPeriod++;
	}

	SRM_VERIFY(firstLiborAlive != m_nbLibors - 1,
		"Can't calibrate to this swaption : last libor reset < swaption expiry", 
		method); 

	// returns index of last libor with reset on or after swapEnd
	// don't change the < back again!
	lastLiborUsed = firstLiborAlive;
	while ( (lastLiborUsed < m_nbLibors - 1) &&
		    (m_resetDates[lastLiborUsed+1] < swapEnd)
		  )
	{
		lastLiborUsed++;
	}

	// get the appropriate weights for the Rebonato formula of the
	// swaption volatility
	weights = swap.getWeights();

	SRM_VERIFY(weights.size() == nbFixPay, 
		"Error : weights array has wrong size (internal error, getModelBlackVolatility)", 
		method);

	// get the necessary factor loadings to compute the correlation
	// we need separate loadings for upper and lower because first Libor
	// reset date for lower will be before, or at, the swap start date
	startLower = m_resetDates[firstLiborAlive];

	if (startLower < TOO_SMALL)
	{
		throw ModelException(method, "Error in Libor reset dates schedule (getModelBlackVolatility)");
	}

	getFactorVols(
		m_facLdArray,
		firstLiborAlive,
		lastLiborUsed,
		0, swapStart);

	liEnd = firstLiborAlive;
	for (i = 0; i < nbFixPay; i++ )
	{
		liStart = liEnd;
		ljEnd = firstLiborAlive; 

		for (j = 0; j < nbFixPay; j++ )
		{
			li = liStart;
			ljStart = ljEnd;

			covLower = 0.;
			while ( (li <= lastLiborUsed) &&
			    	(m_resetDates[li] < endDates[i]) )
			{
				weight_i = (Maths::min(m_payDates[li],endDates[i])
					      - Maths::max(m_resetDates[li],beginDates[i]))
					/ swapStep[i];
				lj = ljStart;

				while ( (lj <= lastLiborUsed ) &&
					    (m_resetDates[lj] < endDates[j]) )
				{
					weight_j = (Maths::min(m_payDates[lj],endDates[j])
						      - Maths::max(m_resetDates[lj],beginDates[j])) / swapStep[j];

					for (factor = 0; factor < m_nbFactors; factor++)
					{
						covLower += weight_i * weight_j * (
						   (m_facLdArray[li][3*factor]   * m_facLdArray[lj][3*factor] +
							m_facLdArray[li][3*factor+1] * m_facLdArray[lj][3*factor+1] +
							m_facLdArray[li][3*factor+2] * m_facLdArray[lj][3*factor+2]) );
					} // for factor
					lj++;
				}
				li++;
			}

			ljEnd = lj;
			if (m_resetDates[ljEnd] > endDates[j])
			{
				if (ljEnd > 0)
				{
					ljEnd--;
				}
			}

			blackVol += weights[i] * weights[j] * covLower / startLower;

			// the weights used here refer to the expression w_{i}(0)L_{i}(0)/Sr(0)
			// in the analytic Rebonato formula for the swaption Black volatility
		} // for j

		liEnd = li;
		// check if there will be a stub in next round
		// and adjust accordingly
		if (m_resetDates[liEnd] > endDates[i])
		{
			if (liEnd > 0)
			{
				liEnd--;
			}
		}

	} // for i

	weights.clear();

	blackVol = sqrt(blackVol);

	return blackVol;
}

//****************************************************************************//
//  end : getModelBlackVolatility                                             //
//****************************************************************************//

/** setModelSwaptionGrid */
void SRMRatesLiborDiffuse::setModelSwaptionGrid()
{
	long row, col, swapIndex;

	m_modelSwaptionVols.resize(m_marketExpiries.size());

	swapIndex = 0;
	for (row = 0; row < (int)m_marketExpiries.size(); row++)
	{
		m_modelSwaptionVols[row].resize(m_marketTenors.size());
		for (col = 0; col < (int)m_marketTenors.size(); col++)
		{
			m_modelSwaptionVols[row][col] = 
				getModelBlackVolatility( *(m_gridSwap[swapIndex]));
			swapIndex++;
		}
	}

	return;
}


//****************************************************************************//
//  end : setModelSwaptionGrid                                                //
//****************************************************************************//

/** setModelCorrelationMatrix */
void SRMRatesLiborDiffuse::setModelCorrelationMatrix(double time1, double time2)
{
	long i = -1,
		j = -1,
		factor = -1;

	double cov_ij = -1,
		var_i  = -1,
		var_j  = -1,
		sqrtVar_i = -1,
		sqrtVar_j = -1;

	const long maxNbCorr = Maths::min(100L,m_nbLibors-1);

	getFactorVols(m_facLdArray, 0, maxNbCorr-1, time1, time2);

	for (i=0; i<maxNbCorr; i++)
	{
		for (j=0; j<maxNbCorr; j++)
		{
			var_i = 0.0;
			var_j = 0.0;
			cov_ij = 0.0;
			for (factor=0; factor<m_nbFactors; factor++)
			{
				cov_ij += (m_facLdArray[i][3*factor] * m_facLdArray[j][3*factor] +
					m_facLdArray[i][3*factor+1] * m_facLdArray[j][3*factor+1] +
					m_facLdArray[i][3*factor+2] * m_facLdArray[j][3*factor+2]);

				var_i += (m_facLdArray[i][3*factor] * m_facLdArray[i][3*factor] +
					m_facLdArray[i][3*factor+1] * m_facLdArray[i][3*factor+1] +
					m_facLdArray[i][3*factor+2] * m_facLdArray[i][3*factor+2]);

				var_j += (m_facLdArray[j][3*factor] * m_facLdArray[j][3*factor] +
					m_facLdArray[j][3*factor+1] * m_facLdArray[j][3*factor+1] +
					m_facLdArray[j][3*factor+2] * m_facLdArray[j][3*factor+2]);
			}

			sqrtVar_i = sqrt(var_i);
			sqrtVar_j = sqrt(var_j);
			if ( ( sqrtVar_i > TOO_SMALL ) && ( sqrtVar_j > TOO_SMALL ) )
			{
				m_modelCorrelations[i][j] = cov_ij/(sqrtVar_i * sqrtVar_j);
			}
			else
			{
				// correlation not defined - devide by zero
				m_modelCorrelations[i][j] = 1000.0;
			}

		}
	}

	return;
}

//****************************************************************************//
//  end : setModelCorrelationMatrix                                           //
//****************************************************************************//							

// get vol error for swaption at row = calibRow and column = calibColumn
double SRMRatesLiborDiffuse::getVolError(long calibRow, long calibColumn)
{
	static const string method = "SRMRatesLiborDiffuse::getVolError";

	int row = -1,
		column = -1,
		swaptionIndex = -1,
		nbVols = 0,
		nbExpiries = m_marketExpiries.size(),
		nbTenors = m_marketTenors.size();

	double error = 0.0,
		marketBlackVol = -1. ,
		modelBlackVol  = -1. ,
		currentError   = -1. ;

	if ( ( calibRow < 0 ) || ( calibRow >= nbExpiries ) )
	{
		throw ModelException(method, "Invalid row index in SRMRatesLiborDiffuse::getVolError");
	}

	if ( ( calibColumn < 0 ) || ( calibColumn >= nbTenors ) )
	{
		throw ModelException(method, "Invalid column index in SRMRatesLiborDiffuse::getVolError");
	}

	// calculate error between market - model vols for current
	// values of the parameters (matrices m_A, m_B, m_D)
	swaptionIndex = calibColumn;

	for (row = 0; row < calibRow; row++)
	{
		// update swaptionIndex
		swaptionIndex += nbTenors;
	}

	marketBlackVol = m_marketSwaptionVols[calibRow][calibColumn];
	modelBlackVol = getModelBlackVolatility(*(m_gridSwap[swaptionIndex]));
	error = (marketBlackVol - modelBlackVol) *
		(marketBlackVol - modelBlackVol);

	return error;
}

/** getSwapRateCorr */
double SRMRatesLiborDiffuse::getSwapRateCorr(
	SRMSwapClass &swap1,
	SRMSwapClass &swap2,
	double start,
	double end)
{
	static const string method = "SRMRatesLiborDiffuse::getSwapRateCorr";

	int i  = -1,
		j  = -1,
		li = -1, // libor index for i
		lj = -1, // libor index for j
		liStart = -1,  
		liEnd = -1,
		ljStart = -1,
		ljEnd = -1,
		factor = -1,
		firstLiborIndex1 = -1,
		firstLiborIndex2 = -1,
		lastLiborIndex1 = -1,
		lastLiborIndex2 = -1,
		firstLiborIndex = -1,
		lastLiborIndex = -1,
		nbFixPay1 = -1,
		nbFixPay2 = -1;

	double 
		covLower = -1.,
		startLower1 = -1.,
		startLower2 = -1. ,
		correlation = 0.,
		var1 = 0.,
		var2 = 0.,
		sqrtVar1 = 0.,
		sqrtVar2 = 0.,
		weight_i = -1.,
		weight_j = -1.;

	vector<double> 
		swapStep1, swapStep2,
		weights1, weights2, 
		swapBeginTime1,
		swapEndTime1,
		swapBeginTime2,
		swapEndTime2;

	const DateTimeArray
		swapResetDates1 = swap1.getResetDates(),
		swapResetDates2 = swap2.getResetDates(),
		swapPayDates1 = swap1.getPayDates(),
		swapPayDates2 = swap2.getPayDates();

	nbFixPay1  = swapPayDates1.size();
	double swapStart1 = today.yearFrac(swap1.getStartDate());
	double swapEnd1   = today.yearFrac(swap1.getMatDate()); 

	nbFixPay2  = swapPayDates2.size();
	double swapStart2 = today.yearFrac(swap2.getStartDate());
	double swapEnd2   = today.yearFrac(swap2.getMatDate()); 

	weights1.resize(nbFixPay1);
	swapStep1.resize(nbFixPay1);
	swapBeginTime1.resize(nbFixPay1);
	swapEndTime1.resize(nbFixPay1);

	weights2.resize(nbFixPay2);
	swapStep2.resize(nbFixPay2);
	swapBeginTime2.resize(nbFixPay2);
	swapEndTime2.resize(nbFixPay2);

	double nextDate = today.yearFrac(swapResetDates1[0]);
	double previousDate;
	for(i = 0; i < nbFixPay1; i++)
	{
		previousDate = nextDate;
		nextDate = today.yearFrac(swapPayDates1[i]);
		swapBeginTime1[i] = previousDate;
		swapEndTime1[i] = nextDate;
		swapStep1[i] = (double) (nextDate - previousDate);
	}

	nextDate = today.yearFrac(swapResetDates2[0]);
	for(i = 0; i < nbFixPay2; i++)
	{
		previousDate = nextDate;
		nextDate = today.yearFrac(swapPayDates2[i]);
		swapBeginTime2[i] = previousDate;
		swapEndTime2[i] = nextDate;
		swapStep2[i] = (double) (nextDate - previousDate);
	}

	// find indices of first and last LIBORS of the swap

	// returns index of first libor with reset on or before swapStart
	firstLiborIndex1 = 0;
	while ( 
		(firstLiborIndex1 < m_nbLibors - 1) &&
		(m_resetDates[firstLiborIndex1+1] <= swapStart1) 
		)
	{
		firstLiborIndex1++;
	}

	SRM_VERIFY(
		firstLiborIndex1 != m_nbLibors - 1,
		"Can't calibrate to this swaption : last libor reset < swaption expiry",
		method);

	// returns index of last libor with reset on or after swapEnd
	// don't change the < back again!
	lastLiborIndex1 = firstLiborIndex1;
	while ( 
		(lastLiborIndex1 < m_nbLibors - 1) &&
		(m_resetDates[lastLiborIndex1+1] < swapEnd1) 
		)
	{
		lastLiborIndex1++;
	}

	// returns index of first libor with reset on or before swapStart
	firstLiborIndex2 = 0;
	while ( 
		(firstLiborIndex2 < m_nbLibors - 1) &&
		(m_resetDates[firstLiborIndex2+1] <= swapStart2) 
		)
	{
		firstLiborIndex2++;
	}

	SRM_VERIFY(
		firstLiborIndex2 != m_nbLibors - 1,
		"Can't calibrate to this swaption : last libor reset < swaption expiry",
		method);

	// returns index of last libor with reset on or after swapEnd
	// don't change the < back again!
	lastLiborIndex2 = firstLiborIndex2;
	while ( 
		(lastLiborIndex2 < m_nbLibors - 1) &&
		(m_resetDates[lastLiborIndex2+1] < swapEnd2) 
		)
	{
		lastLiborIndex2++;
	}

	// get the appropriate weights for the Rebonato formula of the
	// swaption volatility
	weights1 = swap1.getWeights();

	SRM_VERIFY(
		weights1.size() == nbFixPay1,
		"Error : weights array has wrong size (internal error, getModelBlackVolatility)",
		method);

	// get the appropriate weights for the Rebonato formula of the
	// swaption volatility
	weights2 = swap2.getWeights();

	SRM_VERIFY(
		weights2.size() == nbFixPay2,
		"Error : weights array has wrong size (internal error, getModelBlackVolatility)",
		method);

	// get the necessary factor loadings to compute the correlation
	// we need separate loadings for upper and lower because first Libor
	// reset date for lower will be before, or at, the swap start date
	startLower1 = m_resetDates[firstLiborIndex1];
	startLower2 = m_resetDates[firstLiborIndex2];

	if (startLower1 < TOO_SMALL || startLower2 < TOO_SMALL) 
	{
		throw ModelException(method, "Error in Libor reset dates schedule");
	}

	firstLiborIndex = Maths::min(firstLiborIndex1, firstLiborIndex2);
	lastLiborIndex  = Maths::max(lastLiborIndex1, lastLiborIndex2);

	getFactorVols(m_facLdArray, 
		firstLiborIndex, 
		lastLiborIndex, 
		start, end);

	// covariance first
	liEnd = firstLiborIndex1;
	for (i = 0; i < nbFixPay1; i++ )
	{	
		liStart = liEnd;
		ljEnd = firstLiborIndex2;

		for (j = 0; j < nbFixPay2; j++ )
		{
			li = liStart;
			ljStart = ljEnd;

			covLower = 0.;
			while ( (li <= lastLiborIndex1 ) && 
				(m_resetDates[li] < swapEndTime1[i]) )
			{
				weight_i = (Maths::min(m_payDates[li],swapEndTime1[i])
					- Maths::max(m_resetDates[li],swapBeginTime1[i]))
					/ swapStep1[i];
				lj = ljStart;

				while ( (lj <= lastLiborIndex2 ) && 
					(m_resetDates[lj] < swapEndTime2[j]) )
				{
					weight_j = (Maths::min(m_payDates[lj],swapEndTime2[j])
						- Maths::max(m_resetDates[lj],swapBeginTime2[j])) / swapStep2[j];	

					for (factor = 0; factor < m_nbFactors; factor++)
					{
						covLower += weight_i * weight_j *
							((m_facLdArray[li][3*factor] * m_facLdArray[lj][3*factor] +
							m_facLdArray[li][3*factor+1] * m_facLdArray[lj][3*factor+1] +
							m_facLdArray[li][3*factor+2] * m_facLdArray[lj][3*factor+2]) );

					} // for factor

					lj++;
				}
				li++;
			}

			ljEnd = lj;
			if (m_resetDates[ljEnd] > swapEndTime2[j])
			{
				if (ljEnd > 0)
				{
					ljEnd--;
				}
			}

			correlation += weights1[i] * weights2[j] * covLower;

			// the weights used here refer to the expression w_{i}(0)L_{i}(0)/Sr(0)
			// in the analytic Rebonato formula for the swaption Black volatility
		} // for j

		liEnd = li;
		// check if there will be a stub in next round 
		// and adjust accordingly
		if (m_resetDates[liEnd] > swapEndTime1[i])
		{
			if (liEnd > 0)
			{
				liEnd--;
			}
		}

	} // for i

	liEnd = firstLiborIndex1;
	for (i = 0; i < nbFixPay1; i++ )
	{	
		liStart = liEnd;
		ljEnd = firstLiborIndex1;

		for (j = 0; j < nbFixPay1; j++ )
		{
			li = liStart;
			ljStart = ljEnd;

			covLower = 0.;
			while ( (li <= lastLiborIndex1 ) && 
				(m_resetDates[li] < swapEndTime1[i]) )
			{
				weight_i = (Maths::min(m_payDates[li],swapEndTime1[i])
					- Maths::max(m_resetDates[li],swapBeginTime1[i]))
					/ swapStep1[i];
				lj = ljStart;

				while ( (lj <= lastLiborIndex1 ) && 
					(m_resetDates[lj] < swapEndTime1[j]) )
				{
					weight_j = (Maths::min(m_payDates[lj],swapEndTime1[j])
						- Maths::max(m_resetDates[lj],swapBeginTime1[j])) / swapStep1[j];	

					for (factor = 0; factor < m_nbFactors; factor++)
					{
						covLower += weight_i * weight_j *
							((m_facLdArray[li][3*factor]   * m_facLdArray[lj][3*factor] +
							m_facLdArray[li][3*factor+1] * m_facLdArray[lj][3*factor+1] +
							m_facLdArray[li][3*factor+2] * m_facLdArray[lj][3*factor+2]) );

					} // for factor

					lj++;
				}
				li++;
			}

			ljEnd = lj;
			if (m_resetDates[ljEnd] > swapEndTime1[j])
			{
				if (ljEnd > 0)
				{
					ljEnd--;
				}
			}

			var1 += weights1[i] * weights1[j] * covLower;

			// the weights used here refer to the expression w_{i}(0)L_{i}(0)/Sr(0)
			// in the analytic Rebonato formula for the swaption Black volatility
		} // for j

		liEnd = li;
		// check if there will be a stub in next round 
		// and adjust accordingly
		if (m_resetDates[liEnd] > swapEndTime1[i])
		{
			if (liEnd > 0)
			{
				liEnd--;
			}
		}

	} // for i

	liEnd = firstLiborIndex2;
	for (i = 0; i < nbFixPay2; i++ )
	{	
		liStart = liEnd;
		ljEnd = firstLiborIndex2;

		for (j = 0; j < nbFixPay2; j++ )
		{
			li = liStart;
			ljStart = ljEnd;

			covLower = 0.;
			while ( (li <= lastLiborIndex2 ) && 
				(m_resetDates[li] < swapEndTime2[i]) )
			{
				weight_i = (Maths::min(m_payDates[li],swapEndTime2[i])
					- Maths::max(m_resetDates[li],swapBeginTime2[i]))
					/ swapStep2[i];
				lj = ljStart;

				while ( (lj <= lastLiborIndex2 ) && 
					(m_resetDates[lj] < swapEndTime2[j]) )
				{
					weight_j = (Maths::min(m_payDates[lj],swapEndTime2[j])
						- Maths::max(m_resetDates[lj],swapBeginTime2[j])) / swapStep2[j];	

					for (factor = 0; factor < m_nbFactors; factor++)
					{
						covLower += weight_i * weight_j *
							((m_facLdArray[li][3*factor]   * m_facLdArray[lj][3*factor] +
							m_facLdArray[li][3*factor+1] * m_facLdArray[lj][3*factor+1] +
							m_facLdArray[li][3*factor+2] * m_facLdArray[lj][3*factor+2]) );

					} // for factor

					lj++;
				}
				li++;
			}

			ljEnd = lj;
			if (m_resetDates[ljEnd] > swapEndTime2[j])
			{
				if (ljEnd > 0)
				{
					ljEnd--;
				}
			}

			var2 += weights2[i] * weights2[j] * covLower;

			// the weights used here refer to the expression w_{i}(0)L_{i}(0)/Sr(0)
			// in the analytic Rebonato formula for the swaption Black volatility
		} // for j

		liEnd = li;
		// check if there will be a stub in next round 
		// and adjust accordingly
		if (m_resetDates[liEnd] > swapEndTime2[i])
		{
			if (liEnd > 0)
			{
				liEnd--;
			}
		}

	} // for i

	weights1.clear();
	weights2.clear();

	sqrtVar1 = sqrt(var1);
	sqrtVar2 = sqrt(var2);


	if ( ( sqrtVar1 > TOO_SMALL ) && ( sqrtVar2 > TOO_SMALL ) )
	{
		correlation /= ( sqrtVar1 * sqrtVar2 );
	}
	else
	{
		throw ModelException(method, "Error in calibration : swap rate variance too small");
	}


	return correlation;
}

//****************************************************************************//
//  end : getSwapRateCorr                                                     //
//****************************************************************************//

/** getTotalError */
double SRMRatesLiborDiffuse::getTotalError(void)
{
	static const string method = "SRMRatesLiborDiffuse::getTotalError";

	long 
		i = -1,
		swaptionIndex = -1,
		nbCalibSwaptions = m_calibIndices.size();

	double 
		totalError = 0.,
		currentError = 0.,
		marketBlackVol = 0.,
		modelBlackVol = 0.;

	SRM_VERIFY(nbCalibSwaptions >= 1,
		"Error : no swaptions used for calibration",
		method);

	for(i = 0; i < nbCalibSwaptions; i++)
	{
		swaptionIndex = m_calibIndices[i]; 
		marketBlackVol = m_marketCalibVol[i];

		modelBlackVol = getModelBlackVolatility(*(m_gridSwap[swaptionIndex]));
		currentError = 
			(marketBlackVol - modelBlackVol) * 
			(marketBlackVol - modelBlackVol);

		totalError += currentError / nbCalibSwaptions;
	}

	// add the constraints on swap rate correlations

	return (totalError);
}

//*****************************************************************//
//  end : getTotalError                                            //
//*****************************************************************//

/** getCorError */
double SRMRatesLiborDiffuse::getCorError(void)
{
	size_t i;

	size_t 
		nbCorr = m_targetCorr.size();

	double 
		corrStart = -1, 
		corrEnd = -1,
		modelCorr = -1.0,
		result = 0.0;

	SRMSwapClassSP swap1, swap2;

	for (i = 0; i < nbCorr; i++)
	{
		swap1 = m_gridSwap[ m_corrSwapIndices[2*i] ];
		swap2 = m_gridSwap[ m_corrSwapIndices[2*i+1] ];

		corrStart = 0L;
		corrEnd   = today.yearFrac(swap1->getStartDate());
		corrEnd   = Maths::min(corrEnd, 
			today.yearFrac(swap2->getStartDate()) );

		modelCorr = getSwapRateCorr(
			*(swap1.get()),
			*(swap2.get()),
			corrStart,
			corrEnd);

		result += 
			(m_targetCorr[i] - modelCorr) * 
			(m_targetCorr[i] - modelCorr);
	}

	if (nbCorr > 0.)
	{
		return result/nbCorr;
	}
	else
	{
		return 0.;
	}
}

//*****************************************************************//
//  end : getCorError                                              //
//*****************************************************************//

DRLIB_END_NAMESPACE



