//----------------------------------------------------------------------------
//
//   Group       : Credit Hybrids QR
//
//   Filename    : sCIDaffineProcesses.hpp
//
//   Description : CIR, CIRcst, JumpProcesses
//
//   Date        : Nov 2005
//
//
//----------------------------------------------------------------------------

#ifndef SCID_AFFINE_PROCESSES_HPP
#define SCID_AFFINE_PROCESSES_HPP

#include "edginc/Array.hpp"
#include "edginc/Maths.hpp"
#include "edginc/DoubleMatrix.hpp"
#include "edginc/ran2.h"
#include "edginc/LinearInterpolator.hpp"
#include "edginc/LevyProcesses.h"
#include <map>


DRLIB_BEGIN_NAMESPACE




class MARKET_DLL AffineJump
{
public:
	AffineJump(){};
	AffineJump(double frequency, array<double> &SizeOfJumps, array<double> &Impact, double Decay){ setParameters(frequency, SizeOfJumps, Impact, Decay);};
	virtual ~AffineJump(){};

	double freq, decay;
	array<double> jumpSize, impact; 

	AffineJump & operator=(const AffineJump& b)
	{
		if(this == &b) return *this;
		freq = b.freq;
		decay = b.decay;
		impact = b.impact;
		jumpSize = b.jumpSize;
		return *this;
	};

	void setParameters(double frequency, array<double> &SizeOfJumps, array<double> &Impact, double Decay);
	void setParameters(double frequency, double jSize, array<double> &initRate, double Impact, double Decay);

	void survProba(double linear,
				   array<double> &times, 
				   double *SurvProb,     // (O) (name,index time) -> (index time)*m_nbNames+name
   				   double futureTime = 0,
				   double *lambdaFutureTime = 0);	
		// E_t(exp(-\int_t^T linear*\lambda_u du)), where t=futureTime
	    // lambda is equal to lambdaFutureTime at t, and T is times
		// does it for all names, at all values of times and saves it in SurvProb

	void condSurvProba(double linear,
					   array<double> &times, 
					   int nbJumps,
					   Uniform &jumpTimes,   // will provide the time of the jumps
					   double maturity,      // when multiplied by maturity
					   double *SurvProb,         // (O) (name,index time) -> (index time)*m_nbNames+name
	   				   double futureTime = 0,
					   double *lambdaFutureTime = 0);		  
		// same things, but conditional on the values of poisson process
		// gives the values of the number of jump at maturity
		// and simulates inside the time of jump (using jumpTimes)
			
		
	void survProbaMult(
				   double linear,
				   array<double> &times, 
				   double *SurvProb,        // (O) (name,index time) -> (index time)*m_nbNames+name
   				   double futureTime = 0,
				   double *lambdaFutureTime = 0);		  
	double condSurvProbaMult(        // return first Jump ... return 999 if no jump
					double linear,
					array<double> &times, 
				    int nbJumps,
					Uniform &jumpTimes,   // will provide the time of the jumps
  				    double maturity,      // when multiplied by maturity
					double *SurvProb,         // (O) (name,index time) -> (index time)*m_nbNames+name
	   				double futureTime = 0,
					double *lambdaFutureTime = 0);		  
		// the last two functions are as above, but instead of doing SurvProb = something, does SurvProb *= something
		// used a lot in the tranche pricer
		

	void condSpreadAndHazard(
					array<double> &times,
				    int nbJumps,
					Uniform &jumpTimes,
					double maturity,
					Uniform &jumpUniform,    // gives uniform in R^m_nbNames;
					double *intensity,       // (O) (name,index time) -> (index time)*m_nbNames+name
					double *intIntensity);   // (O) (name,index time) -> (index time)*m_nbNames+name
			// used in full MC.. simulates the spreads, conditional on the nb Jumps (simulated somewhere else)
			// simulates the jumpSize, and the jumpTimes
			// output lambda and its integral in intensity and intIntensity


	void testAffineJump(array<double> &times, double *SurvProb, double *SurvProbCond, int nbPaths, long seed, double futureTime, double maturity, double lambdaFutureTime);
	
private:
	//optimize phi
    //inline double phi(double t) { return (1.0 - exp(-decay*t)) / decay; }
	inline double phi(double t);
	inline double getJump(double jumpSize, double impact, double unif) // convert a uniform random variable to the product of a Bernoulli times an exponential
	{
		if (unif>impact) return 0;
		return -jumpSize*log( unif/impact );
	}

	array<double> timeOfJumps;
    std::map<double, double> phiMap;
};



// function which deals with the idiosyncratic part of the spread
// CIR process with time dependent theta!
// store the data, and compute survival proba and conditional one, and help simulate the process

class MARKET_DLL CIR
{
public:
	CIR(){};
	CIR(array<double> &initValue, array<double> &vol, double decay){ setParameters(initValue, vol, decay); };
	CIR(array<double> &initValue, double vol, double decay){setParameters(initValue, vol, decay);};

	double kappa;
	array<double> lambda0, sigma, thetaLimit;
	array<double> theta, thetaDiff;   // size nbNames * timeTheta.size(), nbNames * timeThetaDiff.size()
									 // (name,time)->name*timeTheta+time
	array<double> timeTheta, timeThetaDiff;   // theta is piecewise constant, and the value just afer timeTheta[i] is theta[i]...


	void setParameters(array<double> &initValue, array<double> &vol, double decay)
	{ 
		lambda0 = initValue; 
		sigma = vol; 
		kappa = decay; 
		thetaLimit.resize(lambda0.size());
		for (int m=0; m<thetaLimit.size(); m++)
			thetaLimit[m] = vol[m]*vol[m]/(2*decay);
		alpha.resize(lambda0.size());
		beta.resize(lambda0.size());
	};

	void setParameters(array<double> &initValue, double vol, double decay)
	{ 
		int m;
		lambda0 = initValue; 
		sigma.resize(initValue.size());
		for (m=0; m<sigma.size(); m++) sigma[m] = sqrt(lambda0[m])*vol;
		kappa = decay; 
		thetaLimit.resize(lambda0.size());
		for (m=0; m<thetaLimit.size(); m++)
			thetaLimit[m] = sigma[m]*sigma[m]/(2*decay);
		alpha.resize(lambda0.size());
		beta.resize(lambda0.size());
	};
	inline double thetaShift(int m, int k, double para) { return Maths::max(thetaLimit[m], theta[m*timeTheta.size()+k] + para * lambda0[m]); }
	inline double thetaShiftDiff(int m, int k, double para) { return Maths::max(thetaLimit[m], thetaDiff[m*timeThetaDiff.size()+k] + para * lambda0[m]); }

	CIR & operator=(const CIR& b)
	{
		if(this == &b) return *this;
		lambda0 = b.lambda0;
		sigma = b.sigma;
		kappa = b.kappa;
		theta = b.theta;
		timeTheta = b.timeTheta;
		thetaDiff = b.thetaDiff;
		timeThetaDiff = b.timeThetaDiff;
		return *this;
	};

	void survProbaMult(         // unconditional survival probability, possibly starting in the future.. 
								// multiply this values and put it into SurvProb. 
								// does this for all names
              double linear, 
			  double para,
			  array<double> &timeSurvProb,
			  double *SurvProb,  // name, indexTime -> name*timeSurvProb.size() + indexTime
			  double futureTime = 0,
			  double *lambdaFutureTime = 0);		  

	void survProba(				// same, but store the value in SurvProb, rather than multiply it.
              double linear, 
			  double para,
			  array<double> &timeSurvProb,
			  double *SurvProb,    // name, indexTime -> name*timeSurvProb.size() + indexTime
			  double futureTime = 0,
			  double *lambdaFutureTime = 0);		  
	
	void TimeDepTheta(			// calibrate to single names: given SurvProb, find theta which correspond to it
				  array<double> &time,
				  double *SurvProb);		// name, indexTime -> name*time.size() + indexTime // gives the SurvProb at time

//	void changeTheta(array<double> &newTimeTheta);
	void setThetaDiff(array<double> &newTimeTheta);  // change the time associated to theta, 
													// to adapt the time of simulation of full MC
													// save it in ???diff (timeThetaDiff, thetaDiff)

	void diffuse(double linear, 
				 double para,
				 double *BrownianIncrement, // name, index -> name*(timeTheta.size()-1)+index
				 double *CIRprocess);       // name, index -> name*(timeTheta.size())+index
													// given the Brownian increments, 
													// gives the values of the CIR processes
	void diffuseNoVol(double linear,
					  double para,
					  double *CIRprocess);         // name, index -> name*(timeTheta.size())+index
													// no vol, so it is an ODE instead of SDE..

/// USED ONLY FOR THE TREE DEFINITION OF SCID

	void survProbaMult(
	          double *f_atTimeTheta,
			  array<double> &timeSurvProb,
			  double *f_atTimeSurvProb,
			  double *SurvProb,
			  double h, // time step in the Runge Kutta method of order 4
			  double futureTime = 0,
			  double *lambdaFutureTime = 0);

	void diffuse(double *BrownianIncrement, // name, index -> name*(timeTheta.size()-1)+index
				 double *CIRprocess);       // name, index -> name*(timeTheta.size())+index
											// given the Brownian increments, 
											// gives the values of the CIR processes
	void diffuseNoVol(double *CIRprocess);         // name, index -> name*(timeTheta.size())+index
													// no vol, so it is an ODE instead of SDE..
	
	void testClosedForm(array<double> &ftimeTheta, double h, 
						 double futureTime, array<double> &lambdaFutureTime,
						 array<double> &timeSurvProb, double **SurvProbOut);
	void testSimulation(array<double> &ftimeTheta,
			   long seed, int nbPaths, double timeSteps,
			   double futureTime, array<double> &lambdaFutureTime,
			   array<double> &timeSurvProb, double **SurvProbOut);
private:
	array<double> alpha, beta; // to avoid reallocating memory all the time
	array<double> expDecay;
};




// model the common factor.. only one dimensional, and theta is not time dependent
// so it is simplier, except we need to be able to diffuse from the future
// which is not so hard
// similar structure as above

class MARKET_DLL CIRcst
{
public:
	CIRcst(){};
	CIRcst(double initValue, double vol, double decay){setParameters(initValue, vol, decay);};

	double lambda0, sigma, kappa, theta, thetaLimit;

	void setParameters(double initValue, double vol, double decay);

	CIRcst & operator=(const CIRcst& b)
	{
		if(this == &b) return *this;
		lambda0 = b.lambda0;
		sigma = b.sigma;
		kappa = b.kappa;
		theta = b.theta;
		thetaLimit = b.thetaLimit;
		return *this;
	};

	inline double thetaShift(double para) { return Maths::max(thetaLimit, theta + para*lambda0); };

	void survProbaMult(
			  double ams,
              double linear, 
			  double para,
			  array<double> &timeSurvProb,
			  double *SurvProb,
			  double futureTime = 0,
			  double lambdaFutureTime = 0);		  

	void survProba(
			  double ams,
              double linear, 
			  double para,
			  array<double> &timeSurvProb,
			  double *SurvProb,
			  double futureTime = 0,
			  double lambdaFutureTime = 0);		  



	void FindTheta(double time, double SurvProb);		// find theta such that survProba(time) = SurvProb
	void setTheta(double th){ theta=th; };

	void diffuse(double linear, 
				 double para,
				 double *BrownianIncrement, // of size timesBM.size()-1
				 double *expDecay,          // of size timesBM.size() -1 = exp(-kappa(t_{k+1}-t_k))
		 		 array<double> &timesBM, 
				 double *CIRprocess,
				 double futureTime = 0,
				 double lambdaFutureTime = 0);		 // of size timesBM.size()
	void condSurvProba(
				 double linear, 
				 double para,
				 double *BrownianIncrement, // of size timesBM.size() -1 
				 double *expDecay,          // of size timesBM.size() -1 = exp(-kappa(t_{k+1}-t_k))
	             array<double> &timesBM,      
				 double *CIRprocess,        // of size timesBM.size()
				 array<double> &timeSurvProb,
				 double *SurvProb,
				 double futureTime = 0,
				 double lambdaFutureTime = 0);

/// USED ONLY IN THE TREE DEFINITION OF 
	void survProba(
		  array<double> &t, 
		  double *ft,
		  array<double> &timeSurvProb,
		  array<double> &ftimeSurvProb,
		  double *SurvProb,
		  double futureTime = 0,
		  double lambdaFutureTime = 0,
		  double h = 0.5);

	void survProbaMult(
		  array<double> &t, 
		  double *ft,
		  array<double> &timeSurvProb,
		  array<double> &ftimeSurvProb,
		  double *SurvProb,
		  double futureTime = 0,
		  double lambdaFutureTime = 0,
		  double h = 0.5);

	void diffuse(double *BrownianIncrement, // of size timesBM.size()-1
				 double *expDecay,          // of size timesBM.size() -1 = exp(-kappa(t_{k+1}-t_k))
		 		 array<double> &timesBM, 
				 double *CIRprocess,
				 double futureTime = 0,
				 double lambdaFutureTime = 0);		 // of size timesBM.size()
	void condSurvProba(
				 array<double> &f_at_t_k,
				 double *BrownianIncrement, // of size timesBM.size() -1 
				 double *expDecay,          // of size timesBM.size() -1 = exp(-kappa(t_{k+1}-t_k))
	             array<double> &timesBM,      
				 double *CIRprocess,        // of size timesBM.size()
				 array<double> &timeSurvProb,
				 double *SurvProb,
				 double futureTime = 0,
				 double lambdaFutureTime = 0);

	void testSurvProb(double linear, double para, long seed, double timeSteps, int nbPaths, double futureTime, double lambdaFutureTime, array<double> &timeSurvProb, double *SurvProb);
	void testSurvProb2(array<double> &t, array<double> &ft, 
				 long seed, double timeSteps, int nbPaths, 
				 double futureTime, double lambdaFutureTime, 
				 array<double> &timeSurvProb, double *SurvProb);
};


inline double alphaCIR(double kappa, double theta, double sigma, double t)
{
	double delta = sqrt(kappa*kappa + 2.0 * sigma*sigma);
	if( sigma != 0) 
	{
		double expr = (delta - kappa + (delta + kappa) * exp(delta*t)) / (2 * delta);
		return kappa * theta / (sigma*sigma) * ( (delta + kappa) * t - 2 * log(expr) );
	} 
	else return theta * (1/kappa * (1 - exp(-kappa*t)) - t);
}

inline double betaCIR(double kappa, double sigma, double t) 
{
	double delta = sqrt(kappa*kappa + 2.0 * sigma*sigma);
	if(sigma != 0) return 2 * (1 - exp(delta * t)) / (delta - kappa + (delta + kappa) * exp(delta*t) );
		else return -1/kappa * (1-exp(-kappa*t));
}



inline void piecewiseLinearFunction(array<double> &t, double *f, double *wantedTimes, int wantedTimesSize, double *values)
{
	int k=0;
	int i=0;
	while ((i<wantedTimesSize)&&(wantedTimes[i]<=t[0]))
	{
		values[i] = f[0];
		i++;
	}
	for (; i<wantedTimesSize; i++)
	{
		while ((k+1<t.size())&&(t[k+1]<wantedTimes[i])) k++;
		if (k==t.size()-1)
		{
			for (; i<wantedTimesSize; i++)
				values[i] = f[k];
		}
		else
			values[i]=( f[k]*(t[k+1]-wantedTimes[i])+f[k+1]*(wantedTimes[i]-t[k]) ) / (t[k+1]-t[k]);
	}
}


// solve 
// d(alpha(t)) = d(t) + kappa * alpha(t) - 0.5*sigma^2*alpha(t)^2
// d(beta(t)) = - kappa*theta*alpha(t)
// with d linear between T1 and T2, with values dT1 and dT2 at T1 and T2,
// with alpha and beta known at T2 (with values alpha, beta)
// and we find alpha and beta at T1 (stored in alpha, beta)
// we use a Runge-Kutta method of order 4, with stepsize h
// first function is just the function in the ODE

inline void CIRlinear(double &kappa, double &theta, double &sigma, double &a, double &b, double t, double beta,
				 double &fa, double &fb) // (O)
{
	fb = a+b*t+kappa*beta - sigma*sigma*.5*beta*beta;
	fa = -kappa*theta*beta;
}


inline void RicattiCIRlinear(double kappa, double theta, double sigma,
							 double T2, double dT2,
							 double T1, double dT1,
							 double &alpha, double &beta, // (I/O)
							 double h)
{
	if (T1>=T2) return; 

	double k1a,k2a,k3a,k4a, k1b,k2b,k3b,k4b,hh;
	double tn=T2;
	
	double a = (dT1*T2-dT2*T1)/(T2-T1);
	double b = (dT2-dT1)/(T2-T1);

	while (tn>T1+1e-15)
	{
		hh = min(tn-T1,h);

		CIRlinear(kappa,theta,sigma,a,b,tn,	       beta,              k1a, k1b);
		CIRlinear(kappa,theta,sigma,a,b,tn-0.5*hh, beta - 0.5*hh*k1b, k2a, k2b);
		CIRlinear(kappa,theta,sigma,a,b,tn-0.5*hh, beta - 0.5*hh*k2b, k3a, k3b);
		CIRlinear(kappa,theta,sigma,a,b,tn-hh,	   beta - hh*k3b,	  k4a, k4b);

		alpha -= 0.1666666666666*hh*(k1a+2*k2a+2*k3a+k4a);
		beta  -= 0.1666666666666*hh*(k1b+2*k2b+2*k3b+k4b);
		tn -= hh;
	}
}

// EULER METHOD... JUST FOR TESTING
inline void RicattiCIRlinear2(double kappa, double theta, double sigma,
							 double T2, double dT2,
							 double T1, double dT1,
							 double &alpha, double &beta, // (I/O)
							 double h)
{
	if (T1>=T2) return; 

	double k1a,k1b,hh;
	double tn=T2;
	
	double a = (dT1*T2-dT2*T1)/(T2-T1);
	double b = (dT2-dT1)/(T2-T1);
	while (tn>T1+1e-15)
	{
		hh = min(tn-T1,h);
		CIRlinear(kappa,theta,sigma,a,b,tn,	beta, k1a, k1b);
		alpha -= hh*k1a;
		beta  -= hh*k1b;
		tn -= hh;
	}
}



DRLIB_END_NAMESPACE
#endif
