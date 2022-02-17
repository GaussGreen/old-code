#include "edginc/config.hpp"
#include "edginc/SCIDaffineProcesses.hpp"
#include "edginc/Algorithm.hpp"

DRLIB_BEGIN_NAMESPACE


inline double integrate(double t1, 
						double t2, 
						int &k,    // time[j]<=t1<time[j+1] for j>=k .. modified such that time[k]<=t2
						array<double> &time,
						double *values)  // values piecewise linear function, equal to values[k] at time[k]// values of size time.size()
{
	// we first ensure that we have // time[k]<=t1<time[k+1]
	while ( (k+1<time.size())&&(time[k+1]<=t1) ) k++;
	if (k+1>=time.size())  // max time<t1 .. then we extrapolate values to a constant function
		return (t2-t1)*values[k];
	else
	{
		double t=min(t2,time[k+1]);
		double sum = (t-t1)*values[k]+(values[k+1]-values[k])/(time[k+1]-time[k])*0.5*(Maths::square(t-time[k])-Maths::square(t1-time[k]));
		k++;
		while ((k+1<time.size())&&(time[k+1]<t2))
		{
			sum += 0.5*(time[k+1]-time[k])*(values[k+1]+values[k]);
			k++;
		}
		if (k+1<time.size())  // time[k]<t1 .. then we extrapolate values to a constant function
			return sum+(t2-time[k])*values[k]+(values[k+1]-values[k])*0.5*Maths::square(t2-time[k])/(time[k+1]-time[k]);
		else
			return sum+(t2-time[k])*values[k];
	}
}

inline double integrate(double t1, 
						double t2, 
						int &k,    // time[j]<=t1<time[j+1] for j>=k .. modified such that time[k]<=t2
						array<double> &time,
						double *f, double *g)  // values piecewise linear function, equal to values[k] at time[k]// values of size time.size()
{
	// we first ensure that we have // time[k]<=t1<time[k+1]
	while ( (k+1<time.size())&&(time[k+1]<=t1) ) k++;
	if (k+1>=time.size())  // max time<t1 .. then we extrapolate values to a constant function
		return (t2-t1)*f[k]*g[k];
	else
	{
		double t=min(t2,time[k+1]);
		double sum = (t-t1)*f[k]*g[k]+(f[k+1]*g[k+1]-f[k]*g[k])/(time[k+1]-time[k])*0.5*(Maths::square(t-time[k])-Maths::square(t1-time[k]));
		k++;
		while ((k+1<time.size())&&(time[k+1]<t2))
		{
			sum += 0.5*(time[k+1]-time[k])*(f[k+1]*g[k+1]+f[k]*g[k]);
			k++;
		}
		if (k+1<time.size())  // time[k]<t1 .. then we extrapolate values to a constant function
			return sum+(t2-time[k])*f[k]*g[k]+(f[k+1]*g[k+1]-f[k]*g[k])*0.5*Maths::square(t2-time[k])/(time[k+1]-time[k]);
		else
			return sum+(t2-time[k])*f[k]*g[k];
	}
}










void CIR::survProba(
	          double linear, 
		      double para, 
			  array<double> &timeSurvProb,
			  double *SurvProb,
			  double futureTime,
			  double *lambdaFutureTime)		  
{
	for (int k=0; k<lambda0.size()*timeSurvProb.size(); k++) SurvProb[k]=1;
	survProbaMult(linear,para,timeSurvProb,SurvProb,futureTime,lambdaFutureTime);
};

void CIR::survProbaMult(
	          double linear, 
		      double para, 
			  array<double> &timeSurvProb,
			  double *SurvProb,
			  double futureTime,
			  double *lambdaFutureTime)		  
{
	double t1,t2, sqrtLinear = sqrt(linear);
	int i, j, m, j0=0;
	
	if (futureTime==0) 
	{
		while ( (j0<timeTheta.size()-1) && (timeTheta[j0]<futureTime + 1e-10) ) j0++;
		j0 = Maths::max(0,j0-1);
		for (i=0; i<timeSurvProb.size(); i++)
		{
			if (timeSurvProb[i]>futureTime) 
			{
				for (m=0; m<alpha.size(); m++) alpha[m]=0;
				j=j0;
				while ( (j+1<timeTheta.size()-1) && (timeTheta[j+1]<timeSurvProb[i]) )
				{
					t2 = Maths::max(0.0,timeSurvProb[i]-timeTheta[j]);
					t1 = Maths::max(0.0,timeSurvProb[i]-timeTheta[j+1]);					
					for (m=0; m<lambda0.size(); m++)
						alpha[m] += linear*thetaShift(m,j,para)*(alphaCIR(kappa,1,sqrtLinear*sigma[m], t2) - alphaCIR(kappa,1,sqrtLinear * sigma[m], t1));
					j++;
				}
				if (timeSurvProb[i]>timeTheta[j])
				{
					t2 = Maths::max(0.0,timeSurvProb[i]-timeTheta[j]);
					t1 = 0;
					for (m=0; m<lambda0.size(); m++)
						alpha[m] += linear*thetaShift(m,j,para)*(alphaCIR(kappa,1,sqrtLinear*sigma[m], t2) - alphaCIR(kappa,1,sqrtLinear*sigma[m], t1));
				}
				for (m=0; m<lambda0.size(); m++)
					SurvProb[m*timeSurvProb.size()+i] *= exp( alpha[m] + linear*lambda0[m]*(1+para)*betaCIR(kappa, sqrtLinear*sigma[m], timeSurvProb[i] - futureTime) );
			}
			// else multiply SurvProb[i] by 1, i.e. do nothing
		}
	}
	else
	{
		while ( (j0<timeTheta.size()-1) && (timeTheta[j0]<futureTime + 1e-10) ) j0++;
		j0 = Maths::max(0,j0-1);
		for (i=0; i<timeSurvProb.size(); i++)
		{
			if (timeSurvProb[i]>futureTime) 
			{
				for (m=0; m<alpha.size(); m++) alpha[m]=0;
				j=j0;
				while ( (j+1<timeTheta.size()-1) && (timeTheta[j+1]<timeSurvProb[i]) )
				{
					t2 = Maths::max(0.0,timeSurvProb[i]-Maths::max(futureTime,timeTheta[j]));
					t1 = Maths::max(0.0,timeSurvProb[i]-Maths::max(futureTime,timeTheta[j+1]));
					for (m=0; m<lambda0.size(); m++)
						alpha[m] += linear*thetaShift(m,j,para)*(alphaCIR(kappa,1,sqrtLinear*sigma[m], t2) - alphaCIR(kappa,1,sqrtLinear * sigma[m], t1));
					j++;
				}
				if (timeSurvProb[i]>timeTheta[j])
				{
					t2 = Maths::max(0.0,timeSurvProb[i]-Maths::max(futureTime,timeTheta[j]));
					t1 = Maths::max(0.0,timeSurvProb[i]-Maths::max(futureTime,timeSurvProb[i]));
					for (m=0; m<lambda0.size(); m++)
						alpha[m] += linear*thetaShift(m,j,para)*(alphaCIR(kappa,1,sqrtLinear*sigma[m], t2) - alphaCIR(kappa,1,sqrtLinear*sigma[m], t1));
				}
				for (m=0; m<lambda0.size(); m++)
					SurvProb[m*timeSurvProb.size()+i] *= exp( alpha[m] + lambdaFutureTime[m]*betaCIR(kappa, sqrtLinear*sigma[m], timeSurvProb[i] - futureTime) );
			}
			// else multiply SurvProb[i] by 1, i.e. do nothing
		}
	}
};


void CIR::survProbaMult(
	          double *f_atTimeTheta,
			  array<double> &timeSurvProb,
			  double *ftimeSurvProb,
			  double *SurvProb,
			  double h,
			  double futureTime,
			  double *lambdaFutureTime)
{
	double ffutureTime;
	piecewiseLinearFunction(timeTheta,&f_atTimeTheta[0], &futureTime, 1, &ffutureTime);

	
	int high = timeTheta.size()-1;
	int high2, m;
	for (int i=timeSurvProb.size()-1; i>=0; i--)
	{
		if (timeSurvProb[i] <= futureTime) return; // no change to SurvProb[j], j<=i
		for (m=0; m<alpha.size(); m++) { alpha[m]=0; beta[m]=0; };
		while ((high>=0)&&(timeTheta[high]>=timeSurvProb[i])) high--;
		// now high is the greatest integer such that timeSurvProb[i]>timeTheta[high], unless maybe high==-1
		high2=high;
		if (f_atTimeTheta[high2]<futureTime)
			for (m=0; m<alpha.size(); m++)
				RicattiCIRlinear(kappa, theta[m*timeTheta.size()+high2], sigma[m], timeSurvProb[i], ftimeSurvProb[i], futureTime, ffutureTime, alpha[m],beta[m],h);
		else  
		{
			for (m=0; m<alpha.size(); m++)
				RicattiCIRlinear(kappa, theta[m*timeTheta.size()+high2], sigma[m], timeSurvProb[i], ftimeSurvProb[i], timeTheta[high2], f_atTimeTheta[high2], alpha[m],beta[m],h);
			while ((high2>0)&&(timeTheta[high2-1]>futureTime))
			{
				for (m=0; m<alpha.size(); m++)
					RicattiCIRlinear(kappa, theta[m*timeTheta.size()+high2-1], sigma[m], timeTheta[high2], f_atTimeTheta[high2], timeTheta[high2-1], f_atTimeTheta[high2-1], alpha[m], beta[m], h);
				high2--;
			}
			for (m=0; m<alpha.size(); m++)
				RicattiCIRlinear(kappa, theta[m*timeTheta.size()+high2-1], sigma[m], timeTheta[high2], f_atTimeTheta[high2], futureTime, ffutureTime, alpha[m],beta[m],h);
		}
		if (futureTime==0)
		{
			for (m=0; m<alpha.size(); m++)
				SurvProb[m*timeSurvProb.size()+i] *= exp(alpha[m] + lambda0[m]*beta[m]);
		}
		else
		{
			for (m=0; m<alpha.size(); m++)
				SurvProb[m*timeSurvProb.size()+i] *= exp(alpha[m] + lambdaFutureTime[m]*beta[m]);
		}
	}
};


void CIR::TimeDepTheta(array<double> &time,	double *SurvProb)
{
	double sum;
	theta.resize(lambda0.size()*(time.size()-1));
	timeTheta.resize(time.size()-1);
	int i,j,m;
	for (i=0; i< time.size()-1; i++) timeTheta[i]=time[i];

	
	for (m=0; m<lambda0.size(); m++)
		for (i=0; i< time.size()-1; i++)
		{
			sum = 0;
			for (j=0; j<i; j++)
				sum += theta[m*timeTheta.size()+j]*(alphaCIR(kappa,1,sigma[m],time[i+1]-time[j]) - alphaCIR(kappa,1,sigma[m],time[i+1]-time[j+1]));


			theta[m*timeTheta.size()+i] = ( log(SurvProb[m*time.size()+i+1])-lambda0[m]*betaCIR(kappa,sigma[m],time[i+1]) - sum ) / alphaCIR(kappa, 1, sigma[m], time[i+1]-time[i]);
		
			if (theta[m*timeTheta.size()+i] < thetaLimit[m]) 
			{
				theta[m*timeTheta.size()+i] = thetaLimit[m];
				sum += theta[m*timeTheta.size()+i]*alphaCIR(kappa, 1, sigma[m], time[i+1]-time[i]);
				SurvProb[m*time.size()+i+1]=exp( sum + lambda0[m]*betaCIR(kappa, sigma[m], time[i+1]) );
			}
		}
};

void CIR::setThetaDiff(array<double> &newTimeTheta)
{
	timeThetaDiff.resize(newTimeTheta.size()-1);
	thetaDiff.resize(lambda0.size()*timeThetaDiff.size());  // theta is piecewise constant, and the value just afer timeTheta[i] is theta[i]...
	int k,l,m;
	for (k=0; k<thetaDiff.size(); k++) thetaDiff[k]=0;
	expDecay.resize(timeThetaDiff.size());

	double delta;
	for (k=0; k<newTimeTheta.size()-1; k++)
		for (l=0; l<timeTheta.size(); l++)
		{
			if (l+1==timeTheta.size()) delta = newTimeTheta[k+1] - Maths::max(timeTheta[l],newTimeTheta[k]);
				else delta = Maths::min(timeTheta[l+1],newTimeTheta[k+1])-Maths::max(timeTheta[l],newTimeTheta[k]);
			if (delta>0) for (m=0; m<lambda0.size(); m++) thetaDiff[m*timeThetaDiff.size()+k] += theta[m*timeTheta.size()+l]*delta;
		}
	for (k=0; k<newTimeTheta.size()-1; k++)
	{
		delta = 1/(newTimeTheta[k+1]-newTimeTheta[k]);
		for (m=0; m<lambda0.size(); m++) thetaDiff[m*timeThetaDiff.size()+k] *= delta;
		timeThetaDiff[k] = newTimeTheta[k];
		expDecay[k] = exp(-0.5*kappa*(newTimeTheta[k+1]-newTimeTheta[k]));
	}
}


void CIR::diffuse(double linear, 
				  double para,
				  double *BrownianIncrement, // name, k -> name*(timeThetaDiff.size())+k
				  double *CIRprocess)        // name, k -> name*(timeThetaDiff.size()+1)+k
{
	int k,m;
	for (m=0; m<lambda0.size(); m++) CIRprocess[m*(1+timeThetaDiff.size())] = lambda0[m]*(1+para);
	double thetaPrimeDecayed;
	for (k=0; k<timeThetaDiff.size(); k++)
	{
		for (m=0; m<lambda0.size(); m++) 
		{
			thetaPrimeDecayed = (thetaShiftDiff(m,k,para) - 0.25*sigma[m]*sigma[m]/kappa)*(1-expDecay[k]);
			CIRprocess[m*(1+timeThetaDiff.size())+k+1] = Maths::max(0.0, thetaPrimeDecayed + Maths::square( sqrt(Maths::max(0.0,thetaPrimeDecayed + CIRprocess[m*(1+timeThetaDiff.size())+k] * expDecay[k])) + 0.5 * sigma[m] * BrownianIncrement[m*timeThetaDiff.size()+k]) * expDecay[k]);
		}
	}
	for (k=0; k<timeThetaDiff.size()+1; k++)
		for (m=0; m<lambda0.size(); m++) 
			CIRprocess[m*(1+timeThetaDiff.size())+k] *= linear;
}

void CIR::diffuseNoVol(double linear,
					   double para,
					   double *CIRprocess)
{
	int k,m;
	for (m=0; m<lambda0.size(); m++) CIRprocess[m*(1+timeThetaDiff.size())] = lambda0[m]*(1+para);
	double decayFact;
	for (k=0; k<timeThetaDiff.size(); k++)
	{
		decayFact = expDecay[k]*expDecay[k];
		for (m=0; m<lambda0.size(); m++) 
			CIRprocess[m*(1+timeThetaDiff.size())+k+1] = thetaShiftDiff(m,k,para)*(1-decayFact) + CIRprocess[m*(1+timeThetaDiff.size())+k]*decayFact;
	}
	for (k=0; k<timeThetaDiff.size()+1; k++)
		for (m=0; m<lambda0.size(); m++) 
			CIRprocess[m*(1+timeThetaDiff.size())+k] *= linear;
}

void CIR::diffuse(double *BrownianIncrement, // name, k -> name*(timeThetaDiff.size())+k
				  double *CIRprocess)        // name, k -> name*(timeThetaDiff.size()+1)+k
{
	int k,m;
	for (m=0; m<lambda0.size(); m++) CIRprocess[m*(1+timeThetaDiff.size())] = lambda0[m];
	double thetaPrimeDecayed;
	for (k=0; k<timeThetaDiff.size(); k++)
		for (m=0; m<lambda0.size(); m++) 
		{
			thetaPrimeDecayed = (thetaDiff[m*timeThetaDiff.size()+k] - 0.25*sigma[m]*sigma[m]/kappa)*(1-expDecay[k]);
			CIRprocess[m*(1+timeThetaDiff.size())+k+1] = Maths::max(0.0, thetaPrimeDecayed + Maths::square( sqrt(Maths::max(0.0,thetaPrimeDecayed + CIRprocess[m*(1+timeThetaDiff.size())+k] * expDecay[k])) + 0.5 * sigma[m] * BrownianIncrement[m*timeThetaDiff.size()+k]) * expDecay[k]);
		}
}

void CIR::diffuseNoVol(double *CIRprocess)
{
	int k,m;
	for (m=0; m<lambda0.size(); m++) CIRprocess[m*(1+timeThetaDiff.size())] = lambda0[m];
	double decayFact;
	for (k=0; k<timeThetaDiff.size(); k++)
	{
		decayFact = expDecay[k]*expDecay[k];
		for (m=0; m<lambda0.size(); m++) 
			CIRprocess[m*(1+timeThetaDiff.size())+k+1] = thetaDiff[m*timeThetaDiff.size()+k]*(1-decayFact) + CIRprocess[m*(1+timeThetaDiff.size())+k]*decayFact;
	}
}


void CIR::testClosedForm(array<double> &ftimeTheta, double h, 
						 double futureTime, array<double> &lambdaFutureTime,
						 array<double> &timeSurvProb, double **SurvProbOut)
{
/*	array<double> ftimeSurvProb(timeSurvProb.size());
	piecewiseLinearFunction(timeTheta,&ftimeTheta[0], &timeSurvProb[0], timeSurvProb.size(), &ftimeSurvProb[0]);
	SurvProbOut.resize(timeSurvProb.size(),lambdaFutureTime.size(), 1.0);
	array<double> SurvProbOut2(timeSurvProb.size()*lambdaFutureTime.size(),1.0);
	survProbaMult(
	          &ftimeTheta[0],
			  timeSurvProb,
			  &ftimeSurvProb[0],
			  &SurvProbOut2[0],
			  h,
			  futureTime,
			  &lambdaFutureTime[0]);
	for (int i=0; i<timeSurvProb.size(); i++)
		for (int m=0; m<lambdaFutureTime.size(); m++)
			SurvProbOut[i][m] = SurvProbOut2[m*timeSurvProb.size()+i];

*/
}

void CIR::testSimulation(array<double> &ftimeTheta,
			   long seed, int nbPaths, double timeSteps,
			   double futureTime, array<double> &lambdaFutureTime,
			   array<double> &timeSurvProb, double **SurvProbOut)
{
/*	int k,l;
	SurvProbOut.resize(timeSurvProb.size(),lambdaFutureTime.size());
	SurvProbOut.fill(0.0);
	Matrix<double> SurvProbTemp(timeSurvProb.size(),lambdaFutureTime.size(), 0.0);
	double lastTime = timeSurvProb[timeSurvProb.size()-1];
	int nbSteps = floor(lastTime/timeSteps)+1;

	array<double> timeSimulation(nbSteps+1); 
	for (k=0; k<=nbSteps; k++) timeSimulation[k] = (k*lastTime)/nbSteps;

	setThetaDiff(timeSimulation);

	array<double> f_at_t(timeSimulation.size());
	piecewiseLinearFunction(timeTheta,&ftimeTheta[0], &timeSimulation[0], timeSimulation.size(), &f_at_t[0]);

	BrownianMotion BM(timeSimulation, seed);
	BM.setIpath(0);
	array<double> m_BMincrements(lambdaFutureTime.size()*timeThetaDiff.size());
	array<double> CIRprocess(lambdaFutureTime.size()*timeSimulation.size(), 0.0);

	array<double> expDecay(timeSimulation.size()-1);
	for (k=0; k<expDecay.size(); k++) expDecay[k]=exp(-0.5*kappa*(timeSimulation[k+1]-timeSimulation[k]));

	for (int mc=0; mc<nbPaths; mc++)
	{
		for (int m=0; m<lambdaFutureTime.size(); m++)
			BM.GetBMincrements(&m_BMincrements[m*timeThetaDiff.size()]);
		diffuse(&m_BMincrements[0],&CIRprocess[0]);
	
		SurvProbTemp.fill(0);

		for (m=0; m<lambdaFutureTime.size(); m++)
		{
			for (k=0; (timeSurvProb[k]<1e-10)&&(k<timeSurvProb.size()); k++) SurvProbTemp[k][m] = 0;
			if (k<timeSurvProb.size())
			{
				int interval = 0;
				SurvProbTemp[k][m] = integrate(0,timeSurvProb[k], interval, timeSimulation, &CIRprocess[m*(timeThetaDiff.size())],&f_at_t[0]);
				for (; k<timeSurvProb.size()-1; k++) SurvProbTemp[k+1][m] = SurvProbTemp[k][m] + integrate(timeSurvProb[k], timeSurvProb[k+1], interval, timeSimulation, &CIRprocess[m*(timeThetaDiff.size())], &f_at_t[0]);
			}
			for (k=0; k<timeSurvProb.size(); k++) SurvProbOut[k][m] += exp(-SurvProbTemp[k][m]);
		}
	}
	for (l=0; l<timeSurvProb.size(); l++) 
		for (int m=0; m<lambdaFutureTime.size(); m++)
			SurvProbOut[l][m] /= double(nbPaths);
	l=l;*/
}



//////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////




void CIRcst::setParameters(double initValue, double vol, double decay)
{ 
	lambda0 = initValue; 
	sigma = vol; 
	kappa = decay; 
	thetaLimit = vol*vol/(2*decay);
};


void CIRcst::survProba(
 			  double ams,
              double linear, 
			  double para,
			  array<double> &timeSurvProb,
			  double *SurvProb,
			  double futureTime,
			  double lambdaFutureTime)		  
{
	for (int k=0; k<timeSurvProb.size(); k++) SurvProb[k]=1;
	survProbaMult(ams, linear, para, timeSurvProb, SurvProb, futureTime, lambdaFutureTime);
};

void CIRcst::survProbaMult(
			  double ams,
              double linear, 
			  double para,
			  array<double> &timeSurvProb,
			  double *SurvProb,
			  double futureTime,
			  double lambdaFutureTime)		  
{
	if (futureTime==0) lambdaFutureTime = lambda0*(1+para)*linear;
	double scale = ams*linear;
	double alpha, beta, sqrtScale = sqrt(scale);
	for (int i=0; i<timeSurvProb.size(); i++)
	{
		if (timeSurvProb[i]>futureTime) 
		{
			alpha = alphaCIR(kappa,scale*thetaShift(para),sqrtScale*sigma , timeSurvProb[i] - futureTime);
			beta = betaCIR(kappa, sqrtScale*sigma, timeSurvProb[i] - futureTime);
			SurvProb[i] *= exp(alpha + ams*lambdaFutureTime*beta);
		}
		// else multiply SurvProb[i] by 1, i.e. do nothing
	}
};

void CIRcst::survProba(
		  array<double> &t, 
		  double *ft,
		  array<double> &timeSurvProb,
		  array<double> &ftimeSurvProb,
		  double *SurvProb,
		  double futureTime,
		  double lambdaFutureTime,
		  double h)
{
	for (int k=0; k<timeSurvProb.size(); k++) SurvProb[k]=1;
	survProbaMult(t, ft, timeSurvProb, ftimeSurvProb, SurvProb, futureTime, lambdaFutureTime, h);
}

void CIRcst::survProbaMult(
		  array<double> &t, 
		  double *ft,
		  array<double> &timeSurvProb,
		  array<double> &ftimeSurvProb,
		  double *SurvProb,
		  double futureTime,
		  double lambdaFutureTime,
		  double h)
{
	double ffutureTime;
	if (futureTime<=t[0]) ffutureTime=ft[0];
		else piecewiseLinearFunction(t,&ft[0], &futureTime, 1, &ffutureTime);

	double alpha, beta;
	int high = t.size()-1;
	int high2;
	for (int i=timeSurvProb.size()-1; i>=0; i--)
	{
		if (timeSurvProb[i] <= futureTime) return; // no change to SurvProb[j], j<=i
		alpha=0;
		beta=0;
		while ((high>=0)&&(t[high]>=timeSurvProb[i])) high--;
		// now high is the greatest integer such that timeSurvProb[i]>t[high], unless maybe high==-1
		if (high==-1)
			RicattiCIRlinear(kappa, theta, sigma, timeSurvProb[i], ftimeSurvProb[i], futureTime, ffutureTime, alpha,beta,h);
		else
		{
			high2=high;
			if (t[high2]<futureTime)
				RicattiCIRlinear(kappa, theta, sigma, timeSurvProb[i], ftimeSurvProb[i], futureTime, ffutureTime, alpha, beta, h);
			else
			{
				RicattiCIRlinear(kappa, theta, sigma, timeSurvProb[i], ftimeSurvProb[i], t[high2], ft[high2], alpha, beta, h);
				while ((high2>0)&&(t[high2-1]>futureTime))
				{
					RicattiCIRlinear(kappa, theta, sigma, t[high2], ft[high2], t[high2-1], ft[high2-1], alpha, beta, h);
					high2--;
				}
				RicattiCIRlinear(kappa, theta, sigma, t[high2], ft[high2], futureTime, ffutureTime, alpha, beta, h);
			}
		}
		SurvProb[i] *= exp(alpha + lambda0*beta);
	}
}


void CIRcst::FindTheta(double time, double SurvProb)		// find theta such that survProba(time) = SurvProb
{
	double alpha = alphaCIR(kappa,1,sigma,time);
	double beta = betaCIR(kappa, sigma, time);
	setTheta(Maths::max(( log(SurvProb) - lambda0*beta ) / alpha, sigma*sigma/(2*kappa)));
}


void CIRcst::diffuse(
				  double linear, 	
				  double para,
				  double *BrownianIncrement, // of size timesBM.size()-1
				  double *expDecay,          // of size timesBM.size() -1 = exp(-kappa(t_{k+1}-t_k))
		 		  array<double> &timesBM, 
				  double *CIRprocess,         // of size timesBM.size()
				  double futureTime,
				  double lambdaFutureTime)		
{
	if (futureTime==0) lambdaFutureTime = lambda0*(1+para);
			else lambdaFutureTime /= linear;  
			// linear just multiplies the SDE by linear.. we are diffusing the SDE, and then multiplying everything by linear
			// hence we need to divide first
	if (futureTime<timesBM[0]) futureTime = timesBM[0]; // or throw an error message...

	// no diffusion up to futureTime.. it sets up at least CIRprocess[0]
	int k;
	for (k=0; (timesBM[k]<futureTime+1e-10)&&(k<timesBM.size()); k++) CIRprocess[k] = lambdaFutureTime;

	if (k<timesBM.size())
	{
		double thetaPrimeDecayed;
		double thetaStrato = thetaShift(para) - 0.25*sigma*sigma/kappa;
		for (; k<timesBM.size(); k++)
		{
			thetaPrimeDecayed = thetaStrato*(1-expDecay[k-1]);
			CIRprocess[k] = Maths::max(0.0,thetaPrimeDecayed + Maths::square( sqrt(Maths::max(0.0,thetaPrimeDecayed + CIRprocess[k-1]* expDecay[k-1])) + 0.5 * sigma * BrownianIncrement[k-1]) * expDecay[k-1]);
		}
	}
	for (k=0; k<timesBM.size(); k++) CIRprocess[k] *= linear;
}

void CIRcst::condSurvProba(
				  double linear, 
				  double para,
				  double *BrownianIncrement, // of size timesBM.size() -1 
				  double *expDecay,          // of size timesBM.size() -1 = exp(-kappa(t_{k+1}-t_k))
	              array<double> &timesBM,      
				  double *CIRprocess,        // of size timesBM.size()
				  array<double> &timeSurvProb,
				  double *SurvProb,                // of size timesSurvProb.size() 
				  double futureTime,
				  double lambdaFutureTime)				
{ 
	diffuse(linear, para, BrownianIncrement, expDecay, timesBM, CIRprocess,futureTime, lambdaFutureTime);
	// no diffusion up to futureTime.. set up at least CIRprocess[0]
	int k;
	for (k=0; (timeSurvProb[k]<futureTime+1e-10)&&(k<timeSurvProb.size()); k++) SurvProb[k] = 0;

	if (k<timeSurvProb.size())
	{
		int interval = 0;
		SurvProb[k] = integrate(futureTime,timeSurvProb[k], interval, timesBM, CIRprocess);
		for (; k<timeSurvProb.size()-1; k++) SurvProb[k+1] = SurvProb[k] + integrate(timeSurvProb[k], timeSurvProb[k+1], interval, timesBM, CIRprocess);
	}
	for (k=0; k<timeSurvProb.size(); k++) SurvProb[k]=exp(-SurvProb[k]);
}


void CIRcst::diffuse(double *BrownianIncrement, // of size timesBM.size()-1
					 double *expDecay,          // of size timesBM.size() -1 = exp(-kappa(t_{k+1}-t_k))
		 			 array<double> &timesBM, 
					 double *CIRprocess,
					 double futureTime,
					 double lambdaFutureTime)
{
	if (futureTime==0) lambdaFutureTime = lambda0;
	if (futureTime<timesBM[0]) futureTime = timesBM[0]; // or throw an error message...
	// no diffusion up to futureTime.. it sets up at least CIRprocess[0]
	int k;
	for (k=0; (timesBM[k]<futureTime+1e-10)&&(k<timesBM.size()); k++) CIRprocess[k] = lambdaFutureTime;

	if (k<timesBM.size())
	{
		double thetaPrimeDecayed;
		double thetaStrato = theta - 0.25*sigma*sigma/kappa;
		for (; k<timesBM.size(); k++)
		{
			thetaPrimeDecayed = thetaStrato*(1-expDecay[k-1]);
			CIRprocess[k] = Maths::max(0.0,thetaPrimeDecayed + Maths::square( sqrt(Maths::max(0.0,thetaPrimeDecayed + CIRprocess[k-1]* expDecay[k-1])) + 0.5 * sigma * BrownianIncrement[k-1]) * expDecay[k-1]);
		}
	}
}
// of size timesBM.size()

void CIRcst::condSurvProba(
				 array<double> &f_at_timesBM,
				 double *BrownianIncrement, // of size timesBM.size() -1 
				 double *expDecay,          // of size timesBM.size() -1 = exp(-kappa(t_{k+1}-t_k))
	             array<double> &timesBM,      
				 double *CIRprocess,        // of size timesBM.size()
				 array<double> &timeSurvProb,
				 double *SurvProb,
				 double futureTime,
				 double lambdaFutureTime)
{
	diffuse(BrownianIncrement, expDecay, timesBM, CIRprocess,futureTime, lambdaFutureTime);
	// no diffusion up to futureTime.. set up at least CIRprocess[0]
	int k;
	for (k=0; (timeSurvProb[k]<futureTime+1e-10)&&(k<timeSurvProb.size()); k++) SurvProb[k] = 0;

	if (k<timeSurvProb.size())
	{
		int interval = 0;
		SurvProb[k] = integrate(futureTime,timeSurvProb[k], interval, timesBM, CIRprocess,&f_at_timesBM[0]);
		for (; k<timeSurvProb.size()-1; k++) SurvProb[k+1] = SurvProb[k] + integrate(timeSurvProb[k], timeSurvProb[k+1], interval, timesBM, CIRprocess,&f_at_timesBM[0]);
	}
	for (k=0; k<timeSurvProb.size(); k++) SurvProb[k]=exp(-SurvProb[k]);

}




void CIRcst::testSurvProb(double linear, double para, long seed, double timeSteps, int nbPaths, double futureTime, double lambdaFutureTime, array<double> &timeSurvProb, double *SurvProb)
{
	/*double lastTime = timeSurvProb[timeSurvProb.size()-1];
	int nbSteps = floor(lastTime/timeSteps)+1;

	array<double> timeSimulation(nbSteps+1); 
	int k,l,mc;
	for (k=0; k<=nbSteps; k++) timeSimulation[k] = (k*lastTime)/nbSteps;

	GaussianRandomSequence rng;
	rng.initializeObject(seed, nbSteps*nbPaths); 
	rng.populatearray(0);
	double * BMincrements = rng.getarray();
	
	Matrix<double> m_BMincrements(nbPaths,nbSteps);
	for (k=0; k<nbPaths; k++) 
		for (l=0; l<nbSteps; l++)
			m_BMincrements[k][l] = BMincrements[k*nbSteps+l]*sqrt(timeSimulation[l+1]-timeSimulation[l]);
	array<double> CIRprocess(timeSimulation.size()), SurvProbTemp(timeSurvProb.size());
	for (l=0; l<timeSurvProb.size(); l++) SurvProb[l]=0;
	array<double> expDecay(timeSimulation.size()-1);
	for (k=0; k<expDecay.size(); k++) expDecay[k]=exp(-0.5*kappa*(timeSimulation[k+1]-timeSimulation[k]));
	for (mc=0; mc<nbPaths; mc++)
	{
		condSurvProba(linear, para, &m_BMincrements[mc][0], &expDecay[0], timeSimulation, &CIRprocess[0], timeSurvProb, &SurvProbTemp[0], futureTime, lambdaFutureTime);
		for (l=0; l<timeSurvProb.size(); l++) SurvProb[l]+=SurvProbTemp[l];
	}
	for (l=0; l<timeSurvProb.size(); l++) SurvProb[l] /= double(nbPaths);*/
}

void CIRcst::testSurvProb2(array<double> &t, array<double> &ft, long seed, double timeSteps, int nbPaths, 
				 double futureTime, double lambdaFutureTime, array<double> &timeSurvProb, double *SurvProb)
{
/*	double lastTime = timeSurvProb[timeSurvProb.size()-1];
	int nbSteps = floor(lastTime/timeSteps)+1;

	array<double> timeSimulation(nbSteps+1); 
	int k, l;
	for (k=0; k<=nbSteps; k++) timeSimulation[k] = (k*lastTime)/nbSteps;

	array<double> f_at_t(timeSimulation.size());
	piecewiseLinearFunction(t, &ft[0], &timeSimulation[0], timeSimulation.size(), &f_at_t[0]);

	GaussianRandomSequence rng;
	rng.initializeObject(seed, nbSteps*nbPaths); 
	rng.populatearray(0);
	double * BMincrements = rng.getarray();
	Matrix<double> m_BMincrements(nbPaths,nbSteps);
	for (k=0; k<nbPaths; k++) 
		for (l=0; l<nbSteps; l++)
			m_BMincrements[k][l] = BMincrements[k*nbSteps+l]*sqrt(timeSimulation[l+1]-timeSimulation[l]);

	array<double> CIRprocess(timeSimulation.size()), SurvProbTemp(timeSurvProb.size());

	for (l=0; l<timeSurvProb.size(); l++) SurvProb[l]=0;

	array<double> expDecay(timeSimulation.size()-1);
	for (k=0; k<expDecay.size(); k++) expDecay[k]=exp(-0.5*kappa*(timeSimulation[k+1]-timeSimulation[k]));
	for (int mc=0; mc<nbPaths; mc++)
	{
		condSurvProba(f_at_t, &m_BMincrements[mc][0], &expDecay[0], timeSimulation, &CIRprocess[0], timeSurvProb, &SurvProbTemp[0], futureTime, lambdaFutureTime);
		for (l=0; l<timeSurvProb.size(); l++) SurvProb[l]+=SurvProbTemp[l];
	}
	for (l=0; l<timeSurvProb.size(); l++) SurvProb[l] /= double(nbPaths);*/
}





void AffineJump::setParameters(double frequency, array<double> &SizeOfJumps, array<double> &Impact, double Decay)
{
	freq = frequency;
	jumpSize = SizeOfJumps;
	impact = Impact;
	decay = Decay;
	timeOfJumps.resize(100); //// to be changed by something which depends on the frequency
}

void AffineJump::setParameters(double frequency, double jSize, array<double> &initRates, double Impact, double Decay)
{
	freq = frequency;
	jumpSize = initRates;
	for (int i=0; i<jumpSize.size(); i++) jumpSize[i] *= jSize;
	impact.resize(initRates.size());
	for (int i=0; i<impact.size(); i++) impact[i]=Impact;
	decay = Decay;
	timeOfJumps.resize(100); //// to be changed by something which depends on the frequency

}



void AffineJump::survProba( 
					double linear,
					array<double> &times, 
					double *SurvProb,
   					double futureTime,
					double *lambdaFutureTime)		  
{
	for (int k=0; k<times.size(); k++) SurvProb[k]=1;
	survProbaMult(linear,times,SurvProb,futureTime,lambdaFutureTime);
}


void AffineJump::condSurvProba(
					   double linear,
					   array<double> &times, 
					   int nbJumps,
					   Uniform &jumpTimes,   // will provide the time of the jumps
					   double maturity,		 // when multiplied by maturity
					   double *SurvProb,
	   				   double futureTime,
					   double *lambdaFutureTime)
{
	for (int k=0; k<times.size(); k++) SurvProb[k]=1;
	condSurvProbaMult(linear, times, nbJumps,jumpTimes,maturity,SurvProb,futureTime,lambdaFutureTime);
}

void AffineJump::survProbaMult( 
					double linear,
					array<double> &times, 
					double *SurvProb,      // (O) (name,index time) -> (index time)*m_nbNames+name
   					double futureTime,
					double *lambdaFutureTime)		  
{
	int k,m;
	double alpha, beta;
	if (futureTime == 0) 
	{
		for (k=0; k<times.size(); k++)
		{
			beta = - phi(times[k]);
			for (m=0; m<jumpSize.size(); m++)
			{
				alpha = freq * impact[m] / (linear*jumpSize[m] + decay) * (-linear*jumpSize[m] *(times[k]-futureTime) + log(1-linear*jumpSize[m]*beta));
				SurvProb[m*times.size()+k] *= exp(alpha);
			}
			// else multiply SurvProb[k] by 1, i.e. do nothing...
		}
	}
	else
	{
		for (k=0; k<times.size(); k++)
		{
			if (times[k]>futureTime)
			{
				beta = - phi(times[k]-futureTime);
				for (m=0; m<jumpSize.size(); m++)
				{
					alpha = freq * impact[m] / (linear*jumpSize[m] + decay) * (-linear*jumpSize[m] *(times[k]-futureTime) + log(1-linear*jumpSize[m]*beta));
					SurvProb[m*times.size()+k] *= exp(alpha + beta * linear* lambdaFutureTime[m]);
				}
			} 
			// else multiply SurvProb[k] by 1, i.e. do nothing...
		}
	}

}


double AffineJump::condSurvProbaMult(
					   double linear,
					   array<double> &times, 
					   int nbJumps,          // nbJumps between futureTime and maturity
					   Uniform &jumpTimes,   // will provide the time of the jumps
					   double maturity,		 
					   double *SurvProb,
	   				   double futureTime,
					   double *lambdaFutureTime)
{
	int k,l,m;
	double beta, timeOfJump;
	double firstJump = 999;
    long timeSize = times.size();
	if (futureTime == 0) 
	{
		for (k=0; k<nbJumps; k++)
		{
			timeOfJump = *jumpTimes.getVector()*maturity;
			if (firstJump>timeOfJump) firstJump = timeOfJump;
			for (l=0; l<timeSize; l++)
			{
				if (timeOfJump<times[l])
				{
					beta = phi(times[l]-timeOfJump);
					for (m=0; m<jumpSize.size(); m++) SurvProb[m*timeSize+l] *= 1 - impact[m]*(1 - 1/(1 + linear*jumpSize[m]*beta));
				}
			}
		}
	}
	else
	{
		for (l=0; l<timeSize; l++)
			if (times[l]>futureTime) 
				for (m=0; m<jumpSize.size(); m++)
					SurvProb[m*timeSize+l]  *= exp(-linear*lambdaFutureTime[m]*phi(times[l]-futureTime));
		for (k=0; k<nbJumps; k++)
		{
			timeOfJump = futureTime + (*jumpTimes.getVector())*(maturity - futureTime);
			if (firstJump>timeOfJump) firstJump = timeOfJump;
			for (l=0; l<timeSize; l++)
			{
				if (timeOfJump<times[l])
				{
					beta = phi(times[l]-timeOfJump);
					for (m=0; m<jumpSize.size(); m++) 
						SurvProb[m*timeSize+l] *= 1 - impact[m]*(1 - 1/(1 + linear*jumpSize[m]*beta));
				}
			}
		}
	}
	return firstJump;
}



void AffineJump::condSpreadAndHazard(        // give SurvSpread and hazard rates
					array<double> &times,
				    int nbJumps,
					Uniform &jumpTimes,
					double maturity,
					Uniform &jumpUniform,
					double *intensity,       // (O) (name,index time) -> (name)*times.size()+index time
					double *intIntensity)    // (O) (name,index time) -> (name)*times.size()+index time
{
	int k,l,m;
	double phit, expt, actualJumpSize;
	for (k=0; k<nbJumps; k++) timeOfJumps[k] = *jumpTimes.getVector()*maturity;
	Algorithm::shellSort(timeOfJumps, 0, nbJumps-1);

	for (k=0; k<times.size(); k++)
		for (m=0; m<jumpSize.size(); m++)
		{
			intensity[m*times.size()+k] = 0;
			intIntensity[m*times.size()+k] = 0;
		}
	for (l=0; l<nbJumps; l++)
	{
		double * jumps = jumpUniform.getVector(); 

		for (k=0; k<times.size(); k++)
		{
			if (timeOfJumps[l]<times[k])
			{
				expt = exp(-decay*(times[k] - timeOfJumps[l]));
				phit = (1 - expt) / decay;
				for (m=0; m<jumpSize.size(); m++)
				{
					actualJumpSize = getJump(jumpSize[m], impact[m], jumps[m]);
					intensity[m*times.size()+k] += actualJumpSize*expt;
					intIntensity[m*times.size()+k] += actualJumpSize*phit;
				}
			}
		}
	}
/*
	int k,m;
	nbJumps = MIN(timeOfJumps.size(),nbJumps);
	for (k=0; k<nbJumps; k++) timeOfJumps[k] = *jumpTimes.getVector()*maturity;
	SORT(&timeOfJumps[0], 0, nbJumps);
	int nextJump = 0;
	
	double expt = exp(-decay*times[0]);
	double phit = (1 - expt)/decay;
	double actualJumpSize;
	k=0;
	while (timeOfJumps[nextJump]>times[k]) 
	{
		for (m=0; m<jumpSize.size(); m++)
		{
			intensity[m*times.size()+k] = 0;
			intIntensity[m*times.size()+k] = 0;
		}
		k++;
	}
	for (; (k<times.size())&&(nextJump<nbJumps); k++)
	{
		expt = exp(-decay*(times[k]-times[k-1]));
		phit = (1 - expt)/decay;
		for (m=0; m<jumpSize.size(); m++)
		{
			intensity[m*times.size()+k] = intensity[m*times.size()+k-1]*expt;
			intIntensity[m*times.size()+k] = intIntensity[m*times.size()+k-1]
										+ intensity[m*times.size()+k-1]*phit;
		}
		while (timeOfJumps[nextJump]<times[k]) // jump in (t[k-1],t[k])
		{
			expt = exp(-decay*(times[k]-timeOfJumps[nextJump]));
			phit = (1 - expt)/decay;
			double * jumps = jumpUniform.getVector(); 
			for (m=0; m<jumpSize.size(); m++)
			{
				actualJumpSize = getJump(jumpSize[m], impact[m], jumps[m]);
				intensity[m*times.size()+k] += actualJumpSize*expt;
				intIntensity[m*times.size()+k] += actualJumpSize*phit;
			}
			nextJump++;
		}
	}
	for (; k<times.size(); k++)
	{
		expt = exp(-decay*(times[k]-times[k-1]));
		phit = exp(-decay*times[k-1])*(1 - expt)/decay;
		for (m=0; m<jumpSize.size(); m++)
		{
			intensity[m*times.size()+k] = intensity[m*times.size()+k-1]*expt;
			intIntensity[m*times.size()+k] = intIntensity[m*times.size()+k-1]
										+ intensity[m*times.size()+k-1]*phit;
		}
	}*/
}

void AffineJump::testAffineJump(array<double> &times, double *SurvProb, double *SurvProbCond, int nbPaths, long seed, double futureTime, double maturity, double lambdaFutureTime)
{ // works if AffineJump is set with one name only!
	PoissonProcess Poisson(&freq,1, maturity - futureTime, seed);
	Poisson.setIpath(0);
	Uniform jumpTimes(seed+1,1);

	survProba(1, times, SurvProb, futureTime, &lambdaFutureTime);
	array<double> SurvProbTemp(times.size(),0.0), SurvProbTemp2(times.size(),0.0);
	int nbJumps;
	double probaZeroJump = Poisson.probaNoJump();
	for (int mc=0; mc<nbPaths; mc++)
	{
		Poisson.GetPoissonNbJumps(true, &nbJumps);
		condSurvProba(1, times, nbJumps, jumpTimes, maturity, &SurvProbTemp[0], futureTime, &lambdaFutureTime);
		for (int k=0; k<times.size(); k++) 
		{
			SurvProbTemp2[k]+= SurvProbTemp[k];
		}
	}
	for (int k=0; k<times.size(); k++) 
		if (times[k]>futureTime)
			SurvProbCond[k] = (1-probaZeroJump)*SurvProbTemp2[k] / double(nbPaths)
						+ probaZeroJump*exp(-lambdaFutureTime*phi(times[k]-futureTime));
};

double AffineJump::phi(double t) { 
    std::map<double, double>::iterator pos = phiMap.find(t);
    if (pos != phiMap.end())
        return pos->second;
    else {
        double value = (1.0 - exp(-decay*t)) / decay; 
        phiMap.insert(std::make_pair(t,value));
        return value;
    }
};















DRLIB_END_NAMESPACE
