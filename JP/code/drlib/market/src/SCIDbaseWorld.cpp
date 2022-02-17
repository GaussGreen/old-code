#include "edginc/config.hpp"
#include "edginc/SCIDbaseWorld.hpp"

DRLIB_BEGIN_NAMESPACE

void BaseWorld::survProba(
			  bool idio, bool cf, bool jump,
			  double linear,
			  double para,
			  array<double> &timeSurvProb,
			  double *SurvProb,   // (O) (name,index time) -> name*timeSurvProb.size()+index time
			  double futureTime,
			  double *lambdaFutureTimeIdio,
			  double lambdaFutureTimeCf,
			  DoubleMatrix *lambdaFutureTimeJump)
{
	if (futureTime==0)
	{
		if ((cf)&&(HasCommonFactor))
			for (int m=0; m<m_ams.size(); m++) 
				m_cfCIR.survProba(m_ams[m],linear,para,timeSurvProb,&SurvProb[m*timeSurvProb.size()]);
			else for (int k=0; k<m_ams.size()*timeSurvProb.size(); k++) SurvProb[k]=1;
		if (idio) m_idioCIR.survProbaMult(linear,para,timeSurvProb,SurvProb);
		if (jump) for (unsigned int r=0; r<m_jumps.size(); r++)
			m_jumps[r].survProbaMult(linear,timeSurvProb,SurvProb);
	}
	else
	{
		if ((cf)&&(HasCommonFactor))
		{
			for (int m=0; m<m_ams.size(); m++)
				m_cfCIR.survProba(m_ams[m],linear,para,timeSurvProb,&SurvProb[m*timeSurvProb.size()],futureTime,lambdaFutureTimeCf);
		}
		else 
			for (int k=0; k<m_ams.size()*timeSurvProb.size(); k++) SurvProb[k]=1;
		if (idio) 
			m_idioCIR.survProbaMult(linear,para,timeSurvProb,SurvProb,futureTime,lambdaFutureTimeIdio);
		if (jump) 
		{
			for (unsigned int r=0; r<m_jumps.size(); r++) 
				m_jumps[r].survProbaMult(linear,timeSurvProb,SurvProb,futureTime,&((*lambdaFutureTimeJump)[r][0]));
		}
	}
}


void BaseWorld::CalibrateDiffusion(
					   int nbNames,
					   DateTime today,
					   DateTimeArray &singleNameDates,
					   DoubleMatrix &survProbas,  // size (nbNames,singleNameDates.size())
					   DateTime &cfDate,
					   const double volIdio,
					   const double volCf,
					   const double decay,
					   const double ratio)
{
	int k,m;
	unsigned r;
	DoubleArray singleNameTimesWithToday(singleNameDates.size()+1,0.0);
	for (k=0; k<singleNameDates.size(); k++)
		singleNameTimesWithToday[k+1] = today.yearFrac(singleNameDates[k]);

	int indexMaturityCf=0;
	while (singleNameDates[indexMaturityCf]<cfDate) indexMaturityCf++;
	double maturityCf = singleNameTimesWithToday[indexMaturityCf+1];

	DoubleArray baseRates(nbNames);
	for (m=0; m<nbNames; m++)
		baseRates[m] = - log( survProbas[0][m] ) / singleNameTimesWithToday[1];

// find whether we have jumps or not
	HasIdioVol = (volIdio>1e-10);
/////////// find common factor parameters
	m_ams.resize(nbNames);
	// find theta for Common Factor
	m_cfCIR.setParameters(1,volCf,decay);
	double sumLambda0=0, sumTheta=0;
	for (m=0; m<nbNames; m++) 
	{	
		sumTheta += (log(survProbas[indexMaturityCf][m]) - baseRates[m]*betaCIR(decay,volCf*sqrt(baseRates[m]),maturityCf)) / alphaCIR(decay,1,volCf*sqrt(baseRates[m]),maturityCf);
		sumLambda0 += baseRates[m];
	}
	m_cfCIR.setTheta(sumTheta/sumLambda0);
	// now find the m_ams
	array<double> matCf(1,maturityCf);
	double meanGamma = m_cfCIR.theta*maturityCf 
		             - (m_cfCIR.theta - m_cfCIR.lambda0)*(1-exp(-m_cfCIR.kappa*maturityCf))/m_cfCIR.kappa;
	array<double> SurvProb2(nbNames,1.0);
	for (r=0; r<m_jumps.size(); r++) m_jumps[r].survProbaMult(1,matCf,&SurvProb2[0]);
	double sum=0;
	for (m=0; m<nbNames; m++)
	{
		m_ams[m] = Maths::max(0.0, -ratio*log(survProbas[indexMaturityCf][m]/SurvProb2[m])/meanGamma);
		sum += m_ams[m];
	}
	HasCommonFactor = (sum>1e-10);

	array<double> SurvProbBase(nbNames*singleNameTimesWithToday.size(), 1.0);
	for (r=0; r<m_jumps.size(); r++) 
		m_jumps[r].survProbaMult(1,singleNameTimesWithToday,&SurvProbBase[0]);
	for (m=0; m<nbNames; m++)
	{
		if (HasCommonFactor) m_cfCIR.survProbaMult(m_ams[m],1.0,0.0,singleNameTimesWithToday,&SurvProbBase[m*singleNameTimesWithToday.size()]);
		for (k=1; k<singleNameTimesWithToday.size(); k++) 
			SurvProbBase[m*singleNameTimesWithToday.size()+k] = Maths::min( 1.0, survProbas[k-1][m] / SurvProbBase[m*singleNameTimesWithToday.size()+k] );
	}
		
	// we convert SurvProbBase to SDE parameters...
	array<double> lambda0(nbNames);
	for (m=0; m<nbNames; m++)
		lambda0[m] = Maths::max(0.2*baseRates[m],baseRates[m]-m_ams[m]);
	
	m_idioCIR.setParameters(lambda0, volIdio , decay);
	m_idioCIR.TimeDepTheta(singleNameTimesWithToday, &SurvProbBase[0]);
}


void BaseWorld::setTimes(array<double> &timeSimulation)
{
	m_timeSimulation = timeSimulation;
	cir.resize(m_timeSimulation.size());
	m_expDecayTimes.resize(m_timeSimulation.size()-1);
	for (int k=0; k<m_timeSimulation.size()-1; k++)
		m_expDecayTimes[k] = exp(-0.5*m_cfCIR.kappa*(m_timeSimulation[k+1]-m_timeSimulation[k]));
}


void BaseWorld::InitializeCondSurvProb(
				   double linear,
				   double para,
				   array<double> &timeSurvProb,
				   double futureTime, 
				   double *lambdaFutureTimeIdio)  // initialize memory and set idioSurvProb;
{
	m_idioSurvProb.resize(m_ams.size()*timeSurvProb.size());
	m_condSurvProb.resize(m_ams.size()*timeSurvProb.size());
	SurvProbTemp.resize(timeSurvProb.size()); 
	if (futureTime==0) m_idioCIR.survProba(linear,para,timeSurvProb,&m_idioSurvProb[0]);
		else m_idioCIR.survProba(linear,para,timeSurvProb,&m_idioSurvProb[0],futureTime,lambdaFutureTimeIdio);
}

double BaseWorld::getCondSurvProba(
				   double linear,
				   double para,
				   array<double> &timeSurvProb,
				   double *bmIncrements,
				   int *nbJumps,
				   Uniform *jumpTimes,
				   double maturity,         // provides the jump Times (nbJumps of them)...
   				   double futureTime, 
				   double lambdaFutureTimeCf,
				   DoubleMatrix *lambdaFutureTimeJump)
{
	int k,m;
	unsigned r;
	m_condSurvProb = m_idioSurvProb;
	double firstJump = 999, firstJump2;
	if (futureTime!=0)
	{
		if (HasCommonFactor)
		{
			m_cfCIR.condSurvProba(linear,para,bmIncrements,&m_expDecayTimes[0],m_timeSimulation,&cir[0],timeSurvProb,&SurvProbTemp[0],futureTime,lambdaFutureTimeCf);
			for (m=0; m<m_ams.size(); m++)
				for (k=0; k<timeSurvProb.size(); k++) 
					m_condSurvProb[m*timeSurvProb.size()+k] *= pow(SurvProbTemp[k],m_ams[m]);
		}
		if (nbJumps!=0)
		{
			for (r=0; r<m_jumps.size(); r++)
			{
 				firstJump2 = m_jumps[r].condSurvProbaMult(linear,timeSurvProb,nbJumps[r],*jumpTimes,maturity,&m_condSurvProb[0],futureTime,&((*lambdaFutureTimeJump)[r][0]));
				if (firstJump2<firstJump) firstJump = firstJump2;
			}
		}
		else
		if (lambdaFutureTimeJump!=0)
		{
			for (r=0; r<m_jumps.size(); r++)
			{
				if ((*lambdaFutureTimeJump)[r][0]!=0)
 					m_jumps[r].condSurvProbaMult(linear,timeSurvProb,0,*jumpTimes,maturity,&m_condSurvProb[0],futureTime,&((*lambdaFutureTimeJump)[r][0]));
			}
		}
	}
	else
	{
		if (HasCommonFactor)
		{
			m_cfCIR.condSurvProba(linear,para,bmIncrements,&m_expDecayTimes[0],m_timeSimulation,&cir[0],timeSurvProb,&SurvProbTemp[0]);
			for (m=0; m<m_ams.size(); m++)
				for (k=0; k<timeSurvProb.size(); k++) 
					m_condSurvProb[m*timeSurvProb.size()+k] *= pow(SurvProbTemp[k],m_ams[m]);
		}
		if (nbJumps!=0)
			for (r=0; r<m_jumps.size(); r++) 
			{
				firstJump2 = m_jumps[r].condSurvProbaMult(linear,timeSurvProb,nbJumps[r],*jumpTimes,maturity,&m_condSurvProb[0]);
				if (firstJump2<firstJump) firstJump = firstJump2;
			}
	}
	return firstJump;
}


DRLIB_END_NAMESPACE

