#include "edginc/config.hpp"
#include "edginc/JumpAffine.hpp"
#include "edginc/Algorithm.hpp"

DRLIB_BEGIN_NAMESPACE


void PoissonDiffuseHelp::Initialize(DoubleArray &_avFreq)
{ 
	avFreq = _avFreq;
	nbJumps.resize(avFreq.size());
	jumpTimes.resize(avFreq.size());
	for (int i=0; i<avFreq.size(); ++i)
		jumpTimes[i] = DoubleArraySP(new DoubleArray());
	sumFreq = accumulate(avFreq.begin(),avFreq.end(),0.0); 
};

int PoissonDiffuseHelp::simulatePoisson(double freq, bool cond, double unif)  // to the uniform number unif associate number of jump
{
	double proba, Sproba;
	int nbJumps;
	if (cond==true)
	{
		nbJumps = 1;
		proba = Sproba = freq / (exp(freq) - 1.0);
	}
	else
	{
		nbJumps = 0;
		proba = Sproba = exp(-freq);
	}
	while ((unif>Sproba)&&(Sproba<1-1e-20)) // second test to make sure we do not spend the rest of our life in this loop
	{
		nbJumps++;
		proba *= freq / nbJumps;
		Sproba += proba;
	}; 
	return nbJumps;
}

int PoissonDiffuseHelp::simulateBinomial(double p, int N, double unif)
{
	int n=0;
	double proba, Sproba;
	proba = Sproba = pow(1.0-p,N);
	while ((unif>Sproba)&&(n<N))
	{
		proba *= p*(N-n)/( (1-p)*(n+1) );
		Sproba += proba;
		++n;
	}
	return n;
}

int PoissonDiffuseHelp::GetPoissonNbJumps(double T, 
										  bool cond,
										  double unif,  // unif is the uniform rv which will give the total nb of jumps
										  Sequence *rng, // random number then used to select which Poisson process jumped
										  long& iPath,
										  double t) // if cond false, nb Jumps conditionned on at least one jump...
{

	// we simulate the nb Jumps of the sum of the Poisson processes...
	int totalNbJumps = simulatePoisson(sumFreq*(T-t),cond,unif);

	// now we need to decide which one has it
	double alpha;
	int nbJumpsLeft = totalNbJumps;
	double freqLeft = sumFreq;
	int i;
	for (i=0; (i<avFreq.size())&&(nbJumpsLeft>0); ++i)
	{
		alpha = avFreq[i] / freqLeft;
		if (alpha>1-1e-9) nbJumps[i] = nbJumpsLeft;
		else
		{
			rng->populateVector(iPath);
			iPath++;
			nbJumps[i] = simulateBinomial(alpha,nbJumpsLeft,*(rng->getVector()));
		}
		nbJumpsLeft -= nbJumps[i];
		freqLeft -= avFreq[i];
	}
	for (; i<avFreq.size(); ++i) nbJumps[i]=0;
	firstJump = 1e10;
	return totalNbJumps;
}

void PoissonDiffuseHelp::simulateJumpTimeCondNbJumps(double T, 
													 int index,
													 Sequence *rng, 
													 long& iPath,
													 double t)
{
	if (nbJumps[index]>0)
	{
		if (nbJumps[index]>jumpTimes[index]->size())  jumpTimes[index]->resize(nbJumps[index]);
		for (int i=0; i<nbJumps[index]; ++i)
		{
			rng->populateVector(iPath);
			iPath++;
			double jtime = t + *(rng->getVector()) * (T-t);
			(*jumpTimes[index])[i] = jtime;
			if (jtime<firstJump)
				firstJump = jtime;
		}
	}
}



/////////////////////////////////////////////////////////////////////////////
IObject* JumpAffine::defaultConstructor()
{
	return new JumpAffine();
}

CClassConstSP const JumpAffine::TYPE =
CClass::registerClassLoadMethod("JumpAffine",
								typeid(JumpAffine),
								JumpAffine::load);

void JumpAffine::load(CClassSP& clazz)
{
	clazz->setPublic();
	REGISTER(JumpAffine, clazz);
	SUPERCLASS(CObject);
	IMPLEMENTS(IAffineProcess);
	EMPTY_SHELL_METHOD(defaultConstructor);
}

JumpAffine::JumpAffine():
IAffineProcess(),
CObject(TYPE),
guessPosf(0),
h(0.25),
pieceT(DoubleArraySP(new DoubleArray(1,100)))
{
};

JumpAffine::JumpAffine(double _lambda0, double _kappa, double _freq, IJumpDistributionSP _jump):
IAffineProcess(),
CObject(TYPE),
guessPosf(0),
h(0.25),
pieceT(DoubleArraySP(new DoubleArray(1,100)))
{
	setParameters(_lambda0, _kappa, _freq, _jump) ; 
};



void JumpAffine::setFs(LinearInterpolantSP &_f)
{
	f=_f; 
	guessPosf = 0;
};


void JumpAffine::setParameters(double _lambda0, double _kappa, double _freq, IJumpDistributionSP _jump) 
{ 
	kappa=_kappa; 
	lambda0=_lambda0; 
	freq = _freq;
	jump = _jump;
};


double JumpAffine::laplace(double u, double x) const
{
	int index = lower_bound(pieceT->begin(),pieceT->end(),u)-pieceT->begin();
	if (index==pieceT->size()) --index;
	// t_{i} <= u < t_{i+1}
	return jump->laplace(index, u, x);
}


double JumpAffine::minusLogSurvProba(double bigT, 
									double futureTime, 
									double lambdaFutureTime)
{
	QLIB_VERIFY( (pieceT.get()!=0) && (pieceT->size()>0), "pieceT not constructed");
	if (futureTime==0.0) lambdaFutureTime = lambda0;
	return minusLogSurvProbability(this,bigT,futureTime,lambdaFutureTime); 
};


void VFjumpPartHelper(double K0, double K1, double l0, double l1, double laplace, double rho, 
					  double beta,
					  double &fa, double &fb)
{
	fa = - K0*beta - l0*(laplace-1.0);
	fb = rho - K1*beta - l1*(laplace-1.0);
}

void JumpAffine::ricattiVF(int index, double u, double beta, double &fa, double &fb) const
{
	VFjumpPartHelper(0.0, - kappa, freq, 0.0, jump->laplace(index, u, -beta),
		f->valueWithGuess(u,guessPosf), beta, fa, fb);
}


void JumpAffine::simulateJumpSizes(UniformRandomSequence *rng,
								   long& iPath)
{
	if (nbJumps>jumpSizes.size()) jumpSizes.resize(nbJumps);
	for (int i=0; i<nbJumps; ++i)
		jumpSizes[i] = jump->simulateOneJump((*jumpTimes)[i], rng, iPath);
}

double JumpAffine::readLambdaT(double T,
							   double futureTime,
							   double lambdaFutureTime)
{
	if (futureTime == 0 ) lambdaFutureTime = lambda0;
	double lambda = lambdaFutureTime*exp(-kappa*(T-futureTime));
	for (int i=0; i<nbJumps; ++i)
		if ( (lambdaFutureTime<(*jumpTimes)[i]) && ((*jumpTimes)[i]<=T) )
			lambda += jumpSizes[i]*exp(-kappa*(T-(*jumpTimes)[i]));
	return lambda;
}

double JumpAffine::read_F_LambdaT(double T, 
								  double futureTime, 
								  double lambdaFutureTime)
{
	return f->valueWithGuess(T,guessPosf)*readLambdaT(T,futureTime,lambdaFutureTime);
}

double JumpAffine::integrateExp(double t) // return \int_s^{s+t} exp(-kappa*(u-s))du = \int_0^{t} exp(-kappa*u)du 
{
	if (kappa<1e-10) return t;
	else return (1.0 - exp(-kappa*t))/kappa;
}

double JumpAffine::integrateIntegrateExp(double t) // return \int_s^{s+t} exp(-kappa*(u-s))du = \int_0^{t} exp(-kappa*u)du 
{
	if (kappa<1e-10) return t*t*0.5;
	else return ( exp(-kappa*t) - 1.0 + kappa*t ) / (kappa*kappa);
}


double JumpAffine::integrateExpF(double s, double t)  // return \int_s^t f_u exp(-kappa*(u-s))du
{
	if (s>=t) return 0.0;
	QLIB_VERIFY(pieceT->size()>0, "pieceT not initialized");
	int index = upper_bound(pieceT->begin(),pieceT->end(),s) - pieceT->begin();
	if (index==pieceT->size()) --index;
	double fx, dfx, t1 = s, t2 = s, integrale = 0.0;

	if (kappa<1e-10)
	{
		while (t2<t)
		{
			t2 = min((*pieceT)[index],t);
			fx = f->valueWithGuess(t1,guessPosf);
			dfx = f->valueWithGuess(t1,1,guessPosf);
			// we now add /int_{t1}^{t2} (fx+dfx(u-t1))*exp(-kappa(u-s))du = 
			integrale += (t2-t1)*(fx+.5*dfx*(t2-t1));
			t1 = t2;
			++index;
		} 

	}
	else
	{	
		double expont1 = 1.0, expont2;
		while (t2<t)
		{
			t2 = min((*pieceT)[index],t);
			expont2 = exp(-kappa*(t2-s));
			fx = f->valueWithGuess(t1,guessPosf);
			dfx = f->valueWithGuess(t1,1,guessPosf);
			// we now add /int_{t1}^{t2} (fx+dfx(u-t1))*exp(-kappa(u-s))du = 
			integrale += ( expont1*(dfx+fx*kappa) - expont2*(dfx*(1.0 + kappa*(t2-t1)) + fx*kappa) ) / (kappa*kappa);
			expont1 = expont2;
			t1 = t2;
			++index;
		} 
	}
	return integrale;
}


double JumpAffine::readMinusLogSurvProba(bool weKnowJumpSizes,
										 double T, 
										 double futureTime, 
										 double lambdaFutureTime)
{
	if (futureTime == 0 ) lambdaFutureTime = lambda0;
	double intlambda;
	if (weKnowJumpSizes)
	{
		intlambda = lambdaFutureTime*integrateExp(T-futureTime);
		for (int i=0; i<nbJumps; ++i)
			if ( (lambdaFutureTime<(*jumpTimes)[i]) && ((*jumpTimes)[i]<=T) )
				intlambda += jumpSizes[i]*integrateExp(T-(*jumpTimes)[i]);
	}
	else
	{
		intlambda = lambdaFutureTime*integrateExp(T-futureTime);
		for (int i=0; i<nbJumps; ++i)
			if ( (lambdaFutureTime<(*jumpTimes)[i]) && ((*jumpTimes)[i]<=T) )
				intlambda += -log( laplace((*jumpTimes)[i],integrateExp(T-(*jumpTimes)[i])) );
	}
	return intlambda;
}

double JumpAffine::read_F_MinusLogSurvProba(bool weKnowJumpSizes,
											double T, 
											double futureTime, 
											double lambdaFutureTime)
{
	if (futureTime == 0) lambdaFutureTime = lambda0;
	double intlambda;
	if (weKnowJumpSizes)
	{
		intlambda = lambdaFutureTime*integrateExpF(futureTime,T);
		for (int i=0; i<nbJumps; ++i)
			if ( (lambdaFutureTime<(*jumpTimes)[i]) && ((*jumpTimes)[i]<=T) )
				intlambda += jumpSizes[i]*integrateExpF((*jumpTimes)[i],T);
	}
	else
	{
		intlambda = lambdaFutureTime*integrateExpF(futureTime,T);
		for (int i=0; i<nbJumps; ++i)
			if ( (lambdaFutureTime<(*jumpTimes)[i]) && ((*jumpTimes)[i]<=T) )
				intlambda += -log( laplace((*jumpTimes)[i],integrateExpF((*jumpTimes)[i],T) ) );
	}
	return intlambda;
}


/* external symbol to allow class to be forced to be linked in */
bool JumpAffineLoad(){
	return (JumpAffine::TYPE != 0);
}


class JumpAffinesolver: public Func1D::NoDeriv
{
public:
	JumpAffinesolver(JumpAffine *_jumpAffine,IJumpDistributionSP _jumpDist):
	  jumpAffine(_jumpAffine), 
	  jumpDist(_jumpDist){};

	virtual double operator()(double  x) const
	{
		jumpDist->storeParameters(index);
		jumpDist->multiplyJump(index,x);
		double answer = jumpAffine->minusLogSurvProba(timeCalib) - target;
		jumpDist->restoreParameters(index);
		return answer;
	}
	double target;
	int index;
	mutable JumpAffine *jumpAffine;
	mutable IJumpDistributionSP jumpDist;
	double timeCalib;
};


void JumpAffine::calibrate(DoubleArray &t, double *spt, double ratio, const double MAXMULTJUMP)
{
	QLIB_VERIFY(t.size()>0, "not enough data to calibrate");
	pieceT = DoubleArraySP(new DoubleArray(t));
	jump->setPiecewiseTimes(pieceT);

	JumpAffinesolver evaluation(this,jump);
    
	LinearInterpolantSP newf = LinearInterpolantSP( new LinearInterpolant(DoubleArray(1,0.0),DoubleArray(1,1.0)));
	setFs(newf);

	ZBrent solver(1e-5);
	ZBracPositive bracketer;
	double guess1, guess2, sol = 1.0;
	
	for (int i=0; i<t.size(); ++i)
	{
		evaluation.index = i;
		evaluation.target = -ratio*log(spt[i]);
		evaluation.timeCalib = t[i];
		if ((sol<1e-10)||(sol>MAXMULTJUMP*0.9)) guess2 = 1.0; 		else guess2 = sol;

		if (evaluation(0.0)>0) sol = 0.0;
			else if (evaluation(MAXMULTJUMP)<0) sol = MAXMULTJUMP;
				else
				{
					bracketer.bracket(evaluation,guess1,guess2);
					sol = solver.solve(evaluation,guess1,guess2);
				}
		jump->multiplyJump(i,sol);
	}
}


DRLIB_END_NAMESPACE
