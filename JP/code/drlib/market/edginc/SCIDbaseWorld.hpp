#ifndef ___SCIDBASEWORLD_HPP
#define ___SCIDBASEWORLD_HPP

#include "edginc/SCIDaffineProcesses.hpp"
#include "edginc/DateTime.hpp"

DRLIB_BEGIN_NAMESPACE

class MARKET_DLL BaseWorld
{
public:
	BaseWorld(){};

	CIR m_idioCIR;     // the idiosyncratic parameters (and ability to compute SurvProb and cond SurvProb)
	CIRcst m_cfCIR;		// the common factor
	array<double> m_ams;  // size m_nbNames; // the scaling factor for the common factor
	std::vector<AffineJump> m_jumps; // size m_nbMarkets // the jump parameters (and ability to compute SurvProb and cond SurvProb)
	bool HasCommonFactor, HasIdioVol; // we don't have HasJump as we know whether we have jumps by nbJumps;
	
	array<double> m_idioSurvProb, m_condSurvProb; // size nbNames * times.size() // initialize in setIdioSurvProb
							// storage for some Surv Proba

	void survProba(  // unconditional surv proba, possibly computed in the future
			  bool idio, bool cf, bool jump,
			  double linear,
			  double para,
			  array<double> &timeSurvProb,
			  double *SurvProb,  // (O) (name,index time) -> name*timeSurvProb.size()+index time
			  double futureTime = 0,
			  double *lambdaFutureTimeIdio = 0,
			  double lambdaFutureTimeCf = 0,
			  DoubleMatrix *lambdaFutureTimeJump = 0);
	
	void setTimes(array<double> &timeSimulation);  // time at which we simulate the common factor
	// compute idiosyncratic survival probability, and memory management
	void InitializeCondSurvProb(           
		 		   double linear,
				   double para,
				   array<double> &timeSurvProb,
				   double futureTime = 0,
				   double *lambdaFutureTimeIdio = 0);  // initialize memory and set idioSurvProb;
	// conditional on Poisson processes and common factor, return survival proba of all names at timeSurvProb
	// used only by the fast MC
	double getCondSurvProba(
		 		   double linear,
				   double para,
				   array<double> &timeSurvProb,
				   double *bmIncrements,
				   int *nbJumps = 0,
				   Uniform *jumpTimes = 0,        // provides the jump Times 
				   double maturity = 0 ,         
   				   double futureTime = 0, 
				   double lambdaFutureTimeCf = 0,
				   DoubleMatrix *lambdaFutureTimeJump = 0);  // return the time of first jump



	// calibrate the diffusions, assuming the jump markets are already given
	void CalibrateDiffusion(
					   int nbNames,
					   DateTime today,
					   DateTimeArray &singleNameDates,
					   DoubleMatrix &survProbas,  // size (nbNames,singleNameDates.size())
					   DateTime &cfDate,
					   const double volIdio,
					   const double volCf,
					   const double decay,
					   const double ratio);

						
private:
	array<double> cir, SurvProbTemp; // size 	timeSimulation.size(), size times.size()
	array<double> m_timeSimulation, m_expDecayTimes;
};

class MARKET_DLL Spread // to store spread of one given world in MC simulation
{
public:
	Spread(){};
	~Spread(){};
	void initialize(int nbNames, int nbMarkets, int nbTimes) 
	{ 
		idio.resize(nbNames*nbTimes); 
		cf.resize(nbTimes); 
		baseTotal.resize(nbNames*nbTimes); 
		total.resize(nbNames*nbTimes); 
		integral.resize(nbNames*nbTimes);
		baseIntegral.resize(nbNames*nbTimes);
	};

	vector<double> idio;  // size nbNames * nbTimes
	vector<double> cf;    // size nbTimes
	vector<double> baseTotal;  // size nbNames * nbTimes
	vector<double> total;  // size nbNames * nbTimes
	vector<double> baseIntegral;  // size nbNames * nbTimes
	vector<double> integral;  // size nbNames * nbTimes
};

DRLIB_END_NAMESPACE

#endif
