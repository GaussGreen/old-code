#ifndef __INTEGRATIONTEST_H__
#define __INTEGRATIONTEST_H__

#include "INTIntegrator.h"

//----------------------------------------------------------------------------------------------

typedef enum
{ 
	userdefinedalpha, 
	optimisedalpha, 
	haselgrovealpha, 
	cyclotomicalpha, 
	bakeralpha, 
	niederreiteralpha, 
	primerootalpha 
} 
Alphatype;

static int numberofalphatypes = 7;	

//----------------------------------------------------------------------------------------------

typedef enum 
{ 
	none, 
	stratifiedmontecarlo, 
	randomwalk 
} 
Optimisationmethod;

//----------------------------------------------------------------------------------------------

//double error(double approximate, double exact);

class  CIntegrationtest
{
private:
	CMontecarlointegrator _montecarlointegrator;
	CDiophantineintegrator _diophantineintegrator;
	
	vector< CAlpha* > _alpha;
	vector< vector<CWeight*> > _weight;
	vector< vector<CPeriodization*> > _periodization;

	void init(bool includehaselgrovealpha,
		 	  int smoothness, 
			  bool includecyclotomicalpha,
			  bool includebakeralpha,
			  bool rationalexponents, 
			  double maximalexponent,
			  bool includeniederreiteralpha,
			  bool includeprimerootalpha,

			  bool includearithmeticmeanweight,
			  bool includehaselgroveweight,
			  bool includeniederreiterweight,
			  bool includesugiharamurotaweight,

			  bool includenoperiodization,
			  bool includebakerperiodization,
			  bool includehaselgroveperiodization,
			  bool includealgebraicperiodization,
			  bool includetrigonometricperiodization,
			  bool includenormalperiodization);

public:
	CIntegrationtest(bool includehaselgrovealpha = true,
					 int smoothness = 2, 
					 bool includecyclotomicalpha = true,
					 bool includebakeralpha = true,
					 bool rationalexponents = true, 
					 double maximalexponent = 7.0,
					 bool includeniederreiteralpha = true,
					 bool includeprimerootalpha = true,

					 bool includearithmeticmeanweight = true,
					 bool includehaselgroveweight = true,
					 bool includeniederreiterweight = true,
					 bool includesugiharamurotaweight = true,

					 bool includenoperiodization = true,
					 bool includebakerperiodization = true,
					 bool includehaselgroveperiodization = true,
					 bool includealgebraicperiodization = true,
					 bool includetrigonometricperiodization = true,
					 bool includenormalperiodization = true)
	{
		init(includehaselgrovealpha,
        	 smoothness, 
			 includecyclotomicalpha,
			 includebakeralpha,
			 rationalexponents, 
			 maximalexponent,
			 includeniederreiteralpha,
			 includeprimerootalpha,

			 includearithmeticmeanweight,
        	 includehaselgroveweight,
			 includeniederreiterweight,
			 includesugiharamurotaweight,

			 includenoperiodization,
        	 includebakerperiodization,
			 includehaselgroveperiodization,
			 includealgebraicperiodization,
			 includetrigonometricperiodization,
			 includenormalperiodization);

		_alpha[0] = NULL;
		_alpha[1] = NULL;

	}
	CIntegrationtest(EqSpDoubleVector& userdefinedalpha, 
					 bool includehaselgrovealpha = true,
					 int smoothness = 2, 
					 bool includecyclotomicalpha = true,
					 bool includebakeralpha = true,
					 bool rationalexponents = true, 
					 double maximalexponent = 7.0,
					 bool includeniederreiteralpha = true,
					 bool includeprimerootalpha = true,

					 bool includearithmeticmeanweight = true,
					 bool includehaselgroveweight = true,
					 bool includeniederreiterweight = true,
					 bool includesugiharamurotaweight = true,

					 bool includenoperiodization = true,
					 bool includebakerperiodization = true,
					 bool includehaselgroveperiodization = true,
					 bool includealgebraicperiodization = true,
					 bool includetrigonometricperiodization = true,
					 bool includenormalperiodization = true)
	{
		init(includehaselgrovealpha,
        	 smoothness, 
			 includecyclotomicalpha,
			 includebakeralpha,
			 rationalexponents, 
			 maximalexponent,
			 includeniederreiteralpha,
			 includeprimerootalpha,

			 includearithmeticmeanweight,
        	 includehaselgroveweight,
			 includeniederreiterweight,
			 includesugiharamurotaweight,

			 includenoperiodization,
        	 includebakerperiodization,
			 includehaselgroveperiodization,
			 includealgebraicperiodization,
			 includetrigonometricperiodization,
			 includenormalperiodization);

		_alpha[0] = new CUseralpha(userdefinedalpha);
		_alpha[1] = NULL;
	}
	~CIntegrationtest()
	{
		int i, j; 
		for (i = 0; i < numberofalphatypes; i++) 
		{
			delete _alpha[i];
		}
		for (i = 0; i < numberofweighttypes; i++)
		{
			for (j = 0; j < _weight[i].size(); j++)
			{
				delete _weight[i][j];
			}
		}
		for (i = 0; i < numberofperiodizationtypes; i++)
		{
			for (j = 0; j < _periodization[i].size(); j++)
			{
				delete _periodization[i][j];
			}
		}
	}

	void plotalphadiscrepancies(int k, int N);
	
	void compareerrors(CIntegrand& integrand, CDomain& domain, double analyticintegral, 
					   int Nmin, int Nstep, int Nmax, const char* filename,
					   int desiredrateofconvergence, int desireddegreeofsymmetry);

	void perform(CIntegrand& integrand, CDomain& domain, double analyticintegral, 
				 bool plotintegrand = true, 
				 bool tabulateerrors = true, int N = 5000, 
				 bool ploterrorslowrange = true, bool ploterrorsmediumrange = true, 
				 bool ploterrorshighrange = true, int Nmin = 500, int Nstep = 100, int Nmax = 5000,
				 int desiredrateofconvergence = -1, int desireddegreeofsymmetry = -1);
};

//----------------------------------------------------------------------------------------------

#endif
