#include	"stdafx.h"


//#include <eqspmcint/StdAfx.h>
#include "MCInt\INTSTL.h"

#include "INTIntegrator.h"
#include "MCInt\INTFinanceIntegrands.h"
#include "MCInt\INTRng.h"

#include "MCInt\INTTestIntegrands.h"
#include "MCInt\INTTestUtilities.h"
#include "MCInt\INTTestIntegrands.h"
#include "MCInt\INTTestEngine.h"
#include "MCInt\INTTestIntegrator.h"



/*Included to test if the classes are exported */
#include "MCInt\INTSystematic.h"



// Parameters k: dimension, class of option in {1 = Asian, 2 = Basket}, Q: number of instances

void testsystematicQ(int k, int optionclass, int Q, double corr)
{
  std::cout << "Perform systematic test: " << std::endl;
  std::cout << "-------------------------------------------------------------------" << std::endl << std::endl;

  // Instantiate all alphas
  CAlpha *has2a = new CHaselgrovealpha(2);
  CAlpha *has4a = new CHaselgrovealpha(4);
  CAlpha *cyca = new CCyclotomicalpha();
  double m = 3.0;
  CAlpha *baka = new CBakeralpha(m);
  CAlpha *niea = new CNiederreiteralpha();
  CAlpha *pra = new CPrimerootalpha();

  // Instantiate all weights with r in {1, ..., 4}
  int r;
  CWeights* amw = new CArithmeticmeanweights();
   vector<CWeights*> hasw(5, (CWeights*)NULL), niew(5, (CWeights*)NULL), smw(5, (CWeights*)NULL);  
  //vector<CWeights*> hasw(5, (CWeights*)NULL);


  for (r = 1; r <= 4; r++)
    {
      hasw[r] = new CHaselgroveweights(r);
      niew[r] = new CNiederreiterweights(r);
      smw[r] = new CSugiharaMurotaweights(r);
    }

  // Instantiate all periodizations

  // Zero Type
  CPeriodization *idp = new CIdentityperiodization;
  CPeriodization *bakp = new CBakerperiodization;
  
  // Finite Type with d in {0, ..., 5}
  int d;
  vector<CPeriodization*> polyp(6, (CPeriodization*)NULL), trigp(6, (CPeriodization*)NULL), krep(6, (CPeriodization*)NULL);

  for (d = 0; d <= 5; d++)
    {
      polyp[d] = new CPolynomialperiodization(d);
      trigp[d] = new CTrigonometricperiodization(d);
      krep[d] = new CKressperiodization(d);
    }

  // Infinite Type
  double a = 1.0;
  int p = 1;
  double sigma = 5.0;
  CPeriodization *TANHp = new CTANHperiodization(a, p);
  CPeriodization *DEp = new CDEperiodization(a);
  CPeriodization *gaussp = new CGaussianperiodization(sigma);

  // ----------

  // Domain
  CDomain *D = new CUnitcube(k);

  // ----------

  // Generate Q problem instances at random
  vector<CIntegrand*> f(Q);
  STLDoubleVector I(Q);

  int q = 0;
  while (q < Q)
    {
      // Common parameters, see [AcworthBroadieGlassermann96], pp. 8+9
      int type = 1;

      // Interest rate rho ~ U[0, 10%)
      double rho = randomUab(0.0, 0.1);

      // Expiration date T = 1 year
      double T = 1.0;

      if (optionclass == 1)
        {
          // Asian option, see [AcworthBroadieGlassermann96], p. 5

          // Evenly spaced monitoring dates t_i = i*T/m
          STLDoubleVector t;
          equalt(t, k, T);

          // Initial asset price S0 = 100
          double S0 = 100.0;

          // Strike K ~ U[80, 120)
          double K = randomUab(80.0, 120.0);

          // Annual volatility sigma ~ U[10%, 60%)
          double vol = randomUab(0.1, 0.6);
          
          f[q] = new CAsianintegrand(t, type, S0, K, rho, T, vol);
        }
      else if (optionclass == 2)
        {
          // Basket option, see [AcworthBroadieGlassermann96], pp. 8+9

          // Initial asset prices S0^i = 100
          STLDoubleVector S0(k, 100.0);

          // Annual volatilities sigma_i ~ U[10%, 60%)
          STLDoubleVector vol;
          univol(vol, k, 0.1, 0.6);

          // No correlation
          STLDoubleVectorVector C;
          equalcorr(C, k, corr);

          // Strikes K ~ U[70, 110)
          double K = randomUab(70.0, 110.0);

          f[q] = new CBasketintegrand(S0, vol, C, type, K, rho, T);
        }

      // If the analytic option price is below 0.5, the generated instance is discarded
      I[q] = f[q]->I(k);

      if (I[q] < 0.5) 
        {
          delete f[q];
        }
      else 
        {
          std::cout << (q+1) << ". option: " << std::endl; f[q]->report(); std::cout << std::endl << std::endl;
          std::cout << "Analytic integral over [0, 1)^" << k << " is I = " << I[q] << ". " << std::endl << std::endl;
        
          q++;
        }
    }

  // ----------

  // Numbers of sample points
  STLLongVector N;
  int Nsize = 4;
  CreateVectorClass<long>::createvector(N, Nsize, 1250L, 5000L, 20000L, 80000L);

  // ----------

  // Alphas in {has2a, has4a, cyca, baka, niea, pra}
  vector<CAlpha*> alphas;
  int asize = 4;
  CreateVectorClass<CAlpha*>::createvector(alphas, asize, cyca, baka, niea, pra);
  
  // Weights in {amw, hasw[1..4], niew[1..4], smw[1..4]}
  vector<CWeights*> weights;
  int wsize = 3;
  CreateVectorClass<CWeights*>::createvector(weights, wsize, amw, hasw[2], niew[2]);

  // Periodizations in {idp, bakp, polyp[0..5], trigp[0..5], krep[0..5], TANHp, DEp, gaussp}
  vector<CPeriodization*> periodizations;
  int psize = 4;
  CreateVectorClass<CPeriodization*>::createvector(periodizations, psize, idp, bakp, polyp[2], polyp[3]);

  // ----------

  // Categories
  STLIntegerVectorVector cat;
//cat.resize(2);
//createvector(cat[0], 3, 2,0,1);
//createvector(cat[1], 3, 2,1,0);

  // Plot all levels
  STLIntegerVectorVector levels;

  // ----------

  // Perform systematic test
  systematic(f, *D, I, N, alphas, weights, periodizations, cat, levels);

  // ----------

  // Tidy up

  // Domain
  delete D;
  // Integrands
  for (q = 0; q < Q; q++) delete f[q];
  // Alphas
  delete has2a; delete has4a; delete cyca; delete baka; delete niea; delete pra;
  // Weights
  delete amw; 
  for (r = 1; r <= 4; r++) { delete hasw[r]; delete niew[r]; delete smw[r]; }
  // Periodizations
  delete idp; delete bakp; 
  for (d = 0; d <= 5; d++) { delete polyp[d]; delete trigp[d]; delete krep[d]; }
  delete TANHp; delete DEp; delete gaussp;
}



//----------------------------------------------------------------------------------------------



int main(int argc, char* argv[])
{
// Initialise RNG
  //  srand((unsigned)time(NULL));
  //  ran2init(- abs(rand()));
	ran2init(-12345);

  printf("Integration is running\n\n\n");
  fflush(stdout);




  // Test string operation.

  string str = i2s(10);
  std::cout<<str<<":"<<std::endl;
  
  long llll = 10;
  str = l2s(llll);
  std::cout<<str<<":"<<std::endl;

  str = d2s(10.12345678901234567890);
  std::cout<<str<<":"<<std::endl;






  // Utilities
  //testutilities();

  // Integrands
  //testdomains();
  //testmiscintegrands();
  //testfinanceutilities();
  //testfinanceintegrands();
  
  // Engine
  //testalphas();
  //testweights();
  //testperiodizations();
  //testintegrator();

  
  int Q = 250;

  int asian = 1;

//  testsystematicQ(10, asian, Q, 0.0);

  testsystematic(10, asian, Q);
  testsystematic(50, asian, Q);
  testsystematic(100, asian, Q);

  int basket = 2;
  double corr = 0.0;

  testsystematic(5, basket, Q, corr);
  testsystematic(50, basket, Q, corr);
  testsystematic(100, basket, Q, corr);

  corr = 0.3;

  testsystematic(10, basket, Q, corr);
  testsystematic(50, basket, Q, corr);
  testsystematic(100, basket, Q, corr);
  
  return 0;
  
}
        
//----------------------------------------------------------------------------------------------

