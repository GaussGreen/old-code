#include	"stdafx.h"



#include "MCInt\INTAlpha.h"
#include "MCInt\INTWeights.h"
#include "MCInt\INTPeriodization.h"

#include "MCInt\INTTestEngine.h"

//----------------------------------------------------------------------------------------------

void testalphas(void)
{
  std::cout << "Alphas: " << std::endl;
  std::cout << "-------------------------------------------------------------------" << std::endl << std::endl;
  
  int k = 2;
  STLDoubleVector alpha;
  CreateVectorClass<double>::createvector(alpha, k, 0.234908349, 0.634534342);
  CUseralpha usera(alpha);

  usera.report(k); std::cout << std::endl;
  std::cout << std::endl;
  int N = 100;
  usera.plotdistribution(k, N);

  // ----------

  int s = 2;
  CHaselgrovealpha hasa(2);

  hasa.report(k); std::cout << std::endl;
  std::cout << std::endl;
  hasa.plotdistribution(k, N);

  // ----------

  CCyclotomicalpha cyca;

  cyca.report(k); std::cout << std::endl;
  std::cout << std::endl;
  cyca.plotdistribution(k, N);

  // ----------

  double m = 3.0;
  CBakeralpha baka(m);

  baka.report(k); std::cout << std::endl;
  std::cout << std::endl;
  baka.plotdistribution(k, N);

  // ----------

  CNiederreiteralpha niea;

  niea.report(k); std::cout << std::endl;
  std::cout << std::endl;
  niea.plotdistribution(k, N);

  // ----------

  CPrimerootalpha pra;

  pra.report(k); std::cout << std::endl;
  std::cout << std::endl;
  pra.plotdistribution(k, N);
}

//----------------------------------------------------------------------------------------------

void testweights(void)
{
  std::cout << "Weights: " << std::endl;
  std::cout << "-------------------------------------------------------------------" << std::endl << std::endl;

  CArithmeticmeanweights amw;

  amw.report(); std::cout << std::endl;
  std::cout << std::endl;
  int N = 100;
  amw.plotdistribution(N);

  // ----------

  int r = 3;
  CHaselgroveweights hasw(r);

  hasw.report(); std::cout << std::endl;
  std::cout << std::endl;
  hasw.plotdistribution(N);

  // ----------

  CNiederreiterweights niew(r);

  niew.report(); std::cout << std::endl;
  std::cout << std::endl;
  niew.plotdistribution(N);

  // ----------

  CSugiharaMurotaweights smw(r);

  smw.report(); std::cout << std::endl;
  std::cout << std::endl;
  smw.plotdistribution(N);
}

//----------------------------------------------------------------------------------------------

void testperiodizations(void)
{
  std::cout << "Periodizations: " << std::endl;
  std::cout << "-------------------------------------------------------------------" << std::endl << std::endl;

  // Plot both phi(x) and phi1(x)
  int choice = 2;

  CIdentityperiodization idp;

  idp.report(); std::cout << std::endl;
  std::cout << std::endl;
  idp.plot(choice);

  // ----------

  CBakerperiodization bakp;

  bakp.report(); std::cout << std::endl;
  std::cout << std::endl;
  bakp.plot(choice);

  // ----------

  int d = 3;
  CPolynomialperiodization polyp(d);

  polyp.report(); std::cout << std::endl;
  std::cout << std::endl;
  polyp.plot(choice);

  // ----------

  CTrigonometricperiodization trigp(d);

  trigp.report(); std::cout << std::endl;
  std::cout << std::endl;
  trigp.plot(choice);

  // ----------

  CKressperiodization kp(d);

  kp.report(); std::cout << std::endl;
  std::cout << std::endl;
  kp.plot(choice);

  // ----------

  double a = 1.0;
  int p = 1;
  CTANHperiodization TANHp(a, p);

  TANHp.report(); std::cout << std::endl;
  std::cout << std::endl;
  TANHp.plot(choice);

  // ----------

  a = 1.0;
  CDEperiodization DEp(a);

  DEp.report(); std::cout << std::endl;
  std::cout << std::endl;
  DEp.plot(choice);

  // ----------

  double sigma = 5.0;
  CGaussianperiodization gaussp(sigma);

  gaussp.report(); std::cout << std::endl;
  std::cout << std::endl;
  gaussp.plot(choice);
}

//----------------------------------------------------------------------------------------------

