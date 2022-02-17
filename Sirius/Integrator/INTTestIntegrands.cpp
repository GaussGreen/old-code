#include	"stdafx.h"



#include "MCInt\INTMiscIntegrands.h"
#include "MCInt\INTFinanceIntegrands.h"

#include "MCInt\INTTestIntegrands.h"

//----------------------------------------------------------------------------------------------

void testdomains(void)
{
	int k;
  double res = 10;

  std::cout << "Domains: " << std::endl;
  std::cout << "-------------------------------------------------------------------" << std::endl << std::endl;

  // ----------

  k = 4;
  STLDoubleVectorVector M;
  CreateVectorClass<double>::creatematrix(M, 4, 7.0,2.0,1.0,5.0, 3.0,1.0,2.0,1.0, 1.0,6.0,1.0,1.0, 2.0,1.0,2.0,6.0);
  CParallelepiped P(k, M);
  P.report();
  std::cout << std::endl << std::endl;

  // ----------

  k = 3;
  CGeneralcube GC(k, 0.0,1.0,-0.5,0.5,-0.3,-0.1);
  GC.report();
  std::cout << std::endl << std::endl;

  // ----------

  k = 5;
  CCube C(k, 2,4);
  C.report();
  std::cout << std::endl << std::endl;

  // ----------

  k = 5;
  CUnitcube UC(k);
  UC.report();
  std::cout << std::endl << std::endl;
}

//----------------------------------------------------------------------------------------------

void testfinanceutilities(void)
{
  std::cout << "Finance utilities: " << std::endl;
  std::cout << "-------------------------------------------------------------------" << std::endl << std::endl;

  std::cout << "Black-Scholes formula for options: ";
  std::cout << "BS(call, 100, 100, 0, 20, 0.2) = " << BS(1, 100, 100, 0, 20, 0.2) << std::endl;
  std::cout << std::endl;

  // ----------

  std::cout << "Input utilities: " << std::endl << std::endl;

  // ----------

  std::cout << "Time steps t: 0 = t_0 <= t_1 <= ... <= t_k = T: " << std::endl;

  int k = 5;
  double T = 1.0;
  std::cout << "T = " << T << std::endl << std::endl;

  STLDoubleVector t;
  std::cout << "Value(s) generated uniformly at random: "; 
  std::cout << "t = " << unit(t, k, T) << std::endl;
  std::cout << "Equal value in all dimensions: "; 
  std::cout << "t = " << equalt(t, k, T) << std::endl;
  CreateVectorClass<double>::createvector(t, k-1, 0.6, 0.2, 0.3, 0.4);
  std::cout << "User-specified values: "; 
  std::cout << "t = " << usert(t, k, T) << std::endl;
  std::cout << std::endl;

  // ----------

  std::cout << "Volatilities vol: (vol_1, ..., vol_k) in [vol_min, vol_max)^k: " << std::endl;

  k = 3;
  double vol_min = 0.2;
  double vol_max = 0.5;
  std::cout << "[vol_min, vol_max) = [" << vol_min << ", " << vol_max << ")" << std::endl << std::endl;

  STLDoubleVector vol;
  std::cout << "Value(s) generated uniformly at random: "; 
  std::cout << "vol = " << univol(vol, k, vol_min, vol_max) << std::endl;
  std::cout << "Equal value in all dimensions: "; 
  std::cout << "vol = " << equalvol(vol, k, 0.3) << std::endl;
  CreateVectorClass<double>::createvector(vol, k, 0.6, 0.2, 0.3);
  std::cout << "User-specified values: "; 
  std::cout << "vol = " << uservol(vol, k) << std::endl;
  std::cout << std::endl;

  // ----------

  std::cout << "Correlations C: kxk-matrix with C_ii = 1, C_ij = C_ji in [corr_min, corr_max): " << std::endl;

  k = 4;
  double corr_min = 0.05;
  double corr_max = 0.15;
  std::cout << "[corr_min, corr_max) = [" << corr_min << ", " << corr_max << ")" << std::endl << std::endl;

  STLDoubleVectorVector C;
  std::cout << "Value(s) generated uniformly at random: " << std::endl;
  std::cout << "C = " << unicorr(C, k, corr_min, corr_max) << std::endl;
  std::cout << "Equal value in all dimensions: " << std::endl;
  std::cout << "C = " << equalcorr(C, k, 0.3) << std::endl;
  std::cout << "User-specified values: " << std::endl;
  std::cout << "C = " << usercorr(C, k, 0.2, 0.1, 0.3, 0.0, 0.2, 0.3) << std::endl;
  std::cout << std::endl;

  // ----------

  std::cout << "Stock prices S: (S_1, ..., S_k) in [S_min, S_max)^k: " << std::endl;

  k = 4;
  double S_min = 50.0;
  double S_max = 200.0;
  std::cout << "[S_min, vol_max) = [" << S_min << ", " << S_max << ")" << std::endl << std::endl;

  STLDoubleVector S;
  std::cout << "Value(s) generated uniformly at random: "; 
  std::cout << "S0 = " << uniS(S, k, S_min, S_max) << std::endl;
  std::cout << "Equal value in all dimensions: "; 
  std::cout << "S0 = " << equalS(S, k, 100.0) << std::endl;
  CreateVectorClass<double>::createvector(S, k, 70.0, 150.0, 30.0, 90.0);
  std::cout << "User-specified values: "; 
  std::cout << "S0 = " << userS(S, k) << std::endl;
  std::cout << std::endl;
}

//----------------------------------------------------------------------------------------------

void testfinanceintegrands(void)
{
  std::cout << "Finance integrands: " << std::endl;
  std::cout << "-------------------------------------------------------------------" << std::endl << std::endl;

  int k = 2;
  
  int res = 20;
  int margin = 1;
  CUnitcube unitcube(k);

  int type = 1;
  double K = 110.0, rho = 0.0, T = 1.0;
  double s0 = 100.0, vol = 0.2;

  // ----------

  STLDoubleVector t;
  unit(t, k, T);

  CEuropeanintegrand eoi(t, type, s0, K, rho, T, vol);

  eoi.report(); std::cout << std::endl;
  std::cout << std::endl;
  eoi.plot(unitcube, "", res, margin, true);

  // ----------

  unit(t, k, T);

  CAsianintegrand aoi(t, type, s0, K, rho, T, vol);

  aoi.report(); std::cout << std::endl;
  std::cout << std::endl;
  aoi.plot(unitcube, "", res, margin, true);

  // ----------

  STLDoubleVector S0;
  double S0_min = 50.0, S0_max = 200.0;
  uniS(S0, k, S0_min, S0_max);

  STLDoubleVector Vol;
  double vol_min = 0.1, vol_max = 0.3;
  univol(Vol, k, vol_min, vol_max);

  STLDoubleVectorVector C;
  double corr_min = 0.05, corr_max = 0.1;
  unicorr(C, k, corr_min, corr_max);

  CBasketintegrand boi(S0, Vol, C, type, K, rho, T);
  
  boi.report(); std::cout << std::endl;
  std::cout << std::endl;
  boi.plot(unitcube, "", res, margin, true);
}

//----------------------------------------------------------------------------------------------

