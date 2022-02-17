#include "StdAfx.h"
#include "Rnd.h"
#include "MlEqMaths.h"




void test()
{
		int ndim = 5;
		int npaths = 1000;


		CUnitcube* D = new CUnitcube(ndim);
	
		CArithmeticmeanweights* weights = new CArithmeticmeanweights;// weight
		CAlpha* alpha = new CPrimerootalpha;
		CBakerperiodization* periodization=new CBakerperiodization;
		CNTparameters* par = new CNTparameters(alpha,weights,periodization);
		CNTintegrator* NTintegrator= new CNTintegrator;

		Diophantine dio(ndim,
				npaths,
				alpha ,	
				weights,
				D,
				periodization,
				NTintegrator,
				par);



}



