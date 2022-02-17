// pgcd.cpp: implementation of the pgcd class.
//
//////////////////////////////////////////////////////////////////////

#include "ICMKernel\util\icm_pgcd.h"

#include "ICMKernel\util\icm_utils.h"

// #include "Quantifyer.h"

pgcd::pgcd()
{
}
//-------------------------------------------------------------------------------
pgcd::~pgcd()
{
}
//-------------------------------------------------------------------------------
long
pgcd::solve ( long a, long b)
{
	while (b != 0) 
	{ 
		unsigned long c = a % b; 
		a = b; 
		b = c; 
	} 
	return a; 
}
 
long
pgcd::solve (const std::vector<long>& a)

{
	if (a.size()==0) ICMTHROW(ERR_INVALID_ARGUMENT,"ite_gcd: empty argument"); 
	if (a.size()==1) return a[0]; 
	std::vector<long>::const_iterator it = a.begin(); 
	unsigned long pgcd = solve(a[0],a[1]); 
	int n = a.size() + 2; 
	while (pgcd!=1 && it<a.end())
	{
		pgcd = solve(pgcd,*it); 
		++it; 
	}
	return pgcd ;
}

//-------------------------------------------------------------------------------
/*long
pgcd::solve (const std::vector<long>& a)
{
//EP_ASSERT(a.size()>1)
		if(a.size()==1) return a[0]; 
		if(a.size()==2) return solve(a[0], a[1]);
		else {
			// std::vector<long> b;
			// for(int k=1;k<a.size();k++) b.push_back(solve(a[0], a[k]));
			std::vector<long> b(a.size()-1);
			std::vector<long>::iterator bit=b.begin(); 
			for(int k=1;k<a.size();k++) { *bit=solve(a[0], a[k]) ; ++bit; } 
			return solve(b);
		}
}
*/

double pgcd::round(double value)
{
// 	QUANTIFYER("pgcd::rounded"); 
	double result = 0.;

	double arrondi = floor(value);
	result = arrondi;

	if (fabs(arrondi - value)>0)
	{
		if (fabs(arrondi-value) <  fabs(arrondi+1.-value)) 
		{
			if (fabs(arrondi-value) <  fabs(arrondi-1.-value)) 
				result = arrondi;
			else
				result = arrondi -1.;

		}else if (fabs(arrondi-value) <  fabs(arrondi-1. -value)) 
		{
			if (fabs(arrondi-value) <  fabs(arrondi+1. -value)) 
				result = arrondi;
			else
				result = arrondi +1.;

		}
		else if (fabs(arrondi-1. -value) <  fabs(arrondi+1.-value)) 
		{
			if (fabs(arrondi-1. -value) <  fabs(arrondi-value)) 
				result = arrondi - 1.;
			else
				result = arrondi;
		}
		 else if (fabs(arrondi-1. -value) <  fabs(arrondi-value)) 
		{
			 if (fabs(arrondi-1. -value) <  fabs(arrondi+1.-value)) 
				result = arrondi - 1.;
			 else	
			 	result = arrondi + 1.;
		}
		else if (fabs(arrondi+1. -value) <  fabs(arrondi-value)) 
		{
			if (fabs(arrondi+1. -value) <  fabs(arrondi-1. -value)) 
				result = arrondi + 1.;
			else
				result = arrondi -1.;
		}
		 else if (fabs(arrondi+1. -value) <  fabs(arrondi-1.-value)) 
		{
			 if (fabs(arrondi+1. -value) <  fabs(arrondi-value)) 
				result = arrondi + 1.;
			 else	
			 	result = arrondi;
		}
			
	}

	return (result);
}


//-------------------------------------------------------------------------------

double
pgcd::solve (const ARM_Vector& a )
{
	if (a.size()<=1)
	throw Exception(__LINE__, __FILE__, ERR_UNIMP_METHOD_CALL,
           "solve error !");

	std::vector<long> LR_long(a.size());
	// Determine norm
	double norm = 1.;
	for(int k=0;k<a.size();k++) 
	{
		double tmp = pgcd::round(a.Elt(k));
		LR_long[k] = tmp;
	}
	return solve(LR_long) ;
}  

//------------------------ for Credit_Manager Must be deleted after Credit_Manager
double
pgcd::solve (const std::vector<double>& a )
{
	if (a.size()<=1)
	throw Exception(__LINE__, __FILE__, ERR_UNIMP_METHOD_CALL,
           "solve error !");

	std::vector<long> LR_long(a.size());
	// Determine norm
	double norm = 1.;
	for(int k=0;k<a.size();k++) 
	{
		double tmp = pgcd::round(a[k]);
		LR_long[k] = tmp;
	}
	return solve(LR_long) ;
}  
//-------------------------------------------------------------------------------
bool
pgcd::test(std::string& errStr)
{
	bool ret(true);
	
	errStr = "Testing Euclide Algorithm";
/*
	long res = pgcd::solve (5,15);
	SDMTEST(res==5);
	res = pgcd::solve (5,13);
	SDMTEST(res==1);
	res = pgcd::solve (400,26000);
	SDMTEST(res==400);

	std::vector <long> a(5);
	a[0] = 100; a[1] = 200; a[2] = 50; a[3] = 150; a[4] = 5050;
	res = solve(a);
	SDMTEST(res==50);

	std::vector <double> b(5);
	b[0] = 0.1; b[1] = 0.2; b[2] = 0.05; b[3] = 0.150; b[4] = 5.050;
	double r = solve(b);
	SDMTEST(r==0.05);
*/
	return ret;
}
//-------------------------------------------------------------------------------
