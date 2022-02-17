

#ifndef MLEQMATHSH
#define MLEQMATHSH


#include "cMatrix.h"
#include "solve.h"



class CVector;

class  MlEqMaths
{

protected:

	static double factln(int n);
	static double gammln(double xx);


 public:
  static void rkqs(double y[], double dydx[], int n, double *x, double htry, double eps,
	    double yscal[], double *hdid, double *hnext,
	    void (*derivs)(double, double [], double[], void *),
	    void *vp);

  static double gasdev(long *idum);
  static double gasdev1(long *idum,double lo,double lup);
  static void rkck(double y[], double dydx[], int n, double x, double h,
		   double yout[], double yerr[], 
		   void (*derivs)(double, double[], double[], void *),
		   void *vp);


 static void svdcmp(CMatrix& a, CVector& w, CMatrix& v);
 static void odeint(CVector& ystart,int nvar,double x1,double x2,double eps,double h1,double hmin,int*nok,int* nbad,
									void (*derivs)(double,double[],double[],void *),
									void (*rkqs)(double[],double[],int,double*,double,double,double[],double*,double*,void (*)(double,double[],double[],void *), void *),
									double* xp,double** yp,int& kount,void *vp,const int kmax,double dxsav,double scale = 10.0);

  static double Min  (const double& a, const double& b);
  static double Max  (const double& a, const double& b);
  static int Min  (const int& a, const int& b);
  static int Max  (const int& a, const int& b);

  static bool   equal(const double& a, const double& b);

  // searches array xx to find i such that xx[ii] <= x <= xx[i+1]	
  // return -1 if x < xx[0] returns n-1 if x > xx[n-1]
  static void locate(const CVector& xx,double x,int& i);
  static void locate(const GVector<long>& xx,long x,int& i); // searches array xx to find i such that xx[ii] <= x <= xx[i+1]	
  static void Locate(const GVector<long>& xx,long x,int& i); // searches array xx to find i such that xx[ii] <= x <= xx[i+1]	
  static void locate(int ndim, const CVector& xx,double x,int& i);
  static void locate(const CVector& xx,double x,int n,int& i); // searches array xx to find i such that xx[ii] <= x <= xx[i+1]	takes only n elements out of array


  static double linearInterp(const CVector& x, const CVector& y,double xx,int method=0);
  static double linearInterp(int ndim,CVector& x,CVector& y,double xx,int method);
  static double linearInterp(const CVector& x, const CVector& y,double xx,int method,bool& isEdge);

  
  static double constantInterp(const CVector& x, const CVector& y,double xx);
  static double constantInterp(const CVector& x, const CVector& y,double xx,bool & isEdge);


  static double math_factorial(int x, int y);
  static double binomial(int n,int m);

  static double Hermite(int n,double x,double y);
  static double Hermite(int n,double x,double y,double prevVal,double prevPrevVal);

  static int math_cubic_roots(double a, double b, double c, double d, double* roots);

  static int iround(double x, int how);
  static double dround(double x, int how);
  static int is_integer(double x);

  static double dsgn(double x);
  static double dsgne(double x);
  static double dsqr(double x);
  static double dcube(double x);
  static double dmax(double x, double y);
  static double dmin(double x, double y);
  static int isgn(double x);
  static int isgne(double x);
  static int imax(int x, int y);
  static int imin(int x, int y);
 
  static void Max(double& current_max,int& i,CVector& exclude,CVector& x);
  static void Min(double& current_min,int& i,CVector& exclude,CVector& x);
  static double Min(double x,double y);
  static double Max(double x,double y);

  static void ludcmp(CMatrix& a, GVector<int>& indx, double& d);
  static void lubksb(CMatrix& a, GVector<int>& indx, CVector& b);

//  static void Hpsort(CVector& x,GVector< long >& map);
  static void Hpsort(CVector& x);
  static bool deqz(double x);

  static double MSIGN(double a,double b);
  static double ncchisq(double chisq, double lambda, double v, double rel_accuracy=1e-6);
  static void dGauleg( double x1, double x2, CVector& x, CVector& w, int n, bool normalseWeights = false );
  static void choldc(CMatrix& a,CVector& p);
//  static void choldc(GMatrix < CMatrix >& a,GVector < CVector > & p);

  static double RoundDouble(double doValue, int nPrecision);	
  static const int MATH_ROUND_NONE ;               // 0;//		"No rounding"
  static const int MATH_ROUND_DOWN ;               // 1;//		"Round down towards negative infinity"
  static const int MATH_ROUND_UP;                  // 2;//		"Round up towards positive infinity"
  static const int MATH_ROUND_TO_ZERO;             // 3;//		"Round towards zero"
  static const int MATH_ROUND_FROM_ZERO;	   // 4;//		"Round away from zero"
  static const int MATH_ROUND_NEAREST;	     	   // 5;//		"Round to the nearest integer"

};

bool invertMatrix(CMatrix& inverse,const CMatrix& in);

class fitPolynomial: public LSGRGSolver
{
	CMatrix m_data;// [idata][2]
	
public:

	int m_degree;

	void  initialize(CVector& initialGuess,	CMatrix& data,
					   CMatrix& xValBounds,CMatrix& ObjectiveBounds,
					   double InitialTolerance,double FinalTolerance,
					   double StoppingTolerance,
					   CMatrix& NonZeroPartialDerivativesSpecification,int outputFlag);


//	void  PartialDerivatives(CVector& xVal,CVector& partialDerivatives){};
	void  ObjectiveFcn(double* gVals,double* xVals);

};


// root finding stuff
enum enum_rootfind_method
{
  eROOTFIND_BISECTION=1,
  eROOTFIND_SECANT=2,
  eROOTFIND_FALSE_POSITION=3,
  eROOTFIND_RIDDERS=4,
  eROOTFIND_BRENT=5,
  eROOTFIND_CBLFR=6,
  eROOTFIND_NEWTON=7,
  eROOTFIND_BRENT_GROWRANGE=8,
  eROOTFIND_BRENT_NOGROW=9
};


 int rootfind_solve(int rootfind_flag, bool (*fn)(double x,void* vp,double* f), 
										   double lower_x, double upper_x,
										   double accuracy, int max_tries, 
										   double find_y, void* vp, double* found_x,
										   int (*dfndx)(double x,void* vp,double* f) );




class ARTrackingPursue
{

// Adil Reghai tracking pursue algorithm

	double scalarProduct(CVector&x, CVector& y);
	double scalarProduct(int j, CVector& x);

	void multiplyVec(CVector& x,double f);
	double variance(CVector& x);

	CVector m_targetVec;//[ndim]
	CMatrix m_projections;//[ndim][nvec]

	double m_epsilon;
	double m_stoppingTolerance;
	int m_maxiter;
	
	CVector m_stdev;//[nvec]
	double m_stdevx;

public:

	CVector m_weights;//[nvec]	
	void project();
	void initialize(CVector& targetVec,CMatrix& projections,double epsilon,double varReduction,int maxiter=100);

};

void createCholesky(CMatrix &cholesky,CMatrix& correl);
//void createCholesky(GVector<CMatrix> &cholesky,CMatrix& correl);

class Cholesky
{
	CMatrix m_cholesky;

public:

	int m_numberOfPeriods;
	int m_numberOfFactors;


	void initialize(CMatrix& correl,int numberOfPeriods);
	Cholesky(CMatrix& correl,int numberOfPeriods);
	Cholesky(){};

	double getCholesky(int isclice,int ifactor,int jfactor);
};


void merge( std::vector<long>& outPutSet,std::vector<int>& mapping, const std::vector<long>& inPutArray_a,const std::vector<long>& inPutArray_b );
void merge( std::vector<long>& outPutSet, const std::vector<long>& inPutArray_a,const std::vector<long>& inPutArray_b );
void merge( CVector& outPutSet, CVector& inPutArray_a,CVector& inPutArray_b ,double tolerance); 
void merge( std::vector<double>& outPutSet, const std::vector<double>& inPutArray_a,const std::vector<double>& inPutArray_b ,double tolerance);
void merge( GVector<long>& outPutSet, GVector<long>& inPutArray_a,GVector<long>& inPutArray_b );



class RKI_class {
public:
	void (*func)(double, double[], void *);
	void *data;
};


void Runge_Kutta_Integrate_new_vector(CVector& integral, void (*fcn)(double, double[], void *p), const CVector& limits, double h_min, double eps, double bump, void *vp, double h_init, double scale);
void Runge_Kutta_Integrate_new(double& integral,void (*fcn)(double, double[], void *p),const CVector& limits,double h_min,double eps,double bump, void *vp, double h_init ,double scale);
void Runge_Kutta_Integrate_buckets(CVector& integralBuckets,void (*fcn)(double, double[], void *p),const CVector& limits,double h_min,double eps,double bump, void *vp, double h_init ,double scale);



#endif 
