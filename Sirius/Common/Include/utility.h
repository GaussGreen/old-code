
#ifndef UTILITYH
#define UTILITYH

#include "cMatrix.h"

#include "MlEqHandles.h"
#include "solve.h"

using namespace std;

class MlEqZeroCurve;

enum Event{
	EarlyExercise				= 1,
	CashFlow					= 2,
	BarrierEvent				= 3
	};


void trim(std::string* psz);




// This nice little class is used to set a given variable a value on
// class construction and set it to another value on destruction.
// I use it in MlEqVolatilityStructure::Scenario_ShiftSkew::Apply
// as a recursion blocker.
template<class T> class flip_flop
{
public:
	explicit flip_flop(T& variable, const T& start, const T& end) : m_variable(variable)
	{		
		m_variable = start;
		m_end = end;
	}

	virtual ~flip_flop(void)
	{
		m_variable = m_end;
	}

protected:
	T&	m_variable;
	T	m_end;
};



void ConvertConstantNotionalRamStrike(NotionalTypeEnum* pType, long nStartDate, long nValuationDate, double* pfNotional, double* pfStrike);

class curve
{
	public:

	int m_method;

	CVector m_xVals;
	CVector m_yVals;

	public:

	curve(CVector& xArray,CVector& yArray);
	curve(){};
	int initialize(CVector& xArray,CVector& yArray);
	curve& operator=(const curve& rhs);

	double getYValue(double x);
	double getXValue(int i);
	int getsize();
};


//double Bs(double forward,double vol,double maturity,double strike,double discount_factor,int cp);
double Bs(double forward,double vol,double maturity,double strike,double discount_factor,int cp,double blend=0.0);
double Bs(double forward,double vol,double maturity,double strike,double discount_factor,int cp,double blend,double volvol,int nsig =4,int npoints =12);

void double_knock_out_rebate(double & price,double mat,double spot,double forward,
						  double discount_factor,double vol,double blend,double barrier_up,
						  double barrier_do,double rebate);
void double_knock_out_rebate(double & price,double mat,double spot,double forward,
						  double discount_factor,double vol,double barrier_up,
						  double barrier_do,double rebate);
void double_knock_out_rebate(double & price,double mat,double spot,double forward,
						  double discount_factor,double vol,double blend,double barrier_up,
						  double barrier_do,double rebate,double volvol,int nsig /*=4*/,int npoints /*=12*/);


double MlEqBSImpliedVol(
						double option_val,		// option value
						double forward,			// forward price
						double maturity,		// maturity
						double strike,			// strike
						double discount_factor,	// discount factor
						int cp,					// 1 : call, -1 : put, 0 : forward
						int rootfind_flag=8,		// rootfinder flag 
						double accuracy=1e-5,		// accuracy (optional) 
						double lower_x=0.0,			// lower bound to search (optional) 
						double upper_x=1.5			// upper bound to search (optional) 
);





void STLVectorFromCVector(vector < double > & stlVec,const CVector & in);
void CVectorFromSTLVector(CVector& vec, const vector < double > & in);

void STLVectorVectorFromCMatrix(vector < vector < double > > & stlVec,const CMatrix& in );
void STLVectorVectorTransposeFromCMatrix(vector < vector < double > > & stlVec,const CMatrix& in );
void CMatrixFromSTLVectorVector(CMatrix& mat, const vector < vector < double > > & in );


double bridge(double delta_t,double U,double L,double X_i,double X_f,double vol);

///////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////

class MlEqMonteCarlo;

class product
{	
	protected:

	virtual void fillValues(CMatrix& value,CMatrix& path_array,MlEqMonteCarlo& mc);

	public:

//	CMatrix m_spotfixings;//[iasset][isate]

	GVector<long> m_payoffDates;
	int m_numberPayouts;
	CMatrix m_results;

	virtual void payout(CMatrix& value,CMatrix& pathArray,CVector& discountFactors,MlEqMonteCarlo& m){};
	virtual void setUp(CMatrix& value,MlEqMonteCarlo& mc);
	
	virtual ~product(){};
	product(){};

	virtual void createResult(CMatrix& result,MlEqMonteCarlo& m);
	virtual void accruePayout(CMatrix& result,CMatrix& value);

};	

double getObjectFromCMatrix(CMatrix& objects,int iRow,int iCol,int numberOfSlices);

template<class T> const T& getObjectFromVectorVector(vector <vector < T > >& objects,int iRow,int iCol,int numberOfSlices)
{
	int nsize = objects.size();
	if ( nsize == 1 )
	{
		if ( objects[0].size() <= iCol ){
			throw("error indexing Column array");
		}

		return objects[0][iCol];
	}
	else if ( nsize == numberOfSlices )
	{
		if ( objects[iRow].size() <= iCol ){
			throw("error indexing strikearray");
		}

		return objects[iRow][iCol];
	}
	else{
		throw("indexing error ");
	}

};



double  getObjectFromCVector(CVector&  objects,int index);

template<class T>  T& getObjectFromVector(vector < T > & objects,int index)
{
	int nsize = objects.size();
	if ( nsize == 1 )
	{
/*		if ( objects.size() <= index ){
			throw("error indexing Column array");
		}
*/

		return objects[0];
	}
	else 
	{
		if ( objects.size() <= index ){
			throw("error indexing strikearray");
		}

		return objects[index];
	}

};

template<class T> const T& getObjectFromVector(const GVector < T > & objects,int index)
{
	int nsize = objects.getsize();
	if ( nsize == 1 ){
		return objects[0];
	}
	else {
		if ( objects.getsize() <= index ){
			throw("error indexing strikearray");
		}

		return objects[index];
	}

};





double Hermit(int n, double x, double y );
void HermiteOptions(CVector& result, CVector& strikes,const  CVector& hermitCoeff,double maturity,double forward,int returnVolFlag=1,double elasticity=0,double volvol=0,int ngauss=0);



class hermiteSkew
{

	friend void hermitNormalization(double x, double z[], void *vp);
	friend void hermitOptions(double x, double z[], void *vp);

	double	m_sqrtMat;
	double	m_sqrt_2pi;	
	double	m_normalization;
	CVector m_strikes;
	CVector m_coeff;
	int		m_iter;
	double  m_factor;
	double  m_logNormalization;

	double	newtonUpdate(double hermit1,double forward);

public:

	double m_maturity;
	double m_forward;

	CVector m_hermitCoeff;

	hermiteSkew(){m_iter=0;m_factor=1;};
	void initialize(double maturity,double forward,const CVector& hermitCoeff,bool adjustATM=true);
	void calculateOptions(CVector& result, CVector& strikes,bool returnVolFlag=true);
	double adjustSlope(double targetSlope,double forward,double eps);

	double getStateVar(double x );

};







  

#endif
