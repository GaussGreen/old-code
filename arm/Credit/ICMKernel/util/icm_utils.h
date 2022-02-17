#ifndef __ICM_UTILS_H__
#define __ICM_UTILS_H__

#include "ARMKernel\glob\linalg.h"
#include "ARMKernel\crv\zerocurv.h"
#include "ICMKernel\glob\icm_enums.h"
#include "ICMKernel\util\icm_pgcd.h"
#include <vector>
#include <limits>
#include <math.h>


class ARM_VolCurve; 
class ARM_VolLInterpol;
class ORDER_
{
	public :
	double id ;

	public :
		ORDER_(double d) : id(d) {}

		bool operator < (const ORDER_ & rhs) const 
		{ return ((id < rhs.id) && (fabs(id - rhs.id)>DB_TOL));	}
} ;


// void FreePointerTabChar(char**& Tab,const int& size);	//Free an array of char** with size size
// char** CopyTabString(char** Tab,const int& size);		//Allocate a char** of dim "size" and copy Tab in this new array 
// int FindRowInTab(char* chaine,  char** Tab,const int& size);		//Return the index of chaine in Tab, if not found return CREDIT_DEFAULT_VALUE

// rounded of a double
// static inline int round(double x) // { return (int) (x+0.5); }
// { return pgcd::round(x); }
static inline double round(double x) { 
	return	pgcd::round(x);
}

// rounding to N digits. 
//		round(2.15,1) = 2.2 
//		round(2.149,1) = 2.1
//		round(-1.475,2) = -1.48	
static inline double round(double x,unsigned dig)
{
	int t=1; for(int i=0;i<dig;i++) t*=10; 
	if (x<0) return ceil(x*t-0.5)/t ;
	return floor(x*t+0.5)/t ;
}

//Merge Two std:vector
void MergeVector(const vector<double>& v1,
				 const vector<double>& v2,
				 vector<double>& vout);	


class ARM_IRIndex;
/** void ComputeFwdDates(ARM_Vector* flowStartDates, 
					 ARM_Vector* flowEndDates,
                     ARM_Vector* paymentDates, 
					 ARM_Vector* resetDates,
					 ARM_IRIndex* irIndex,
					 ARM_Vector*& FwdStartDates,
					 ARM_Vector*& FwdEndDates);

**/ 
//	------------------------------------------------------------------------
//	CountYears using const Arguments 
static inline double CountYears(int DayCount,const ARM_Date& date1,const ARM_Date&date2)
{ 
	return CountYears(DayCount,const_cast<ARM_Date&>(date1),const_cast<ARM_Date&>(date2)); 
}
//	------------------------------------------------------------------------
//	CountYears using const Arguments 
static inline double CountYears(int DayCount,int nbDays, const ARM_Date& date1,const ARM_Date&date2)
{ 
	switch (DayCount) 
    {
        case KNOBASE: return 1.; 
        case K30_360 : case KACTUAL_360 : return nbDays/360.; 
        case KACTUAL_365 : return nbDays/365.;  
	} ;
	// otherwise: call the initial ARM method (recomputes nbDays).
	return CountYears(DayCount,date1,date2); 
	// ICMTHROW(ERR_INVALID_ARGUMENT,"CountYears not available for this base ="<<DayCount); 
}
//	------------------------------------------------------------------------
//	DaysBetweenDates using  const Arguments 
static inline double DaysBetweenDates(int DayCount,const ARM_Date& date1,const ARM_Date&date2)
{ 
	return DaysBetweenDates(DayCount,unconst(date1),unconst(date2)); 
}

//	securing clone
//	
//		bypass the "NON CONST" clone for ARM_Object
//		downcast to same type as initial object
//
//		ARM_DefaultCurve* item = dyn_clone(myConcreteDefCurvePointer ) ; 
//
template <class T> static inline  T* dyn_clone(const T* item) 
{
	if (!item) ICMTHROW(ERR_INVALID_ARGUMENT,"Can't clone: null argument"); 
	return dynamic_cast<T*>( const_cast<T*>(item) ->Clone() ); 
}
//	------------------------------------------------------------------------
//
//	Comparing numeric values. See BOOST lib for reference.
//
//	Function names comes from Fortran operators:
//
//	LEQ =	Lower or Equal	
//	LT  =	Lower Than
//	GEQ =	Greater or Equal
//	GT	=	Greater Than
//	EQ	=	Equal
//	NEQ	=	Not Equal
//
//	All operators are deduced from the LT ( < ) defined as follows :
//
//		A<B		ssi		A-B < -machine epsilon
//				ssi		A-B < -4* ABS(A) * epsilon
//						or 
//						A-B < -4* ABS(B) * epsilon
//
//		4 is choosen
//			* 1 rouding in A-B
//			* 1 u_rounding fabs 
//			* 2 roudings in product	
// 
//		OR could be replaced by AND
//		
template<class NUM> static inline NUM NUM_abs(NUM x) { return x< 0 ? -x: x; }
template <class NUM> static inline bool lt(NUM left,NUM right,NUM eps = std::numeric_limits<NUM>::epsilon()) 
{ 
    // NUM d1   = NUM_div( left - right , NUM_abs( right ) );
	// NUM d2   = NUM_div( left - right, NUM_abs( left ) );
	return ( 
		( (left - right) < -4.*eps*NUM_abs(left) ) 
		|| 
		( (left - right) < -4.*eps*NUM_abs(right) )
		);
}
template <class NUM> static inline bool gt(NUM x,NUM y,NUM eps = std::numeric_limits<NUM>::epsilon()) 
{	
	return lt(y,x,eps); 
}
template <class NUM> static inline bool geq(NUM x,NUM y,NUM eps = std::numeric_limits<NUM>::epsilon()) 
{	
	return !lt(x,y,eps); 
}
template <class NUM> static inline bool leq(NUM x,NUM y,NUM eps = std::numeric_limits<NUM>::epsilon()) 
{	
	return !lt(y,x,eps); 
}
template <class NUM> static inline bool eq(NUM x,NUM y,NUM eps = std::numeric_limits<NUM>::epsilon()) 
{	
	return  leq(x,y,eps) && leq(y,x,eps) ; 
}
template <class NUM> static inline bool neq(NUM x,NUM y,NUM eps = std::numeric_limits<NUM>::epsilon()) 
{	
	return  !eq(x,y,eps) ; 
}

// -- local Helper for comparing T with a given tolerance
template <class T> class lt_comp
{
private: T itsTol ; 
public: lt_comp(const T& tol=std::numeric_limits<T>::epsilon()) : itsTol(tol) {}
public: bool operator()(const T&x,const T&y) const
		{ 
			return lt(x-y,-itsTol) ; 
		} 
}; 
// -- local Helper for comparing doubles with a TOL tolerance
class tol_comp
{
private: double itsTol ; 
public: tol_comp(const double& tol) : itsTol(tol) {}
public: bool operator()(const double&x,const double&y) 
		{ 
			return lt(x-y,-itsTol) ; 
		} 
}; 

// typedef std::set<double, lt_comp<double,std::numeric_limits<double>::epsilon() > > double_eps_set; 

//	JLA 
//
//	Outputting to streams.... 
std::ostream& operator<<(std::ostream&,const ARM_Vector&); 
//std::ostream& operator<<(std::ostream&,const ARM_Date&); 

// Onin
double Mean(const vector<double>& V); //Average of std::vector
double Sigma(const vector<double>& V, double Mean = -999.0); 
double VectorInterpol(const vector<double>& X,const vector<double>& Y, const double& x, int Method = K_LINEAR);

double LinearVectorInterpol(const vector<double>& X,const vector<double>& Y, const double& x);
double FlatVectorInterpol(const vector<double>& X,const vector<double>& Y, const double& x,bool valueinf= true);
double FlatVectorInterpol(const ARM_Vector& X,const ARM_Vector& Y, const double& x, bool valueinf= true);

void Bornes(const vector<double>& X,const double& x, int& _inf, int& _sup,bool &equal);
// double vector_max(const vector<double>& V);
// double vector_min(const vector<double>& V);
// double vector_sum(const vector<double>& V);

//CC
std::vector<std::string> SortVectorString(const std::vector<std::string>& vToSort, std::vector<int>& vPermut,const  int  incOrdec);
ARM_Vector SortWithConstraints(const ARM_Vector& vToSort, const std::vector<std::string>& vConstraintes, ARM_Vector& vPermut, const int intOrder);
bool Inf2Contraintes( const double& dref, const double& dtest, const std::string& sref, const std::string& stest);
// concat the line vector (of the matrix) into a string 
std::vector<std::string> concatStringVector(const std::vector<std::vector<std::string> >& VVString);
int FindIndexInVector(ARM_Vector& V,double value);

//gives nthLine,nthCol for ARM_VolLInterpol info is "strike|maturity"
void SearchCoordForVolatility(ARM_VolLInterpol* volatility,const std::string&  info,int& nthLine,int& nthCol); 

//	info is of the following form : 
//
//		[K:1]				Strike Value
//		[M:d:20060612]		for Maturity
//		[M:yf:0.25687]		for Year fraction	
//		[M:t:1Y]			for Tenor
//
//	output
//		line,col>0
//		line=0	: entire line
//		col=0	: entire col
//
//		otherwise Exception. 
//
//	also supported: 
//		[K:-]				returns 0 for Col
//		[M:d:-]				returns 0 for Row
//		
void SearchCoordForVolatility(const ARM_VolCurve&vol,const std::string& info,int& foundLine,int& foundCol); 

//Generate a schedule in days starting at start until end with a step of 'step'
ARM_Vector* GenerateIntSch(const double& start,const double& end, const int& step = 1);

//Generate a schedule in days starting at start until end with a step of 'step'
// with condition (genday-sch(i))>step else genday = sch(i)
ARM_Vector* GenerateIntSchInludingOtherSchedule(const double& start,
												const double& end, 
												ARM_Vector* sch,
												const int& step = 1);


//	------------------------------------------------------------------------
//		Returns the roll date in the specified convention
//
ARM_Date DateRoll(const ARM_Date& asof, qCDS_ADJ adj) ; 


//	------------------------------------------------------------------------
//		SHIFT an IR Curve to a specific Date (Theta...)
//
ARM_ZeroCurve*	GenerateIRCurveMovedInTime(ARM_ZeroCurve* IRCurve, ARM_Date TheDate);

void GenerateScheduleYF(ARM_Date& AsOf, ARM_Date& startdate, ARM_Date& enddate,long frequency,
						ARM_Vector& YFstartdates,ARM_Vector& YFenddates);

double ComputeEL_CDO_LHP(const double& tranche_down, const double& tranche_up, const double& beta_down,
							   const double& beta_up,  const double& Pdef,  double recovery);

double ComputeEL_LHP(const double& x, const double& beta,  const double& Pdef,  double recovery);
void ExtractVectorBeetwen2values(ARM_Vector* v,double start,double end, vector<double>& vout,int paytiming=K_ARREARS);
int getDefLegInterestRule(int feelegIntRule); 
void CptPeriodDates(const double& start,const double& end,
					const int& resetgap,const int& resettiming,const int& paytiming,ARM_Currency* ccy,
					double yfmaturity,qCDS_ADJ fwdadj,int paygap,double ForcePayDate,
					int resetweekday, int resetoccur,
					double& paydate,double& reset,double& fwdstart,double& fwdend);
void CptResetDateWeekly(const double& start,const double& end,int resetweekday,double yfmaturity, ARM_Currency* ccy,qCDS_ADJ fwdadj,ARM_Vector& resetDates, ARM_Vector& fwdstart, ARM_Vector& fwdend);

class ICM_Correlation;
void DeduceYearTermsForTSR(ICM_Correlation* correlation, 
						   double yf,
						   double& yt1_corr,
						   double& yt2_corr,
						   double& yt1_pdef,
						   double& yt2_pdef);

// Function will return the index name 
//		given by GetIndexType() for ARM index
//		given by First element of the underlying for Credit Index 
//
std::string GetIndexName(const ARM_IRIndex&index) ;
void vc2stdv(const int size,double* V,vector<double>& out);


double QuickCDSpv(ARM_Vector& schedule, const int& periodmatu,  double rate,double recovery,
					  const vector<double>& el,double notfeeleg,double notdefleg, 
					  ARM_ZeroCurve* zc, double& feepv,double& defpv);


double QuickFwdSpread(ARM_Vector& curveschedule,
					  ARM_Vector& defaultproba,
					  double yfstart,
					  double yfend,
					  ARM_ZeroCurve* zc,
					  double recovery,
					  double step = 0.25);

double QuickCdsDuration(ARM_Vector& curveschedule,
					  ARM_Vector& defaultproba,
					  double yfstart,
					  double yfend,
					  ARM_ZeroCurve* zc,
					  double recovery,
					  double step = 0.25);

#endif
