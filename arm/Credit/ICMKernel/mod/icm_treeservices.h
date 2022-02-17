#ifndef	_ICM_Tree_Services_
#define _ICM_Tree_Services_

#include "ICMKernel/mod/icm_binomialtree.h"

//	--	Forward
class ICM_DefaultCurve ;

//	-------------------------------------------------------------
//		This is the "standard" binomial tree.
//			cf. Fakher docs
//	-------------------------------------------------------------
class	BinomialData
{
public:
	double	h;		// h(i,j)	intensity
	double	p;		// p(i,j)	proba to be at (i,j) alive
public:
	BinomialData() : h(0),p(0) {}	//	init done to 0 
} ;

//	-------------------------------------------------------------
static inline 
std::ostream& operator<<(std::ostream&o,const BinomialData&ref)
{ o<<"[h="<<ref.h<<" p="<<ref.p<<"]"; return o; }

typedef BinomialTree<BinomialData>	StdBinomialTree ;


//	------------------------------------------------------------
//		"Financial" Services are located here
//	------------------------------------------------------------
class TreeServices
{
public:
	//
	//	--	Calibration Method 
	//
	//		This will calibrate computes (h,p) for each node 
	//		( assuming transitions probas are 0. 5)
	//
	static void calibrate(StdBinomialTree& tree,
		ICM_DefaultCurve& defCurve,
		double sigma) ;
	//
	//	--	Calibration Method 
	//
	//		This will calibrate computes (h,p) for each node 
	//		( assuming transitions probas are 0. 5)
	//
	static void calibrate2(StdBinomialTree& tree,
		ICM_DefaultCurve& defCurve,
		double sigma) ;
	//
	//	--	Pricing Method 
	//
	//		Backward pricing of a set of known flows occuring a 
	//		given yearfractions . Backward is done until yearfraction 0
	//		(first step of the tree)
	//
	//		flows is a 2D matrix , col1=yearfractions, col2=flows.
	//		yearfractions should be a subset of the tree discretisation
	//
	//		If "default", losses is payed at the next slice . This is a rough approximate
	//
	//		The output is a single double value. 
	//
	static double backwardprice(const StdBinomialTree& tree,const ICM_QMatrix<double>&flows,ARM_ZeroCurve&discCurve) ; 
	static double backwardprice_default(const StdBinomialTree& tree,
		double losses,
		ARM_ZeroCurve&discCurve) ; 
	//
	//	--	Pricing Method 
	//
	//		Backward pricing of a set of known flows occuring a 
	//		given yearfractions. 
	//		
	//		Backwarding is done from one startSlicePosition  
	//		to an endSlicePosition. startSlicePosition>=endSlicePosition < tree.depth()
	//
	//		flows is a 2D matrix , col1=yearfractions, col2=flows.
	//		yearfractions should be a subset of the tree discretisation
	//		flows outside [startSlicePosition,endSlicePosition] are ignored. 
	//
	//		If "default", losses is payed at the next slice . This is a rough approximate
	//		The output is a ICM_QMatrix<double> that will be resized 
	//		to fit tree.slice(endSlicePosition).size() . 
	//
	
	static void backwardprice(const StdBinomialTree& tree,
							  unsigned long backwardFrom,
							  unsigned long backwardTo,
							  const ICM_QMatrix<double>&flows,
							  ARM_ZeroCurve&discCurve,
							  const ICM_QMatrix<double>&statesFrom, 
							  ICM_QMatrix<double>&ret,
							  ICM_QMatrix<double>&ret_Greeks) ; 
	
	static void  backwardprice_default(const StdBinomialTree& tree,
									   unsigned long backwardFrom,
									   unsigned long backwardTo,
									   double losses,
									   ARM_ZeroCurve&discCurve,
									   ICM_QMatrix<double>&ret) ; 
	
	//
	//		flows is an array describing the contigent payments
	//			col0	yf of start observation date 
	//			col1	yf of end observation date
	//			col2	yf of payment date
	//			col3	payment value
	//
	//		It is assumed that the period defined by ]col0,col1] do not overlap
	//		and are ordered (this is not checked at run time)
	//
	//		When backwarding to slice at T(i), we condider the line of flows where
	//
	//			yfStartObservation <  T(i) <= yfEndObservation 
	//
	static void  backwardprice_default(const StdBinomialTree& tree,
									   unsigned long backwardFrom,
									   unsigned long backwardTo,
									   const ICM_QMatrix<double> &flows,
									   ARM_ZeroCurve&discCurve,
									   const ICM_QMatrix<double> &statesFrom,
									   ICM_QMatrix<double>&ret);

	static void  backwardprice_default3(const StdBinomialTree& tree,
									    unsigned long backwardFrom,
										unsigned long backwardTo,
										const ICM_QMatrix<double> &flows,
										ARM_ZeroCurve&discCurve,
										const ICM_QMatrix<double> &statesFrom,
										ICM_QMatrix<double>&ret,
										ICM_QMatrix<double>&ret_Greeks);

	static void backwardprice_BN(const StdBinomialTree& tree,
							     unsigned long backwardFrom,
								 unsigned long backwardTo,
								 const StdBinomialTree& flows1,
								 const ICM_QMatrix<double>&PayDates1,
								 const ICM_QMatrix<double>&flows2,
								 ARM_ZeroCurve&discCurve,
								 const ICM_QMatrix<double>&statesFrom,								 
								 ICM_QMatrix<double>&ret,
								 ICM_QMatrix<double>&ret_Greeks); 

	//	--	Helpers 
	static void dumpProbas(std::ostream& o,const StdBinomialTree& tree) ;
	static void dumpIntensities(std::ostream& o,const StdBinomialTree& tree) ;
	
	// Dichotomic research of indexes

	static void dichotomic_research(const StdBinomialTree& tree,
								   double yearterm,
								   unsigned long minindex, 
								   unsigned long maxindex,
								   ICM_QMatrix<double>&ret);
} ;
#endif // _ICM_Tree_Services_