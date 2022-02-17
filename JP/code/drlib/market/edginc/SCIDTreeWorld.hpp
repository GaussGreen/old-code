#ifndef _SCIDtreeWorld_H
#define _SCIDtreeWorld_H

#include <vector>
#include "edginc/DoubleMatrix.hpp"
#include "edginc/Array.hpp"
#include "edginc/LinearInterpolator.hpp"
#include "edginc/YieldCurve.hpp"


DRLIB_BEGIN_NAMESPACE


class MARKET_DLL SimpleTree: public virtual VirtualDestructorBase
{
public:
	SimpleTree();
	SimpleTree(SimpleTree& tree, smartPtr<SimpleTree> &parent);
	SimpleTree(smartPtr<SimpleTree> &_parent, double _t, double _ft, double w);
	SimpleTree(smartPtr<SimpleTree> &_parent, double _s, double _t, double _fs, double _ft, double _w, size_t _generation);

	~SimpleTree(){};                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                              

	void set();
	void set(SimpleTree& _tree, smartPtr<SimpleTree> &parent);
	void set(smartPtr<SimpleTree> &_parent, double _t, double _ft, double _w); 
	void set(smartPtr<SimpleTree> &_parent, double _s, double _t, double _fs, double _ft, double _w, size_t _generation);

	bool check(); // for debugging purposes.. check the sum of weights is correct and that children[i]->parent is this

	void modifyFt(double _ft);
	void modifyWeight(double _weight);
	void feed(double *times, size_t size, size_t gen, vector<double> &slopes);

	void getLeavesUpToGeneration(vector < smartPtr<SimpleTree> > &leaves, size_t generation);
	void getWorldsFromTree(DoubleArray &x, DoubleArray &y, LinearInterpolantSP *f=0);
	void getWorldsFromTree(LinearInterpolantSP &f);

	void deleteZeroChildren();
	void CreateNearbyLeaves(vector< smartPtr< SimpleTree > > &newLeaves, DoubleArray &epsilons, double minbound, double maxbound, bool add);
	void CreateNearbyLeavesMultipleGens(vector< vector< smartPtr< SimpleTree > > > &newLeaves, size_t maxGen, DoubleArray &epsilons, double minbound, double maxbound, bool add);
	void reconcile(double epsilon, bool add);
	void changePositions(double x, double minbound, double maxbound);
	void changeWeights(double newWeight);
	void resetWeights(size_t gen);
	void addUpWeights(vector< smartPtr<SimpleTree> > &leaves);

	double weight, s, t, fs, ft;  // weight of this node, node deals with f between s and t (s<t)
								  // and f is equal to fs at s, and ft at t. Linear in between
	size_t generation;  // how far away from the creator of the tree!!

	vector< smartPtr<SimpleTree> > children;
	smartPtr<SimpleTree> parent;


// data a_sociated to this node
	CDoubleMatrixSP TEL;
	DoubleArraySP RA,DL;
	double PEL;
};	

typedef smartPtr<SimpleTree> SimpleTreeSP;

class MARKET_DLL SCIDtreeWorld: public virtual VirtualDestructorBase
								
{
public:
	SCIDtreeWorld(){};

	static double calibrateQuadProg(DoubleArraySP &newWeights,
		vector<DoubleArray> &MTMs,
		vector<DoubleArray> &PELconstraints,
		DoubleMatrix &weightConstraints,
		DoubleArray &sumWeightsConstraints,
		double stressPEL,
		DoubleArray& smoothing);

	static void convertTELtoLegs(const DoubleArray &kmin,
		const DoubleArray &kmax,
		const DoubleMatrix &TEL,
		const DateTime &today,
		const DateTimeArray &TELdates,
		const CashFlowArray &cf,
		const DateTime &maturity,
		const YieldCurveSP &discount,
		double *DL,   // (O)
		double *RA);  // (O)  must both be of size kmin.size()


	bool	check() {return tree->check();}
	void	globalFeed(DoubleMatrix &a_t, DoubleArray &w);
	void	feed(size_t gen, DoubleArrayArray &newPos);
	void    particularFeed(vector<SimpleTreeSP> &leaves, const DoubleArray &ratios, double minbound, double maxbound);
	
	void	getWorldsFromTree(vector<LinearInterpolantSP> &f, vector<SimpleTreeSP> &leaves); 
	void	getWorldsFromTree(DoubleArray &x, vector<DoubleArray> &y, vector<SimpleTreeSP> &leaves); 

	void	getInverseDistanceBetweenLeaves(vector < smartPtr<SimpleTree> > &leaves, DoubleArray &invDist);
	void	getLeaves(vector < SimpleTreeSP > &leaves, int gen = -1) { if (gen>=0) tree->getLeavesUpToGeneration(leaves,gen); 
																			   else tree->getLeavesUpToGeneration(leaves, treeTimes.size()); }

	void	readTree(DoubleMatrix &valuesAtTreeTime); // for debugging

	void	readTEL(SimpleTreeSP &leave, DoubleMatrix & TEL);
	void	readPEL(SimpleTreeSP &leave, DoubleArray & PEL);
	void	readLegs(SimpleTreeSP &leave, DoubleMatrix& DL, DoubleMatrix& RA);

	size_t  getMaxNbGeneration() { return treeTimes.size(); }
	
	void    setTreeTimes(DoubleArray &_treeTimes, DateTime _today);
	void	setTreeDates(DateTimeArray &_treeDates, DateTime _today);
	void	setTELtimes(DoubleArray &times);
	void	setTrancheData( DoubleArray &_kmin, 
							DoubleArray &_kmax, 
							int _coupon,
							YieldCurveSP _discount,
							DayCountConventionSP _dcc);
	void	fromTELtoLegsOneGeneration(vector<SimpleTreeSP> &leaves);
	void	fromTELtoLegsNGeneration(vector<SimpleTreeSP> &leaves, int lastGen);

	double  findWeights(size_t genLow, size_t genHigh,
						const DoubleMatrix &parSpread,
						const DoubleArray &PEL,
						double stressPEL,
						DoubleArray& smoothing,
						bool deleteZeroWorlds);        // par spread for all tranches and PEL At all MTM Dates	
													// but we calibrate just for one generation
	DateTime	   today;	
	DateTimeArray  treeDates, TELdates;
	DoubleArray    treeTimes, TELtimes;
	IntArray       TELindexes;
	DoubleArray	   kmin, kmax;
	int			   coupon;
	YieldCurveSP   discount;
	DayCountConventionSP dcc;

private:
	void TELtoLegs( int lastTELdate,
					DoubleMatrix &TEL,
					DateTime maturity,
					DoubleArray &DL, 
					DoubleArray &RA);



	SimpleTreeSP tree;
};

typedef smartPtr<SCIDtreeWorld> SCIDtreeWorldSP;


DRLIB_END_NAMESPACE

#endif
