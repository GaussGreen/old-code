
#include "firsttoinc.h"
#include "icm_treeservices.h"
#include "ICMKernel/crv/icm_defaultcurve.h"
#include "ICMKernel/util/icm_utils.h"
#include <nag.h>
#include "ICMKernel\util\icm_RootFinder1D.h"

//	-----------------------------------------------------------------
//		Local functor used in the calibration (determining h i,0)
//	-----------------------------------------------------------------
class	localExpFunction
{
public : 
	localExpFunction(
		StdBinomialTree::slice_t::iterator begin,
		StdBinomialTree::slice_t::iterator end,
		double prevDT,double nextDT,double sigma,double survival,double q) :
	   _begin(begin),_end(end),_prevDT(prevDT),_nextDT(nextDT),_sigma(sigma),_survival(survival),_q(q){}
	double operator()(double x) const 
	{
		double res=0; 
		unsigned int j; 
		StdBinomialTree::slice_t::iterator it; 
		for(j=0, it=_begin; it!=_end; ++it,++j) 
			res += it->data().p*exp(-x*_nextDT*exp(j*_sigma*sqrt(_prevDT/(_q*(1.-_q))))); 
		return res - _survival ;
	}
private:
	StdBinomialTree::slice_t::iterator _begin,_end; 
	double _prevDT,_nextDT,_sigma,_survival,_q; 
} ;
//	-----------------------------------------------------------------
//		Local functor used in the calibration (determining q (down proba))
//	-----------------------------------------------------------------
class	localFunction
{
public : 
	localFunction(
		double DT,double sigma) :
	   _DT(DT),_sigma(sigma) {}
	double operator()(double x) const 
	{
		double tmp = _sigma*sqrt(_DT/(x*(1.-x))); 
		return log(x + (1.-x)*exp(tmp))-(1.-x)*tmp-0.5*_sigma*_sigma*_DT ;
	}
private:
	double _DT,_sigma; 
} ;
//	-----------------------------------------------------------------
//static 
void
TreeServices::calibrate(StdBinomialTree& aTree,
						ICM_DefaultCurve& defCurve,
						double sigma)
{
	//	--	Can't work on less than 2 level trees
	if (aTree.depth()<=2) 
		ICMTHROW(ERR_INVALID_ARGUMENT,"TreeServices::calibrate: Can't work on less that two level trees"); 

	//	--	For each node both proba are set to .
	StdBinomialTree::iterator it=aTree.begin(); 

	ICMLOG("TreeServices::calibrate: calibrating a tree of depth "<<aTree.depth()); 
	// ICMLOG("... assigning .5 probas"); 
	while (it!=aTree.end()) 
	{
		StdBinomialTree::node_elt_t::child_iterator child = it->begin(); 
		child->first= 0.5; 
		++child; 
		child->first= 0.5; 
		++it; 
	}

	//	--	Initialisation of data().p=0 for each node is not required
	//		due to the data constructor

	//	--	For each slice 
	//
	unsigned int i = 0 ; 

	for(i=0;i<aTree.depth()-1; i++) 
	{
		// ICMLOG("... processing slice "<<i); 
		StdBinomialTree::slice_t aSlice=aTree.slice(i); 
		// double DT = aTree.DT(i) ;
		StdBinomialTree::slice_t::iterator it = aSlice.begin(); 
		StdBinomialTree::node_data_t& baseData= it->data(); 
		//	--	We solve the equation giving h(i,0) 
		//		for i=0 we use explicit formula,
		//		otherwise this is a Dichotomic search

		if (i==0) 
		{
			//	--	For the first node, 
			//		h(0,0) = -ln ( Survival(T1) ) / T1 
			//		p(0,0) = 1. 
			baseData.h= - log( defCurve.SurvivalProba(aTree.Time(i+1)) ) / aTree.DT(i); 
			baseData.p= 1.; 
			// ICMLOG("... initial node processed "<<baseData); 
		}
		else 
		{
			// double q = aTree.node(i-1,0).begin()->first	; 
			//	search f(x)=0 
			localExpFunction f(
				aSlice.begin(),
				aSlice.end(),
				aTree.DT(i-1),
				aTree.DT(i),
				sigma,
				defCurve.SurvivalProba(aTree.Time(i+1)),
				0.5
				); 
			double result = RootFinder1D(f).Dichotomy(0,1,100, 1.e-12,1.e-12); 
			baseData.h = result; 
			// ICMLOG("... h value "<<result); 
			//	--	All h(i,j) are obtained from h(i,0) by 
			//		multiplying with exp(2*sigma*sqrt(DTi)) 
			//		
			// ICMLOG("... propagating h to current Slice"); 
			double coeff = exp(2.*sigma*sqrt(aTree.DT(i-1))); 
			// double coeff = exp(sigma*sqrt(aTree.DT(i-1)/(q*(1.-q))); 
			for (unsigned int j=1;j<aSlice.size(); j++) 
			{
				aSlice.data(j).h=aSlice.data(j-1).h*coeff; 
				// ICMLOG("... h(.,"<<j<<")="<<aSlice.data(j).h); 
			}
		}
		

		//	--	For all existing children of this slice : 
		//		the proba is updated 
		// ICMLOG("... computing p values for next slice"); 
		for (it=aSlice.begin();it!=aSlice.end();++it) 
		{
			StdBinomialTree::node_elt_t::child_iterator child ; 
			for (child=it->begin(); child!=it->end(); ++child) 
			{
				if (child->second)
				{
					child->second->data().p += child->first * exp( - it->data().h * aTree.DT(i) ) * it->data().p ; 
				}
			}
		}
	}	
}

//	-----------------------------------------------------------------
//static 
void
TreeServices::calibrate2(StdBinomialTree& aTree,
						ICM_DefaultCurve& defCurve,
						double sigma)
{
	//	--	Can't work on less than 2 level trees
	if (aTree.depth()<=2) 
		ICMTHROW(ERR_INVALID_ARGUMENT,"TreeServices::calibrate: Can't work on less that two level trees"); 

	//	--	For each node both proba are set to .
	StdBinomialTree::iterator it=aTree.begin(); 

	ICMLOG("TreeServices::calibrate: calibrating a tree of depth "<<aTree.depth()); 

	//	Compute the q (down) probability for each slice, and set it on the tree. 
	unsigned int i = 0 ; 
	for(i=0;i<aTree.depth()-1;i++)
	{
		localFunction f(aTree.DT(i),sigma); 
		double result = RootFinder1D(f).Dichotomy(0,1,100, 1.e-12,1.e-12); 
//		ICMLOG("Transition proba (down) at slice "<<i<<" = "<<result); 
		StdBinomialTree::slice_t::iterator it = aTree.slice(i).begin(); 
		while (it!=aTree.slice(i).end()) 
		{
			it->begin()->first = result; 
			(it->begin()+1)->first = 1.-result; 
			++it; 
		}
	}

	//	--	Initialisation of data().p=0 for each node is not required
	//		due to the data constructor

	for(i=0;i<aTree.depth()-1; i++) 
	{
		StdBinomialTree::slice_t aSlice=aTree.slice(i); 
		StdBinomialTree::slice_t::iterator it = aSlice.begin(); 
		
		for (it=aSlice.begin();it!=aSlice.end();++it) 
		{
			StdBinomialTree::node_elt_t::child_iterator child ; 
			for (child=it->begin(); child!=it->end(); ++child) 
			{
				if (child->second)
				{
					child->second->data().p = 0.; 
				}
			}
		}
	}

	//	--	For each slice 
	//
	
	for(i=0;i<aTree.depth()-1; i++) 
	{
		// ICMLOG("... processing slice "<<i); 
		StdBinomialTree::slice_t aSlice=aTree.slice(i); 
		// double DT = aTree.DT(i) ;
		StdBinomialTree::slice_t::iterator it = aSlice.begin(); 
		StdBinomialTree::node_data_t& baseData= it->data(); 
		//	--	We solve the equation giving h(i,0) 
		//		for i=0 we use explicit formula,
		//		otherwise this is a Dichotomic search

		if (i==0) 
		{
			//	--	For the first node, 
			//		h(0,0) = -ln ( Survival(T1) ) / T1 
			//		p(0,0) = 1. 
			baseData.h= - log( defCurve.SurvivalProba(aTree.Time(i+1)) ) / aTree.DT(i); 
			baseData.p= 1.; 
			// ICMLOG("... initial node processed "<<baseData); 
		}
		else 
		{
			double q = aTree.node(i-1,0).begin()->first	; 

			//	search f(x)=0 
			// localExpFunction f(aSlice.begin(),aSlice.end(),aTree.DT(i),sigma,defCurve.SurvivalProba(aTree.Time(i+1))); 
			localExpFunction f(
				aSlice.begin(),
				aSlice.end(),
				aTree.DT(i-1),
				aTree.DT(i),
				sigma,
				defCurve.SurvivalProba(aTree.Time(i+1)),
				q); 
			double result = RootFinder1D(f).Dichotomy(0,1,100, 1.e-12,1.e-12); 
			baseData.h = result; 
			// ICMLOG("... h value "<<result); 
			//	--	All h(i,j) are obtained from h(i,0) by 
			//		multiplying with exp(2*sigma*sqrt(DTi)) 
			//		
			// ICMLOG("... propagating h to current Slice"); 
			// double coeff = exp(2.*sigma*sqrt(aTree.DT(i-1))); 
			double coeff = exp(sigma*sqrt(aTree.DT(i-1)/(q*(1.-q)))); 
			for (unsigned int j=1;j<aSlice.size(); j++) 
			{
				aSlice.data(j).h=aSlice.data(j-1).h*coeff; 
				//ICMLOG("... h(.,"<<j<<")="<<aSlice.data(j).h);
			}
		}
		

		//	--	For all existing children of this slice : 
		//		the proba is updated 
		// ICMLOG("... computing p values for next slice"); 
		for (it=aSlice.begin();it!=aSlice.end();++it) 
		{
			StdBinomialTree::node_elt_t::child_iterator child ; 
			for (child=it->begin(); child!=it->end(); ++child) 
			{
				if (child->second)
				{
					child->second->data().p += child->first * exp( - it->data().h * aTree.DT(i) ) * it->data().p ; 
				}
			}
		}
	}	
}

//	-----------------------------------------------------------------
double 
TreeServices::backwardprice(const StdBinomialTree& tree,const ICM_QMatrix<double>&flows,
							ARM_ZeroCurve& discCurve)
{

	// 
	long backwardFrom(tree.depth()-1); 
	long backwardTo(0); 
	ICM_QMatrix<double> statesFrom(tree.slice(backwardFrom).size(),1); 
	ICM_QMatrix<double> ret ; 
	ICM_QMatrix<double> ret_greeks ; 
	backwardprice(tree,backwardFrom,backwardTo,flows,discCurve,statesFrom,ret,ret_greeks); 
	return ret(0,0) ; 
	// 
		 
}
//	-----------------------------------------------------------------
//static 
void TreeServices::backwardprice(const StdBinomialTree& tree,
								 unsigned long backwardFrom,
								 unsigned long backwardTo,
								 const ICM_QMatrix<double>&flows,
								 ARM_ZeroCurve&discCurve,
								 const ICM_QMatrix<double>&statesFrom, 
								 ICM_QMatrix<double>&ret,
								 ICM_QMatrix<double>&ret_greeks) 
{
	//	--	Can't compute	
	if (tree.depth()<=1) 
		ICMTHROW(ERR_INVALID_ARGUMENT,"TreeServices::backwardprice: degenerated tree (depth="<<tree.depth()<<")"); 
	if (backwardTo>backwardFrom) 
		ICMTHROW(ERR_INVALID_ARGUMENT,"TreeServices::backwardprice: can't backword from slice "<<backwardFrom<<" to slice "<<backwardTo); 
	if (backwardFrom>=tree.depth()) 
		ICMTHROW(ERR_INVALID_ARGUMENT,"TreeServices::backwardprice: startSlicePosition "<<backwardFrom<<" out of bounds (" << tree.depth() ); 
	if (statesFrom.Getnbrows() < tree.slice(backwardTo).size() ) 
		ICMTHROW(ERR_INVALID_ARGUMENT,"TreeServices::backwardprice: statesFrom dimension incorrect"); 

	//	--	Resize the output 
	ret.Resize(tree.slice(backwardTo).size(),1); 
	ret_greeks.Resize(6,1) ;

	//	--	No flows values
	if (flows.Getnbcols()<=1) 
		ICMTHROW(ERR_INVALID_ARGUMENT,"TreeServices::backwardprice: bad dims flows (2 cols required)") ; 

	//	--	From now, we have at least two slices in the tree, 
	//		and at least one flow.
 
	
	//	--	Initialisation
	//		
	std::vector<double> tmp1(tree.slice(backwardFrom).size()); 
	std::vector<double>	tmp2(tree.slice(backwardFrom).size()); 
	std::vector<double>::iterator root1=tmp1.begin(); 
	std::vector<double>::iterator root2=tmp2.begin(); 
	
	//	The states are setted to their initial value 
	for(unsigned int i=0;i<tree.slice(backwardFrom).size();i++) 
		tmp1[i]=statesFrom(i,0); 

	//	Ini
	int iTree ;		// for the tree
	long iFlow ;	// for the flows
	iFlow=flows.Getnbrows()-1; 

	//	We look for the flow happening at backwardFrom time. 
	iTree=backwardFrom; 
	double T=tree.Time(iTree); 
	while (iFlow>=0 && gt(flows(iFlow,0),T))iFlow--; 

	//	if found, we add it to tmp1 
	if (iFlow>=0 && eq(flows(iFlow,0),T))	
	{
//		ICMLOG("(init) Slice= "<<iTree<<" treeTime="<<T<<" matched flowTime="<<flows(iFlow,0)<<" added Flow="<<flows(iFlow,1)); 
		for(unsigned int j=0;j<tree.slice(iTree).size();j++)
			tmp1[j] += flows(iFlow,1) ; 
	}

	//	--	Iteration 
	long signedFrom(backwardFrom); 
	long signedTo(backwardTo); 
	unsigned int i0 = 0 ;

	for(iTree=signedFrom-1L; iTree>=signedTo+1L;iTree--) 
	{
		//	--	Bacwarding c1 to c2
		std::vector<double>::iterator c1(root1); 
		std::vector<double>::iterator c2(root2); 

		StdBinomialTree::const_slice_t aSlice(tree.slice(iTree)); 
		StdBinomialTree::const_slice_t::const_iterator it = aSlice.begin(); 
		
		while (it!=aSlice.end()) 
		{
			StdBinomialTree::node_elt_t::const_child_iterator child1(it->begin()); 
			StdBinomialTree::node_elt_t::const_child_iterator child2(it->begin()+1); 
			*c2 =     
					child1->first *  (*c1) 
				+	child2->first *  (*(c1+1)) ; 
				
			*c2 *= exp( -it->data().h * tree.DT(iTree) ); 
			*c2 *= discCurve.DiscountPrice(tree.Time(iTree+1))/discCurve.DiscountPrice(tree.Time(iTree)) ;
			//ICMLOG("Slice= "<<iTree<<" treeTime="<<tree.Time(iTree)<<" matched value"<<*c2 );
			if(iTree<=2)
			{
				ret_greeks(i0,0) = *c2 ;
				i0++ ;
			}			

			++c2; 
			++c1; 
			++it; 
		}

		//	-- adding to c2 the flow at iTree time. 
		c2=root2; 
		double T=tree.Time(iTree); 
		while (iFlow>=0 && gt(flows(iFlow,0),T)) iFlow--; 
		//	if found, we add it to tmp1 
		if (iFlow>=0 && eq(flows(iFlow,0),T))	
		{
			// ICMLOG("(backw) Slice= "<<iTree<<" treeTime="<<T<<" matched flowTime="<<flows(iFlow,0)<<" added Flow="<<flows(iFlow,1)); 
			for(unsigned int j=0;j<tree.slice(iTree).size();j++) { *c2 += flows(iFlow,1) ; ++c2; } 
		}

		//	--	At the end of the loop, we swap root1 and root2
		//		so that root1 reprensents once again the latest computed values
		std::swap(root1,root2);
	}

	//	--	We finally compute the backwardTo slice, without considering any flow that could exist.
	{
		iTree=backwardTo; 
		std::vector<double>::iterator c1(root1); 
		std::vector<double>::iterator c2(root2); 

		StdBinomialTree::const_slice_t aSlice(tree.slice(iTree)); 
		StdBinomialTree::const_slice_t::const_iterator it = aSlice.begin(); 
		
		while (it!=aSlice.end()) 
		{
			StdBinomialTree::node_elt_t::const_child_iterator child1(it->begin()); 
			StdBinomialTree::node_elt_t::const_child_iterator child2(it->begin()+1); 
 			*c2 =     
 					child1->first *  (*c1) 
 				+	child2->first *  (*(c1+1)) ; 


			*c2 *= exp( -it->data().h * tree.DT(iTree) ); 
			*c2 *= discCurve.DiscountPrice(tree.Time(iTree+1))/discCurve.DiscountPrice(tree.Time(iTree)) ;
			//ICMLOG("Slice= "<<iTree<<" treeTime="<<tree.Time(iTree)<<" matched value"<<*c2 );
			ret_greeks(i0,0) = *c2 ;
			++c2; 
			++c1; 
			++it; 
		}
	}


	//	--	Desired values are in c2
	//		
	std::vector<double>::iterator c2(root2);
	for(unsigned int j=0;j<tree.slice(backwardTo).size();j++) 
		ret(j,0)= *(c2++) ;
}

//	-----------------------------------------------------------------
double 
TreeServices::backwardprice_default(const StdBinomialTree& tree,
									double losses,
									ARM_ZeroCurve& discCurve)
{

	//
	long backwardFrom=tree.depth()-1; 
	long backwardTo=0; 
	ICM_QMatrix<double> ret; 
	backwardprice_default(tree,backwardFrom,backwardTo,losses,discCurve,ret); 
	return ret(0,0); 

}
//	-----------------------------------------------------------------
// static 
void 
TreeServices::backwardprice_default(const StdBinomialTree& tree,
						  unsigned long backwardFrom,
						  unsigned long backwardTo,
						  double losses,
						  ARM_ZeroCurve&discCurve,
						  ICM_QMatrix<double>&ret) 
{
	//	--	Can't compute	
	if (tree.depth()<=1) 
		ICMTHROW(ERR_INVALID_ARGUMENT,"TreeServices::backwardprice_default: degenerated tree (depth="<<tree.depth()<<")"); 
	if (backwardTo>backwardFrom) 
		ICMTHROW(ERR_INVALID_ARGUMENT,"TreeServices::backwardprice_default: can't backword from slice "<<backwardFrom<<" to slice "<<backwardTo); 
	if (backwardFrom>=tree.depth()) 
		ICMTHROW(ERR_INVALID_ARGUMENT,"TreeServices::backwardprice_default: startSlicePosition "<<backwardFrom<<" out of bounds (" << tree.depth() ); 
	
	//	--	Resize the output 
	ret.Resize(tree.slice(backwardTo).size(),1); 

	
	//	--	Initialisation
	//		
	std::vector<double> tmp1(tree.slice(backwardFrom).size()); 
	std::vector<double>	tmp2(tree.slice(backwardFrom).size()); 
	std::vector<double>::iterator root1=tmp1.begin(); 
	std::vector<double>::iterator root2=tmp2.begin(); 

	//	We populate c1 with 0 values. 
	std::fill(tmp1.begin(),tmp1.end(),0.); 

	//	--	Iteration 
	int iTree ;	// for the tree

	long signedTo(backwardTo); 
	long signedFrom(backwardFrom); 
	for(iTree=signedFrom-1L; iTree>=signedTo;iTree--) 
	{
		//	--	When entering this loop, we compute root2 from root1

		//	--	Adding the external flow to root1
		double T=tree.Time(iTree); 
		//	--	Backward computation of root2
		StdBinomialTree::const_slice_t aSlice(tree.slice(iTree)); 
		StdBinomialTree::const_slice_t::const_iterator it = aSlice.begin(); 
		std::vector<double>::iterator c1(root1); 
		std::vector<double>::iterator c2(root2); 
		while (it!=aSlice.end()) 
		{
			StdBinomialTree::node_elt_t::const_child_iterator child1(it->begin()); 
			StdBinomialTree::node_elt_t::const_child_iterator child2(it->begin()+1); 
	
			double local = exp( -it->data().h * tree.DT(iTree) ); 
			
 			*c2 =	losses * (1.-local) 
				 +	local* (
							child1->first *  (*c1)
						+	child2->first * (*(c1+1))
					) ; 


			*c2 *=  discCurve.DiscountPrice(tree.Time(iTree+1))/discCurve.DiscountPrice(tree.Time(iTree)) ;
			++c2; 
			++c1; 
			++it; 
		}
		//	--	At the end of the loop, we swap root1 and root2
		//		so that root1 reprensents once again the latest computed values
		std::swap(root1,root2); 
	}


	//	--	After the loop, root1 first elements is the computed value
	//		
	std::vector<double>::iterator c1(root1); 
	for(unsigned int j=0;j<tree.slice(backwardTo).size();j++) 
		ret(j,0)= *(c1++) ;
}

//	-----------------------------------------------------------------
// static 
void 
TreeServices::backwardprice_default(const StdBinomialTree& tree,
									unsigned long backwardFrom,
									unsigned long backwardTo,
									const ICM_QMatrix<double>& flows,
									ARM_ZeroCurve&discCurve,
									const ICM_QMatrix<double> &statesFrom,
									ICM_QMatrix<double>&ret) 
{
	//	--	Can't compute	
	if (tree.depth()<=1) 
		ICMTHROW(ERR_INVALID_ARGUMENT,"TreeServices::backwardprice_default: degenerated tree (depth="<<tree.depth()<<")"); 
	if (backwardTo>backwardFrom) 
		ICMTHROW(ERR_INVALID_ARGUMENT,"TreeServices::backwardprice_default: can't backword from slice "<<backwardFrom<<" to slice "<<backwardTo); 
	if (backwardFrom>=tree.depth()) 
		ICMTHROW(ERR_INVALID_ARGUMENT,"TreeServices::backwardprice_default: startSlicePosition "<<backwardFrom<<" out of bounds (" << tree.depth() ); 
	if (statesFrom.Getnbrows() < tree.slice(backwardTo).size() ) 
		ICMTHROW(ERR_INVALID_ARGUMENT,"TreeServices::backwardprice_default: statesFrom dimension incorrect"); 
	
	//	--	Resize the output 
	ret.Resize(tree.slice(backwardTo).size(),1); 

	//	--	No flows values
	if (flows.Getnbcols()<4) 
		ICMTHROW(ERR_INVALID_ARGUMENT,"TreeServices::backwardprice_default: bad dims flows (4 cols required)") ; 
	
	//	--	Initialisation
	//		
	std::vector<double> tmp1(tree.slice(backwardFrom).size()); 
	std::vector<double>	tmp2(tree.slice(backwardFrom).size()); 
	std::vector<double>::iterator root1=tmp1.begin(); 
	std::vector<double>::iterator root2=tmp2.begin(); 

	//	The states are setted to their initial value 
	for(unsigned int i=0;i<tree.slice(backwardFrom).size();i++) 
		tmp1[i]=statesFrom(i,0); 

	//	--	Iteration 
	int iTree ;	// for the tree
	int iFlow ; // for the flows
	
	long signedTo(backwardTo);		
	long signedFrom(backwardFrom); 
	iFlow = flows.Getnbrows()-1;
	for(iTree=signedFrom-1L; iTree>=signedTo;iTree--) 
	{
		//	--	When entering this loop, we compute root2 from root1
		double T=tree.Time(iTree) ; 

		//	--	We lool for the payment description that will apply at this level
		//	
		//		yfStartObservation <= T < yfEndObservation 
		//
		while 
			(
				iFlow>=0 
			&& 
				!(  leq(flows(iFlow,0),T) &&  lt(T,flows(iFlow,1)) ) 
			) iFlow--; 
		if (iFlow<0) 
			ICMTHROW(ERR_INVALID_ARGUMENT,"TreeServices::backwardprice_default: Can't find a contingent payment description for slice "<<iTree<<" at T="<<T); 
//		ICMLOG("Slice "<<iTree<<" at "<<T<<": obsStart="<<flows(iFlow,0)<<" obsEnd="<<flows(iFlow,1)<<" yfPay="<<flows(iFlow,2)<<" payValue="<<flows(iFlow,3)); 
		
		//	--	Backward computation of root2
		StdBinomialTree::const_slice_t aSlice(tree.slice(iTree)); 
		StdBinomialTree::const_slice_t::const_iterator it = aSlice.begin(); 
		std::vector<double>::iterator c1(root1); 
		std::vector<double>::iterator c2(root2); 
		while (it!=aSlice.end()) 
		{
			StdBinomialTree::node_elt_t::const_child_iterator child1(it->begin()); 
			StdBinomialTree::node_elt_t::const_child_iterator child2(it->begin()+1); 
	
			double local = exp( -it->data().h * tree.DT(iTree) ); 
			
			double ctgPart =
					(1.-local)
				*	flows(iFlow,3) 
				*	discCurve.DiscountPrice(flows(iFlow,2))/discCurve.DiscountPrice(tree.Time(iTree)) ; 

			double nonctgPart =
					local
				*	( child1->first *  (*c1) +	child2->first * (*(c1+1)) )
				*	discCurve.DiscountPrice(tree.Time(iTree+1))/discCurve.DiscountPrice(tree.Time(iTree)) ; 

			*c2 = ctgPart + nonctgPart ;

			++c2; 
			++c1; 
			++it; 
		}
		//	--	At the end of the loop, we swap root1 and root2
		//		so that root1 reprensents once again the latest computed values
		std::swap(root1,root2); 
	}


	//	--	After the loop, root1 first elements is the computed value
	//		
	std::vector<double>::iterator c1(root1); 
	for(unsigned int j=0;j<tree.slice(backwardTo).size();j++) 
		ret(j,0)= *(c1++) ;
}

//	-----------------------------------------------------------------
//static 
void
TreeServices::dumpProbas(std::ostream& o,const StdBinomialTree& aTree)
{
	o<<"-- Dumping Tree : Probas "<<std::endl ;
	for(unsigned int i=0;i<aTree.depth();i++)
	{
		for(unsigned j=0;j<=i;j++) 
			o<<aTree.data(i,j).p<<" "; 
		o<<std::endl; 
	}
	o<<"-- End of Dumping Tree : Probas "<<std::endl ;
}
//	-----------------------------------------------------------------
//static 
void
TreeServices::dumpIntensities(std::ostream& o,const StdBinomialTree& aTree)
{
	o<<"-- Dumping Tree : Intensities "<<std::endl ;
	for(unsigned int i=0;i<aTree.depth();i++)
	{
		for(unsigned j=0;j<=i;j++) 
			o<<aTree.data(i,j).h<<" "; 
		o<<std::endl; 
	}
	o<<"-- End of Dumping Tree : Intensities "<<std::endl ;
}

//	-----------------------------------------------------------------
// static 
void 
TreeServices::backwardprice_default3(const StdBinomialTree& tree,
									unsigned long backwardFrom,
									unsigned long backwardTo,
									const ICM_QMatrix<double>& flows,
									ARM_ZeroCurve&discCurve,
									const ICM_QMatrix<double> &statesFrom,
									ICM_QMatrix<double>&ret,
									ICM_QMatrix<double>&ret_greeks) 
{
	//	--	Can't compute	
	if (tree.depth()<=1) 
		ICMTHROW(ERR_INVALID_ARGUMENT,"TreeServices::backwardprice_default3: degenerated tree (depth="<<tree.depth()<<")"); 
	if (backwardTo>backwardFrom) 
		ICMTHROW(ERR_INVALID_ARGUMENT,"TreeServices::backwardprice_default3: can't backword from slice "<<backwardFrom<<" to slice "<<backwardTo); 
	if (backwardFrom>=tree.depth()) 
		ICMTHROW(ERR_INVALID_ARGUMENT,"TreeServices::backwardprice_default3: startSlicePosition "<<backwardFrom<<" out of bounds (" << tree.depth() ); 
	if (statesFrom.Getnbrows() < tree.slice(backwardTo).size() ) 
		ICMTHROW(ERR_INVALID_ARGUMENT,"TreeServices::backwardprice_default3: statesFrom dimension incorrect"); 
	
	//	--	Resize the output 
	ret.Resize(tree.slice(backwardTo).size(),1); 
	ret_greeks.Resize(6,1); 

	//	--	No flows values
	if (flows.Getnbcols()<2) 
		ICMTHROW(ERR_INVALID_ARGUMENT,"TreeServices::backwardprice_default3: bad dims flows (4 cols required)") ; 
	
	//	--	Initialisation
	//		
	std::vector<double> tmp1(tree.slice(backwardFrom).size()); 
	std::vector<double>	tmp2(tree.slice(backwardFrom).size()); 
	std::vector<double>::iterator root1=tmp1.begin(); 
	std::vector<double>::iterator root2=tmp2.begin(); 

	//	The states are setted to their initial value 
	for(unsigned int i=0;i<tree.slice(backwardFrom).size();i++) 
		tmp1[i]=statesFrom(i,0); 

	//	--	Iteration 
	int iTree ;	// for the tree
	//int iFlow ; // for the flows
	
	long signedTo(backwardTo);		
	long signedFrom(backwardFrom);

	int i0 = 0 ;
	
	for(iTree=signedFrom-1L; iTree>=signedTo;iTree--) 
	{
		//	--	Backward computation of root2
		StdBinomialTree::const_slice_t aSlice(tree.slice(iTree)); 
		StdBinomialTree::const_slice_t::const_iterator it = aSlice.begin(); 
		std::vector<double>::iterator c1(root1); 
		std::vector<double>::iterator c2(root2); 
		while (it!=aSlice.end()) 
		{
			StdBinomialTree::node_elt_t::const_child_iterator child1(it->begin()); 
			StdBinomialTree::node_elt_t::const_child_iterator child2(it->begin()+1); 
	
			double local = exp( -it->data().h * tree.DT(iTree) ); 
			
			double ctgPart =
					(1.-local)
				*	flows(iTree,1) 
				*	discCurve.DiscountPrice(flows(iTree,0))/discCurve.DiscountPrice(tree.Time(iTree)) ; 

			double nonctgPart =
					local
				*	( child1->first *  (*c1) +	child2->first * (*(c1+1)) )
				*	discCurve.DiscountPrice(tree.Time(iTree+1))/discCurve.DiscountPrice(tree.Time(iTree)) ;

			*c2 = ctgPart + nonctgPart ;

			if(iTree<=2)
			{
				ret_greeks(i0,0) = *c2 ;
				i0++ ;
			}			

			++c2; 
			++c1; 
			++it; 
		}
		//	--	At the end of the loop, we swap root1 and root2
		//		so that root1 reprensents once again the latest computed values
		std::swap(root1,root2); 
	}


	//	--	After the loop, root1 first elements is the computed value
	//		
	std::vector<double>::iterator c1(root1); 
	for(unsigned int j=0;j<tree.slice(backwardTo).size();j++) 
		ret(j,0)= *(c1++) ;
}

//	-----------------------------------------------------------------
// static 
void 
TreeServices::backwardprice_BN(const StdBinomialTree& tree,
							   unsigned long backwardFrom,
							   unsigned long backwardTo,
							   const StdBinomialTree& flows1,
							   const ICM_QMatrix<double>&PayDates1,
							   const ICM_QMatrix<double>&flows2,
							   ARM_ZeroCurve&discCurve,
							   const ICM_QMatrix<double>&statesFrom,								 
							   ICM_QMatrix<double>&ret,
							   ICM_QMatrix<double>&ret_greeks)
{

	//	--	Can't compute	
	if (tree.depth()<=1) 
		ICMTHROW(ERR_INVALID_ARGUMENT,"TreeServices::backwardprice_BN: degenerated tree (depth="<<tree.depth()<<")"); 
	if (backwardTo>backwardFrom) 
		ICMTHROW(ERR_INVALID_ARGUMENT,"TreeServices::backwardprice_BN: can't backword from slice "<<backwardFrom<<" to slice "<<backwardTo); 
	if (backwardFrom>=tree.depth()) 
		ICMTHROW(ERR_INVALID_ARGUMENT,"TreeServices::backwardprice_BN: startSlicePosition "<<backwardFrom<<" out of bounds (" << tree.depth() ); 
	if (statesFrom.Getnbrows() < tree.slice(backwardTo).size() ) 
		ICMTHROW(ERR_INVALID_ARGUMENT,"TreeServices::backwardprice_BN: statesFrom dimension incorrect"); 
	
	//	--	Resize the output 
	ret.Resize(tree.slice(backwardTo).size(),1); 
	ret_greeks.Resize(6,1) ;

	//	--	No flows values
	if (flows2.Getnbcols()<=1) 
		ICMTHROW(ERR_INVALID_ARGUMENT,"TreeServices::backwardprice_BN: bad dims flows (2 cols required)") ; 
	
	//	--	Initialisation
	//		
	std::vector<double> tmp1(flows1.slice(backwardFrom).size()); 
	std::vector<double>	tmp2(flows1.slice(backwardFrom).size()); 
	std::vector<double>::iterator root1=tmp1.begin(); 
	std::vector<double>::iterator root2=tmp2.begin(); 

	//	The states are setted to their initial value 
	for(unsigned int i=0;i<flows1.slice(backwardFrom).size();i++) 
		tmp1[i]=statesFrom(i,0); 

	//	--	Iteration 
	int iTree ;	// for the tree
	//int iFlow ; // for the flows
	
	long signedTo(backwardTo);		
	long signedFrom(backwardFrom);

	unsigned int i0 = 0;

	for(iTree=signedFrom-1L; iTree>=signedTo;iTree--) 
	{
		//	--	Backward computation of root2
		StdBinomialTree::const_slice_t aFlowSlice(flows1.slice(iTree)); 
		StdBinomialTree::const_slice_t aSlice(tree.slice(iTree)); 
		StdBinomialTree::const_slice_t::const_iterator it = aSlice.begin(); 
		StdBinomialTree::const_slice_t::const_iterator Flowit = aFlowSlice.begin(); 
		std::vector<double>::iterator c1(root1); 
		std::vector<double>::iterator c2(root2); 
		while (it!=aSlice.end() && Flowit!=aFlowSlice.end()) 
		{
			StdBinomialTree::node_elt_t::const_child_iterator child1(it->begin()); 
			StdBinomialTree::node_elt_t::const_child_iterator child2(it->begin()+1); 
	
			double local = exp( -it->data().h * tree.DT(iTree) ); 
			double ExerciceValue = Flowit->data().h ;
			double ContiuationValue =
						local
						* (child1->first *  (*c1) +	child2->first * (*(c1+1)))
						* discCurve.DiscountPrice(tree.Time(iTree+1))/discCurve.DiscountPrice(tree.Time(iTree))
						+ 
						(1-local)
						*flows2(iTree,0)
						* discCurve.DiscountPrice(flows2(iTree,1))/discCurve.DiscountPrice(tree.Time(iTree));
									  
			*c2 = std::_cpp_max(ExerciceValue, ContiuationValue);
			
			if(iTree<=2)
			{
				ret_greeks(i0,0) = *c2 ;
				i0++ ;
			}			
			
			++c2; 
			++c1; 
			++it;
			++Flowit; 
		}
		//	--	At the end of the loop, we swap root1 and root2
		//		so that root1 reprensents once again the latest computed values
		std::swap(root1,root2); 
	}


	//	--	After the loop, root1 first elements is the computed value
	//		
	std::vector<double>::iterator c1(root1); 
	for(unsigned int j=0;j<tree.slice(backwardTo).size();j++) 
		ret(j,0)= *(c1++) ;
}

void 
TreeServices::dichotomic_research(const StdBinomialTree& tree,
								 double yearterm,
								 unsigned long minindex, 
								 unsigned long maxindex,
								 ICM_QMatrix<double>&ret)
{

	//	--	Can't compute	
	if (tree.depth()<=1) 
		ICMTHROW(ERR_INVALID_ARGUMENT,"TreeServices::dichotomic_research: degenerated tree (depth="<<tree.depth()<<")"); 
	if (minindex>maxindex) 
		ICMTHROW(ERR_INVALID_ARGUMENT,"TreeServices::dichotomic_research: "<<minindex<<" is greater than "<<maxindex); 
	if (maxindex>=tree.depth()) 
		ICMTHROW(ERR_INVALID_ARGUMENT,"TreeServices::dichotomic_research: "<<maxindex<<" out of bounds (" << tree.depth() ); 
	if (yearterm > tree.Time(maxindex) || yearterm < tree.Time(minindex))
		ICMTHROW(ERR_INVALID_ARGUMENT,"TreeServices::dichotomic_research: "<<yearterm<<" out of bounds (" << tree.depth() ); 

	//	--	Resize the output 
	ret.Resize(2,1); 

	// Reserach of the 2 slices which enclose yearterm

	if(yearterm > tree.Time(minindex) && yearterm < tree.Time(maxindex) && maxindex-minindex == 1 )
	{
		ret(0,0) = minindex ;
		ret(1,0) = maxindex ;
	}
	else
	{
		if (yearterm == tree.Time(maxindex))
		{
			ret(0,0)= maxindex-1 ; // slice just before yearterm
			ret(1,0)= maxindex	 ; // slice after or equal to yearterm
		}
		else
		{
			if (yearterm == tree.Time(minindex))
			{
				ret(0,0) = std::_cpp_max((double)(minindex - 1),0.)	;	// slice just before yearterm
				ret(1,0) = minindex 	;						// slice after or equal to yearterm
			}
			else
			{
				unsigned int mid = (minindex + maxindex)/2 ;
				if(yearterm >= tree.Time(mid))
					dichotomic_research(tree, yearterm, mid, maxindex, ret) ;
				else
				{
					if(yearterm < tree.Time(mid))
						dichotomic_research(tree, yearterm, minindex, mid, ret) ;
				}
			}	
		}
	}
}