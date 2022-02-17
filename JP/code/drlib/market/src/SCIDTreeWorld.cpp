#include "edginc/config.hpp"
#include "edginc/SCIDTreeWorld.hpp"
#include "edginc/EffectiveCurve.hpp"
#include "edginc/SwapTool.hpp"
#include "edginc/QuadraticProg.hpp"

DRLIB_BEGIN_NAMESPACE

void SCIDtreeWorld::getInverseDistanceBetweenLeaves(vector < smartPtr<SimpleTree> > &leaves, DoubleArray &invDist)
{
	invDist.resize(leaves.size());
	fill(invDist.begin(),invDist.end(),0.0);
	for (int i=0; i<leaves.size()-1; ++i)
	{
		if (leaves[i+1]->parent==leaves[i]->parent)
		{
			double a = leaves[i]->ft;
			double b = leaves[i+1]->ft;
			QLIB_VERIFY(a<b, "worlds not in order");
			invDist[i] = 1.0/(b-a);
		}
	}
}

void SCIDtreeWorld::globalFeed(DoubleMatrix &a_t, DoubleArray &w)
{
	tree = SimpleTreeSP (new SimpleTree());
	SimpleTreeSP toBeFed;
	int nbWorlds = a_t.numRows();
	int nbMaturities = a_t.numCols();
	QLIB_VERIFY(nbMaturities>=treeDates.size(), "wrong number of maturities");
	QLIB_VERIFY(nbWorlds == w.size(), "wrong number of worlds");
	for (int i=0; i<nbWorlds; ++i)
	{
		toBeFed = tree;
        for (int j=0; j<treeDates.size(); ++j)
		{
			double nextPos = a_t[j][i];
			size_t k=0;
			while ( (k<toBeFed->children.size()) && (abs(toBeFed->children[k]->ft-nextPos)>1e-10) ) ++k;
			if (k==toBeFed->children.size())
			{
				if (j==0)
					toBeFed->children.push_back(SimpleTreeSP (new SimpleTree(toBeFed, 0, treeTimes[j], nextPos, nextPos, w[i],j+1)));
				else
					toBeFed->children.push_back(SimpleTreeSP (new SimpleTree(toBeFed, treeTimes[j-1], treeTimes[j], toBeFed->ft, nextPos, w[i],j+1)));
				toBeFed = toBeFed->children.back();
			}
			else
			{
				toBeFed = toBeFed->children[k];
				toBeFed->weight += w[i];
			}
		}
	}
}

void SCIDtreeWorld::feed(size_t gen, DoubleArrayArray &newPos)
{
	QLIB_VERIFY(gen>=0 && int(gen)<treeTimes.size(), "wrong generation in sequentialFeed");
	if (gen==0)
	{
		QLIB_VERIFY(newPos.size()>0, "nothing in newPos");
		int size = newPos[0].size();
		QLIB_VERIFY(size>0, "nothing in newPos[0]");
		tree = SimpleTreeSP (new SimpleTree());
		tree->children.resize(size);
		double w = tree->weight/double(size);
		for (int i=0; i<size; i++)
			tree->children[i] = SimpleTreeSP (new SimpleTree(tree, 0, treeTimes[0],newPos[0][i],newPos[0][i], w,1));
	}
	else
	{
		vector<SimpleTreeSP> leaves;
		getLeaves(leaves,gen);
		QLIB_VERIFY(leaves.size()>0, "nothing to feed in sequentialFeed");
		QLIB_VERIFY(leaves.size()==newPos.size(), "newPos of the wrong size");
		for (size_t i=0; i<leaves.size(); ++i)
		{
			double w = leaves[i]->weight / newPos[i].size(); 
			for (int j=0; j<newPos[i].size(); ++j )
				leaves[i]->children.push_back(SimpleTreeSP (new SimpleTree( leaves[i],
																			treeTimes[gen-1],
																			treeTimes[gen],
																			leaves[i]->ft,
																			newPos[i][j],
																			w,
																			gen+1)));
		}
	}
}

void SCIDtreeWorld::particularFeed(vector<SimpleTreeSP> &leaves, const DoubleArray &ratios, double minbound, double maxbound)
{
	QLIB_VERIFY(leaves.size()==ratios.size(), "not corresponding number of leaves and ratios");
	double t,nextFt,ft = minbound;
	maxbound *= 0.8;
	for (size_t i=0; i<leaves.size(); ++i)
	{
		t = treeTimes[leaves[i]->generation];

		if (i<leaves.size()-1) nextFt = Maths::min(leaves[i]->ft*ratios[i],maxbound);
						  else nextFt = maxbound;
		ft = Maths::max(ft*0.8,Maths::min(leaves[i]->ft*ratios[i],nextFt/0.8));

		double weight = leaves[i]->weight;
		if (leaves[i]->children.size()!=1) 
		{
			leaves[i]->children.clear();
			leaves[i]->children.push_back(SimpleTreeSP(new SimpleTree(leaves[i], t , ft, weight)));
		}
		else leaves[i]->children[0]->set(leaves[i], t , ft, weight);
		leaves[i] = leaves[i]->children[0];
	}
}



void SCIDtreeWorld::readTree(DoubleMatrix &valuesAtTreeTime)
{
	vector< SimpleTreeSP > leaves;
	getLeaves(leaves);
	valuesAtTreeTime.resize(leaves.size(), treeTimes.size()+1);
	SimpleTreeSP upTree;
	for (size_t i=0; i<leaves.size(); i++)
	{
		upTree = leaves[i];
		valuesAtTreeTime[i][0] = upTree->weight;
		for (int j=treeTimes.size()-1; j>=0; j--)
		{
			valuesAtTreeTime[i][j+1] = upTree->ft;
			upTree = upTree->parent;
		}
	}
}

void SCIDtreeWorld::setTreeTimes(DoubleArray &_treeTimes, DateTime _today)
{ 
	today = _today;
	treeTimes = _treeTimes; 
	treeDates.resize(treeTimes.size());
	for (int i=0; i<treeDates.size(); i++)
		treeDates[i] = today.rollDate(floor(treeTimes[i]*365+0.5));
}

void SCIDtreeWorld::setTreeDates(DateTimeArray &_treeDates, DateTime _today)
{ 
	today = _today;
	treeDates = _treeDates; 
	treeTimes.resize(treeDates.size());
	for (int i=0; i<treeDates.size(); i++)
		treeTimes[i] = today.yearFrac(treeDates[i]);
}


void SCIDtreeWorld::setTELtimes(DoubleArray &times)
{
	size_t maxNbGeneration = treeTimes.size();
	QLIB_VERIFY(maxNbGeneration>0, "Empty tree!!");
	TELtimes = times;   // TELtimes[TELindexes[gen-1]+1 to TELindexes[gen] included] are in generation gen
	TELdates.resize(times.size());
	for (int i=0; i<times.size(); i++)
		TELdates[i] = today.rollDate( floor(365.0*times[i]+0.5) );
	TELindexes.resize(maxNbGeneration+1);
	int j=0;
	TELindexes[0]=-1;
	for (size_t i=1; i<=maxNbGeneration; i++)
	{
		while ( (j<times.size()) && (times[j]<=treeTimes[i-1]) ) ++j;
		TELindexes[i]=j-1;
	}
}

void SCIDtreeWorld::setTrancheData( DoubleArray &_kmin, 
									DoubleArray &_kmax, 
									int _coupon,
									YieldCurveSP _discount,
									DayCountConventionSP _dcc)
{
	size_t maxNbGeneration = treeTimes.size();
	QLIB_VERIFY(maxNbGeneration>0, "Empty tree!!");
	kmin = _kmin;
	kmax = _kmax;
	discount = _discount;
	coupon = _coupon;
	dcc = _dcc;
}


void SCIDtreeWorld::readTEL(SimpleTreeSP &leave, DoubleMatrix & TEL)
{
	size_t gen = leave->generation, index;
	if (gen==0) return;
	int nbTranches = leave->TEL->numRows();
	TEL.resize(TELindexes[gen]+1, nbTranches);
	SimpleTreeSP copyLeave = leave;
	while (copyLeave->generation!=0)
	{
		gen = copyLeave->generation;
		index = TELindexes[gen-1]+1;
		int col = copyLeave->TEL->numCols();
		int row = copyLeave->TEL->numRows();
		QLIB_VERIFY( (col==TELindexes[gen] - TELindexes[gen-1] ) && (row = nbTranches) , "TEL of wrong size in readTEL");
		for (int j=index; j<=TELindexes[gen]; j++)
			for (int k=0; k<nbTranches; k++)
				TEL[j][k] = (*(copyLeave->TEL))[j-index][k];
		copyLeave = copyLeave->parent;
	}
}


void SCIDtreeWorld::readPEL(SimpleTreeSP &leave, DoubleArray & PEL)
{
	size_t gen = leave->generation;
	if (gen==0) return;
	PEL.resize(gen);
	SimpleTreeSP copyLeave = leave;
	while (copyLeave->generation!=0)
	{
		--gen;
		PEL[gen] = copyLeave->PEL;
		copyLeave = copyLeave->parent;
	}
}

void SCIDtreeWorld::readLegs(SimpleTreeSP &leave, DoubleMatrix& DL, DoubleMatrix& RA)
{
	size_t gen = leave->generation;
	if (gen==0) return;
	DL.resize(kmin.size(), gen);
	RA.resize(kmin.size(), gen);
	SimpleTreeSP copyLeave = leave;
	while (copyLeave->generation!=0)
	{
		--gen;
		QLIB_VERIFY(copyLeave->DL->size()==kmin.size(), "wrong number of tranche in the tree");
		for (int t=0; t<kmin.size(); ++t)
		{
			DL[t][gen] = (*copyLeave->DL)[t];
			RA[t][gen] = (*copyLeave->RA)[t];
		}
		copyLeave = copyLeave->parent;
	}
}


void SCIDtreeWorld::TELtoLegs(int lastTELdate,
				  DoubleMatrix &TEL,
				  DateTime maturity,
				  DoubleArray &DL, 
				  DoubleArray &RA)
{
	QLIB_VERIFY(TEL.numRows()==kmin.size(), "wrong number of tranche expected losses");
	QLIB_VERIFY(TEL.numCols()>lastTELdate, "not enough TEL in TELtoLegs");
	QLIB_VERIFY(lastTELdate<TELdates.size(), "lastTELdates out of range");
	DL.resize(kmin.size());
	RA.resize(kmin.size());

	CashFlowArray cf;
	cf = SwapTool::cashflows(today, maturity, false, 1.0, coupon, "M", &(*dcc));
	cf.back().amount -= 1.0; // ugly
	DateTimeArray dates(TELdates.begin(), TELdates.begin()+lastTELdate+1);
	convertTELtoLegs(kmin, kmax, TEL, today, dates, cf, maturity, discount, &DL[0], &RA[0]);
	

}


void SCIDtreeWorld::fromTELtoLegsOneGeneration(vector<SimpleTreeSP> &leaves )
{
	QLIB_VERIFY(leaves.size()>0, "nothing to do in here");
	int gen = leaves[0]->generation;
	DoubleMatrix TEL;
	for (size_t l=0; l<leaves.size(); l++)
	{
		readTEL(leaves[l], TEL);
		TELtoLegs(TELindexes[gen],TEL,treeDates[gen-1],*(leaves[l]->DL),*(leaves[l]->RA));
	}
}

void SCIDtreeWorld::fromTELtoLegsNGeneration(vector<SimpleTreeSP> &leaves, int lastGen )
{
	int gen = leaves[0]->generation;
	if (gen<=lastGen) fromTELtoLegsOneGeneration(leaves);
	if (gen<lastGen) 
		for (size_t i=0; i<leaves.size(); ++i)
			fromTELtoLegsNGeneration(leaves[i]->children, lastGen);
}



double SCIDtreeWorld::findWeights(size_t genLow, size_t genHigh,
								const DoubleMatrix &parSpread,
								const DoubleArray &PEL,
								double stressPEL,
								DoubleArray& smoothing,
								bool deleteZeroWorlds)        // par spread for all tranches and PEL At all MTM Dates	
														  // but we calibrate just for one generation
{
	QLIB_VERIFY(genLow>0 && genHigh<=getMaxNbGeneration() && genLow<=genHigh, "bad generation number");
	vector < SimpleTreeSP > prevLeaves;
	tree->getLeavesUpToGeneration(prevLeaves, genLow-1);
	size_t nbSumWeightsConstraints = prevLeaves.size();
	vector< vector< SimpleTreeSP > > leaves(prevLeaves.size());
	int nbWeights = 0;
	for (size_t l=0; l<nbSumWeightsConstraints; ++l)
	{
		prevLeaves[l]->getLeavesUpToGeneration(leaves[l], genHigh);
		nbWeights += leaves[l].size();
		QLIB_VERIFY(leaves[l].size()>0,  "no descendant for this leave");
	}
	int nbTranches = parSpread.numCols();

	QLIB_VERIFY(nbWeights>0, "no worlds!!");

	vector<DoubleArray> MTMs(nbWeights);
	vector<DoubleArray> PELs(nbWeights);
	DoubleMatrix  weightConstraints(nbSumWeightsConstraints, nbWeights);
	weightConstraints.fill(0.0);
	DoubleArray sumWeightsConstraints(nbSumWeightsConstraints);

	// collect Default Legs and Risky Annuity and PEL
	int countLeaf = 0;
	for (size_t l=0; l<nbSumWeightsConstraints; ++l)
		for (size_t m=0; m<leaves[l].size(); ++m)
		{
			SimpleTreeSP leaf=leaves[l][m];
			while (leaf->generation>=genLow)
			{
				// dealing with MTM constraints
				for (int tr=0; tr<nbTranches; ++tr)
				{
					double ps = parSpread[tr][leaf->generation-1];
					double upf = 0.0; // to modified when supported
					double DL = (*(leaf->DL))[tr];
					double RA = (*(leaf->RA))[tr];
					MTMs[countLeaf].push_back( DL - ps*RA - upf );
				}
				// dealing with PEL constraints
				PELs[countLeaf].push_back(leaf->PEL - PEL[leaf->generation-1] );

				leaf = leaf->parent;
			}
			++countLeaf;
		}

	// starts with the equality constraints
	countLeaf = 0;
	for (size_t i=0; i<nbSumWeightsConstraints; i++)
	{
		sumWeightsConstraints[i] = prevLeaves[i]->weight;
		for (size_t j=0; j<leaves[i].size(); j++)
		{
			weightConstraints[i][countLeaf] = 1.0;
			++countLeaf;
		}
	} // sum of weights equality constraints
	DoubleArray g(nbWeights, 0.0);
	DoubleArraySP newWeights;
	double sol = calibrateQuadProg(newWeights, MTMs, PELs, weightConstraints, sumWeightsConstraints, stressPEL, smoothing);

	int count = 0;
	for (size_t l=0; l<nbSumWeightsConstraints; ++l)
	{
		for (size_t m=0; m<leaves[l].size(); ++m)
		{
			leaves[l][m]->changeWeights((*newWeights)[count]);
			++count;
		}
		prevLeaves[l]->addUpWeights(leaves[l]);
	}
	if (deleteZeroWorlds)
	{
		for (size_t i=0; i<nbSumWeightsConstraints; i++)
			prevLeaves[i]->deleteZeroChildren();
	}

	return sol;
}

///////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////
SimpleTree::SimpleTree()
{
	set();
}

SimpleTree::SimpleTree(SimpleTree& _tree, SimpleTreeSP &parent)
{
	set(_tree,parent);
}

SimpleTree::SimpleTree(smartPtr<SimpleTree> &_parent, double _t, double _ft, double w)
{ 
	set(_parent, _t, _ft, w); 
}
SimpleTree::SimpleTree(smartPtr<SimpleTree> &_parent, double _s, double _t, double _fs, double _ft, double _w, size_t _generation) 
{ 
	set(_parent, _s, _t, _fs, _ft, _w, _generation); 
}

void SimpleTree::set() 
{ 
	weight=1.0; generation=0;  
};

void SimpleTree::set(SimpleTree& _tree, SimpleTreeSP &parent)
{
	set(parent, _tree.s, _tree.t, _tree.fs, _tree.ft, _tree.weight, _tree.generation); 
	for (size_t i=0; i<_tree.children.size(); ++i)
	{
		SimpleTreeSP voidtree = SimpleTreeSP(new SimpleTree()); children.push_back(voidtree);
		SimpleTreeSP thistree = SimpleTreeSP(this); children.back()->set(*_tree.children[i],thistree);
	}
}

void SimpleTree::set(smartPtr<SimpleTree> &_parent, double _t, double _ft, double _w) 
{ 
	set(_parent,_parent->t,_t,_parent->ft,_ft,_w, _parent->generation+1); 
};

void SimpleTree::set(SimpleTreeSP &_parent, double _s, double _t, double _fs, double _ft, double _w, size_t _generation) 
{ 
	parent=_parent; 
	s=_s; t=_t; 
	fs=_fs; ft=_ft; 
	weight=_w; 
	generation=_generation; 
	TEL = CDoubleMatrixSP (new DoubleMatrix());
	RA = DoubleArraySP (new DoubleArray());
	DL = DoubleArraySP (new DoubleArray());
}

void SimpleTree::resetWeights(size_t gen)
{
	if (generation<gen) 
	{
		weight=0;
		for (size_t i=0; i<children.size(); ++i)
			children[i]->resetWeights(gen);
	}
}

void SimpleTree::modifyFt(double _ft)
{
	if (generation==1) fs=_ft;
	ft = _ft;
	for (size_t i=0; i<children.size(); ++i)
		children[i]->fs=fs;
}

void SimpleTree::modifyWeight(double _weight)
{
	double ratio = _weight/weight;
	weight = _weight;
	for (size_t i=0; i<children.size(); ++i)
		children[i]->modifyWeight(ratio*children[i]->weight);
}

void SimpleTree::addUpWeights(vector<SimpleTreeSP> &leaves)
{
	QLIB_VERIFY(leaves.size()>0, "no tree elements in AddUpWeight");
	size_t gen = leaves[0]->generation;
	for (size_t i=1; i<leaves.size(); ++i) 
		QLIB_VERIFY(gen==leaves[i]->generation, "wrong generation in AddUpWeight");
	// first reset the weights between tree and leaves
	for (size_t i=0; i<children.size(); ++i)
		children[i]->resetWeights(gen);
	// then add up the weight
	for (size_t i=0; i<leaves.size(); ++i)
	{
		SimpleTreeSP leaf = leaves[i];
		double w = leaf->weight;
		leaf = leaf->parent;
		while (leaf->generation>generation)
		{
			leaf->weight += w;
			leaf = leaf->parent;
		}
		QLIB_VERIFY(leaf.get() == this, "leaves do not descent from tree in AddUpWeight");
	}
}


void SimpleTree::changeWeights(double newWeight)
{
	weight = newWeight;
	size_t nbChildren = children.size();
	for (size_t i=0; i<nbChildren; i++)
		children[i]->changeWeights(newWeight/nbChildren);
}

void SimpleTree::reconcile(double epsilon, bool add)
{
	for (size_t i=0; i<children.size(); ++i)
	{
		IntArray close;
		double fti = children[i]->ft;
		double weighti = children[i]->weight;
		for (size_t j=i+1; j<children.size(); ++j)
		{
			bool nearby;
			if (add) nearby = abs(children[j]->ft-fti)<=epsilon+1e-15;
				else nearby = abs(children[j]->ft/fti-1)<=epsilon+1e-15;
			if ( (children[j]->weight>1e-12) && nearby)
				close.push_back(j);
		}
		if (close.size()>0)
		{
			double pos = weighti*fti;
			double sumweight = weighti;
			for (int k=0; k<close.size(); ++k)
			{
				sumweight += children[close[k]]->weight;
				pos += children[close[k]]->weight*children[close[k]]->ft;
				children[close[k]]->changeWeights(0.0);
			}
			children[i]->changeWeights(sumweight);
			QLIB_VERIFY(sumweight>0, "zero weight -> division by zero!");
			pos /= sumweight;
			children[i]->fs=pos;
			children[i]->changePositions(pos,pos,pos);
			deleteZeroChildren();
		}
	}
}


double add_mult(double x, double y, double add, double lowbound, double highbound)
{
	double z;
	if (add) z = x+y;
	else z = x*(1+y);
	return Maths::max(lowbound,Maths::min(highbound,z));
}

void SimpleTree::CreateNearbyLeaves(vector<SimpleTreeSP> &newLeaves, DoubleArray &epsilons, double minbound, double maxbound, bool add)
{
	DoubleArray newft(epsilons.size());
	QLIB_VERIFY(generation>0, "we can't deal with the common ancester of the tree");
	DoubleArray oldft(parent->children.size());
	for (int j=0; j<oldft.size(); ++j) oldft[j]= parent->children[j]->ft;
	IntArray whichEpsilons;
	for (int j=0; j<epsilons.size(); ++j)
	{
		newft[j] = add_mult(ft,epsilons[j],add,minbound,maxbound);
		bool found = false; 
		for (int i=0; i<oldft.size() && !found; ++i)
			if (abs(newft[j]-oldft[i])<1e-10) found=true; 
		if (!found) 
			whichEpsilons.push_back(j);
	};
	modifyWeight(weight/(1.0 + whichEpsilons.size()));

	for (int j=0; j<whichEpsilons.size(); ++j)
	{
		parent->children.push_back( SimpleTreeSP(new SimpleTree(*this,parent)) );
		parent->children.back()->modifyFt(newft[whichEpsilons[j]]);
		newLeaves.push_back(parent->children.back());
	}
}

void SimpleTree::CreateNearbyLeavesMultipleGens(vector< vector<SimpleTreeSP> > &newLeaves, size_t maxGen, DoubleArray &epsilons, double minbound, double maxbound, bool add)
{
	if (generation<=maxGen)
	{	
		int newLeavesOldSize = newLeaves[generation].size();
		CreateNearbyLeaves(newLeaves[generation], epsilons, minbound, maxbound, add);
		size_t newLeavesNewSize = newLeaves[generation].size();
		size_t nbChildren = children.size();
		for (size_t j=0; j<nbChildren; ++j)
			children[j]->CreateNearbyLeavesMultipleGens(newLeaves, maxGen, epsilons, minbound, maxbound, add);
		for (size_t i=newLeavesOldSize; i<newLeavesNewSize; ++i)
			for (size_t j=0; j<newLeaves[generation][i]->children.size(); ++j)
				newLeaves[generation][i]->children[j]->CreateNearbyLeavesMultipleGens(newLeaves, maxGen, epsilons, minbound, maxbound, add);
	}
}

void SimpleTree::changePositions(double x, double minbound, double maxbound)
{
	if (generation==1)
	{
		double fst = Maths::min(maxbound, Maths::max(minbound, x));
		fs = fst;
		ft = fst;
		for (size_t j=0; j<children.size(); ++j) children[j]->fs = fst;
	}
	else
	{
		ft = Maths::min(maxbound, Maths::max(minbound, fs*x));
		for (size_t j=0; j<children.size(); ++j)
			children[j]->fs = ft;
	}
}


void SimpleTree::deleteZeroChildren()
{
	if (children.size()==0) return;

	vector<SimpleTreeSP> newChildren;
	for (size_t i=0; i<children.size(); i++)
	{
		if (children[i]->weight>1e-12) 
		{
			children[i]->deleteZeroChildren();
			newChildren.push_back(children[i]);
		}
	}
	children = newChildren;
}

void SimpleTree::feed(double *times, size_t size, size_t gen, vector<double> &slopes)
{
	if (size>0)
	{
		size_t i=0;
		while ((i<size)&&(times[i]<=t)) i++;
		if (i<size)
		{
			children.resize(slopes.size());
			for (size_t j=0; j<slopes.size(); j++)
			{
				SimpleTreeSP thisTree = SimpleTreeSP(this);
				children[j] = SimpleTreeSP (new SimpleTree( thisTree,t,times[i],ft,ft*slopes[j], weight/double(slopes.size()), gen));
				children[j]->feed(&times[i+1], size-i-1, gen+1,slopes);
			}
		}
	}
}


void SimpleTree::getWorldsFromTree(DoubleArray &x, DoubleArray &y, LinearInterpolantSP *f)
{
	QLIB_VERIFY(generation>0, "bad generation in getWorldsFromTree");
	x.resize(generation);
	y.resize(generation);
	size_t i = x.size()-1;
	y[i] = ft;
	x[i] = t;
	SimpleTreeSP tempTree = parent;
	while (tempTree->generation > 0 )
	{
		--i;
		y[i] = tempTree->ft;
		x[i] = tempTree->t;
		tempTree = tempTree->parent;
	}
	if (f!=0)
		*f = LinearInterpolantSP(new LinearInterpolant(x,y));
}

void SimpleTree::getWorldsFromTree(LinearInterpolantSP &f)
{
	QLIB_VERIFY(generation>0, "bad generation in getWorldsFromTree");
	DoubleArray x(generation), y(generation);
	size_t i = generation-1;
	y[i] = ft;
	x[i] = t;
	SimpleTreeSP tempTree = parent;
	while (tempTree->generation > 0 )
	{
		--i;
		y[i] = tempTree->ft;
		x[i] = tempTree->t;
		tempTree = tempTree->parent;
	}
	f = LinearInterpolantSP(new LinearInterpolant(x,y));
}


void SCIDtreeWorld::getWorldsFromTree(vector<LinearInterpolantSP> &f, vector<SimpleTreeSP> &leaves) 
{ 
	DoubleArray x, y; 
	f.resize(leaves.size());
	for (size_t i=0; i<leaves.size(); ++i) 
		leaves[i]->getWorldsFromTree(x,y,&f[i]); 
};

void SCIDtreeWorld::getWorldsFromTree(DoubleArray &x, vector<DoubleArray> &y, vector<SimpleTreeSP> &leaves) 
{ 
	DoubleArray xx;
	y.resize(leaves.size());
	if (leaves.size()>0) leaves[0]->getWorldsFromTree(x,y[0]); 
	xx=x;
	for (size_t i=0; i<leaves.size(); ++i) 
	{
		leaves[i]->getWorldsFromTree(x,y[i]); 
		QLIB_VERIFY(x==xx, "leaves not of the same generation in getWorldsFromTree");
	}
};


void SimpleTree::getLeavesUpToGeneration(vector < SimpleTreeSP > &leaves, size_t gen)
{
	if (gen==generation) leaves.push_back(SimpleTreeSP(this));
	else if (gen>generation) 
		for (size_t i=0; i<children.size(); i++) 
			children[i]->getLeavesUpToGeneration(leaves,gen);
}

bool SimpleTree::check()
{
	if (children.size()==0) return true;
	double w=0.0;
	for (size_t i=0; i<children.size(); ++i)
	{
		w += children[i]->weight;
		if (children[i]->parent.get() != this) return false;
		if (!children[i]->check()) return false;
	}
	return (abs(w-weight)<1e-10);
}

void SCIDtreeWorld::convertTELtoLegs ( const DoubleArray &kmin,
						const DoubleArray &kmax,
						const DoubleMatrix &TEL,
						const DateTime &today,
						const DateTimeArray &TELdates,
						const CashFlowArray &cf,
						const DateTime &maturity,
						const YieldCurveSP &discount,
						double *DL,   // (O)
						double *RA)  // (O)  must both be of size kmin.size()
{
	QLIB_VERIFY(TEL.numRows()==kmin.size(), "wrong number of tranche expected losses");
	QLIB_VERIFY(TEL.numCols()>=TELdates.size(), "not enough TEL in TELtoLegs");
	bool firstDateIsToday = (TELdates[0]==today);
	DateTimeArray dates;
	if (!firstDateIsToday) dates.push_back(today);
	int lastDate=0;
	while ( (lastDate<TELdates.size()) && (TELdates[lastDate]<=maturity) ) ++lastDate;

	dates.insert(dates.end(), TELdates.begin(), TELdates.begin()+lastDate);

	DoubleArray riskyDiscount(dates.size());
	riskyDiscount[0] = 1.0;
	for (int i=0; i<kmin.size(); ++i)
	{
		int index=0;
		if (!firstDateIsToday) ++index;
		double scale = 1.0 / (kmax[i] - kmin[i]);
		for (int j=index; j<dates.size(); ++j)
			riskyDiscount[j] = Maths::max(1e-12,1.0 - scale*TEL[j-index][i]);

		EffectiveCurve trancheCurve(today, discount, dates, riskyDiscount, EffectiveCurve::FLAT_FORWARD);
		DL[i] = trancheCurve.protectionPV(today,today,maturity,IDiscountCurveRisky::RECOVER_1);
		RA[i] = trancheCurve.annuityPV(cf,today,IDiscountCurveRisky::RECOVER_1);
	}
}

double  SCIDtreeWorld::calibrateQuadProg ( DoubleArraySP &newWeights,
		  				    vector<DoubleArray> &MTMs,
							vector<DoubleArray> &PELconstraints,
							DoubleMatrix &weightConstraints,
							DoubleArray &sumWeightsConstraints,
							double stressPEL,
							DoubleArray& smoothing)
{
	size_t nbWeights = MTMs.size();
	QLIB_VERIFY(PELconstraints.size()==nbWeights && nbWeights>0, "inconsistent number of weights in calibrateQuadProg");
	size_t nbSumWeightsConstraints = sumWeightsConstraints.size();
	QLIB_VERIFY(weightConstraints.numCols() == nbSumWeightsConstraints 
		&& weightConstraints.numRows() == nbWeights, "wrong number of constraints in calibrateQuadProg");
	int nbTranches = MTMs[0].size();
	for (size_t i=1; i<nbWeights; ++i) QLIB_VERIFY(MTMs[i].size()==nbTranches, "problem with the number of tranches in calibrateQuadProg" );
	int nbPEL = PELconstraints[0].size();
	for (size_t i=1; i<nbWeights; ++i) QLIB_VERIFY(PELconstraints[i].size()==nbPEL, "problem with the number of portfolio expected losses in calibrateQuadProg" );

	double norm = 0.0;
	DoubleArray Hessian(nbWeights*nbWeights, 0.0);
	for (size_t i=0; i<nbWeights; ++i)
		for (size_t j=0; j<nbWeights; ++j)
		{
			for (int tr=0; tr<nbTranches; ++tr)
				Hessian[i*nbWeights+j] += MTMs[i][tr]*MTMs[j][tr];
			norm += Maths::square(Hessian[i*nbWeights+j]);
		}
	norm = sqrt(norm);
	QLIB_VERIFY(smoothing.size()==nbWeights, "wrong smoothing coefficient");
	Hessian[0] += smoothing[0];
	for (int i=0; i<nbWeights-1; ++i)
	{
		Hessian[i*nbWeights+i+1] += -smoothing[i];
		Hessian[(i+1)*nbWeights+i] += -smoothing[i];
		Hessian[(i+1)*nbWeights+(i+1)] += smoothing[i]+smoothing[i+1];
	}

	if (stressPEL>0)
	{
		for (size_t i=0; i<nbWeights; ++i)
			for (size_t j=0; j<nbWeights; ++j)
				for (int p=0; p<nbPEL; ++p)
					Hessian[i*nbWeights+j] += stressPEL*PELconstraints[i][p]*PELconstraints[j][p];
	}


	int meq = nbSumWeightsConstraints; // nbMaturities will have to be changed to nbPELmaturities of something like this
	if (stressPEL<0) meq += nbPEL;
	int m = meq + nbWeights;  // number of equality plus linear constraints 
	DoubleArray constraints(m*nbWeights, 0.0);
	DoubleArray constraintsEqualTo(m, 0.0);

	// starts with the equality constraints
	for (size_t i=0; i<nbSumWeightsConstraints; i++)
	{
		constraintsEqualTo[i] = sumWeightsConstraints[i];
		for (size_t j=0; j<nbWeights; ++j)
			constraints[i*nbWeights+j] = weightConstraints[i][j];
	} 
	if (stressPEL<0) 
	{
		for (int i=0; i<nbPEL; ++i)
			for (size_t j=0; j<nbWeights; ++j)
				constraints[(nbSumWeightsConstraints+i)*nbWeights+j] = PELconstraints[j][i];
	}
	// then do the inequality constraints, i.e. that the weights are all positive
	for (size_t i=0; i<nbWeights; i++)
		constraints[(meq+i)*nbWeights+i] = 1.0;

	DoubleArray g(nbWeights, 0.0);
	double sol;
	newWeights = DoubleArraySP (new DoubleArray(nbWeights));

	imsl_d_quadratic_prog(
		m, nbWeights, meq, &constraints[0], &constraintsEqualTo[0], &g[0], &Hessian[0],
		IMSL_RETURN_USER, &(*newWeights)[0],
		IMSL_OBJ, &sol,
		0);

	// throw IMSLException if quadratic prog fails
	IMSLException::throwIMSLExceptionIfError();
	//	return QuadraticProg::minimize(m, nbWeights, meq, constraints, constraintsEqualTo, g , Hessian);
	return sol;
}



DRLIB_END_NAMESPACE
