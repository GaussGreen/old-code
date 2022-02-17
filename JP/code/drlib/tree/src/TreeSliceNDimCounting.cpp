#include "edginc/config.hpp"
#include "edginc/TreeSlice.hpp"
#include "edginc/CountingTree.hpp"

DRLIB_BEGIN_NAMESPACE

/******************************** TreeSliceNDimCounting ********************************/

TreeSliceNDimCounting::TreeSliceNDimCounting(const Range& range, const string& curveToDEV) 
: range(range), dim(-1), values(0), deleteWhenDestroyed(true), curveToDEV(curveToDEV)
{}

TreeSliceNDimCounting::~TreeSliceNDimCounting() {resize(-1);}

TreeSliceSP TreeSliceNDimCounting::clone( bool copyValues) const {
    TreeSliceSP slice(new TreeSliceNDimCounting(range, curveToDEV));
    if (copyValues) 
        *slice = *this;
    return slice;
}

/* Resize the slice's array and initialize values to 0. */
void TreeSliceNDimCounting::resize(int n) {
    if (values) {
        double * toDelete = values;
        values = 0;
        if (deleteWhenDestroyed)
            delete[] toDelete;
    }
    if (n>=0) {
        values = new double[n];
        ::memset( values, 0, n * sizeof(double) );
    }
    theValues = values; // for debug purposes
}

void TreeSliceNDimCounting::takeFrom(TreeSliceNDimCounting &s) {
    if (values != s.values) resize(-1); // delete values, if appropriate
    values = s.values;
    deleteWhenDestroyed = s.deleteWhenDestroyed;
    s.deleteWhenDestroyed = false;

    dim = s.dim;
    copyTreeSlice(s);
    theValues = values; // for debug purposes
}

void TreeSliceNDimCounting::expand(int newDim) {
    try {
        if (newDim == dim) {
            return;
        }
        if (newDim < dim) {
            throw ModelException(
				"Unable to expand a " + Format::toString(dim) + 
                " dimension slice into a smaller " + Format::toString(newDim) +
                " dimension slice");
        }
        if (dim < 0) {
            throw ModelException(
				"Cannot expand undefined and uninitialised slice");
        }
        TreeSliceNDimCounting newSlice(range);
        newSlice.copyTreeSlice(*this);
        newSlice.allocDim(newDim);
        int size = newSlice.valuesSize();
        for (int i=0; i<size; ++i) {
            newSlice.values[i] = 0;
        }
        newSlice += *this;
        takeFrom(newSlice);
    }
    catch (exception& e) {
        throw ModelException(e, "TreeSliceNDimCounting::expand, name=\""+name+"\"");
    }
}

/* How big should the values array be? */
int TreeSliceNDimCounting::valuesSize() const {
    int size = 1;
    /* First compute size due to diffusion state-variables */
    const vector<int>& sz = range.dSize;
    switch (dim) {
		case 5: size = sz[4];
		case 4: size *= sz[3];
        case 3: size *= sz[2];
        case 2: size *= sz[1];
        case 1: size *= sz[0]; 
        case 0: break;
        default:
			{
				int k=-1;
				while (k<dim && k++) {
					size *= sz[k];
				}
			}
    }
    /* Now inflate size due to counting-state variables */
    if (range.hasCStates) {
        const vector<int>& csz = range.cSize;
        int k=-1;
        while (k<dim && k++) {
            if (csz[k]>1) size *= csz[k];
        }
    }
    return size;
}

void TreeSliceNDimCounting::allocDim(int newDim) {
    dim = newDim;
    resize(valuesSize());
}

// assign a constant value to this slice - that makes it zero dimensional
TreeSlice& TreeSliceNDimCounting::operator=(double constant) {
    if (dim!=0)
        allocDim(0);
    values[0]=constant;
    return *this;
}

TreeSlice& TreeSliceNDimCounting::operator=(const TreeSlice& s) {
    const TreeSliceNDimCounting *r = dynamic_cast<const TreeSliceNDimCounting*>(&s);
    if (!r)
        throw ModelException("TreeSliceNDimCounting::operator=(const TreeSlice& s)",
        "Cannot assign slices of different types");
    return operator=(*r);
}

void TreeSliceNDimCounting::testTreeStep() const {
    // since I always take the tree-step directly from the tree,
    // this seems to be redundant.
    //if (treeStep != range.treeStep) {
    //    if (dim==0 && *values==0.)
    //        return; // tolerate constant 0. everywhere
    //    throw ModelException(
    //        "Slice \""+name+"\" has an internal treeStep (" +
    //        Format::toString(treeStep) + 
    //        ") which does not match the current time step (" +
    //        Format::toString(range.treeStep) +
    //        "). You probably forgot to register this slice for DEV.");
    //}
}

TreeSliceNDimCounting& TreeSliceNDimCounting::operator=(const TreeSliceNDimCounting& s) {
    try {
        if (dim!=s.dim) allocDim(s.dim);
        ::memcpy(values, s.values, valuesSize() * sizeof(double));
        return *this;

    } catch (exception& e) {
        throw ModelException(e, "TreeSliceNDimCounting::operator=(const TreeSliceRates& s)");
    }
}

// this assumes _zero_ counting states
int TreeSliceNDimCounting::offset( const vector<int>& d) const {
	int i=0;
    const vector<int>& globalMinLimits = range.dLGlobalLimits;
    const vector<int>& sz = range.dSize;
    const vector<int>& globalMinCLimits = range.cLGlobalLimits;
    const vector<int>& csz = range.cSize;
    bool hasC = range.hasCStates;
	int os = d[0] - globalMinLimits[0];
	while (i++ && i<dim) {
		os  = os*sz[i-1] + d[i] - globalMinLimits[i];
        if (hasC && csz[i-1]>1) {
            os = os*csz[i-1] + 0 - globalMinCLimits[i];
        }
	}
	return os;
}

/* Returns an offset for a given (c,d) state. Note that the c and d states
   are interleaved so that the representation of 1D, 2D, 3D, etc. slices is
   consistent. The innermost values coordinate is always a d-state. */
int TreeSliceNDimCounting::offset(const vector<int>& c, const vector<int>& d) const {   
    const vector<int>& globalMinCLimits = range.cLGlobalLimits;
    const vector<int>& csz = range.cSize;
    const vector<int>& globalMinDLimits = range.dLGlobalLimits;
    const vector<int>& dsz = range.dSize;
	int os = d[0] - globalMinDLimits[0];
    bool hasC = range.hasCStates;
    int i = 0;
	while (i++ && i<dim) {
		os  = os*dsz[i-1] + d[i] - globalMinDLimits[i];
        if (hasC && csz[i-1]>1) {
            os = os*csz[i-1] + c[i] - globalMinCLimits[i];
        }
	}
	return os;
}

int TreeSliceNDimCounting::offset( int i, int j, int k ) const
{
    
#if DEBUG
    /* By the way, this only works with null counting-states! */
    if (range.hasCStates) {
        throw ModelException("You may not use this method if there are counting states");
    }
    if (dim!=3) throw ModelException(
        "Calling offset(i,j,k) on slice \""+name+"\" of dim "+Format::toString(dim));
#endif
    int x3 = k - range.dLGlobalLimits[2]; // - ( range.limits.bot3[ i ][ j ] + range.limits.top3[ i ][ j ] ) / 2;
    int x2 = j - range.dLGlobalLimits[1]; // - ( range.limits.bot2[ i ]      + range.limits.top2[ i ]      ) / 2;
    int x1 = i - range.dLGlobalLimits[0];

    return range.dSize[2] * ( range.dSize[1] * x1 + x2 ) + x3;
}

// calc the offet to add to "values" to access an element in a 2D slice
int TreeSliceNDimCounting::offset( int i, int j ) const
{
    /* By the way, this only works with null counting-states! */
#if DEBUG
    /* By the way, this only works with null counting-states! */
    if (range.hasCStates) {
        throw ModelException("You may not use this method if there are counting states");
    }
    if (dim!=2) throw ModelException(
        "Calling offset(i,j) on slice \""+name+"\" of dim "+Format::toString(dim));
#endif
    int x2 = j - range.dLGlobalLimits[1]; // - ( range.limits.bot2[ i ]      + range.limits.top2[ i ]      ) / 2;
    int x1 = i - range.dLGlobalLimits[0];

    return range.dSize[1] * x1 + x2;
}

// calc the offet to add to "values" to access an element in a 1D slice
int TreeSliceNDimCounting::offset( int i ) const
{
    /* By the way, this only works with null counting-states! */
#if DEBUG
    /* By the way, this only works with null counting-states! */
    if (range.hasCStates) {
        throw ModelException("You may not use this method if there are counting states");
    }
    if (dim!=1) throw ModelException(
        "Calling offset(i) on slice \""+name+"\" of dim "+Format::toString(dim));
#endif
    return      i - range.dLGlobalLimits[0];
}

double TreeSliceNDimCounting::getCentre() const {
    testTreeStep();

    switch (dim) {
        case 3: return values[ offset( 0, 0, 0 ) ];
        case 2: return values[ offset( 0, 0 ) ];
        case 1: return values[ offset( 0 ) ];
        case 0: return values[0];
        case -1:
            throw ModelException("TreeSliceNDimCounting::getCentre", 
                                 "Tree slice (name = " + name +
                                 ") has not been initialized");
        default:
			{
				const static vector<int> zAry(dim, 0.0);
				return values[ offset(zAry) ]; // NB this assumes 0 c-states!
			}
    }
}

static const char *errorSliceNotInitialized = "Slice operation on a non initialized slice: ";

// see usage of iterSeek(...) in calc() below
void TreeSliceNDimCounting::iterSeek(int i) const { // get ready to loop on dimension #0 (1st)
    iter = values;
    switch (dim) {
        case 1: 
            iter += offset(i);
            break;
        case 0: 
            break;
        case -1: 
            throw ModelException(errorSliceNotInitialized + name);
        default: ASSERT(0);
    }
}

void TreeSliceNDimCounting::iterSeek(int i, int j) const { // get ready to loop on dimension #1 (2nd)
    iter = values;
    switch (dim) {
        case 2: 
            iter += offset(i,j); 
            break;
        case 1: 
            iter += offset(i);
            break;
        case 0: 
            break;
        case -1: throw ModelException(errorSliceNotInitialized + name);
        default: ASSERT(0);
    }
}

void TreeSliceNDimCounting::iterSeek(int i, int j, int k) const { // get ready to loop on dimension #2 (3rd)
    iter = values;
    switch (dim) {
        case 3: 
            iter += offset(i,j,k); 
            break;
        case 2: 
            iter += offset(i,j);
            break;
        case 1: 
            iter += offset(i);
            break;
        case 0: 
            break;
        case -1: 
            throw ModelException(errorSliceNotInitialized + name);
        default: ASSERT(0);
    }
}

void TreeSliceNDimCounting::iterSeek(const vector<int>& oAry) const {
    iter = values;
    switch (dim) {
		case 3: 
            iter += offset(oAry[0],oAry[1],oAry[2]); 
            break;
        case 2: 
            iter += offset(oAry[0],oAry[1]);
            break;
        case 1: 
            iter += offset(oAry[0]);
            break;
        case 0: 
            break;
        case -1: 
            throw ModelException(errorSliceNotInitialized + name);
        default:
            throw ModelException("TreeSliceNDimCounting::iterSeek I don't think you should be here!");
            // this can go now?
			//iter += offset(oAry);
			//if (oAry.size()==dim) {
			//} else if (oAry.size()>0 && oAry.size()<dim) {
			//} else {
			//	throw ModelException("Attempt to iterSeek with dimension ("+
			//		Format::toString(oAry.size())+") higher than slice dimension ("+
		//			Format::toString(dim)+").");
		//	}
    }
}

void TreeSliceNDimCounting::iterSeek(const vector<int>& oAry, const vector<int>& dAry) const {
    throw ModelException("TreeSliceNDimCounting::iterSeek This doesn't work!");
    iter = values;
    switch (dim) {
		case 3: 
            iter += offset(oAry[0],oAry[1],oAry[2]); 
            break;
        case 2: 
            iter += offset(oAry[0],oAry[1]);
            break;
        case 1: 
            iter += offset(oAry[0]);
            break;
        case 0: 
            break;
        case -1: 
            throw ModelException(errorSliceNotInitialized + name);
        default:
            throw ModelException("TreeSliceNDimCounting::iterSeek I don't think you should be here!");
            // this can go now?
			//iter += offset(oAry);
			//if (oAry.size()==dim) {
			//} else if (oAry.size()>0 && oAry.size()<dim) {
			//} else {
			//	throw ModelException("Attempt to iterSeek with dimension ("+
			//		Format::toString(oAry.size())+") higher than slice dimension ("+
			//		Format::toString(dim)+").");
			//}
    }
}

// this is not efficient given that I have changed the limits to be rectangular
// rather than elliptical. You should optimise this later. The current version is just
// a hacked presto-chango version of TreeSliceRates.
// Note also that I haven't extended this yet to dimensions higher than 3.
// [Charles Morcom 6/27/06]
TreeSliceSP TreeSliceNDimCounting::calcSmoothStep() const {
    TreeSliceNDimCountingSP r(new TreeSliceNDimCounting(range));
    r->name = "calcSmoothStep("+name+")";
    r->allocDim(dim);
    double *rValues = r->values;
    throw ModelException("Have not implemented step smoothing yet!");
    switch (dim) {
        case 1: {
            int iStart = 0;//range.minLimits[0];
            int iStop  = 0;//range.maxLimits[1];

            for( int i = iStart; i <= iStop; ++i ) {
                double m=0, v=val(i);
                if (i != iStart) 
                    m = max(m, fabs( val(i-1) - v ));
                if (i != iStop)  
                    m = max(m, fabs( val(i+1) - v ));
                rValues[offset(i)] = m;
            }
            break;
        }
        case 2:
        {
            int iStart = 0;//range.minLimits[0];
            int iStop  = 0;//range.maxLimits[0];
            for( int i = iStart; i <= iStop; ++i ) {
                int jStart = 0;//range.minLimits[1];
                int jStop  = 0;//range.maxLimits[1];

                for( int j = jStart; j <= jStop; ++j ) {
                    double m=0, v=val(i,j);

                    if (j != jStart) 
                    {
                        m = max(m, fabs( val(i,j-1) - v ));
                    }
                    if (j != jStop)  
                    {
                        m = max(m, fabs( val(i,j+1) - v ));
                    }
                    if(i != iStart
                    && j >= 0 //range.minLimits[1]
                    && j <= 0) //range.maxLimits[1]) 
                    {
                        m = max(m, fabs( val(i-1,j) - v ));
                    }
                    if(i != iStop
                    && j >= 0 //range.minLimits[1]
                    && j <= 0) //range.maxLimits[1]) 
                    {
                        m = max(m, fabs( val(i+1,j) - v ));
                    }

                    rValues[offset(i,j)] = m;
                }
            }
            break;
        }
        case 3: {
            int iStart = 0;//range.minLimits[0];
            int iStop  = 0;//range.maxLimits[0];
            for( int i = iStart; i <= iStop; ++i ) {
                int jStart = 0;//range.minLimits[1];
                int jStop  = 0;//range.maxLimits[1];

                for( int j = jStart; j <= jStop; ++j ) {
                    int kStart = 0;//range.minLimits[2];
                    int kStop  = 0;//range.maxLimits[2];

                    for( int k = kStart; k <= kStop; ++k ) 
                    {
                        double m=0, v=val(i,j,k);

                        if (k != kStart) 
                        {
                            m = max(m, fabs( val(i,j,k-1) - v ));
                        }
                        if (k != kStop)  
                        {
                            m = max(m, fabs( val(i,j,k+1) - v ));
                        }
                        if(j != jStart
                        && k >= 0 //range.minLimits[2]
                        && k <= 0) //range.maxLimits[2])
                        {
                            m = max(m, fabs( val(i,j-1,k) - v ));
                        }
                        if(j != jStop
                        && k >= 0 //range.minLimits[2]
                        && k <= 0) //range.maxLimits[2])
                        {
                            m = max(m, fabs( val(i,j+1,k) - v ));
                        }
                        if(i != iStart
                        && j >= 0 //range.minLimits[1]
                        && j <= 0 //range.maxLimits[1]
                        && k >= 0 //range.minLimits[2]
                        && k <= 0) //range.maxLimits[2])
                        {
                            m = max(m, fabs( val(i-1,j,k) - v ));
                        }
                        if(i != iStop
                        && j >= 0 //range.minLimits[1]
                        && j <= 0 //range.maxLimits[1]
                        && k >= 0 //range.minLimits[2]
                        && k <= 0) //range.maxLimits[2])
                        {
                            m = max(m, fabs( val(i+1,j,k) - v ));
                        }
                        rValues[offset(i,j,k)] = m;
                    }
                }
            }
            break;
        }
        default:
            throw ModelException("TreeSliceRates::calcSmoothStep", "unsupported slice dimension");
    }
    return r;
}


DRLIB_END_NAMESPACE
