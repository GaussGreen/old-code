#include "edginc/config.hpp"
#include "edginc/TreeSlice.hpp"

DRLIB_BEGIN_NAMESPACE

/******************************** TreeSliceRates::Range ********************************/

void TreeSliceRates::Range::init(int bot1, int top1, // 1st dimension range
            int bot2, int top2, // 2st dimension range
            int bot3, int top3, // 3st dimension range
            int, int)           
{

    bot[0] = bot1;
    top[0] = top1;
    size[0] = top1 - bot1 + 1;
    if( size[0] <= 0 )
        throw ModelException("top1 < bot1 during construction");

    bot[1] = bot2;
    top[1] = top2;
    size[1] = top2 - bot2 + 1;
    if( size[1] <= 0 )
        throw ModelException("top2 < bot2 during construction");

    bot[2] = bot3;
    top[2] = top3;
    size[2] = top3 - bot3 + 1;
    if( size[2] <= 0 )
        throw ModelException("top3 < bot3 during construction");

    int d = SLICE_RATES_MAX_DIMS;
    treeStep = -1;

    while( --d && size[ d ] <= 1 ) 
    {}
    nDim = d + 1;
}

/******************************** TreeSliceRates ********************************/

TreeSliceRates::TreeSliceRates(const RatesSliceRange &range, const string& curveToDEV, int curveToDEVIdx) :
    curveToDEV(curveToDEV), curveToDEVIdx(curveToDEVIdx),
    dim(-1), values(0), deleteWhenDestroyed(true), range(dynamic_cast<const Range&>(range))
{}

TreeSliceRates::~TreeSliceRates() {allocDim(-1);}

TreeSliceSP TreeSliceRates::clone( bool copyValues) const {
    TreeSliceRatesSP slice( new TreeSliceRates(range, curveToDEV, curveToDEVIdx) );
    if (copyValues) {
        slice->allocDim(dim);
        memcpy(slice->values, values, valuesSize()*sizeof(double));
    }
    slice->treeStep = treeStep;
    return slice;
}

void TreeSliceRates::resize(int n) {
    if (values) {
        double * toDelete = values;
        values = 0;
        if (deleteWhenDestroyed)
            delete[] toDelete;
    }
    if (dim>=0 && n>=0) {
        values = new double[n];
        ::memset( values, 0, n * sizeof(double) );
    }
    else values=0;
    iter=0; // so that the printDebug does not fail
    theValues = values; // for debug purposes
}

void TreeSliceRates::takeFrom(TreeSliceRates &s) {
    if (values != s.values) 
        resize(-1);
    values = s.values;
    deleteWhenDestroyed = s.deleteWhenDestroyed;
    s.deleteWhenDestroyed = false;

    dim = s.dim;
    copyTreeSlice(s);
    theValues = values; // for debug purposes
}

void TreeSliceRates::expand(int newDim) {
    try {
        if (newDim == dim) {
            return;
        }
        if (newDim < dim) {
            throw ModelException("Unable to expand a " + Format::toString(dim) + 
                " dimension slice into a smaller " + Format::toString(newDim) +
                " dimension slice");
        }
        if (dim < 0) {
            throw ModelException("Cannot expand undefined and uninitialised slice");
        }
        TreeSliceRates newSlice(range, curveToDEV, curveToDEVIdx);
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
        throw ModelException(e, "TreeSliceRates::expand, name=\""+name+"\"");
    }
}


int TreeSliceRates::valuesSize() const {
    int size = 1;
    switch (dim) {
        case 3: size = range.size[2];
        case 2: size *= range.size[1];
        case 1: size *= range.size[0]; 
        case 0: break;
        default: return 0;
    }
    return size;
}

void TreeSliceRates::allocDim(int newDim) {
    dim = newDim;
    resize(valuesSize());
}

TreeSlice& TreeSliceRates::operator=(double constant) {
    if (dim!=0)
        allocDim(0);
    values[0]=constant;
    treeStep = range.treeStep;  // set slice to current tree time step
    return *this;
}

TreeSlice& TreeSliceRates::operator=(const TreeSlice& s) {
    const TreeSliceRates *r = dynamic_cast<const TreeSliceRates*>(&s);
    if (!r)
        throw ModelException("TreeSliceRates::operator=(const TreeSlice& s)",
        "Cannot assign slices of different types");
    return operator=(*r);
}

void TreeSliceRates::testTreeStep() const {
    if (treeStep != range.treeStep) {
        // skip check for slices that are specifically set treeStep < 0
        // or if not yet in the time loop
        if (treeStep < 0 || range.treeStep < 0)
            return;

        throw ModelException(
            "Slice \""+name+"\" has an internal treeStep (" +
            Format::toString(treeStep) + 
            ") which does not match the current time step (" +
            Format::toString(range.treeStep) +
            "). You probably forgot to register this slice for DEV.");
    }
}

TreeSliceRates& TreeSliceRates::operator=(const TreeSliceRates& s) {
    try {
        if (dim!=s.dim)
            allocDim(s.dim);

        ::memcpy(values, s.values, valuesSize() * sizeof(double));
        treeStep = s.treeStep;
        testTreeStep();
        return *this;
    }
    catch (exception& e) {
        throw ModelException(e, "TreeSliceRates::operator=(const TreeSliceRates& s)");
    }
}

int TreeSliceRates::offset( int i, int j, int k ) const
{
#if DEBUG
    if (dim!=3) throw ModelException(
        "Calling offset(i,j,k) on slice \""+name+"\" of dim "+Format::toString(dim));
#endif
    int x3 = k - range.bot[2] - ( range.limits.bot3[ i ][ j ] + range.limits.top3[ i ][ j ] ) / 2;
    int x2 = j - range.bot[1] - ( range.limits.bot2[ i ]      + range.limits.top2[ i ]      ) / 2;
    int x1 = i - range.bot[0];

    return range.size[2] * ( range.size[1] * x1 + x2 ) + x3;
}

// calc the offet to add to "values" to access an element in a 2D slice
int TreeSliceRates::offset( int i, int j ) const
{
#if DEBUG
    if (dim!=2) throw ModelException(
        "Calling offset(i,j) on slice \""+name+"\" of dim "+Format::toString(dim));
#endif
    int x2 = j - range.bot[1] - ( range.limits.bot2[ i ]      + range.limits.top2[ i ]      ) / 2;
    int x1 = i - range.bot[0];

    return range.size[1] * x1 + x2;
}

// calc the offet to add to "values" to access an element in a 1D slice
int TreeSliceRates::offset( int i ) const
{
#if DEBUG
    if (dim!=1) throw ModelException(
        "Calling offset(i) on slice \""+name+"\" of dim "+Format::toString(dim));
#endif
    return      i - range.bot[0];
}

double TreeSliceRates::getCentre() const {
    testTreeStep();

    switch (dim) {
        case 3: return values[ offset( 0, 0, 0 ) ];
        case 2: return values[ offset( 0, 0 ) ];
        case 1: return values[ offset( 0 ) ];
        case 0: return values[0];
        case -1:
            throw ModelException("TreeSliceRates::getCentre", 
                                 "Tree slice (name = " + name +
                                 ") has not been initialized");
        default:
            throw ModelException("TreeSliceRates::getCentre", 
                                 "Invalid dimension (" + Format::toString(dim) + 
                                 ") for tree slice, name = " + name);
    }
}

static const char *errorSliceNotInitialized = "Slice operation on a non initialized slice: ";

// see usage of iterSeek(...) in calc() below
void TreeSliceRates::iterSeek(int i) const { // get ready to loop on dimension #0 (1st)
    iter = values;
    incInner=0;
    switch (dim) {
        case 1: 
            iter += offset(i);
            incInner=1; 
            break;
        case 0: 
            incInner=0;
            break;
        case -1: 
            throw ModelException(errorSliceNotInitialized + name);
        default: ASSERT(0);
    }
}

void TreeSliceRates::iterSeek(int i, int j) const { // get ready to loop on dimension #1 (2nd)
    iter = values;
    switch (dim) {
        case 2: 
            iter += offset(i,j); 
            incInner=1; 
            break;
        case 1: 
            iter += offset(i);
            incInner=0; 
            break;
        case 0: 
            incInner=0; 
            break;
        case -1: throw ModelException(errorSliceNotInitialized + name);
        default: ASSERT(0);
    }
}

void TreeSliceRates::iterSeek(int i, int j, int k) const { // get ready to loop on dimension #2 (3rd)
    iter = values;
    switch (dim) {
        case 3: 
            iter += offset(i,j,k); 
            incInner=1; 
            break;
        case 2: 
            iter += offset(i,j);
            incInner=0; 
            break;
        case 1: 
            iter += offset(i);
            incInner=0; 
            break;
        case 0: 
            incInner=0; 
            break;
        case -1: 
            throw ModelException(errorSliceNotInitialized + name);
        default: ASSERT(0);
    }
}

TreeSliceSP TreeSliceRates::calcSmoothStep() const {
    TreeSliceRatesSP r(new TreeSliceRates(range, curveToDEV, curveToDEVIdx));
    r->name = "calcSmoothStep("+name+")";
    r->allocDim(dim);
    double *rValues = r->values;

    switch (dim) {
        case 0: {
            rValues[0] = 0.;
            break;
        }
        case 1: {
            int iStart = range.limits.bot1;
            int iStop  = range.limits.top1;

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
            int iStart = range.limits.bot1;
            int iStop  = range.limits.top1;
            for( int i = iStart; i <= iStop; ++i ) {
                int jStart = range.limits.bot2[ i ];
                int jStop  = range.limits.top2[ i ];

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
                    && j >= range.limits.bot2[ i-1 ]
                    && j <= range.limits.top2[ i-1 ]) 
                    {
                        m = max(m, fabs( val(i-1,j) - v ));
                    }
                    if(i != iStop
                    && j >= range.limits.bot2[ i+1 ]
                    && j <= range.limits.top2[ i+1 ]) 
                    {
                        m = max(m, fabs( val(i+1,j) - v ));
                    }

                    rValues[offset(i,j)] = m;
                }
            }
            break;
        }
        case 3: {
            int iStart = range.limits.bot1;
            int iStop  = range.limits.top1;
            for( int i = iStart; i <= iStop; ++i ) {
                int jStart = range.limits.bot2[ i ];
                int jStop  = range.limits.top2[ i ];

                for( int j = jStart; j <= jStop; ++j ) {
                    int kStart = range.limits.bot3[ i ][ j ];
                    int kStop  = range.limits.top3[ i ][ j ];

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
                        && k >= range.limits.bot3[ i ][ j-1 ]
                        && k <= range.limits.top3[ i ][ j-1 ])
                        {
                            m = max(m, fabs( val(i,j-1,k) - v ));
                        }
                        if(j != jStop
                        && k >= range.limits.bot3[ i ][ j+1 ]
                        && k <= range.limits.top3[ i ][ j+1 ])
                        {
                            m = max(m, fabs( val(i,j+1,k) - v ));
                        }
                        if(i != iStart
                        && j >= range.limits.bot2[ i-1 ]
                        && j <= range.limits.top2[ i-1 ]
                        && k >= range.limits.bot3[ i-1 ][ j ]
                        && k <= range.limits.top3[ i-1 ][ j ])
                        {
                            m = max(m, fabs( val(i-1,j,k) - v ));
                        }
                        if(i != iStop
                        && j >= range.limits.bot2[ i+1 ]
                        && j <= range.limits.top2[ i+1 ]
                        && k >= range.limits.bot3[ i+1 ][ j ]
                        && k <= range.limits.top3[ i+1 ][ j ])
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
    r->treeStep = range.treeStep;
    return r;
}

void TreeSliceRates::prepareSlicesForEval(int nbSlices, TreeSliceRates **s, TreeSliceRates &newSlice) {
    int newDim = 0; // dim of the resulting slice

    for (int i=0; i<nbSlices; ++i) {
        // check that slices are valid
        if (!s[i]) {
            throw ModelException(
                "Expression contains a non TreeSliceRates slice type "
                "in position "+Format::toString(i));
        }
        if (s[i]->dim<0) {
            throw ModelException("Uninitialized slice (position "+Format::toString(i)+")");
        }
        // checks that treeSteps are consistent
        s[i]->testTreeStep(); 
        // calc max of dims
        newDim = Maths::max(newDim, s[i]->dim); 
    }
    // initialize newSlice
    if (newDim != dim) {
        // with a new buffer
        newSlice.allocDim(newDim);
        newSlice.copyTreeSlice(*this);
    } else {
        // use the same buffer than *this
        newSlice.takeFrom(*this);
    }
    s[nbSlices]=&newSlice;
}

void TreeSliceRates::getCalcRange( int & bot, int & top ) const
{
    if( dim == 1 )
    {
        bot = range.limits.bot1;
        top = range.limits.top1;
    }
    else
    {
        bot = 0;
        top = valuesSize() - 1;
    }
}

double * TreeSliceRates::getValues() const
{
    return dim == 1 ? values - range.limits.bot1 : values;
}



void TreeSliceRates::printDetails(char *s) const {
    static char buf[30];
    TreeSlice::printDetails(s);
    strcat(s, " dim=");
    sprintf(buf,"%d", dim);
    strcat(s, buf);
    if (dim>=0) {
        strcat(s, " centre=");
        sprintf(buf, "%lf", getCentre());
        strcat(s, buf);
    }
}

DRLIB_END_NAMESPACE
