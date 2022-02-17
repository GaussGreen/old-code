#include "edginc/config.hpp"
#include "edginc/TreeSlice.hpp"

DRLIB_BEGIN_NAMESPACE

/******************************** TreeSliceRatesCompact::Range ********************************/

void TreeSliceRatesCompact::Range::init(int fullsize1,
                                        int fullsize2,
                                        int fullsize3,
                                        int fullsize4)
{
    fullsize[0] = fullsize1;
    fullsize[1] = fullsize2;
    fullsize[2] = fullsize3;
    fullsize[3] = fullsize4;

    if( fullsize[0] < 0 )
        throw ModelException("size 0 negative during construction");
    if( fullsize[1] < 0 )
        throw ModelException("size 1 negative during construction");
    if( fullsize[2] < 0 )
        throw ModelException("size 2 negative during construction");
    if( fullsize[3] < 0 ) 
        throw ModelException("size 3 negative during construction");

    int d = SLICE_RATES_MAX_DIMS;
    treeStep = -1;

    while( --d && fullsize[ d ] <= 1 ) 
    {}
    nDim = d + 1;
}

/******************************** TreeSliceRatesCompact ********************************/

TreeSliceRatesCompact::TreeSliceRatesCompact(const RatesSliceRange &range, const string& curveToDEV, int curveToDEVIdx) :
    curveToDEV(curveToDEV), curveToDEVIdx(curveToDEVIdx),
    dim(-1), values(0), deleteWhenDestroyed(true), range(dynamic_cast<const Range&>(range))
{}

TreeSliceRatesCompact::~TreeSliceRatesCompact() {allocDim(-1);}

TreeSliceSP TreeSliceRatesCompact::clone( bool copyValues) const {
    TreeSliceRatesCompactSP slice( new TreeSliceRatesCompact(range, curveToDEV, curveToDEVIdx) );
    if (copyValues) {
        slice->allocDim(dim);
        memcpy(slice->values, values, valuesSize()*sizeof(double));
    }
    slice->treeStep = treeStep;
    return slice;
}

void TreeSliceRatesCompact::resize(int n) {
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

void TreeSliceRatesCompact::takeFrom(TreeSliceRatesCompact &s) {
    if (values != s.values) 
        resize(-1);
    values = s.values;
    deleteWhenDestroyed = s.deleteWhenDestroyed;
    s.deleteWhenDestroyed = false;

    dim = s.dim;
    copyTreeSlice(s);
    theValues = values; // for debug purposes
}

void TreeSliceRatesCompact::expand(int newDim) {
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
        TreeSliceRatesCompact newSlice(range, curveToDEV, curveToDEVIdx);
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
        throw ModelException(e, "TreeSliceRatesCompact::expand, name=\""+name+"\"");
    }
}


int TreeSliceRatesCompact::valuesSize() const {
    return range.fullsize[dim];
}

void TreeSliceRatesCompact::allocDim(int newDim) {
    dim = newDim;
    resize(valuesSize());
}

TreeSlice& TreeSliceRatesCompact::operator=(double constant) {
    if (dim!=0)
        allocDim(0);
    values[0]=constant;
    treeStep = range.treeStep;  // set slice to current tree time step
    return *this;
}

TreeSlice& TreeSliceRatesCompact::operator=(const TreeSlice& s) {
    const TreeSliceRatesCompact *r = dynamic_cast<const TreeSliceRatesCompact*>(&s);
    if (!r)
        throw ModelException("TreeSliceRatesCompact::operator=(const TreeSlice& s)",
        "Cannot assign slices of different types");
    return operator=(*r);
}

void TreeSliceRatesCompact::testTreeStep() const {
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

TreeSliceRatesCompact& TreeSliceRatesCompact::operator=(const TreeSliceRatesCompact& s) {
    try {
        if (dim!=s.dim)
            allocDim(s.dim);

        ::memcpy(values, s.values, valuesSize() * sizeof(double));
        treeStep = s.treeStep;
        testTreeStep();
        return *this;
    }
    catch (exception& e) {
        throw ModelException(e, "TreeSliceRatesCompact::operator=(const TreeSliceRatesCompact& s)");
    }
}

int TreeSliceRatesCompact::offset( int i, int j, int k , int L) const
{
#if DEBUG
    if (dim!=4) throw ModelException(
        "Calling offset(i,j,k,L) on slice \""+name+"\" of dim "+Format::toString(dim));
#endif

    return range.offsets.offset4[i][j][k] + L;
}

int TreeSliceRatesCompact::offset( int i, int j, int k ) const
{
#if DEBUG
    if (dim!=3) throw ModelException(
        "Calling offset(i,j,k) on slice \""+name+"\" of dim "+Format::toString(dim));
#endif

    return range.offsets.offset3[i][j] + k;
}

// calc the offet to add to "values" to access an element in a 2D slice
int TreeSliceRatesCompact::offset( int i, int j ) const
{
#if DEBUG
    if (dim!=2) throw ModelException(
        "Calling offset(i,j) on slice \""+name+"\" of dim "+Format::toString(dim));
#endif

    return range.offsets.offset2[i] + j;
}

// calc the offet to add to "values" to access an element in a 1D slice
int TreeSliceRatesCompact::offset( int i ) const
{
#if DEBUG
    if (dim!=1) throw ModelException(
        "Calling offset(i) on slice \""+name+"\" of dim "+Format::toString(dim));
#endif

    return range.offsets.offset1 + i;
}

double TreeSliceRatesCompact::getCentre() const {
    testTreeStep();

    switch (dim) {
        case 4: return values[ offset( 0, 0, 0, 0 ) ];
        case 3: return values[ offset( 0, 0, 0 ) ];
        case 2: return values[ offset( 0, 0 ) ];
        case 1: return values[ offset( 0 ) ];
        case 0: return values[0];
        case -1:
            throw ModelException("TreeSliceRatesCompact::getCentre", 
                                 "Tree slice (name = " + name +
                                 ") has not been initialized");
        default:
            throw ModelException("TreeSliceRatesCompact::getCentre", 
                                 "Invalid dimension (" + Format::toString(dim) + 
                                 ") for tree slice, name = " + name);
    }
}

static const char *errorSliceNotInitialized = "Slice operation on a non initialized slice: ";

// see usage of iterSeek(...) in calc() below
void TreeSliceRatesCompact::iterSeek(int i) const { // get ready to loop on dimension #0 (1st)
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

void TreeSliceRatesCompact::iterSeek(int i, int j) const { // get ready to loop on dimension #1 (2nd)
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

void TreeSliceRatesCompact::iterSeek(int i, int j, int k) const { // get ready to loop on dimension #2 (3rd)
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

void TreeSliceRatesCompact::iterSeek(int i, int j, int k, int L ) const { // get ready to loop on dimension #2 (3rd)
    iter = values;
    switch (dim) {
        case 4: 
            iter += offset(i,j,k,L); 
            incInner=1; 
            break;
        case 3: 
            iter += offset(i,j,k); 
            incInner=0; 
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

TreeSliceSP TreeSliceRatesCompact::calcSmoothStep() const {
    TreeSliceRatesCompactSP r(new TreeSliceRatesCompact(range, curveToDEV, curveToDEVIdx));
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
        case 4:
        {
            int iStart = range.limits.bot1;
            int iStop  = range.limits.top1;

            for( int i = iStart; i <= iStop; ++i )
            {
                int jStart = range.limits.bot2[ i ];
                int jStop  = range.limits.top2[ i ];

                for( int j = jStart; j <= jStop; ++j )
                {
                    int kStart = range.limits.bot3[ i ][ j ];
                    int kStop  = range.limits.top3[ i ][ j ];

                    for( int k = kStart; k <= kStop; ++k ) 
                    {
                        int LStart = range.limits.bot4[ i ][ j ][ k ];
                        int LStop  = range.limits.top4[ i ][ j ][ k ];

                        for( int L = LStart; L <= LStop; ++L ) 
                        {
                            double m=0, v=val(i,j,k,L);

                            if (L != LStart)
                            {
                                m = max(m, fabs( val(i,j,k,L-1) - v));
                            }

                            if (L != LStop)
                            {
                                m = max(m, fabs( val(i,j,k,L+1) - v));
                            }

                            if ((k != kStart)                             &&
                                (L >= range.limits.bot4[ i ][ j ][ k-1 ]) &&
                                (L <= range.limits.top4[ i ][ j ][ k-1 ]))
                            {
                                m = max(m, fabs( val(i,j,k-1,L) - v));
                            }

                            if ((k != kStop)                              &&
                                (L >= range.limits.bot4[ i ][ j ][ k+1 ]) &&
                                (L <= range.limits.top4[ i ][ j ][ k+1 ]))
                            {
                                m = max(m, fabs( val(i,j,k+1,L) - v));
                            }

                            if ((j != jStart)                             &&
                                (k >= range.limits.bot3[ i ][ j-1 ])      &&
                                (k <= range.limits.top3[ i ][ j-1 ])      &&
                                (L >= range.limits.bot4[ i ][ j-1 ][ k ]) &&
                                (L <= range.limits.top4[ i ][ j-1 ][ k ]))
                            {
                                m = max(m, fabs( val(i,j-1,k,L) - v));
                            }

                            if ((j != jStop)                              &&
                                (k >= range.limits.bot3[ i ][ j+1 ])      &&
                                (k <= range.limits.top3[ i ][ j+1 ])      &&
                                (L >= range.limits.bot4[ i ][ j+1 ][ k ]) &&
                                (L <= range.limits.top4[ i ][ j+1 ][ k ]))
                            {
                                m = max(m, fabs( val(i,j+1,k,L) - v));
                            }

                            if ((i != iStart)                             &&
                                (j >= range.limits.bot2[ i-1 ])           &&
                                (j <= range.limits.top2[ i-1 ])           &&
                                (k >= range.limits.bot3[ i-1 ][ j ])      &&
                                (k <= range.limits.top3[ i-1 ][ j ])      &&
                                (L >= range.limits.bot4[ i-1 ][ j ][ k ]) &&
                                (L <= range.limits.top4[ i-1 ][ j ][ k ]))
                            {
                                m = max(m, fabs( val(i-1,j,k,L) - v));
                            }

                            if ((i != iStop)                              &&
                                (j >= range.limits.bot2[ i+1 ])           &&
                                (j <= range.limits.top2[ i+1 ])           &&
                                (k >= range.limits.bot3[ i+1 ][ j ])      &&
                                (k <= range.limits.top3[ i+1 ][ j ])      &&
                                (L >= range.limits.bot4[ i+1 ][ j ][ k ]) &&
                                (L <= range.limits.top4[ i+1 ][ j ][ k ]))
                            {
                                m = max(m, fabs( val(i+1,j,k,L) - v));
                            }
                        }
                    }
                }
            }

            break;
        }                    
        default:
            throw ModelException("TreeSliceRatesCompact::calcSmoothStep", "unsupported slice dimension");
    }
    r->treeStep = range.treeStep;
    return r;
}

void TreeSliceRatesCompact::prepareSlicesForEval(int nbSlices, TreeSliceRatesCompact **s, TreeSliceRatesCompact &newSlice) {
    int newDim = 0; // dim of the resulting slice

    for (int i=0; i<nbSlices; ++i) {
        // check that slices are valid
        if (!s[i]) {
            throw ModelException(
                "Expression contains a non TreeSliceRatesCompact slice type "
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

void TreeSliceRatesCompact::getCalcRange( int & bot, int & top ) const
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

double * TreeSliceRatesCompact::getValues() const
{
    return dim == 1 ? values - range.limits.bot1 : values;
}



void TreeSliceRatesCompact::printDetails(char *s) const {
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
