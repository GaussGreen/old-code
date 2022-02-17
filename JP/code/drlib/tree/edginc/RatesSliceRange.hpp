#ifndef RATES_SLICE_RANGE_HPP
#define RATES_SLICE_RANGE_HPP

#include "edginc/DECLARE.hpp"

DRLIB_BEGIN_NAMESPACE

#define SLICE_RATES_MAX_DIMS 4

/*************************** types ***************************/
class TREE_DLL RatesSliceRange {
public:
    struct Limits // dynamic size for each dimension depending on outter dimensions
    {
        int     bot1,     top1;
        int   * bot2,   * top2;
        int  ** bot3,  ** top3;
        int *** bot4, *** top4;
    };

    struct Offsets
    {
        int    offset1;
        int   *offset2;
        int  **offset3;
        int ***offset4;
    };

    int bot     [ SLICE_RATES_MAX_DIMS ]; // for each dimension, starting index from outside
    int top     [ SLICE_RATES_MAX_DIMS ]; // for each dimension, ending index from outside
    int size    [ SLICE_RATES_MAX_DIMS ]; // for each dimension, top - bot + 1;
    int fullsize[ SLICE_RATES_MAX_DIMS ]; // for each dimension, full slice size;
    int nDim;                             // number of dimensions of the tree
    Limits  limits;                       // current range limits
    Offsets offsets;                      // current node offsets
    int     treeStep;                     // current time step tree is set at

private:

    virtual void init(int bot1,     int top1,              // 1st dimension range
                      int bot2,     int top2,              // 2st dimension range
                      int bot3,     int top3,              // 3st dimension range
                      int bot4 = 0, int top4 = 0)          // 4th dimension range
    {
        throw ModelException("RatesSliceRange::init()", "incorrect slice init");
    }

    virtual void init(int fullsize1,
                      int fullsize2,
                      int fullsize3,
                      int fullsize4)
    {
        throw ModelException("RatesSliceRange::init()", "incorrect slice init");
    }
};

typedef refCountPtr<RatesSliceRange> RatesSliceRangeSP;

DRLIB_END_NAMESPACE

#endif
