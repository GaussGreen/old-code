//----------------------------------------------------------------------------
//
//   Group       : EDG Derivatives Research
//
//   Filename    : LocalVolGrid.hpp
//
//   Description : local vol grid support interpolations
//
//   Author      : Ning Shen
//
//   Date        : 6 Jan 2003
//
//
//----------------------------------------------------------------------------

#ifndef LOCALVOLGRID_HPP
#define LOCALVOLGRID_HPP

#include "edginc/AtomicArray.hpp"
#include "edginc/VolProcessedDVF.hpp"
#include "edginc/ILocalVolGrid.hpp"

DRLIB_BEGIN_NAMESPACE

/***************************************** 
LocalVolGrid class
a grid of local vols that support linear or spline interpolations
*****************************************/
FORWARD_DECLARE_REF_COUNT(LocalVolGrid);

class MCARLO_DLL LocalVolGrid: virtual public ILocalVolGrid {
public:
    /** Generates a grid of local vols for interpolation at each simulation date */
    LocalVolGrid(const DoubleArray& fwds, 
                 CVolProcessedDVFConstSP vol, 
                 const DateTimeArray& futurePathDates,
				 const DateTimeArray& cumDates,
                 int numVolGrid, 
                 double stdevGridRange);

    LocalVolGrid(){};
    
    /** Static method wrapper around constructor */
    static ILocalVolGridSP createLocalVolGrid(const DoubleArray& fwds, 
                                              CVolProcessedDVFConstSP vol, 
                                              const DateTimeArray& futurePathDates,
                                              const DateTimeArray& cumDates,
                                              int numVolGrid, 
                                              double stdevGridRange);

    /** calculate local volatility on a lattice of strikes */
    virtual void CalcLocVol(const CSliceDouble& strikes, int iDate, CSliceDouble& locVol);
    
    /** Linear interpolation of local vol from grid */ 
    virtual double interpLocVar(int step, double logSpot) const;

    virtual double maxDrift() const;

    /** Rolls LocalVolGrid forward */
    virtual void rollTime(const DateTime& valueDate);

private:
    DoubleArrayArray    spotGrid;       //!< Spot grid at which LV is computed
    DoubleArrayArray    lvGrid;         //!< The LocalVol grid
    DoubleArrayArray    lvY2;           //!< Used in spline computations

	DoubleArray         gridStart;      //!< Log of spotGrid[0][iStep]
	DoubleArray         gridSize;       //!< Increment in log spotGrid space

    DateTimeArray       gridDates;      //!< Dates to which the grid refers
    int                 dateOffset;     //!< Date offset for interpolation - used in Theta tweaks

    double              _maxDrift;      //!< Used for QuickGreeks
};


DRLIB_END_NAMESPACE

#endif
