//----------------------------------------------------------------------------
//
//   Group       : Quantitative Research
//
//   Filename    : BetaSkewGrid.cpp
//
//   Description : Class to be able to aggregate BetaSkewGridResultArrays.
//                 Contains the actual result data array plus a method to 
//                 aggregate other objects of the same type
//
//   Author      : Jose Hilera
//
//   Date        : July 2005
//
//----------------------------------------------------------------------------

#include "edginc/config.hpp"
#include "edginc/BetaSkewGrid.hpp"
#include "edginc/Maths.hpp"


DRLIB_BEGIN_NAMESPACE

BetaSkewGrid::BetaSkewGrid(int size) : CObject(TYPE){
    gridResult = BetaSkewGridResultArraySP(new BetaSkewGridResultArray(size));
}

BetaSkewGrid::~BetaSkewGrid() {}


void BetaSkewGrid::setPoint(int index, BetaSkewGridResultSP bs) {
    (*gridResult)[index] = bs;
}


/** scale by factor x (implementation of CombinableResult) */
void BetaSkewGrid::scale(double scaleFactor) {
    for (int i=0; i < gridResult->size(); ++i) {
        (*gridResult)[i]->scale(scaleFactor);
    }
}


/** add another BetaSkewGrid object to this result
 * (implementation of CombinableResult) */
void BetaSkewGrid::add(const CombinableResult& x, 
                                            double scaleFactor) 
{
    bool found;

    // gcc bug: force to IObject before dynamic cast
    const BetaSkewGrid& newCBS = dynamic_cast<const BetaSkewGrid&>(static_cast<const IObject&>(x));

    for (int i=0; i < newCBS.gridResult->size(); ++i) 
    {
        found = false;

        // if this grid result is already present in the current array,
        // add its content (after scaling) - otherwise, concatenate this
        // new grid result to the end of the current array
        for (int j=0; j < gridResult->size() && !found; ++j) {
            if ((*gridResult)[j]->getGridPoint().getMaturity().equals(
                              (*newCBS.gridResult)[i]->getGridPoint().getMaturity(), false) &&
                Maths::equals((*gridResult)[j]->getGridPoint().getStrike(),
                              (*newCBS.gridResult)[i]->getGridPoint().getStrike()))
            {
                // Point already present -> add its contents (scaled)
                (*gridResult)[j]->addToResult(scaleFactor * (*newCBS.gridResult)[i]->getResult());
                found = true;
            }
        }

        if (!found) {
            // Add the new point (scaled) to the current array.
            // Need to create a new gridResult point because we cannot scale
            // the original point.
            gridResult->push_back(
                BetaSkewGridResultSP (
                  new BetaSkewGridResult((*newCBS.gridResult)[i]->getGridPoint(),
                                         scaleFactor * (*newCBS.gridResult)[i]->getResult())));
        }
    }
}


/** write object out in 'output' format - ie suitable for comparing
 * regression files with */
void BetaSkewGrid::outputWrite(const string& linePrefix,
                                                    const string& prefix, 
                                                    ostream& stream) const 
{
    gridResult->outputWrite(linePrefix, prefix, stream);
}


/** Invoked when class is 'loaded' */
void BetaSkewGrid::load(CClassSP& clazz){
    clazz->setPublic(); // make visible to EAS/spreadsheet
    REGISTER(BetaSkewGrid, clazz);
    SUPERCLASS(CObject);
    IMPLEMENTS(CombinableResult);
    EMPTY_SHELL_METHOD(defaultBetaSkewGrid);
    FIELD(gridResult, "The grid result");
}


IObject* BetaSkewGrid::defaultBetaSkewGrid(){
    return new BetaSkewGrid();
}

CClassConstSP const BetaSkewGrid::TYPE = 
    CClass::registerClassLoadMethod("BetaSkewGrid", 
                                    typeid(BetaSkewGrid), BetaSkewGrid::load);

DRLIB_END_NAMESPACE
