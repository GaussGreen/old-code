//----------------------------------------------------------------------------
//
//   Group       : Quantitative Research
//
//   Filename    : IRGrid.cpp
//
//   Description : Class to be able to aggregate IRGridResultArrays.
//                 Contains the actual result data array plus a method to 
//                 aggregate other objects of the same type
//
//   Author      : Jose Hilera
//
//   Date        : July 2005
//
//----------------------------------------------------------------------------

#include "edginc/config.hpp"
#include "edginc/IRGrid.hpp"

DRLIB_BEGIN_NAMESPACE

IRGrid::IRGrid(int size): CObject(TYPE) {
    gridResult = IRGridResultArraySP(new IRGridResultArray(size));
}

IRGrid::~IRGrid() {}


void IRGrid::setPoint(int index, IRGridResult irg) {
    (*gridResult)[index] = irg;
}


/** scale by factor x (implementation of CombinableResult) */
void IRGrid::scale(double scaleFactor) {
    for (int i=0; i < gridResult->size(); ++i) {
        (*gridResult)[i].scale(scaleFactor);
    }
}


/** add another IRGrid object to this result
 * (implementation of CombinableResult) */
void IRGrid::add(const CombinableResult& x, 
                                      double scaleFactor) 
{
    bool found;
    int i, j;

    // gcc bug: force to IObject before dynamic cast
    const IRGrid& newCIR = dynamic_cast<const IRGrid&>(static_cast<const IObject&>(x));

    for (i=0; i < newCIR.gridResult->size(); ++i) 
    {
        found = false;

        // if this grid result is already present in the current array,
        // add its content (after scaling) - otherwise, concatenate this
        // new grid result to the end of the current array
        for (j=0; j < gridResult->size() && !found; ++j) {
            if ((*gridResult)[j].getGridPoint()->getExpiry()->equals(
                      (*newCIR.gridResult)[i].getGridPoint()->getExpiry().get()) &&
                (*gridResult)[j].getGridPoint()->getTenor()->equals(
                      (*newCIR.gridResult)[i].getGridPoint()->getTenor().get()))
            {
                // Point already present -> add its contents (scaled)
                (*gridResult)[j].addToResult(scaleFactor * (*newCIR.gridResult)[i].getResult());
                found = true;
            }
        }

        if (!found) {
            // Add the new point (scaled) to the current array.
            // Need to create a new gridResult point because we cannot scale
            // the original point.
            gridResult->push_back(
                  IRGridResult((*newCIR.gridResult)[i].getGridPoint(),
                               scaleFactor * (*newCIR.gridResult)[i].getResult()));
        }
    }
}


/** write object out in 'output' format - ie suitable for comparing
 * regression files with */
void IRGrid::outputWrite(const string& linePrefix,
                                                    const string& prefix, 
                                                    ostream& stream) const 
{
    gridResult->outputWrite(linePrefix, prefix, stream);
}


/** Invoked when class is 'loaded' */
void IRGrid::load(CClassSP& clazz){
    clazz->setPublic(); // make visible to EAS/spreadsheet
    REGISTER(IRGrid, clazz);
    SUPERCLASS(CObject);
    IMPLEMENTS(CombinableResult);
    EMPTY_SHELL_METHOD(defaultIRGrid);
    FIELD(gridResult, "The grid result");
}


IObject* IRGrid::defaultIRGrid(){
    return new IRGrid();
}

CClassConstSP const IRGrid::TYPE = 
CClass::registerClassLoadMethod("IRGrid", typeid(IRGrid), IRGrid::load);

DRLIB_END_NAMESPACE
