//----------------------------------------------------------------------------
//
//   Group       : EDG Derivatives Research
//
//   Filename    : SPIBondFloor.hpp
//
//   Description : SPI bond floor handling
//
//   Date        : Dec 2005
//
//
//----------------------------------------------------------------------------
#ifndef EDR_SPI_BF_HPP
#define EDR_SPI_BF_HPP

#include "edginc/config.hpp"
#include "edginc/DateTime.hpp"
#include "edginc/SVGenDiscFactor.hpp"

DRLIB_BEGIN_NAMESPACE

// forward declaration
class SPIRunTime;

// Increased sophistication of bond floor means it's time for a class...(or two)
class PRODUCTS_DLL ISPIBondFloor {
public:
    // BL (locked-in level varies during simulation
    // The rest is (currently) static
    virtual const double getLevel(double BL,
                                    int iStep) const = 0;

    virtual const double getLevelToday() const = 0;
    
    virtual void refresh(int iFirstFutureStep) = 0;

    virtual ~ISPIBondFloor();
};
typedef refCountPtr<ISPIBondFloor> ISPIBondFloorSP;

// a factory
class PRODUCTS_DLL SPIBondFloor {
public:
    static ISPIBondFloor* make(SPIRunTime*    dynBaskRT,
                               double         Basis,
                               const DateTime today);
};

DRLIB_END_NAMESPACE
#endif
