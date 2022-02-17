/*!
 *
 * Copyright (c) IXIS CM January 2005 Paris
 *
 *	\file pricingcontext.h
 *
 *  \brief pricing context save information computed during a price
 *
 *	\author  JM Prie
 *	\version 1.0
 *	\date July 2005
 */


#ifndef _INGPINFRA_PRICINGCONTEXT_H
#define _INGPINFRA_PRICINGCONTEXT_H

#include "gpbase/removeidentifiedwarning.h"
#include "gpbase/port.h"
#include "gpbase/typedef.h"

#define UNIMPLEMENTED_CONTEXT_FUNCTION  { ARM_THROW( ERR_INVALID_ARGUMENT, "unimplemented pricing function in ARM_PricingContext"); }

CC_BEGIN_NAMESPACE( ARM )

class ARM_PricingCallContext;
class ARM_PricingDigitalContext;

//////////////////////////////////////////////
/// \struct ARM_PricingContext
/// \brief
/// ARM_PricingContext defines a general
/// container for intermediate pricing datas
//////////////////////////////////////////////
struct ARM_PricingContext
{
    virtual ARM_PricingCallContext* ToCallContext() UNIMPLEMENTED_CONTEXT_FUNCTION;
    virtual ARM_PricingDigitalContext* ToDigitalContext() UNIMPLEMENTED_CONTEXT_FUNCTION;
};

//////////////////////////////////////////////
/// \class ARM_PricingCallContext
/// \brief
/// ARM_PricingCallContext defines a container
/// for intermediate pricing datas for a call
/// computation
//////////////////////////////////////////////
class ARM_PricingCallContext : public ARM_PricingContext
{
private:
    double itsVol;
    double itsShift;
    ARM_GP_VectorPtr itsForward;

public: 
	/// Downcasting
    ARM_PricingCallContext()  : itsVol(0.0), itsShift(0.0), itsForward(NULL) {};
	ARM_PricingCallContext( const ARM_PricingCallContext& rhs )
        : itsVol(rhs.itsVol), itsShift(rhs.itsShift),itsForward(rhs.itsForward) {}
	ARM_PricingCallContext& operator=( const ARM_PricingCallContext& rhs )
    {
        if(&rhs != this)
        {
            this->~ARM_PricingCallContext();
            new (this) ARM_PricingCallContext (rhs);
        }
        return *this;
    }

    virtual ~ARM_PricingCallContext() {}

    virtual ARM_PricingCallContext* ToCallContext() { return this; }

    /// Accessors
    double GetVol() const { return itsVol; }
    void SetVol(double vol) { itsVol=vol; }
    double GetShift() const { return itsShift; }
    void SetShift(double shift) { itsShift=shift; }
    const ARM_GP_VectorPtr& GetForward() const { return itsForward; }
    void SetForward(const ARM_GP_VectorPtr& forward) { itsForward=forward; }
};


//////////////////////////////////////////////
/// \class ARM_PricingDigitalContext
/// \brief
/// ARM_PricingDigitalContext defines a container
/// for intermediate pricing datas for a call
/// computation
//////////////////////////////////////////////
class ARM_PricingDigitalContext : public ARM_PricingContext
{
private:
    std::vector<double> itsVols;
    std::vector<double> itsShifts;
    ARM_VectorPtrVector itsForwards;

public: 
	/// Downcasting
    ARM_PricingDigitalContext() : itsVols(0), itsShifts(0),itsForwards(0) {};
	ARM_PricingDigitalContext( const ARM_PricingDigitalContext& rhs )
        : itsVols(rhs.itsVols),itsShifts(rhs.itsShifts), itsForwards(rhs.itsForwards) {}
	ARM_PricingDigitalContext& operator=( const ARM_PricingDigitalContext& rhs )
    {
        if(&rhs != this)
        {
            this->~ARM_PricingDigitalContext();
            new (this) ARM_PricingDigitalContext (rhs);
        }
        return *this;
    }

    virtual ~ARM_PricingDigitalContext() {}

    virtual ARM_PricingDigitalContext* ToDigitalContext() { return this; }

    /// Accessors
    const std::vector<double>& GetVols() const { return itsVols; }
    void SetVols(std::vector<double> vols) { itsVols=vols; }
    const std::vector<double>& GetShifts() const { return itsShifts; }
    void SetShifts(std::vector<double> shifts) { itsShifts=shifts; }
    const ARM_VectorPtrVector& GetForward() const { return itsForwards; }
    void SetForward(const ARM_VectorPtrVector& forwards) { itsForwards=forwards; }

    void InsertDatas(double vol, double shift, const ARM_GP_VectorPtr& forward)
    { itsVols.push_back(vol); itsShifts.push_back(shift); itsForwards.push_back(forward); }
};



CC_END_NAMESPACE()

#endif

/*----------------------------------------------------------------------------*/
/*---- End of file ----*/

