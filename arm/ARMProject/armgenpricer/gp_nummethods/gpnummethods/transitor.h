/*!
 *
 * Copyright (c) CDC IXIS CM July 2003 Paris
 *
 *	\file transitor.h
 *
 *  \brief
 *
 *	\author  JM Prie, E. Benhamou
 *	\version 1.0
 *	\date December 2004
 */


#ifndef _INGPNUMMETHODS_TRANSITOR_H
#define _INGPNUMMETHODS_TRANSITOR_H

#include "gpbase/port.h"
#include "gpbase/env.h"
#include "gpbase/gpvector.h"
#include "typedef.h"

CC_BEGIN_NAMESPACE( ARM )

/// forward declaration
struct ARM_SliceBase;

class ARM_TransitionStates
{
private:
    ARM_IntVector itsStates;
    std::vector<double> itsProbas;
    size_t itsUsedSize;
public:
    ARM_TransitionStates(size_t maxSize=0)
    : itsStates(maxSize), itsProbas(maxSize),itsUsedSize(0) {}
    virtual ~ARM_TransitionStates() {};

    /// Accessors
    const ARM_IntVector& GetStates() const { return itsStates; }
    const std::vector<double>& GetProbas() const { return itsProbas; }
    size_t GetUsedSize() const { return itsUsedSize; }
    void SetUsedSize(size_t size) { itsUsedSize = size; }

    void SetTransition(size_t transIdx, size_t stateIdx, double proba)
    { itsStates[transIdx] = stateIdx; itsProbas[transIdx] = proba; }
};


///---------------------------
/// ARM_TransitionBase
///---------------------------
struct ARM_TransitorBase : public ARM_RootObject
{
    virtual void ComputeTransitions(const ARM_SliceBase* tmpSlice, size_t stateIdx, const ARM_SliceBase* nextSlice, ARM_TransitionStates& transStates) const = 0;
};

///---------------------------
/// ARM_Transition2D
///---------------------------
struct ARM_Transitor2D : public ARM_TransitorBase
{
    virtual ARM_Object* Clone() const { return new ARM_Transitor2D(*this); }
    virtual void ComputeTransitions(const ARM_SliceBase* tmpSlice, size_t stateIdx, const ARM_SliceBase* nextSlice, ARM_TransitionStates& transStates) const;
    virtual string toString(const string& indent="", const string& nextIndent="") const { return "ARM_Transitor2D";}
};

///---------------------------
/// ARM_Transition3D
///---------------------------
struct ARM_Transitor3D : public ARM_TransitorBase
{
    virtual ARM_Object* Clone() const { return new ARM_Transitor3D(*this); }
    virtual void ComputeTransitions(const ARM_SliceBase* tmpSlice, size_t stateIdx, const ARM_SliceBase* nextSlice, ARM_TransitionStates& transStates) const;
    virtual string toString(const string& indent="", const string& nextIndent="") const { return "ARM_Transitor3D";}
};


///---------------------------
/// ARM_TransitionND
///---------------------------
struct ARM_TransitorND : public ARM_TransitorBase
{
private:
	size_t itsDim;
public:
	ARM_TransitorND( size_t dim ): ARM_TransitorBase(), itsDim(dim) {}

    virtual ARM_Object* Clone() const { return new ARM_TransitorND(*this); }

    virtual void ComputeTransitions(const ARM_SliceBase* tmpSlice, size_t stateIdx, const ARM_SliceBase* nextSlice, ARM_TransitionStates& transStates) const;

    virtual string toString(const string& indent="", const string& nextIndent="") const { return "ARM_TransitorND";}
};

struct ARM_TransitorFactory
{
	static ARM_TransitorBase* CreateTransitor( size_t dim )
	{
		switch( dim )
		{
		case 2 :
			return new ARM_Transitor2D;
		case 3:
			return new ARM_Transitor3D;
		default:
			return new ARM_TransitorND(dim);
		}
	}
};


CC_END_NAMESPACE()

#endif

/*----------------------------------------------------------------------------*/
/*---- End of file ----*/

