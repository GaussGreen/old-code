//----------------------------------------------------------------------------
//
//   Group       : EDG Derivatives Research
//
//   Filename    : HaveNonEqSpread.hpp
//
//   Description : Interface class for getNonEqSpread() method for equity dependent spread
//
//   Author      : Qing Hou
//
//   Date        : 23 Dec 2003
//
//
//----------------------------------------------------------------------------

#ifndef HAVE_NON_EQ_SPREAD_HPP
#define HAVE_NON_EQ_SPREAD_HPP

DRLIB_BEGIN_NAMESPACE

/** abstract interface class  */
class TREE_DLL IHaveNonEqSpread {
public:
	// return 0 if there is no non-eq-related-spread.
	virtual const CleanSpreadCurve *getNonEqSpreads() const=0;
};

// typedef smartPtr<IHaveNonEqSpread> IHaveNonEqSpreadSP;

DRLIB_END_NAMESPACE

#endif
