//----------------------------------------------------------------------------
//
//   Group       : Energy Exotics Derivatives Research
//
//   Filename    : ForwardParallel.hpp
//
//   Description : Sensitivity to a parallel shift in the forward curve
//
//   Author      : Simon Creeger
//
//   Date        : 12 September 2006
//
//
//----------------------------------------------------------------------------

#ifndef FORWARDPARALLEL_HPP
#define FORWARDPARALLEL_HPP
#include "edginc/ScalarShift.hpp"
#include "edginc/smartPtr.hpp"
#include "edginc/Additive.hpp"

DRLIB_BEGIN_NAMESPACE

/** Sens Control for ForwardParallel - a scalar shift where the derivative is
calculated via a one sided tweak operation */
class RISKMGR_DLL ForwardParallel: public ScalarShift,
	public virtual Additive {
public:
	friend class ForwardParallelHelper;
	static CClassConstSP const TYPE;
	const static string NAME;
	const static double DEFAULT_SHIFT;

	/** What an object must implement to be tweakable for FORWARD_PARALLEL */
	class RISKMGR_DLL Shift{
	public:
		friend class ForwardParallelHelper;
		static CClassConstSP const TYPE;
		virtual ~Shift();
		/** Returns the name of the yield curve - used to determine
		whether to tweak the object */
		virtual string sensName(ForwardParallel* shift) const = 0;
		/** Shifts the object using given shift. Return true to make the
		infrastructure keep tweaking the components within the object
		which implements this interface */
		virtual bool sensShift(ForwardParallel* shift) = 0;
	};
	typedef Shift IShift;

	/** What an object must implement to be able to perform a restorable
	tweak for FORWARD_PARALLEL. This is used by {@link SensMgr SensMgr} to
	determine whether it should try and tweak a particular object. */
	class RISKMGR_DLL RestorableShift: public virtual Shift{
	public:
		friend class ForwardParallelHelper;
		static CClassConstSP const TYPE;
		virtual ~RestorableShift();
		/** Restores the object to its original form */
		virtual void sensRestore(ForwardParallel* shift) = 0;
	};
	typedef RestorableShift IRestorableShift;

	/** constructor with explicit shift size */
	ForwardParallel(double shiftSize);

	/** constructor with explicit shift size */
	ForwardParallel(double     shiftSize,
		IModel*    model,
		Control*   control);

	/** Once used to make a shift, this reports the appropriate divisor
	for this sensitivity */
	double divisor() const;

	/** returns the interface identifying what an object has to do in order
	to be support the tweak that this object represents */
	CClassConstSP shiftInterface() const;

	/** returns the interface identifying what an object has to do in order
	to be support the tweak that this object represents which is also
	restorable */
	CClassConstSP restorableShiftInterface() const;


	// methods to shift given objects

	/** Returns true if the supplied object matches the supplied name
	for this sensitivity. The object must implement this sensitivity's
	Shift interface */
	virtual bool nameMatches(const OutputName&         name,
		IObjectConstSP            obj);

	/** Appends the name(s) of the supplied object with respect to
	this sensitivity to the supplied list. The object must
	implement this sensitivity's Shift interface */
	virtual void appendName(OutputNameArray&          namesList,
		IObjectConstSP            obj);

	/**
	* @param obj The object to shift. The object must implement the
	ForwardParallel.Shift interface. The return value indicates whether
	or not components of this object need to be tweaked ie true:
	infrastructure should continue to recurse through components
	tweaking them; false: the infrastructure shouldn't touch any
	components within this object */
	virtual bool shift(IObjectSP obj);


	/**
	*
	* @param obj The object to shift. The
	object must implement the ForwardParallel.RestorableShift interface
	*/
	virtual void restore(IObjectSP obj);

private:
	/** for reflection */
	ForwardParallel();
	ForwardParallel(const ForwardParallel &rhs);
	ForwardParallel& operator=(const ForwardParallel& rhs);
	};

typedef smartConstPtr<ForwardParallel> ForwardParallelConstSP;
typedef smartPtr<ForwardParallel> ForwardParallelSP;

DRLIB_END_NAMESPACE

#endif
