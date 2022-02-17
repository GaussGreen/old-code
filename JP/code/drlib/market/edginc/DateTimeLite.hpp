//----------------------------------------------------------------------------
//
// Group       : CH Quantitative Research
//
// Filename    : DateTimeLite.hpp
//
// Description : Representing a light version of the DateTime class
//				 Used in MC. Enables performance.
//
// Date        : Aug 2006
//
//----------------------------------------------------------------------------


#ifndef DATETIMELITE_HPP
#define DATETIMELITE_HPP

#include "edginc/VirtualDestructorBase.hpp"

DRLIB_BEGIN_NAMESPACE

class MARKET_DLL DateTimeLite : public virtual VirtualDestructorBase
{
public:
    /** Copy Constructor */
	DateTimeLite(const DateTimeLite& rhs);

	/** Constructor to convert from a DateTime */
	DateTimeLite(const DateTime& dateTime = DateTime());

	/** virtual destructor */
	virtual ~DateTimeLite() {}

    /** equal to operator */
	DateTimeLite& operator=(const DateTimeLite& rhs);

    /** equal to operator - rhs is a DateTime */
	DateTimeLite& operator=(const DateTime& rhs);

	/** less than operator */
	bool operator<(const DateTimeLite& rhs) const;

	/** less than operator - rhs is a DateTime */
	bool operator<(const DateTime& rhs) const;

	/** greater than operator */
	bool operator>(const DateTimeLite& rhs) const;

	/** greater than operator - rhs is a DateTime */
	bool operator>(const DateTime& rhs) const;

	/** greater than equal to operator */
	bool operator>=(const DateTimeLite& rhs) const;

	/** greater than equal to operator - rhs is a DateTime */
	bool operator>=(const DateTime& rhs) const;

	/** less than equal to operator */
	bool operator<=(const DateTime& rhs) const;

	/** less than equal to operator - rhs is a DateTime */
	bool operator<=(const DateTimeLite& rhs) const;

	/** less than operator */
	bool operator==(const DateTime& rhs) const;

	/** less than operator - rhs is a DateTime */
	bool operator==(const DateTimeLite& rhs) const;

	/** minus operator */
	int operator-(const DateTimeLite& rhs) const;

	/** minus operator - rhs is a DateTime */
	int operator-(const DateTime& rhs) const;

	/** casting into a DateTime */
	operator DateTime() const;

	/** equivalent to calling getDate() on a DateTime object */
	int getTime() const;

private:
// ## Member variables follow

    /** this should correspond to the internal representation inside the DateTime object */
	int time;
};

typedef vector<DateTimeLite> DateTimeLiteVector;
typedef refCountPtr<DateTimeLiteVector> DateTimeLiteVectorSP;
typedef refCountPtr<DateTimeLiteVector> DateTimeLiteVectorConstSP;

DRLIB_END_NAMESPACE
#endif
