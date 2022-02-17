//----------------------------------------------------------------------------
//
//   Group       : Equity Derivatives Research
//
//   Filename    : GDRDate.hpp
//
//   Description : Wrapper for the date in DR Interface
//
//   Author      : Mark A Robson
//
//   Date        : 4 Dec 2003
//
//
//----------------------------------------------------------------------------

#ifndef EDR_GDRDATE_HPP
#define EDR_GDRDATE_HPP
#include "edginc/DRObject.hpp"
#include "edginc/DRAnalyticsInterface.h"

DRLIB_BEGIN_NAMESPACE
/** Wrapper for the date in DR Interface. This is possibly the worst
    representation of a date ever invented. Please do not use. */
class TOOLKIT_DLL GDRDate: public CObject,
               public virtual IDRObject{
public:
    static CClassConstSP const TYPE;

    virtual ~GDRDate();

    /** Are these objects equal (ie contain the same data) */
    virtual bool equals(IDRObjectConstSP drObj) const;

    /** write object out to writer */
    virtual void write(const string& tag, Writer* writer) const;
    /** populate an empty object from description */
    virtual void import(Reader::Node* elem, Reader* reader);
    /** write object out in 'output' format - ie suitable for comparing
        regression files with */
    virtual void outputWrite(const string& linePrefix,
                             const string& prefix, ostream& stream) const;
    // just returns this
    virtual IObject* clone() const;

    /** Returns a string representation of this date */
    string toString() const;

    /** Create a new date offset from this date */
    GDRDate* rollDate(int offset);

    /** Returns the date being wrapped */
    DRDate dateValue() const;
    
    /** Returns ALIB style date */
    int tDate() const;

    //// constructors
    static GDRDate* create(DRDate date);
    static GDRDate* create(const char* date);
    static GDRDate* create(int tDate); // ALIB style date
private:
    class AddinFnc;
    friend class AddinFnc;
    static DRDate convertDate(const string& date);
    GDRDate(DRDate date);
    GDRDate(const GDRDate& rhs); // don't use
    GDRDate& operator=(const GDRDate& rhs); // don't use
    static void load(CClassSP& clazz);
    GDRDate();
    static IObject* defaultConstructor();
    //// fields
    const DRDate date; // $unregistered. Immutable
};
DRLIB_END_NAMESPACE
#endif
