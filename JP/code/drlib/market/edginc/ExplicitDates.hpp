//----------------------------------------------------------------------------
//
//   Group       : EDG Derivatives Research
//
//   Filename    : ExplicitDates.hpp
//
//   Description : A type of IDateBuilder that provides the convenience of
//                 wrapping a DateTimeArray up behind the IDateBuilder
//                 interface and providing type conversion between DateTimeArray
//                 and IDateBuilder.
//                 This provides an intermediate step before the full
//                 polymorphic IDateBuilder interface is supported (by IMS, for
//                 example).
//
//   Author      : Ben Papworth
//
//   Date        : 24 May 2006
//
//
//----------------------------------------------------------------------------

#ifndef EDR_EXPLICITDATES_HPP
#define EDR_EXPLICITDATES_HPP

#include "edginc/config.hpp"
#include "edginc/DateTime.hpp"
#include "edginc/DateBuilder.hpp"
#include "edginc/Format.hpp"

DRLIB_BEGIN_NAMESPACE

/** Timeline constructed from a DateTimeArray */
class MARKET_DLL ExplicitDates: public CObject, virtual public ISPIDateBuilder {
public:
    static CClassConstSP const TYPE;

    // Construct an ExplicitDates from a DateTimeArray
    ExplicitDates( const DateTimeArray& newDates );

    /** ISPIDateBuilder implementation
    * constructs dates - nothing to do here */
    virtual void constructDates(const IMultiFactors*          assets,
        const ObservationSourceArray& sources,
        const Holiday*                ccyHols,
        bool                          excludeAssetHols);

        /** IDateBuilder implementation
        * Builds explicitly defined sample line */
        virtual DateTimeArraySP dates() const;


        /** Returns the first date in the list */
        virtual const DateTime start() const;

        /** Returns the first date in the list */
        virtual const DateTime end() const;

        /** Returns the number of dates */
        virtual int size() const;

        /** Returns a particular date from the final array of dates */
        virtual const DateTime date(int index) const;

        /** Set the time in all the generated dates */
        virtual void setTime(int time);

        /** which time is currently being used */
        virtual int time() const;

        virtual void validatePop2Object();

        // Conversion function from DateTimeArray to IDateBuilder
        static IObjectSP fromDateTimeArray( const IObjectSP& object,
            CClassConstSP    requiredType) {
                DateTimeArray& newDates = dynamic_cast<DateTimeArray&>(*object);
                return IObjectSP(new ExplicitDates(newDates));
            }

private:
    DateTimeArray theDates;     //!< Array of dates representing date sequence

    static void load(CClassSP& clazz) {
        REGISTER(ExplicitDates, clazz);
        SUPERCLASS(CObject);
        IMPLEMENTS(ISPIDateBuilder);
        EMPTY_SHELL_METHOD(defaultExplicitDates);
        FIELD(theDates, "Dates");
        clazz->setPublic(); // make visible to EAS/spreadsheet

        // Register conversion magic to allow a raw DateTimeArray to be
        // transformed into an <<IDateBuilder>> ExplicitDates. The *interface*
        // type is used to make this type conversion available to methods
        // requiring an IDateBuilder (not an ExplicitDates).
        registerObjectFromArrayMethod( DateTimeArray::TYPE,
            IDateBuilder::TYPE,
            &fromDateTimeArray );

        // Register conversion magic to allow a raw DateTimeArray to be
        // transformed into an <<IDateBuilder>> ExplicitDates. The *interface*
        // type is used to make this type conversion available to methods
        // requiring an IDateBuilder (not an ExplicitDates).
        registerObjectFromArrayMethod( DateTimeArray::TYPE,
            ISPIDateBuilder::TYPE,
            &fromDateTimeArray );
    }

    static IObject* defaultExplicitDates() {
        return new ExplicitDates();
    }

    // Default constructor
    ExplicitDates(): CObject(TYPE) {}

};

DRLIB_END_NAMESPACE
#endif
