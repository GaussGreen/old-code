// ----------------------------------------------------------------------------
//
//   Group       : GCCG QR&D
//
//   Filename    : EnergyContractLabel.hpp
//
//   Description : Labels representing energy futures contract identifiers.
//                 e.g. Jan02 representing January, 2002
//
//   Author      : Sean Chen
//
//   Date        : May 11 2005
//
// ----------------------------------------------------------------------------// ----------------------------------------------------------------------------

#ifndef _EnergyContractLabel_
#define _EnergyContractLabel_

#include "edginc/config.hpp"
#include "edginc/Expiry.hpp"

#include <string>
using namespace std;

DRLIB_BEGIN_NAMESPACE


class TOOLKIT_DLL EnergyContractLabel : public Expiry
{
    public:    
        static CClassConstSP const TYPE;
        friend class EnergyContractLabelHelper;

        EnergyContractLabel();
        EnergyContractLabel(int theMonth, int theYear);
        EnergyContractLabel(const string& theLabelString);
        EnergyContractLabel(const EnergyContractLabel& v);
        virtual ~EnergyContractLabel();

        virtual IObject* clone() const;

        // operators
        EnergyContractLabel& operator=(const EnergyContractLabel& v);
        bool operator==(const EnergyContractLabel& v) const;
        bool operator!=(const EnergyContractLabel& v) const;
        bool operator<=(const EnergyContractLabel& v) const;
        bool operator>=(const EnergyContractLabel& v) const;
        bool operator<(const EnergyContractLabel& v) const;
        bool operator>(const EnergyContractLabel&v) const;

        void setLabel(int theMonth, int theYear);
        void setLabel(const string& theLabelString);
        string asString() const;
        void convertDate(const string& benchMark);
        void addMonths(int numMonths);
        EnergyContractLabel nextContractLabel() const;
        int getMonth() const;
        int getYear() const;
        void setMonth(int theMonth);
        void setYear(int theYear);

        bool isValid() const;
        bool isNull() const;

        string toString() const;
        DateTime toDate(const DateTime& date) const;
        bool equals(const Expiry* expiry) const;

    private:

        int month;
        int year;
        int contractIncrement;  // set as one month and fixed
};


typedef smartConstPtr<EnergyContractLabel> EnergyContractLabelConstSP;
typedef smartPtr<EnergyContractLabel> EnergyContractLabelSP;

typedef array<EnergyContractLabelSP, EnergyContractLabel> EnergyContractLabelArray;
typedef smartPtr<EnergyContractLabelArray> EnergyContractLabelArraySP;


DRLIB_END_NAMESPACE

#endif
