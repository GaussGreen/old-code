// ----------------------------------------------------------------------------
//
//   Group       : GCCG QR&D
//
//   Filename    : EnergyContractLabel.cpp
//
//   Description : Labels representing energy futures contract identifiers.
//                 e.g. Jan02 representing January, 2002
//
//   Author      : Sean Chen
//
//   Date        : May 11 2005
//
// ----------------------------------------------------------------------------

#include "edginc/config.hpp"
#include "edginc/EnergyContractLabel.hpp"
#include "edginc/ModelException.hpp"
#include "edginc/DateTime.hpp"
#include "edginc/Format.hpp"
#include "edginc/Atomic.hpp"

DRLIB_BEGIN_NAMESPACE

const static char* months[] = {"Jan","Feb","Mar","Apr","May","Jun",
                               "Jul","Aug","Sep","Oct","Nov","Dec"};
#define MONTHS_PER_YEAR     12

EnergyContractLabel::EnergyContractLabel() :
              Expiry(TYPE), month(0), year(0), contractIncrement(1)
{
}


EnergyContractLabel::EnergyContractLabel(int theMonth, int theYear) 
: Expiry(TYPE), month(theMonth), year(theYear), contractIncrement(1)
{
}

// only works if label is month and Year format for the mean time
EnergyContractLabel::EnergyContractLabel(const string& labelString)
: Expiry(TYPE), contractIncrement(1)
{
    // pull out delivery month from label - borrowing IMM label conversion
    // labels are in the format (mmmyy or mmmyyyy)
    
    try 
    {
        convertDate(labelString);
    }
    catch (exception& e)
    {
                throw ModelException(e);
    }
}

EnergyContractLabel::EnergyContractLabel(const EnergyContractLabel& v): Expiry(TYPE),
     month(v.month), year(v.year), contractIncrement(v.contractIncrement)
{
}


/** Explicit implementation for performance */
IObject* EnergyContractLabel::clone() const
{
    int& count = getRefCount();
    if (count == 0){
        return new EnergyContractLabel(*this);
    } else {
        return  const_cast<IObject*>((const IObject*)this);
    }
}

EnergyContractLabel::~EnergyContractLabel()
{
}

EnergyContractLabel& EnergyContractLabel::operator=(const EnergyContractLabel &v)
{
    if (&v != this)
    {
        month = v.month;
        year = v.year;
        contractIncrement = v.contractIncrement;
    }
    return *this;
}

void EnergyContractLabel::setLabel(int theMonth, int theYear)
{
    month = theMonth;
    year = theYear;
}

void EnergyContractLabel::setLabel(const string& labelString)
{
    // pull out delivery month from label - borrowing IMM label conversion
    // labels are in the format (mmmyy or mmmyyyy)
    try
    {
        convertDate(labelString);
    }
    catch (exception& e)
    {
            throw ModelException(e);
    }
}

bool EnergyContractLabel::operator==(const EnergyContractLabel& v) const
{
    return     month == v.month && 
               year == v.year &&
               contractIncrement == v.contractIncrement;
}

bool EnergyContractLabel::operator!=(const EnergyContractLabel& v) const
{
    return !(*this == v);
}

bool EnergyContractLabel::operator<=(const EnergyContractLabel& v) const
{
    if (year > v.year) return false;
    else if (year < v.year) return true;
    else if (year == v.year && month <= v.month) return true;
    else return false;
}

bool EnergyContractLabel::operator>=(const EnergyContractLabel& v) const
{
    return !(*this < v);
}

bool EnergyContractLabel::operator<(const EnergyContractLabel& v) const
{
    if (*this == v) return false;
    else if (year > v.year) return false;
    else if (year < v.year) return true;
    else if (year == v.year && month < v.month) return true;
    else return false;
}

bool EnergyContractLabel::operator>(const EnergyContractLabel&v) const
{    
    return !(*this <= v);
}

void EnergyContractLabel::addMonths(int numMonths)
{
    DateTime::MonthDayYear tempDMY(1, month, year);
    tempDMY.month += contractIncrement * numMonths;
    DateTime tempDate = tempDMY.toDateTime();
    month = tempDate.toMDY().month;
    year = tempDate.toMDY().year;
}

EnergyContractLabel EnergyContractLabel::nextContractLabel() const
{
    // slow but useful
    EnergyContractLabel aNewLabel(*this);
    aNewLabel.addMonths(1);
    return aNewLabel;
}

string EnergyContractLabel::asString() const
{
    return Format::toString("%s%.2d", months[month-1], (year<100)?year:year-2000);
}

int EnergyContractLabel::getMonth() const
{
    return month;
}

int EnergyContractLabel::getYear() const
{
    return year;
}

void EnergyContractLabel::setMonth(int theMonth)
{
    month = theMonth;
}

void EnergyContractLabel::setYear(int theYear)
{
    year = theYear;
}

bool EnergyContractLabel::isValid() const
{
    return (month && year);
}

bool EnergyContractLabel::isNull() const
{
    return (month == 0 && year == 0);
}

string EnergyContractLabel::toString() const
{
    return asString();
}

DateTime EnergyContractLabel::toDate(const DateTime& date) const
{
    return date;
}

bool EnergyContractLabel::equals(const Expiry* expiry) const
{
    return *this == *(dynamic_cast<const EnergyContractLabel*>(expiry));
}


//// convert IMM mmm[yy]yy into int month and year
void EnergyContractLabel::convertDate(const string& dateString) {
    
    try {
        /* see if string is remotely in date format */
    int len = dateString.length();
        if (len != 5 && len != 7 )
    {
        string msg="date string (" + dateString + ") not in mmmyy or mmmyyyy format";
        throw ModelException("EnergyContractLabel::convertDate()", msg);
    }

    // chop it into dd, mmm & yyyy
      
        
    string y;
    string m(dateString.substr(0, 3));
    if (len == 5)
        y = dateString.substr(3, 2); 
    else
        y = dateString.substr(3, 4);
        
        // day & year are easy
    int aYear;
        
    sscanf(y.c_str(), "%d", &aYear);
    if(len == 5)
    {
        if (aYear>=0 && aYear<70)
        {
            aYear += 2000;
        }
        else
        {
            aYear += 1900;    
        }
    }
    year = aYear;
    month = 0;

    bool found = false;
    while (month < MONTHS_PER_YEAR && !found)
    {
        if (CString::equalsIgnoreCase(m, months[month]))
        {
                found = true;
        }
        month++; // if found, this converts from 0-1 to 1-12
    }
    if (!found)
    {
        throw ModelException("EnergyContractLabel::convertDate()", "couldn't find month ("+m+")");
    }
    
    }
    catch (exception& e) {
        throw ModelException(e);
    }
}

class EnergyContractLabelHelper
{
public:
    
    static void load(CClassSP& clazz){
        clazz->setPublic(); // make visible to EAS/spreadsheet
        REGISTER(EnergyContractLabel, clazz);
        SUPERCLASS(Expiry);
        EMPTY_SHELL_METHOD(defaultEnergyContractLabel);
        FIELD(month, "MMM of label");
        FIELD(year, "[yy]yy of label");
        FIELD(contractIncrement, "increment btw labels");
        
    }

    static IObject* defaultEnergyContractLabel(){
        return new EnergyContractLabel();
    }
};

CClassConstSP const EnergyContractLabel::TYPE = CClass::registerClassLoadMethod(
    "EnergyContractLabel", typeid(EnergyContractLabel), EnergyContractLabelHelper::load);

DEFINE_TEMPLATE_TYPE(EnergyContractLabelArray);

DRLIB_END_NAMESPACE
