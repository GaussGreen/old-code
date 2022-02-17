//----------------------------------------------------------------------------
//
//   Group       : EDG Derivatives Research
//
//   Filename    : CompoundBasis.hpp
//
//   Description : How rates can be compounded
//
//   Author      : Andrew J Swain
//
//   Date        : 30 January 2001
//
//
//----------------------------------------------------------------------------

#ifndef COMPOUNDBASIS_HPP
#define COMPOUNDBASIS_HPP

#include <string>

using namespace std;  // string

DRLIB_BEGIN_NAMESPACE

class MARKET_DLL CompoundBasis {
public:

    static const int ANNUAL;
    static const int SEMI_ANNUAL;
    static const int QUARTERLY;
    static const int MONTHLY;
    static const int WEEKLY;
    static const int DAILY;
    static const int CONTINUOUS;
    static const int SIMPLE;

    static string toString(int basis);
    static int toInt(const std::string& basis_str);

private:
    CompoundBasis();
};

DRLIB_END_NAMESPACE

#endif

