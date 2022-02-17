//----------------------------------------------------------------------------
//
//   Group       : EDG Derivatives Research
//
//   Filename    : SpreadSheetMode.hpp
//
//   Description : sets up environment used by spreadsheets
//
//   Author      : Andrew J Swain
//
//   Date        : 21 August 2001
//
//
//----------------------------------------------------------------------------

#ifndef SPREADSHEETMODE_HPP
#define SPREADSHEETMODE_HPP

DRLIB_BEGIN_NAMESPACE

/** sets up environment used by spreadsheets */
class TOOLKIT_DLL SpreadSheetMode {
public:
    /** turns it on */
    static void on();

    /** turns it off */
    static void off();

    /** is it on or off ? */
    static bool isOn();

    /** Has the user requested the calculation to abort. Note returns false
        if spreadsheet mode is off. Do not call this function too frequently
        as it has to call back to Excel - which may hit performance */
    static bool abort();

    typedef bool (SpreadSheetAbort)();
    /** Set the method to call to see if the calculation should be aborted.
        This is for use by, for example, the XL Kit. General users should 
        not use it */
    static void setAbortMethod(SpreadSheetAbort* abortMethod);

private:
    static SpreadSheetAbort* abortMethod;
    static bool IS_ON;
};

DRLIB_END_NAMESPACE

#endif




