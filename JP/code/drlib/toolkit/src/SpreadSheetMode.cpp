//----------------------------------------------------------------------------
//
//   Group       : EDG Derivatives Research
//
//   Filename    : SpreadSheetMode.cpp
//
//   Description : sets up environment used by spreadsheets
//
//   Author      : Andrew J Swain
//
//   Date        : 21 August 2001
//
//
//----------------------------------------------------------------------------

#include "edginc/config.hpp"
#include "edginc/SpreadSheetMode.hpp"

DRLIB_BEGIN_NAMESPACE

bool   SpreadSheetMode::IS_ON = false;
SpreadSheetMode::SpreadSheetAbort* SpreadSheetMode::abortMethod = 0;

/** Set the method to call to see if the calculation should be aborted.
    This is for use by, for example, the XL Kit. General users should 
    not use it */
void SpreadSheetMode::setAbortMethod(SpreadSheetAbort* abortMethod){
    SpreadSheetMode::abortMethod = abortMethod;
}

/** Has the user requested the calculation to abort. Note returns false
    if spreadsheet mode is off. Do not call this function too frequently
    as it has to call back to Excel - which may hit performance */
bool SpreadSheetMode::abort(){
    if (IS_ON && abortMethod){
        return abortMethod();
    }
    return false;
}

void SpreadSheetMode::on() {
    IS_ON = true;
}

void SpreadSheetMode::off() {
    IS_ON = false;
}

bool SpreadSheetMode::isOn() {
    return IS_ON;
}

DRLIB_END_NAMESPACE

