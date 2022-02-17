/**
   This file contains various functions that are useful to be able to call
   within a debugger. Originally the functions were in the relevant .cpp file
   but Developer Studio can't cope with this when the application is built as
   a series of dlls. To fix this we include this file once per dll by
   including it each XXXXLib.cpp where XXXX is the directory name.
   Note that these functions are only for use within the debugger - which
   is why they're not inside a class etc.

   DO NOT INCLUDE THIS FILE FROM ANY OTHER FILE APART FROM THE XXXXLib.cpp FILES

   If you want to put some other functions in here which require access to
   other directories, then create eg MarketDebug.hpp and then include this file
   from there. Don't forget to replace the QLIB_TOOLKIT below with eg
   QLIB_MARKET. Then change eg MarketLib.cpp to include MarketDebug.hpp rather
   than this file.
 */
#ifndef QLIB_TOOLKITDEBUG_HPP
#define QLIB_TOOLKITDEBUG_HPP
#if defined (DEBUG) && (defined(QLIB_BUILD_DLL) || defined(QLIB_TOOLKIT))
// try to ensure that this only gets pulled in once per executable

#include "edginc/DateTime.hpp"
#include "edginc/XMLWriter.hpp"
//// handy for debugging. DO NOT USE ANYWHERE ELSE. Nasty memory model
//// Only use in the debugger. Function not on DateTime class as otherwise
//// you have to do drlib::DateTime::pda(dates) to use - which is a bit tedious
char**  pdac(const NAMESPACE::DateTimeArray& dates){
    if (&dates == 0){
        return 0; // avoid crash
    }
    try{
        char** dtArray = new char*[dates.size()+1];
        for (int i = 0; i < dates.size(); i++){
            const string& str = dates[i].toString();
            dtArray[i] = new char[str.size()+1];
            strcpy(dtArray[i], str.c_str());
        }
        dtArray[dates.size()] = 0;
        return dtArray;
    } catch (exception&){
        return 0;
    }        
}

char** pda(NAMESPACE::DateTimeArray& dates){
    // pdac seems fine under the debugger but making pda call pdac doesn't seem
    // to work(?)
    static char dtArray[100][25];
    static char* dtPtrArray[100];
    if (&dates == 0){
        return 0; // avoid crash
    }
    try{
        int numDatesToDo =dates.size() < 99? dates.size(): 99;
        for (int i = 0; i < 99; i++){
            strcpy(dtArray[i], i < numDatesToDo?
                   dates[i].toString().c_str(): "");
            dtPtrArray[i] = dtArray[i];
        }
        dtPtrArray[99] = dtArray[99];
        strcpy(dtArray[99], dates.size() < 100? "": "List truncated");
        return dtPtrArray;
    } catch (exception&){
        return 0;
    }        
}

//// Put in as Vis Studio can't seem to see the p() method outside of the
//// dll
const char* pd(const NAMESPACE::DateTime& date) {
    return date.p();
}

const char* pdc(NAMESPACE::DateTime& date) {
    return date.p();
}

//// handy way of dumping objects to file from within the debugger
void xmlDebug(NAMESPACE::IObject* object, const char* filename){
    NAMESPACE::XMLWriter xml(filename);
    object->write("OBJECT", &xml);
}

//// handy for debugging. DO NOT USE ANYWHERE ELSE. 
//// You can use these to jump from a virtual parent to the actual object
//// itself. Just manually cast the return object in the debugger
void* cast(NAMESPACE::IObject* obj){
    return (dynamic_cast<void*>(obj));
}

const void* castc(const NAMESPACE::IObject* obj){
    return (dynamic_cast<const void*>(obj));
}


#endif
#endif
