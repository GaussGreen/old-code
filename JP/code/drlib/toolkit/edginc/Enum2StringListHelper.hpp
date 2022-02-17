//----------------------------------------------------------------------------
//
//   Group       : EDG Derivatives Research
//
//   Filename    : Enum2StringListHelper.hpp
//
//   Author      : Regis Guichard
//
//   Date        : 30 May 03
//
//
//----------------------------------------------------------------------------

#ifndef EDR_ENUM2STRING_LIST_HELPER_HPP
#define EDR_ENUM2STRING_LIST_HELPER_HPP

#include <string>

DRLIB_BEGIN_NAMESPACE

/** Returns the name of the corresponding type in the form of a string, ie
    does what CClass::getName does for IObject's.
    Should be overloaded appropriately - if not, the empty string will be 
    returned.
    Will be used in Enum2StringListHelper to return the appropriate 
    error message */
template<class T>
string nameForType(T* notused){
    return "";  // empty by default
}

/** Helper class that facilitates the handling of string-type user inputs. 
    Enum2StringListHelper helps convert an enumeration list into a string
    list (and vice versa). 
    The EnumList type should contain an emumeration list starting at 0 
    and terminated by 'NB_ENUMS', as follows:
            struct MyEnumListContainer{
                enum {
                    ENUM1 = 0, 
                    ENUM2,
                    ...,
                    ENUMN,
                    NB_ENUMS
                };
            }; 
    If, in addition, MyEnumListContainer defines the integral type 'defaultIndex'
            struct MyEnumListContainer{
                // as above
                enum {
                    defaultIndex = ENUM2    // for instance
                };
            }; 
    then Enum2StringListHelper's 'getDefaultName' method may be used.
    The user should initialized the private string array 'names', as follows
        typedef Enum2StringListHelper<MyEnumListContainer> MyEnumListContainerHelper
        string MyEnumListContainerHelper::names[
            MyEnumListContainerHelper::EnumList::NB_ENUMS] = {
                "ENUM1",
                "ENUM2",
                ...,
                "ENUMN"
        }; */
template<class EList>
struct Enum2StringListHelper{
    typedef EList EnumList;

    /** Returns the name corresponding to the default index
        as defined in EnumList */
    static string getDefaultName(){
        return names[EnumList::defaultIndex];
    }

    /** Returns the list of names in the form of a comma-separated string */
    static string getNameList(){
        string buffer;
        int iType = 0;
        for (; iType < EnumList::NB_ENUMS - 1; ++iType){
            buffer += "'";
            buffer += names[iType];
            buffer += "', ";
        }
        buffer += "'";
        buffer += names[EnumList::NB_ENUMS - 1];
        buffer += "'";
        return buffer;
    }

    /** Given a name, returns the correpsonding index in the enum list.
        Throws an exception if not found */
    static int getIndex(const string& name){
        static const string method = "Enum2StringListHelper<"
                                     + nameForType<EnumList>(0)
                                     + ">::getIndex";
        int iType = 0;
        for (; iType < EnumList::NB_ENUMS; ++iType){
            if (name == names[iType]){
                return iType;
            }
        }
        throw ModelException(method,
                             "The string '" 
                             + name 
                             + "' could not be found in the following list:\n"
                             + getNameList());
    }

    /** Given an index, returns the matching name.
        Throws an exception if error bound */
    static string getName(int index){
        static const string method = "Enum2StringListHelper<"
                                     + nameForType<EnumList>(0)
                                     + ">::getName";
        if (index < 0
            || index >= EnumList::NB_ENUMS){
            throw ModelException(method,
                                 "Error bound");
        }
        return names[index];
    }

private:
    static string names[];
};

DRLIB_END_NAMESPACE
#endif 
