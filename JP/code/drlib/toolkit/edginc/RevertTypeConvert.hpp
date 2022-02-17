#ifndef EDR_REVERT_HPP
#define EDR_REVERT_HPP

DRLIB_BEGIN_NAMESPACE

/** Defines interface for the ability for one type to convert itself
    into another type.
    This interface is useful for the classes that implements the ITypeConvert.
    This interface provides the 'revert' method that convert a object into an object
    of a different type depending on the interface type (i.e. PYRAMID, KAPITAL,...)

    example: a DependenceMakerGauss object is reverted into a DependenceMaker object 
*/
class TOOLKIT_DLL IRevertTypeConvert {
public:
    static CClassConstSP const TYPE; 
    static string const PYRAMID;

    // given a client, convert myself into the representation that client
    // would have created
    // e.g. for something that uses a "wrapper" to bypass IMS polymorphic 
    // representations, output from this would be a wrapper class with the
    // type string set appropriately and the relevant type instance
    // populated
    // any conversion that contains the origianl should use a clone
    virtual IObjectSP revert(const string& client) const = 0;    

    virtual ~IRevertTypeConvert(); 

protected:
    IRevertTypeConvert(); 

private:

    IRevertTypeConvert(const IRevertTypeConvert& rhs); 
    IRevertTypeConvert& operator=(const IRevertTypeConvert& rhs); 
};


DRLIB_END_NAMESPACE

#endif
