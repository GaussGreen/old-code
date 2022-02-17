
#ifndef EDG_PRODUCTS_LIB_HPP
#define EDG_PRODUCTS_LIB_HPP

DRLIB_BEGIN_NAMESPACE

class PRODUCTS_DLL CProductsLib{
public:
    /** The sole purpose of this method is to ensure that the linker 
        includes all symbols out of the riskmr directory. Many symbols
        are automatically linked because they are used by other classes
        which are already included.
        
        An example of symbols that could be dropped would be an entire class
        representing a product which was referenced by no other classes */
    static void linkInClasses();
private:
    CProductsLib();
};

DRLIB_END_NAMESPACE

#endif
