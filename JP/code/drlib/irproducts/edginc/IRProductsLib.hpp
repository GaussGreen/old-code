#ifndef QLIB_IR_PRODUCTS_LIB_HPP
#define QLIB_IR_PRODUCTS_LIB_HPP

DRLIB_BEGIN_NAMESPACE

class IRPRODUCTS_DLL IRProductsLib{
public:
    /** The sole purpose of this method is to ensure that the linker 
        includes all symbols out of the riskmr directory. Many symbols
        are automatically linked because they are used by other classes
        which are already included.
        
        An example of symbols that could be dropped would be an entire class
        representing a product which was referenced by no other classes */
    static void linkInClasses();
private:
    IRProductsLib();
};

DRLIB_END_NAMESPACE

#endif
