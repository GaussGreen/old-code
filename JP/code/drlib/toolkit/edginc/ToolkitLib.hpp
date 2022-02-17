#ifndef EDG_TOOLKIT_LIB_H
#define EDG_TOOLKIT_LIB_H

DRLIB_BEGIN_NAMESPACE

/** class insures all classes in the toolkit directory are linked in. This is
    probably a temporary measure - will review once the addin is working */
class TOOLKIT_DLL CToolkitLib{
public:
    /** calling this function ensures that the toolkit classes will be linked
        in by the compiler */
    static void linkInClasses();
private:
    // cannot be instantiated
    CToolkitLib();
};
        
DRLIB_END_NAMESPACE

#endif
