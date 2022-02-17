/**
 * @file ExpiryWindow.hpp
 */

#ifndef QLIB_ExpiryWindow_H
#define QLIB_ExpiryWindow_H

#include "edginc/FORWARD_DECLARE.hpp"

DRLIB_BEGIN_NAMESPACE

FORWARD_DECLARE(ExpiryWindow)
FORWARD_DECLARE(Expiry)

#ifndef QLIB_EXPIRYWINDOW_CPP
EXTERN_TEMPLATE(class TOOLKIT_DLL_SP smartConstPtr<ExpiryWindow>);
EXTERN_TEMPLATE(class TOOLKIT_DLL_SP smartPtr<ExpiryWindow>);
EXTERN_TEMPLATE(class TOOLKIT_DLL array<ExpiryWindowSP _COMMA_ ExpiryWindow>);
#else
INSTANTIATE_TEMPLATE(class TOOLKIT_DLL smartConstPtr<ExpiryWindow>);
INSTANTIATE_TEMPLATE(class TOOLKIT_DLL smartPtr<ExpiryWindow>);
INSTANTIATE_TEMPLATE(class TOOLKIT_DLL array<ExpiryWindowSP _COMMA_ ExpiryWindow>);
#endif

/**
 * A time window around an Expiry, defining the period over which values
 * defined at the Expiry may have some influence due to interpolation.
 *
 * Typically previous and next are expiry's predecessor and successor in the
 * sequence (term structure) from which it's drawn.
 */

class TOOLKIT_DLL ExpiryWindow: public CObject {

    static void load(CClassSP& clazz);

public:

    static CClassConstSP const TYPE;

    /**
     * Predecessor in sequence from which expiry is drawn --- beginning of the
     * window.
     *
     * If NULL, the window extends back to the beginning of time.
     */

    ExpiryConstSP previous;

    /**
     * Expiry around which the window is defined.
     *
     * Non-NULL.
     */

    ExpiryConstSP expiry;

    /**
     * Successor in sequence from which expiry is drawn --- end of the window.
     *
     * If NULL, the window extends out to the end of time.
     */

    ExpiryConstSP next;

    /**
     * Constructor.
     *
     * See also SP(), around().
     */

    ExpiryWindow(ExpiryConstSP previous, ExpiryConstSP expiry,
                 ExpiryConstSP next);

    /**
     * Constructor returning a smartPtr.
     */

    static ExpiryWindowSP SP(ExpiryConstSP previous, ExpiryConstSP expiry,
                             ExpiryConstSP next);

    /**
     * An ExpiryWindow around an Expiry in a given ExpiryArray.
     *
     * @a expiry must be a member of @a expiries.  Returns a window around @a
     * expiry bounded by its predecessor (if any) and successor (if any).
     */

    static ExpiryWindowSP around(ExpiryArrayConstSP expiries,
                                 ExpiryConstSP expiry);

    static ExpiryWindowSP find(ExpiryWindowArrayConstSP windows,
                               ExpiryConstSP expiry);

    static ExpiryWindowArraySP series(ExpiryArrayConstSP expiries);

    static ExpiryArrayConstSP expiries(ExpiryWindowArrayConstSP windows);

    ~ExpiryWindow();

    string toString() const;
};

DRLIB_END_NAMESPACE

#endif
