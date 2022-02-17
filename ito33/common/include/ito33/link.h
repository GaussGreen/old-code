/////////////////////////////////////////////////////////////////////////////
// Name:        ito33/link.h
// Purpose:     ITO33_FORCE_LINK_MODULE() & related macros
// Author:      Vadim Zeitlin
// Created:     2004-05-12
// RCS-ID:      $Id: link.h,v 1.2 2004/10/05 09:13:35 pedro Exp $
// Copyright:   (c) 2004 Trilemma LLP
/////////////////////////////////////////////////////////////////////////////

/**
    @file   ito33/link.h
    @brief  Provides macros for forcing linking in the program modules.

    First, an important fact: some, and maybe even most, linkers (notably MSVC
    one) don't link in the object files if they don't contain any symbols
    referenced from other object files (which are linked in).

    Second, the problem: sometimes we don't want to export anything from a
    module but still want it to be linked in because it contains a (static)
    object whose ctor has a non trivial global effect. In fact, this is more
    useful than it appears at first glance as many initializations are done
    like this.

    Finally, the solution: introduce an artificial symbol which we are going to
    reference from some object file which is definitely linked in (e.g. the one
    containing main()) and define it in the modules which we want to link in.
    This is achieved by the macros below.
 */
#ifndef _ITO33_LINK_H_
#define _ITO33_LINK_H_

/**
    Use this macro in a module which should be always linked in.

    @param module_name any valid C++ identifier
 */
#define ITO33_FORCE_LINK_THIS_MODULE(module_name)                             \
                int _ito33_link_dummy_func_##module_name ();                  \
                int _ito33_link_dummy_func_##module_name ()                   \
                {                                                             \
                    return 1;                                                 \
                }                                                             \
                /* just to force a semicolon after macro */ struct Dummy

/**
    Use this macro in the module containing main() function to enforce linking
    in the module which used ITO33_FORCE_LINK_THIS_MODULE() with the same name.

    Of course, this macro doesn't have to be used in the same module which
    contains main() but this is usually the simplest solution.

    @param module_name any valid C++ identifier
 */
#define ITO33_FORCE_LINK_MODULE(module_name)                                  \
                extern int _ito33_link_dummy_func_##module_name ();           \
                static int _ito33_link_dummy_var_##module_name =              \
                               _ito33_link_dummy_func_##module_name ()

#endif // _ITO33_LINK_H_

