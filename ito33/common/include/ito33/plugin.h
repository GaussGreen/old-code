/////////////////////////////////////////////////////////////////////////////
// Name:        ito33/plugin.h
// Purpose:     Plugins interface
// Author:      Vaclav Slavik
// Created:     2005-01-06
// RCS-ID:      $Id: plugin.h,v 1.10 2006/07/18 10:42:32 zeitlin Exp $
// Copyright:   (c) 2004 Trilemma LLP
/////////////////////////////////////////////////////////////////////////////

#ifndef _ITO33_PLUGIN_H_
#define _ITO33_PLUGIN_H_

#include <sstream>

/**
    @file   ito33/plugin.h
    @brief  Simple plugins-loading system.
 */

/**
    Base class for dynamically loadable plugins. Any plugin loaded using
    PluginLoader<T> class must derive from Plugin. A plugin is implemented
    in DLL and made visible using PLUGIN_DEFINE_ENTRY_POINT.

    @see PluginLoader, PLUGIN_DEFINE_ENTRY_POINT
 */
class Plugin
{
public:
  virtual ~Plugin() {}
};

#ifdef _WIN32
  #define PLUGIN_DLLEXPORT __declspec(dllexport)
#else // !_WIN32
  #define PLUGIN_DLLEXPORT
#endif


// helper macros:
#define __PLUGIN_CONCAT(a,b)    __PLUGIN_CONCAT0(a,b)
#define __PLUGIN_CONCAT0(a,b)   a ## _ ## b

/**
    Expands to name of factory function for given @a category, as literal.
  */
#define PLUGIN_FACTORY_FUNC_NAME(category)  \
  __PLUGIN_CONCAT(ITO33_PluginFactory, category)

/// Same as PLUGIN_FACTORY_FUNC_NAME, but expands into a string
#define PLUGIN_FACTORY_FUNC_NAME_STR(category) \
  "ITO33_PluginFactory_" #category


/**
    Expands to C++ code of function's signature. Can be used to either
    declare plugin factory function or at the beginning of its definition.

    See PluginLoader documentation for more information on this function's
    signature and return value.

    For example:
    @code

    // declaration:
    PLUGIN_FACTORY_FUNC_SIGNATURE(DatastoreReader, 1);

    // definition:
    PLUGIN_FACTORY_FUNC_SIGNATURE(DatastoreReader, 1)
    {
      if ( version != 1 )
      {
        *error = "unsupported version";
        return NULL;
      }
      return new MyDatastoreReader;
    }
    @endcode
 */
#define PLUGIN_FACTORY_FUNC_SIGNATURE(category)                           \
  extern "C" PLUGIN_DLLEXPORT                                             \
  Plugin* PLUGIN_FACTORY_FUNC_NAME(category)                              \
          (                                                               \
            /* if an error occurs, its description is stored there */     \
            std::string *error,                                           \
            /* requested interface version */                             \
            int version,                                                  \
            /* further arguments needed for creation or NULL (unused) */  \
            void *parameters                                              \
          )

/**
    Defines plugin factory function for given plugins category. This macro
    should be used in plugin DLL to define entry point(s) neccessary for
    PluginLoader::Load.

    @param category   Plugin category. The meaning is same as in PluginLoader
                      and GetPluginCategory. Note that this argument is _not_
                      a string!
    @param apiversion Version of the API for @a category plugins. Plugins with
                      version different from the one requested can't be loaded.
    @param create     Code snippet used to create a new instance. This may be
                      a function call, a simple new statement, or anything else.

    Example:
      @code
      PLUGIN_DEFINE_ENTRY_POINT(MyPlugin, new MyPlugin())
      @endcode

    @see PluginLoader::Load
 */
#define PLUGIN_DEFINE_ENTRY_POINT(category, apiversion, create)           \
  PLUGIN_FACTORY_FUNC_SIGNATURE(category)                                 \
  {                                                                       \
    if ( apiversion != version )                                          \
    {                                                                     \
      if ( error )                                                        \
      {                                                                   \
        std::stringstream ss;                                             \
        ss << "unsupported version " << version << " requested "          \
           << "(this plugin is version " << apiversion << ")";            \
        *error = ss.str();                                                \
      }                                                                   \
      return NULL;                                                        \
    }                                                                     \
    if ( parameters != NULL )                                             \
    {                                                                     \
      if ( error )                                                        \
        *error = "unexpected parameters passed to the plugin";            \
      return NULL;                                                        \
    }                                                                     \
    return create;                                                        \
  }

#endif // _ITO33_PLUGIN_H_
