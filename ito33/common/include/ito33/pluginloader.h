/////////////////////////////////////////////////////////////////////////////
// Name:        ito33/pluginloader.h
// Purpose:     Plugins interface
// Author:      Vaclav Slavik
// Created:     2005-01-06
// RCS-ID:      $Id: pluginloader.h,v 1.10 2006/05/14 13:37:27 zeitlin Exp $
// Copyright:   (c) 2004 Trilemma LLP
/////////////////////////////////////////////////////////////////////////////

#ifndef _ITO33_PLUGINLOADER_H_
#define _ITO33_PLUGINLOADER_H_

/**
    @file   ito33/pluginloader.h
    @brief  Simple plugins-loading system.
 */

#include "ito33/string.h"
#include "ito33/dynlib.h"

#include "ito33/plugin.h"

/**
    Base class for dynamic plugins loader, PluginLoader<T>.
 */
class PluginLoaderBase
{
public:
  /**
      Adds a path to load plugins from.

      This path (which can be a semicolon or colon-separated list of paths too)
      is inserted in the beginning of the plugins load path.

      @param path a non-NULL and non-empty path string
   */
  static void AddSearchPath(const char *path);

protected:
  /**
      Type of factory function in the DLL.

      Factory function's return value is either pointer to newly created Plugin
      object or NULL in case of failure.

      First argument is pointer to std::string variable; this variable is
      filled with description of the error if the factory fails. It always
      points to valid std::string object and is never NULL.

      Second argument is the requested version of plugin interface.

      Third argument is a pointer to structure with additional parameters
      needed by the factory. It is reserved for future use and not used
      currently.
   */
  typedef Plugin* (*FactoryFuncType)(std::string*, int, void*);


  /**
      Loads the plugin and returns created instance; throws on error.

      @param prefix   DLL name prefix.
      @param name     Name of the plugin to load.
      @param category Plugin's category.
      @param apiver   Version of the API.
      @param handle   if non NULL, filled with the handle of plugin DLL

      @return Instance of the plugin.
   */
  static Plugin *DoLoad(const char *prefix,
                        const char *name,
                        const char *category,
                        int apiver,
                        DllHandle_t *handle = NULL);

  /**
      Returns the full plugin file name.

      This is still just the base name, not the full file name.

      @param prefix   DLL name prefix.
      @param name     Name of the plugin to load.
   */
  static std::string GetPluginFileName(const char *prefix, const char *name);
};

/**
    Template class for dynamic loading of plugins that implement interface T.

    Type T must derive from Plugin and must implement the following static
    method:

      @code
      static const char *GetPluginCategory();
      @endcode

    T::GetPluginCategory returns category of the plugin (e.g. "DatastoreReader")
    and the same value must be used as argument to PLUGIN_DEFINE_ENTRY_POINT
    for plugins derived from T.

    See the Load method for more details.

    @see PluginLoader::Load, PLUGIN_DEFINE_ENTRY_POINT
 */
template<class T> class PluginLoader : public PluginLoaderBase
{
public:
  /**
      Loads plugin with given name from DLL. The name of the DLL from which
      to load the plugin is derived from @a prefix and @a name arguments and
      from the return value of T::GetPluginCategory function. The DLL must be
      in the same directory as application's executable and its name must
      have the following form: plg_&lt;prefix>&lt;name>.dll. @a Name is used to
      identify particular plugin within a category of plugins, @a prefix can
      be used to distinguish DLLs for plugins of different categories.

      The DLL must export factory function that instantiates the plugin; this
      function takes three arguments and returns Plugin pointer (see
      PluginLoaderBase::FactoryFuncType). The factory function must be
      undecorated (extern "C") and named according to the following scheme:

        ITO33_PluginFactory_<Category>

      The factory function takes three arguments:

        std::string *error         - pointer to variable where error message
                                     will be stored in case of error
                                     (never NULL)
        int version                - requested interface version
        void *parameters           - additional parameters for the factory
                                     (unused)

      where 'version' is version of API for this category. The API version is
      same as passed to PLUGIN_DEFINE_ENTRY_POINT and as the value of
      T::PLUGIN_API_VERSION constant, which must be defined.
      It is recommended to use convenience macro PLUGIN_DEFINE_ENTRY_POINT to
      define the factory in plugin DLL.

      On success, the factory function creates new instance of the plugin and
      returns it. On error, it returns NULL and sets the string pointed to by
      'error' to explanatory error message.

      One DLL can contain factories for more than one category.

      The DLL is unloaded only when program terminates unless you retrieve its
      handle using the last parameter of this method and unload it manually
      (this can't be done automatically, unfortunately, as the only moment when
      it would be possible is in plugin dtor, but unloading the DLL from there
      would mean that we return to already unloaded code -- and crash, of course)

      Usage:
        @code
        AutoPtr<MyPlugin> p(PluginLoader<MyPlugin>::Load("my", name));
        @endcode

      @param  prefix  'prefix' part of DLL name. May be empty.
      @param  name    Name of the plugin to load, same as the 'name' part
                      of DLL name.
      @param  handle  if non NULL, filled with the handle of plugin DLL
                      to be unloaded by the caller (otherwise DLL is only
                      unloaded at program shutdown)

      @return Loaded instance of the plugin, to be freed by caller. It is never
              NULL as this function throws on error.
   */
  static T *Load(const std::string& prefix,
                 const std::string& name,
                 DllHandle_t *handle = NULL)
  {
    return static_cast<T*>(DoLoad(prefix.c_str(),
                                  name.c_str(),
                                  T::GetPluginCategory(),
                                  T::PLUGIN_API_VERSION,
                                  handle));
  }
};

#endif // _ITO33_PLUGINLOADER_H_
