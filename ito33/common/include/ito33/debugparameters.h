/////////////////////////////////////////////////////////////////////////////
// Name:        ito33/debugparameters.h
// Purpose:     Struct holding some debug parameters
// Created:     2005/02/01
// RCS-ID:      $Id: debugparameters.h,v 1.3 2006/05/27 20:13:16 zhang Exp $
// Copyright:   (c) 2005 -  Trilemma LLP
/////////////////////////////////////////////////////////////////////////////

/**
    @file   ito33/debugparameters.h
    @brief  Struct holding some debug parameters

    The environment variable ITO33InstDir will be used if defined. Otherwise,
    a default directory is used, which is respectively c:\\ito33 for win32
    and /usr/local/ito33 for unix.
 */

#ifndef _ITO33_DEBUGPARAMETERS_H_
#define _ITO33_DEBUGPARAMETERS_H_

#include "ito33/string.h"

namespace ito33
{


/**
   The debug output parameters

   @noexport
 */
struct DebugParameters
{
  /// Default ctor disables debug
  DebugParameters() { isEnabled = false; }

  /**
     Gets the output file.

     @param defaultName A default name will be used if user doesn't specify one
   */
  std::string GetDebugOutputFile(const char* defaultName) const;
  
  /// Boolean indicates if the debug output is enabled
  bool isEnabled;

  /// The debug output name given by the user
  std::string filename;
  
  /**
      Gets the directory for outputting official debug output files.
      Actually it is "$(ITO33InstDir)/output/" (or "$(ITO33InstDir)\output\"
      for windows) if environment variable 'ITO33InstDir' is defined. Otherwise,
      an exception will be thrown.
   */
  static std::string GetDebugOutputDir();

  /**
      Gets the directory for the developpers to output their debug files.

      It is the value of the environment variable "ITO33DebugDir" if defined.
      Otherwise, it is the value returned by GetDebugOutputDir().

      The return value ends with the file name seperator according to the
      system.
   */
  static std::string GetDebugDir();
};


} // namespace ito33

#endif // #ifndef _ITO33_DEBUGPARAMETERS_H_
