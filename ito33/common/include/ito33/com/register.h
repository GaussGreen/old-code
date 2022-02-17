/////////////////////////////////////////////////////////////////////////////
// Name:        ito33/com/register.h
// Purpose:     Functions for registering and unregistering COM servers
// Author:      Vadim Zeitlin
// Created:     Jan 20, 2004 (extracted from dllmain_impl.cpp and regkey.cpp)
// RCS-ID:      $Id: register.h,v 1.2 2004/10/05 09:13:36 pedro Exp $
// Copyright:   (c) 2003-2004 Trilemma LLP
/////////////////////////////////////////////////////////////////////////////

/**
    @file   ito33/com/register.h
    @brief  Functions for registering and unregistering COM servers
 */

#ifndef _ITO33_OOM_REGISTER_H_
#define _ITO33_OOM_REGISTER_H_

namespace ito33
{

namespace COM
{

/**
    Flags for the server kind.

    The keys created in the registry are different for in and out of process
    servers so we must distinguish between them.
 */
enum ServerKind
{
  Server_InProcess,
  Server_OutOfProcess
};

/**
    Registers a COM server in the registry.

    Each COM server must have an entry (and even several of them) in the
    registry to be known to COM clients. This function creates the absolute
    minimum required.

    If an error occurs during registration, an exception is thrown.

    @param clsid the CLSID of the server as string
    @param progid the version independent progid of the server
    @param version version of the COM interface
    @param description the human readable description of the COM server
    @param kind the kind of the server we're registering
 */
void RegisterCOMServer(const std::string& clsid,
                       const std::string& progid,
                       int version,
                       const std::string& description,
                       ServerKind kind);

/**
    Unregisters an in-process COM server from the registry.

    This method removes the keys and values created by RegisterCOMServer()
    from the registry.

    If an error occurs during unregistration, an exception is thrown.

    @param clsid the CLSID of the server as string
    @param progid the version independent progid of the server
    @param version version of the COM interface
 */
void UnregisterCOMServer(const std::string& clsid,
                         const std::string& progid,
                         int version);

/**
    Register servers for all coclasses in CoClassInfo chain.

    Throws an exception if the registration fails.

    @sa CoClassInfo

    @param kind the kind of the server we're registering
 */
void RegisterAll(ServerKind kind);

/**
    Undo the work of RegisterAll().

    This function doesn't throw an exception as failure to unregister a server
    is usually not catastrophic.

    @return true if unregistered successfuly, false otherwise
 */
bool UnregisterAll();

} // namespace COM

} // namespace ito33

#endif // _ITO33_OOM_REGISTER_H_

