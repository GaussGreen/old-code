//----------------------------------------------------------------------------
//
//   Group       : EDG Derivatives Research
//
//   Filename    : ClientRunnable.cpp
//
//   Description : Interface for client runnable methods via EdrAction
//
//   Author      : Andrew J Swain
//
//   Date        : 14 June 2001
//
//
//----------------------------------------------------------------------------

#include "edginc/config.hpp"
#include "edginc/ClientRunnable.hpp"
#include "edginc/Class.hpp"

DRLIB_BEGIN_NAMESPACE

ClientRunnable::ClientRunnable() {}

ClientRunnable::~ClientRunnable() {}

CClassConstSP const ClientRunnable::TYPE = CClass::registerInterfaceLoadMethod(
    "ClientRunnable", typeid(ClientRunnable), 0);

DRLIB_END_NAMESPACE
