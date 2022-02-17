/////////////////////////////////////////////////////////////////////////////
// Name:        src/XL/sheetobject.cpp
// Purpose:     implementation of SheetObjectManager class
// Author:      Vadim Zeitlin
// Created:     2006-04-14
// RCS-ID:      $Id: sheetobject.cpp,v 1.6 2006/08/21 15:17:54 willy Exp $
// Copyright:   (c) 2006 Trilemma LLP
/////////////////////////////////////////////////////////////////////////////

// ============================================================================
// declarations
// ============================================================================

// ----------------------------------------------------------------------------
// headers
// ----------------------------------------------------------------------------

#include "ito33/XL/oper.h"
#include "ito33/XL/call.h"
#include "ito33/XL/sheetobject.h"

namespace ito33
{

namespace XL
{

// ============================================================================
// SheetObjectManagerBase implementation
// ============================================================================

void *SheetObjectManagerBase::DoGetForThisCell()
{
  Oper xloRef;
  Call(xlfCaller, ByRef(xloRef));
    
  unsigned x, y;
  if ( !xloRef.GetAsSRef(x, y) )
  {
    throw EXCEPTION_MSG(xlerrRef, "No caller cell");
  }

  // check if we already have an object for this cell
  const unsigned long key = MAKELONG(x, y);

  Oper xloSheetRef;
  
  // At this time there is no mean to know what sheet contains the caller cell.
  // We only know the active sheet. The following call to Excel4v works only
  // if the function are in the active sheet, this is why it doesn't work on 
  // Excel recalculation.
  Call(xlSheetId,ByRef(xloSheetRef));

  unsigned long SheetId;
  
  if ( !xloSheetRef.GetSheetId( SheetId ) )
  {
    throw EXCEPTION_MSG(xlerrRef, "Failed to retrieve the active sheet id");
  }  
  
  const SheetsObjects::iterator iterObjects = m_sheetsobjects.find(SheetId);
  
  void *ptr;

  if ( iterObjects != m_sheetsobjects.end() )
  {
    Objects& objects = iterObjects->second;
    const Objects::iterator i = objects.lower_bound(key);
    if ( i == objects.end() || i->first != key )
    {
      ptr = DoCreate();

      objects.insert(i, Objects::value_type(key, ptr));
    }
    else // we have an existing entry for this cell
    {
      ptr = objects[key];
    }
  }
  else
  {
    Objects obj;
    ptr = DoCreate();
    obj.insert(obj.end() ,Objects::value_type(key, ptr ));
    m_sheetsobjects.insert( m_sheetsobjects.end() , SheetsObjects::value_type( SheetId,obj ));
  }

  // and also remember or update the current cell formula
  Oper xloFormula;
  Call(xlfGetFormula, ByRef(xloFormula), Oper(x, y));
  std::string formula;
  xloFormula.As(&formula);
  m_formulae[ptr] = formula;

  return ptr;
}

void SheetObjectManagerBase::PruneInactive()
{
  std::string formula;

  Oper xloFormula;
  ByRef rc(xloFormula);

  const SheetsObjects::iterator sheetsobjend = m_sheetsobjects.end();
  
  for ( SheetsObjects::iterator i= m_sheetsobjects.begin();i != sheetsobjend;++i )
  {
    const Objects::iterator objend = i->second.end();
    bool isSheetEmpty = true;
    for ( Objects::iterator j = i->second.begin(); j != objend; )
    {
      isSheetEmpty = false;
      const unsigned long coord = j->first;
      const Formulae::iterator k = m_formulae.find(j->second);
      ASSERT_MSG( k != m_formulae.end(), "object without corresponding formula?" );

      try
      {
        Ref ref(static_cast<WORD> (HIWORD(coord)),
                static_cast<BYTE> ( LOWORD(coord) ) );

        ref.SetSheetId(i->first);
        
        Call(xlfGetFormula, rc, Oper(ref) );

        // we consider that if the formula changed, then the object is not needed
        // any more -- this doesn't always work, of course, but it's very simple
        // and works often enough
        if ( xloFormula.As(&formula) && formula == k->second )
        {
          // skip the code destroying the object below
          ++j;
          continue;
        }
      }
      catch ( XL::Exception& )
      {
        // maybe the cell doesn't exist any more? in this case we surely don't
        // need the associated object neither
      }

      m_formulae.erase(k);
      DoDelete(j->second);
      j = i->second.erase(j);
    }
    //Remove the element from the map if empty 
    if ( isSheetEmpty )
      i = m_sheetsobjects.erase(i);
  }
}

void *SheetObjectManagerBase::DoFromHandle(HandleArg handle) const
{
  // the stirng must be in the form "#<ptr>"
  unsigned long num = 0;
  const char *p = handle.c_str();
  if ( *p != '#' || !String::ToULong(std::string(p + 1), &num, 16) )
    return NULL;

  void *ptr = reinterpret_cast<void *>(static_cast<intptr_t>(num));
  return DoIsValid(ptr) ? ptr : NULL;
}

void *SheetObjectManagerBase::DoFromHandleOrThrow(HandleArg handle) const
{
  void *p = DoFromHandle(handle);
  if ( !p )
  {
    throw EXCEPTION_MSG(xlerrValue, "Invalid handle");
  }

  return p;
}

XLOPER *SheetObjectManagerBase::DoHandleOf(const void *p) const
{
  // we currently don't check if the object is valid but we might in the future
  return Oper::Return(String::Printf("#%08x", reinterpret_cast<intptr_t>(p)));
}

void SheetObjectManagerBase::DoClear()
{
  const SheetsObjects::iterator sheetsobjend = m_sheetsobjects.end();

  for ( SheetsObjects::iterator i= m_sheetsobjects.begin();i != sheetsobjend;++i )
  {
    const Objects::iterator objend = i->second.end();
    for ( Objects::iterator j = i->second.begin(); j != objend; ++j )
    {
      DoDelete(j->second);
    }
  }

  // no need to empty m_sheetsobjects nor touch m_formulae as we're called from dtor
  // anyhow
}

} // namespace XL

} // namespace ito33
