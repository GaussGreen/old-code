/////////////////////////////////////////////////////////////////////////////
// Name:        utils/gettext.cpp
// Purpose:     implementation of GetTranslation() function
// Author:      Vadim Zeitlin
// Created:     18.12.02
// RCS-ID:      $Id: gettext.cpp,v 1.11 2006/05/27 20:05:06 zhang Exp $
// Copyright:   (c) 2002 Trilemma LLP
/////////////////////////////////////////////////////////////////////////////

#if ( defined(_MSC_VER) && _MSC_VER >= 1400 )
  // 'fopen' was declared deprecated in VC8. But we still can use it in this
  // file
  #define _CRT_SECURE_NO_DEPRECATE 1
#endif

// ============================================================================
// declarations
// ============================================================================

// ----------------------------------------------------------------------------
// headers
// ----------------------------------------------------------------------------

#include "ito33/common.h"
#include "ito33/debug.h"

#include "ito33/gettext.h"
#include "ito33/thread.h"

#include <stdio.h>

#include "ito33/beforestd.h"

#include <string>
#include <map>

#include "ito33/afterstd.h"

#ifdef _WIN32
  // we need this for GetModuleFileName()
  #define WIN32_LEAN_AND_MEAN
  #include <windows.h>
#endif

using namespace ito33;

// ----------------------------------------------------------------------------
// types
// ----------------------------------------------------------------------------

// we need to have types having exactly 8 and 32 bits respectively to
// correspond to the on disk format
typedef unsigned char size_t8;
typedef unsigned int size_t32;

ASSERT_BITSIZE( size_t8, 8 );
ASSERT_BITSIZE( size_t32, 32 );

// this maps has original strings as keys and translated strings as values
typedef std::map<std::string, std::string> TranslationMap;

// ----------------------------------------------------------------------------
// File: small helper class which automatically closes the file in its dtor
// ----------------------------------------------------------------------------

class File
{
public:
  File(FILE *fp) : m_fp(fp) { }
  ~File() { if ( m_fp ) fclose(m_fp); }

  operator FILE *() const { return m_fp; }

private:
  // not copyable
  File(const File&);
  File& operator=(const File&);

  FILE *m_fp;
};

// ----------------------------------------------------------------------------
// MessageCatalog: represents a single (binary) message catalog
// ----------------------------------------------------------------------------

class MessageCatalog
{
public:
  // trivial ctor & dtor
  MessageCatalog();
  ~MessageCatalog();

  // load the catalog from the given disk file, returns true if loaded ok
  bool Load(const char *filename);

  // add the translation from this catalog to the given translation map
  void AddTranslations(TranslationMap& map);

private:
  // this implementation is binary compatible with GNU gettext() version 0.10

  // magic number identifying the .mo format file
  static const size_t32 MSGCATALOG_MAGIC;
  static const size_t32 MSGCATALOG_MAGIC_SW;

  // an entry in the string table
  struct MsgTableEntry
  {
   size_t32   nLen;           // length of the string
   size_t32   ofsString;      // pointer to the string
  };

  // header of a .mo file
  struct MsgCatalogHeader
  {
   size_t32  magic,          // offset +00:  magic id
        revision,       //        +04:  revision
        numStrings;     //        +08:  number of strings in the file
   size_t32  ofsOrigTable,   //        +0C:  start of original string table
        ofsTransTable;  //        +10:  start of translated string table
   size_t32  nHashSize,      //        +14:  hash table size
        ofsHashTable;   //        +18:  offset of hash table start
  };

  // all data is stored here, NULL if no data loaded
  size_t8 *m_pData;

  // data description
  size_t32          m_numStrings;   // number of strings in this domain
  MsgTableEntry    *m_pOrigTable,   // pointer to original   strings
          *m_pTransTable;  //            translated

  // true if need to change endianness on the fly
  bool m_bSwapped;


  // swaps between big and little endian if needed
  size_t32 Swap(size_t32 ui) const
  {
    return m_bSwapped ? (ui << 24) | ((ui & 0xff00) << 8) |
              ((ui >> 8) & 0xff00) | (ui >> 24)
             : ui;
  }

  // return the pointer to the string with the given index
  const char *StringAtOfs(MsgTableEntry *pTable, size_t32 index) const
   { return (const char *)(m_pData + Swap(pTable[index].ofsString)); }
};

// ----------------------------------------------------------------------------
// global data
// ----------------------------------------------------------------------------

static struct TransData
{
  // this critical section protects the entire g_transData
  CriticalSection cs;

  // all loaded translations
  TranslationMap map;

  // if true, try to load translations on demand
  bool tryLoad;

  // ctor is only there to initialize tryLoad to true
  TransData() { tryLoad = true; }

private:
  TransData(const TransData&);
  TransData& operator=(const TransData&);
} g_transData;

// ============================================================================
// MessageCatalog implementation
// ============================================================================

const size_t32 MessageCatalog::MSGCATALOG_MAGIC    = 0x950412de;
const size_t32 MessageCatalog::MSGCATALOG_MAGIC_SW = 0xde120495;

MessageCatalog::MessageCatalog()
{
  m_numStrings = 0;
  m_pData = NULL;
}

MessageCatalog::~MessageCatalog()
{
  delete [] m_pData;
}

bool MessageCatalog::Load(const char *filename)
{
  // open the file
  File f(fopen(filename, "r"));
  if ( !f )
    return false;

  // get the file size
  long len;
  if ( fseek(f, 0, SEEK_END) == -1 || (len = ftell(f)) == -1 )
  {
    // can't seek on the file?
    return false;
  }

  // slurp it in memory
  const size_t nLen = (size_t)len;
  m_pData = new size_t8[nLen];

  rewind(f);
  if ( fread(m_pData, sizeof(size_t8), nLen, f) != nLen )
  {
    // read error
    return false;
  }

  // examine its header
  bool bValid = nLen > sizeof(MsgCatalogHeader);

  MsgCatalogHeader *pHeader = reinterpret_cast<MsgCatalogHeader *>(m_pData);
  if ( bValid )
  {
    // we'll have to swap all the integers if it's true
    m_bSwapped = pHeader->magic == MSGCATALOG_MAGIC_SW;

    // check the magic number
    bValid = m_bSwapped || pHeader->magic == MSGCATALOG_MAGIC;
  }

  if ( !bValid )
  {
    // it's either too short or has incorrect magic number
    return false;
  }

  // initialize the tables
  m_numStrings  = Swap(pHeader->numStrings);
  m_pOrigTable  = (MsgTableEntry *)(m_pData + Swap(pHeader->ofsOrigTable));
  m_pTransTable = (MsgTableEntry *)(m_pData + Swap(pHeader->ofsTransTable));

  // and return success
  return true;
}

void MessageCatalog::AddTranslations(TranslationMap& map)
{
  for ( size_t32 i = 0; i < m_numStrings; i++ )
  {
    std::string key(StringAtOfs(m_pOrigTable, i));
    map[key] = StringAtOfs(m_pTransTable, i);
  }
}

// ============================================================================
// our public API
// ============================================================================

const char *ito33::GetTranslation(const char *str)
{
  if ( !str )
    return NULL;

  // serialize all accesses to the global data
  Lock<CriticalSection> lock(g_transData.cs);

  // if the message catalog hadn't been loaded yet, do it now
  if ( g_transData.tryLoad )
  {
    g_transData.tryLoad = false;

    // what is the name of the catalog we should load? normally it is just
    // the name of the program but we don't have it here except under
    // Windows where we can recover it
#ifdef _WIN32
    char filename[4096];

    // -4 because we want to have enough space to append ".mo"
    if ( ::GetModuleFileName(NULL, filename, SIZEOF(filename) - 4) )
    {
      // find the extension and replace it with ".mo"
      char *pDot = strrchr(filename, '.');
      if ( !pDot )
        pDot = filename + strlen(filename);

      // TODO: add the language code here to allow using different
      //       languages

      *pDot++ = '.';
      *pDot++ = 'm';
      *pDot++ = 'o';
      *pDot++ = '\0';

      MessageCatalog catalog;
      if ( catalog.Load(filename) )
      {
        catalog.AddTranslations(g_transData.map);
      }
      //else: loading msg catalog failed
    }
    //else: GetModuleFileName() failed?
#else
    //#error "Need argv[0] for on demand message catalog loading!"
#endif
  }

  // try to find the translation
  TranslationMap::const_iterator it = g_transData.map.find(str);

  // as the strings in the translation map are immutable, it is ok to return
  // the pointer to the string data directly
  return it == g_transData.map.end() ? str : it->second.c_str();
}


