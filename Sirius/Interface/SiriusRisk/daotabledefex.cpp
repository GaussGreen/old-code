//	daotabledefex.cpp : implementation file
//
//////////////////////////////////////////////////////////////////////

#include "stdafx.h"
#include "daotabledefex.h"

#ifdef _DEBUG
#define new DEBUG_NEW
#undef THIS_FILE
static char THIS_FILE[] = __FILE__;
#endif

IMPLEMENT_DYNAMIC(CDaoTableDefEx, CDaoTableDef)

CDaoTableDefEx::CDaoTableDefEx(CDaoDatabase* pdb) : CDaoTableDef(pdb)
{	
}
