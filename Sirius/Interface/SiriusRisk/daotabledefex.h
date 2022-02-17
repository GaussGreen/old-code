// daotabledefex.h: enhancement to CDaoTableDef to support a 
//			        default constructor with no arguments
//
//////////////////////////////////////////////////////////////////////

#ifndef _DAOTABLEDEFEX_H
#define _DAOTABLEDEFEX_H

#pragma once

class CDaoTableDefEx : public CDaoTableDef
{
public:
	CDaoTableDefEx(CDaoDatabase* pDatabase = NULL);	
	DECLARE_DYNAMIC(CDaoTableDefEx)
};

//{{AFX_INSERT_LOCATION}}

#endif
