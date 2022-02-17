// DRSymbolTable.h: interface for the DRSymbolTable class.
//
//////////////////////////////////////////////////////////////////////

#if !defined(AFX_DRSYMBOLTABLE_H__305831F8_6DBC_11D2_97ED_00C04FD8EB9A__INCLUDED_)
/**@#-*/
#define AFX_DRSYMBOLTABLE_H__305831F8_6DBC_11D2_97ED_00C04FD8EB9A__INCLUDED_
/**@#+*/

#include "kplatdep.h"
#include <map>
#include "kstring.h"

/** Simple symbol table 

  Stores a set of (symbol, value) pairs where both the symbol and value are strings
  */
class KSymbolTable  : public KMap(string, string)
{
public:
	/** Construct an empty table */
	KSymbolTable() {}
	/**@#-*/
	virtual ~KSymbolTable() {}
	/**@#+*/

	/** Gets a symbol from the table.  Replaces all instances of $SymbolName with
		the corresponding value.
		*/
	string get(string);
};

/** Global symbol table */
extern KSymbolTable theSymbolTable;

/**@#-*/
#define PATH_VARIABLES theSymbolTable;
/**@#+*/

#endif // !defined(AFX_DRSYMBOLTABLE_H__305831F8_6DBC_11D2_97ED_00C04FD8EB9A__INCLUDED_)
