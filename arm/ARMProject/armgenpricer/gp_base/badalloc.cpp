#include "gpbase/badalloc.h"

#include <new>
#include <new.h>


_PNH _old_new_handler;
bool hasBeenChanged=false;

/// bad_alloc handler
int bad_alloc_handler(size_t) {
  throw std::bad_alloc();
  return 0;
}

// Switch on the bad allocation handler with exception
void ARM_BadAlloc::SwitchOn()
{
	if (!hasBeenChanged)
	{
		_old_new_handler = _set_new_handler(bad_alloc_handler);
		hasBeenChanged = true;
	}
}

// Switch off the bad allocation handler with exception
void ARM_BadAlloc::SwitchOff()
{
	if (hasBeenChanged)
	{
		_set_new_handler(_old_new_handler);
		hasBeenChanged = false;
	}
}