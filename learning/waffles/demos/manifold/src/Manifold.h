// -------------------------------------------------------------
// The contents of this file may be distributed under the CC0
// license (http://creativecommons.org/publicdomain/zero/1.0/),
// or any compatible license, including (but not limited to) all
// OSI-approved licenses (http://www.opensource.org/licenses).
// -------------------------------------------------------------

#ifndef __MANIFOLD_H__
#define __MANIFOLD_H__

#include "Gui.h"


class ManifoldMenuController : public ControllerBase
{
protected:

public:
	ManifoldMenuController();
	virtual ~ManifoldMenuController();

	void RunModal();
};

#endif // __MANIFOLD_H__
