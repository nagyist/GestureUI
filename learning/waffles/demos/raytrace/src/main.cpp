// -------------------------------------------------------------
// The contents of this file may be distributed under the CC0
// license (http://creativecommons.org/publicdomain/zero/1.0/),
// or any compatible license, including (but not limited to) all
// OSI-approved licenses (http://www.opensource.org/licenses).
// -------------------------------------------------------------

#include <stdio.h>
#include <stdlib.h>
#include <GClasses/GApp.h>
#include <GClasses/GMacros.h>
#include "RayTrace.h"
#ifdef WIN32
#	include <direct.h>
#endif
#include <exception>

const char* g_szAppPath = NULL;

int main(int argc, char *argv[])
{
#ifdef WIN32
	if(_chdir("../../../bin/demos/rib") != 0)
#else
	if(chdir("../../../bin/demos/rib") != 0)
#endif
	{
	}
#ifdef WIN32
	if(_chdir("rib") != 0)
#else
	if(chdir("rib") != 0)
#endif
	{
	}

	int nRet = 0;
	try
	{
		RayTraceController c;
		c.RunModal();
	}
	catch(const std::exception& e)
	{
		fprintf(stderr, "%s\n", e.what());
		nRet = 1;
	}
	catch(...)
	{
		fprintf(stderr, "An exception of unknown type was thrown!");
		nRet = 1;
	}

	return nRet;
}

