/*
	Copyright (C) 2010, Mike Gashler

	This library is free software; you can redistribute it and/or
	modify it under the terms of the GNU Lesser General Public
	License as published by the Free Software Foundation; either
	version 2.1 of the License, or (at your option) any later version.

	see http://www.gnu.org/copyleft/lesser.html
*/

#include "GActivation.h"
#include "GMath.h"
#include "GTwt.h"

namespace GClasses {

GTwtNode* GActivationFunction::toTwt(GTwtDoc* pDoc)
{
	return pDoc->newString(name());
}

// static
GActivationFunction* GActivationFunction::fromTwt(GTwtNode* pNode)
{
	const char* szName = pNode->asString();
	if(strcmp(szName, "logistic") == 0)
		return new GActivationLogistic();
	else if(strcmp(szName, "arctan") == 0)
		return new GActivationArcTan();
	else if(strcmp(szName, "tanh") == 0)
		return new GActivationTanH();
	else if(strcmp(szName, "algebraic") == 0)
		return new GActivationAlgebraic();
	else if(strcmp(szName, "identity") == 0)
		return new GActivationIdentity();
	else if(strcmp(szName, "bend") == 0)
		return new GActivationBend();
	else if(strcmp(szName, "piecewise") == 0)
		return new GActivationPiecewise();
	else
		ThrowError("Unrecognized activation function: ", szName);
	return NULL;
}




// virtual
double GActivationArcTan::halfRange()
{
	return M_PI / 2;
}




// virtual
double GActivationPiecewise::squash(double x)
{
	double d = floor(log(MAX(1.0, ABS(x))) * 2.46630346237643166848294787835207968);
	double a = pow(1.5, d);
	double b = pow(1.5, d + 1);
	if(b >= 700.0)
		return (x >= 0 ? 1.0 : 0.0);
	double t = (ABS(x) - a) / (b - a);
	double v = (1.0 - t) / (exp(-a) + 1.0) + t / (exp(-b) + 1.0);
	return (x >= 0 ? v : 1.0 - v);
}




} // namespace GClasses

