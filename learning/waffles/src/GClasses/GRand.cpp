/*
	Copyright (C) 2006, Mike Gashler

	This library is free software; you can redistribute it and/or
	modify it under the terms of the GNU Lesser General Public
	License as published by the Free Software Foundation; either
	version 2.1 of the License, or (at your option) any later version.

	see http://www.gnu.org/copyleft/lesser.html
*/

#include <string.h>
#include "GRand.h"
#include <math.h>
#include "GMacros.h"
#include <stdlib.h>
#include "sha1.h"
#include "sha2.h"
#include "GHistogram.h"
#include "GTime.h"
#include "GMath.h"
#include "GVec.h"

namespace GClasses {

using std::vector;

COMPILER_ASSERT(sizeof(uint64) == 8);

GRand::GRand(uint64 seed)
{
	setSeed(seed);
}

GRand::~GRand()
{
}

void GRand::setSeed(uint64 seed)
{
	m_b = 0xCA535ACA9535ACB2ull + seed;
	m_a = 0x6CCF6660A66C35E7ull + (seed << 24);
}

uint64 GRand::next(uint64 range)
{
	// Use rejection to find a random value in a range that is a multiple of "range"
	uint64 n = (0xffffffffffffffffull % range) + 1;
	uint64 x;
	do
	{
		x = next();
	} while((x + n) < n);

	// Use modulus to return the final value
	return x % range;
}

double GRand::uniform()
{
	// use 52 random bits for the mantissa
	return (double)(next() & 0xfffffffffffffull) / 4503599627370496.0;
}

double GRand::normal()
{
	double x, y, mag;
	do
	{
		x = uniform() * 2 - 1;
		y = uniform() * 2 - 1;
		mag = x * x + y * y;
	} while(mag >= 1.0 || mag == 0);
	return y * sqrt(-2.0 * log(mag) / mag); // the Box-Muller transform	
}

int GRand::categorical(vector<double>& probabilities)
{
	double d = uniform();
	int i = 0;
	for(vector<double>::iterator it = probabilities.begin(); it != probabilities.end(); it++)
	{
		d -= *it;
		if(d < 0)
			return i;
		i++;
	}
	GAssert(false); // the probabilities are not normalized
	return probabilities.size() - 1;
}

double GRand::exponential()
{
	return -log(uniform());
}

double GRand::cauchy()
{
	return normal() / normal();
}

int GRand::poisson(double mu)
{
	if(mu <= 0)
		ThrowError("invalid parameter");
	double p = 1.0;
	int n = 0;
	if(mu < 30)
	{
		mu = exp(-mu);
		do {
			p *= uniform();
			n++;
		} while(p >= mu);
		return n - 1;
	}
	else
	{
		double u1, u2, x, y;
		double c = 0.767-3.36 / mu;
		double b = M_PI / sqrt(3.0 * mu);
		double a = b * mu;
		if(c <= 0)
			ThrowError("Error generating Poisson deviate");
		double k = log(c) - mu - log(b);
		double ck1 = 0.0;
		double ck2;
		do {
			ck2=0.;
			do {
				u1 = uniform();
				x = (a - log(0.1e-18 + (1.0 - u1) / u1)) / b;
				if(x > -0.5)
					ck2=1.0;
			} while (ck2<0.5);
			n = (int)(x + 0.5);
			u2 = uniform();
			y = 1 + exp(a - b * x);
			ck1 = a - b * x + log(.1e-18 + u2/(y * y));
#ifdef WIN32
			ck2 = k + n * log(.1e-18 + mu) - GMath::logGamma(n + 1.0);
#else
			ck2 = k + n * log(.1e-18 + mu) - lgamma(n + 1.0);
#endif
			if(ck1 <= ck2)
				ck1 = 1.0;
		} while (ck1 < 0.5);
		return n;
	}
}

double GRand::gamma(double alpha)
{
	double x;
	if(alpha <= 0)
		ThrowError("invalid parameter");
	if(alpha == 1)
		return exponential();
	else if(alpha < 1)
	{
		double aa = (alpha + M_E) / M_E;
		double r1, r2;
		do {
			r1 = uniform();
			r2 = uniform();
			if(r1 > 1.0 / aa)
			{
				x = -log(aa * (1.0 - r1) / alpha);
				if(r2 < pow(x, (alpha - 1.0)))
					return x;
			}
			else
			{
				x = pow((aa * r1), (1.0 / alpha));
				if(r2 < exp(-x))
					return x;
			}
		} while(r2 < 2);
	}
	else
	{
		double c1 = alpha-1;
		double c2 = (alpha - 1.0 / (6.0 * alpha)) / c1;
		double c3 = 2.0 / c1;
		double c4 = c3 + 2.0;
		double c5 = 1.0 / sqrt(alpha);
		double r1, r2;
		do {
			do {
				r1 = uniform();
				r2 = uniform();
				if(alpha > 2.5)
					r1 = r2 + c5 * (1.0 - 1.86 * r1);
			} while(r1 <= 0 || r1 >= 1);
			double w = c2 * r2 / r1;
			if((c3 * r1) + w + (1.0 / w) <= c4)
				return c1 * w;
			if((c3 * log(r1)) - log(w) + w < 1)
				return c1 * w;
		} while(r2 < 2);
	}
	ThrowError("Error making random gamma");
	return 0;
}

double GRand::chiSquare(double t)
{
	return gamma(t / 2.0) * 2.0;
}

int GRand::binomial(int n, double p)
{
	int c = 0;
	for(int i = 0; i < n; i++)
	{
		if(uniform() < p)
			c++;
	}
	return c;
}

double GRand::softImpulse(double s)
{
	double y = uniform();
	return 1.0 / (1.0 + pow(1.0 / y - 1.0, 1.0 / s));
}

double GRand::weibull(double gamma)
{
	if(gamma <= 0)
		ThrowError("invalid parameter");
	return pow(exponential(), (1.0 / gamma));
}

void GRand::dirichlet(double* pOutVec, const double* pParams, int dims)
{
	double* pOut = pOutVec;
	const double* pIn = pParams;
	for(int i = 0; i < dims; i++)
		*(pOut++) = gamma(*(pIn++));
	GVec::sumToOne(pOutVec, dims);
}

double GRand::student(double t)
{
	if(t <= 0)
		ThrowError("invalid parameter");
	return normal() / sqrt(chiSquare(t) / t);
}

int GRand::geometric(double p)
{
	if(p < 0 || p > 1)
		ThrowError("invalid parameter");
	return (int)floor(-exponential() / log(1.0 - p));
}

double GRand::f(double t, double u)
{
	if(t <= 0 || u <= 0)
		ThrowError("invalid parameters");
	return chiSquare(t) * u / (t * chiSquare(u));
}

double GRand::logistic()
{
	double y = uniform();
	return log(y) - log(1.0 - y);
}

double GRand::logNormal(double mean, double dev)
{
	return exp(normal() * dev + mean);
}

double GRand::beta(double alpha, double beta)
{
	if(alpha <= 0 || beta <= 0)
		ThrowError("invalid parameters");
	double r = gamma(alpha);
	return r / (r + gamma(beta));
}

void GRand::spherical(double* pOutVec, int dims)
{
	double* pEl = pOutVec;
	for(int i = 0; i < dims; i++)
		*(pEl++) = normal();
	GVec::safeNormalize(pOutVec, dims, this);
}

void GRand::spherical_volume(double* pOutVec, int dims)
{
	spherical(pOutVec, dims);
	GVec::multiply(pOutVec, pow(uniform(), 1.0 / dims), dims);
}

void GRand::cubical(double* pOutVec, int dims)
{
	double* pEl = pOutVec;
	for(int i = 0; i < dims; i++)
		*(pEl++) = uniform();
}

#ifndef NO_TEST_CODE
#define TEST_BIT_HIST_ITERS 100000
void GRand_testBitHistogram()
{
	GRand prng(0);
	size_t counts[64];
	for(size_t i = 0; i < 64; i++)
		counts[i] = 0;
	for(size_t i = 0; i < TEST_BIT_HIST_ITERS; i++)
	{
		unsigned long long n = prng.next();
		for(size_t j = 0; j < 64; j++)
		{
			if(n & (1ull << j))
				counts[j]++;
		}
	}
	for(size_t i = 0; i < 64; i++)
	{
		double d = (double)counts[i] / TEST_BIT_HIST_ITERS;
		if(ABS(d - 0.5) > 0.01)
			ThrowError("Poor bit-wise histogram");
	}
}

#define GRANDUINT_TEST_PRELUDE_SIZE 10000
#define GRANDUINT_TEST_PERIOD_SIZE 100000

void GRand_testSpeed()
{
	// Compare speed with rand(). (Be sure to build optimized, or else the results aren't very meaningful.)
	int i;
	uint64 z;
	double t1,t2,t3;
	t1 = GTime::seconds();
	for(i = 0; i < 100000000; i++)
		z = rand();
	t2 = GTime::seconds();
	GRand gr(0);
	for(i = 0; i < 100000000; i++)
		z = gr.next();
	t3 = GTime::seconds();
	double randtime = t2 - t1;
	double grandtime = t3 - t2;
	if(randtime < grandtime)
		ThrowError("rand is faster than GRand");
}

void GRand_testRange()
{
	// Make sure random doubles are within range
	GRand r(0);
	double min = 0.5;
	double max = 0.5;
	for(int n = 0; n < 100000; n++)
	{
		double d = r.uniform();
		min = MIN(min, d);
		max = MAX(max, d);
	}
	if(min < 0.0 || max > 1.0)
		ThrowError("Out of range");
	if(ABS(min - 0.0) > 0.001)
		ThrowError("poor min");
	if(ABS(max - 1.0) > 0.001)
		ThrowError("poor max");
}

// static
void GRand::test()
{
	GRand_testBitHistogram();
	
	// Test cycle length
	int n;
	uint64 rnd;
	uint64 prev;
	for(n = 0; n < 100; n++)
	{
		GRand r(n);
		for(uint64 j = 0; j < GRANDUINT_TEST_PRELUDE_SIZE; j++)
			r.next();
		uint64 startA = r.m_a;
		uint64 startB = r.m_b;
		prev = r.next();
		for(uint64 j = 0; j < GRANDUINT_TEST_PERIOD_SIZE; j++)
		{
			if(r.m_a == startA || r.m_b == startB)
				ThrowError("Loop too small");
			rnd = r.next();
		}
	}

	GRand_testRange();
	//GRand_testSpeed();
	// todo: add a test for correlations
}
#endif // !NO_TEST_CODE

} // namespace GClasses

