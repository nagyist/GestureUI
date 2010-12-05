// ----------------------------------------------------------------
// The contents of this file are distributed under the CC0 license.
// See http://creativecommons.org/publicdomain/zero/1.0/
// ----------------------------------------------------------------

#include <stdio.h>
#include <stdlib.h>
#include "../GClasses/GApp.h"
#include "../GClasses/GBits.h"
#include "../GClasses/GMacros.h"
#include "../GClasses/GData.h"
#include "../GClasses/GImage.h"
#include "../GClasses/GRand.h"
#include "../GClasses/GFile.h"
#include "../GClasses/GHtml.h"
#include "../GClasses/GTransform.h"
#include "../GClasses/GVec.h"
#include "../GClasses/GText.h"
#include "../GSup/GDirList.h"
#include "../GClasses/GSparseMatrix.h"
#include "../GClasses/GHashTable.h"
#include "../GClasses/GHillClimber.h"
#include "../GClasses/GManifold.h"
#include "../GClasses/GMath.h"
#include "../GClasses/GHeap.h"
#include "../GSup/GRayTrace.h"
#include "../GClasses/GRect.h"
#include "../GClasses/GSystemLearner.h"
#include "../GClasses/GTwt.h"
#include "../wizard/usage.h"
#include <time.h>
#include <iostream>
#ifdef WIN32
#	include <direct.h>
#	include <process.h>
#endif
#include <exception>
#include <string>

using namespace GClasses;
using std::cout;
using std::vector;
using std::string;

void cube(GArgReader& args)
{
	int side = args.pop_uint() - 1;
	GData data(3);
	double* pRow;
	for(int y = 0; y < side; y++)
	{
		double b = ((double)y/* + 0.5*/) / side;
		for(int x = 0; x < side; x++)
		{
			double a = ((double)x/* + 0.5*/) / side;
			pRow = data.newRow(); pRow[0] = 0.0; pRow[1] = 1.0 - a; pRow[2] = b;
			pRow = data.newRow(); pRow[0] = 1.0; pRow[1] = a; pRow[2] = 1.0 - b;
			pRow = data.newRow(); pRow[0] = b; pRow[1] = 0.0; pRow[2] = 1.0 - a;
			pRow = data.newRow(); pRow[0] = 1.0 - b; pRow[1] = 1.0; pRow[2] = a;
			pRow = data.newRow(); pRow[0] = 1.0 - a; pRow[1] = b; pRow[2] = 0.0;
			pRow = data.newRow(); pRow[0] = a; pRow[1] = 1.0 - b; pRow[2] = 1.0;
		}
	}
	pRow = data.newRow(); pRow[0] = 0.0; pRow[1] = 0.0; pRow[2] = 0.0;
	pRow = data.newRow(); pRow[0] = 1.0; pRow[1] = 1.0; pRow[2] = 1.0;
	data.print(cout);
}

void fishBowl(GArgReader& args)
{
	int points = args.pop_uint();

	// Parse options
	unsigned int seed = getpid() * (unsigned int)time(NULL);
	double opening = 0.25;
	while(args.next_is_flag())
	{
		if(args.if_pop("-seed"))
			seed = args.pop_uint();
		else if(args.if_pop("-opening"))
			opening = args.pop_double();
	}
	if(opening >= 2.0)
		ThrowError("opening too big--consumes entire fish bowl");

	// Make the data
	GRand prng(seed);
	GData data(3);
	for(int i = 0; i < points; i++)
	{
		double* pRow = data.newRow();
		int j;
		for(j = 100; j > 0; j--)
		{
			prng.spherical(pRow, 3);
			if(pRow[1] < 1.0 - opening)
				break;
		}
		if(j == 0)
			ThrowError("Failed to find a point on the fish bowl");
	}
	data.sort(1);
	data.print(cout);
}

void Noise(GArgReader& args)
{
	int pats = args.pop_uint();

	// Parse options
	unsigned int nSeed = getpid() * (unsigned int)time(NULL);
	string dist = "gaussian";
	double p1 = 0.0;
	double p2 = 1.0;
	vector<double> probs;
	int dims = 1;
	int vals = 0;
	while(args.next_is_flag())
	{
		if(args.if_pop("-seed"))
			nSeed = args.pop_uint();
		else if(args.if_pop("-dist"))
		{
			dist = args.pop_string();
			if(dist.compare("gaussian") == 0 || dist.compare("normal") == 0)
			{
				p1 = args.pop_double();
				p2 = args.pop_double();
			}
			else if(dist.compare("uniform") == 0)
			{
				p1 = args.pop_double();
				p2 = args.pop_double();
			}
			else if(dist.compare("beta") == 0)
			{
				p1 = args.pop_double();
				p2 = args.pop_double();
			}
			else if(dist.compare("binomial") == 0)
			{
				p1 = (double)args.pop_uint();
				p2 = args.pop_double();
			}
			else if(dist.compare("categorical") == 0)
			{
				vals = args.pop_uint();
				for(int i = 0; i < vals; i++)
					probs.push_back(args.pop_double());
			}
			else if(dist.compare("cauchy") == 0)
			{
				p1 = args.pop_double();
				p2 = args.pop_double();
			}
			else if(dist.compare("chisquare") == 0)
				p1 = args.pop_double();
			else if(dist.compare("exponential") == 0)
				p1 = args.pop_double();
			else if(dist.compare("f") == 0)
			{
				p1 = args.pop_double();
				p2 = args.pop_double();
			}
			else if(dist.compare("gamma") == 0)
			{
				p1 = args.pop_double();
				p2 = args.pop_double();
			}
			else if(dist.compare("geometric") == 0)
				p1 = args.pop_double();
			else if(dist.compare("logistic") == 0)
			{
				p1 = args.pop_double();
				p2 = args.pop_double();
			}
			else if(dist.compare("lognormal") == 0)
			{
				p1 = args.pop_double();
				p2 = args.pop_double();
			}
			else if(dist.compare("poisson") == 0)
				p1 = args.pop_double();
			else if(dist.compare("softimpulse") == 0)
				p1 = args.pop_double();
			else if(dist.compare("spherical") == 0)
			{
				dims = args.pop_uint();
				p1 = args.pop_double();
			}
			else if(dist.compare("student") == 0)
				p1 = args.pop_double();
			else if(dist.compare("weibull") == 0)
				p1 = args.pop_double();
		}
		else
			ThrowError("Invalid option: ", args.peek());
	}

	// Make the data
	GRand prng(nSeed);
	GUniformRelation* pRelation = new GUniformRelation(dims, vals);
	sp_relation pRel = pRelation;
	GData data(pRel);
	if(dist.compare("gaussian") == 0)
	{
		for(int i = 0; i < pats; i++)
			data.newRow()[0] = prng.normal() * p2 + p1;
	}
	else if(dist.compare("uniform") == 0)
	{
		for(int i = 0; i < pats; i++)
			data.newRow()[0] = prng.uniform() * (p2 - p1) + p1;
	}
	else if(dist.compare("beta") == 0)
	{
		for(int i = 0; i < pats; i++)
			data.newRow()[0] = prng.beta(p1, p2);
	}
	else if(dist.compare("binomial") == 0)
	{
		for(int i = 0; i < pats; i++)
			data.newRow()[0] = prng.binomial((int)p1, p2);
	}
	else if(dist.compare("categorical") == 0)
	{
		for(int i = 0; i < pats; i++)
			data.newRow()[0] = prng.categorical(probs);
	}
	else if(dist.compare("cauchy") == 0)
	{
		for(int i = 0; i < pats; i++)
			data.newRow()[0] = prng.cauchy() * p2 + p1;
	}
	else if(dist.compare("chisquare") == 0)
	{
		for(int i = 0; i < pats; i++)
			data.newRow()[0] = prng.chiSquare(p1);
	}
	else if(dist.compare("exponential") == 0)
	{
		for(int i = 0; i < pats; i++)
			data.newRow()[0] = prng.exponential() * p1;
	}
	else if(dist.compare("f") == 0)
	{
		for(int i = 0; i < pats; i++)
			data.newRow()[0] = prng.f(p1, p2);
	}
	else if(dist.compare("gamma") == 0)
	{
		for(int i = 0; i < pats; i++)
			data.newRow()[0] = prng.gamma(p1) / p2;
	}
	else if(dist.compare("geometric") == 0)
	{
		for(int i = 0; i < pats; i++)
			data.newRow()[0] = (double)prng.geometric(p1);
	}
	else if(dist.compare("logistic") == 0)
	{
		for(int i = 0; i < pats; i++)
			data.newRow()[0] = prng.logistic() * p2 + p1;
	}
	else if(dist.compare("lognormal") == 0)
	{
		for(int i = 0; i < pats; i++)
			data.newRow()[0] = prng.logNormal(p1, p2);
	}
	else if(dist.compare("poisson") == 0)
	{
		for(int i = 0; i < pats; i++)
			data.newRow()[0] = (double)prng.poisson(p1);
	}
	else if(dist.compare("softimpulse") == 0)
	{
		for(int i = 0; i < pats; i++)
			data.newRow()[0] = prng.softImpulse(p1);
	}
	else if(dist.compare("spherical") == 0)
	{
		for(int i = 0; i < pats; i++)
		{
			double* pRow = data.newRow();
			prng.spherical(pRow, dims);
			GVec::multiply(pRow, p1, dims);
		}
	}
	else if(dist.compare("student") == 0)
	{
		for(int i = 0; i < pats; i++)
			data.newRow()[0] = prng.student(p1);
	}
	else if(dist.compare("weibull") == 0)
	{
		for(int i = 0; i < pats; i++)
			data.newRow()[0] = prng.weibull(p1);
	}
	else
		ThrowError("Unrecognized distribution: ", dist.c_str());

	data.print(cout);
}

bool isInsideUnitStar(double x, double y)
{
	int count = 0;
	for(int n = 0; n < 5; n++)
	{
		double r = n * 2 * M_PI / 5;
		double c = cos(r);
		double s = sin(r);
		x += s;
		y -= c;
		if((x * s) - (y * c) >= 0)
			count++;
		x -= s;
		y += c;
	}
	return (count >= 4);
}

double LengthOfSwissRoll(double x)
{
#ifdef WIN32
	GAssert(false); // not implemented yet for Win32
	return 0;
#else
	return (x * sqrt(x * x + 1) + asinh(x)) / 2;
#endif
}

void SwissRoll(GArgReader& args)
{
	int points = args.pop_uint();

	// Parse options
	unsigned int nSeed = getpid() * (unsigned int)time(NULL);
	bool cutOutStar = false;
	bool reduced = false;
	while(args.next_is_flag())
	{
		if(args.if_pop("-seed"))
			nSeed = args.pop_uint();
		else if(args.if_pop("-cutoutstar"))
			cutOutStar = true;
		else if(args.if_pop("-reduced"))
			reduced = true;
		else
			ThrowError("Invalid option: ", args.peek());
	}

	// Generate the data
	double unrolledWidth = 48.790603865663797;
	double starCenterX = 43;
	double starSize = 1.5;
	GRand prng(nSeed);
	GArffRelation rel;
	rel.addAttribute("x", 0, NULL);
	rel.addAttribute("y", 0, NULL);
	if(!reduced)
		rel.addAttribute("z", 0, NULL);
	GData data(3);
	for(int n = 0; n < points; n++)
	{
		double t = ((double)n * 8) / points;
		double* pVector = data.newRow();
		pVector[0] = ((t + 2) * sin(t));
		pVector[2] = ((t + 2) * cos(t));
		if(cutOutStar)
		{
			int i;
			for(i = 0; i < 1000; i++)
			{
				pVector[1] = prng.uniform() * 12;
				if(!isInsideUnitStar((n * unrolledWidth / points - starCenterX) / starSize, (6 - pVector[1]) / starSize))
					break;
			}
			if(i >= 1000)
				ThrowError("The star is too big. It severs the manifold.");
		}
		else
			pVector[1] = prng.uniform() * 12;
		if(reduced)
		{
			pVector[0] = pVector[1];
			pVector[1] = LengthOfSwissRoll(t + 2);
			pVector[2] = 0;
		}
	}

	// Print the data
	data.print(cout);
}

double LengthOfSineFunc(void* pThis, double x)
{
	double d = cos(x);
	return sqrt(d * d + 1.0);
}

void SCurve(GArgReader& args)
{
	int points = args.pop_uint();

	// Parse options
	unsigned int nSeed = getpid() * (unsigned int)time(NULL);
	bool reduced = false;
	while(args.next_is_flag())
	{
		if(args.if_pop("-seed"))
			nSeed = args.pop_uint();
		else if(args.if_pop("-reduced"))
			reduced = true;
		else
			ThrowError("Invalid option: ", args.peek());
	}

	// Generate the data
	GRand prng(nSeed);
	GArffRelation rel;
	rel.addAttribute("x", 0, NULL);
	rel.addAttribute("y", 0, NULL);
	if(!reduced)
		rel.addAttribute("z", 0, NULL);
	GData data(3);
	for(int n = 0; n < points; n++)
	{
		double t = ((double)n * 2.2 * M_PI - .1 * M_PI) / points;
		double* pVector = data.newRow();
		pVector[0] = 1.0 - sin(t);
		pVector[1] = t;
		pVector[2] = prng.uniform() * 2;
		if(reduced)
		{
			pVector[0] = pVector[2];
			pVector[1] = (n > 0 ? GMath::integrate(LengthOfSineFunc, 0, t, n + 30, NULL) : 0);
			pVector[2] = 0;
		}
	}

	// Print the data
	data.print(cout);
}

void EntwinedSpirals(GArgReader& args)
{
	int points = args.pop_uint();

	// Parse options
	unsigned int nSeed = getpid() * (unsigned int)time(NULL);
	bool reduced = false;
	while(args.next_is_flag())
	{
		if(args.if_pop("-seed"))
			nSeed = args.pop_uint();
		else if(args.if_pop("-reduced"))
			reduced = true;
		else
			ThrowError("Invalid option: ", args.peek());
	}

	// Generate the data
	GRand prng(nSeed);
	GArffRelation rel;
	rel.addAttribute("x", 0, NULL);
	if(!reduced)
	{
		rel.addAttribute("y", 0, NULL);
		rel.addAttribute("z", 0, NULL);
	}
	GData data(3);

	double dHeight = 3;
	double dWraps = 1.5;
	double dSpiralLength = sqrt((dWraps * 2.0 * M_PI) * (dWraps * 2.0 * M_PI) + dHeight * dHeight);
	double dTotalLength = 2.0 * (dSpiralLength + 1); // radius = 1
	double d;
	for(int n = 0; n < points; n++)
	{
		double t = ((double)n * dTotalLength) / points;
		double* pVector = data.newRow();
		if(reduced)
		{
			pVector[0] = t;
			pVector[1] = 0;
			pVector[2] = 0;
		}
		else
		{
			if(t < dSpiralLength)
			{
				d = (dSpiralLength - t) * dWraps * 2 * M_PI / dSpiralLength; // d = radians
				pVector[0] = -cos(d);
				pVector[1] = dHeight * t / dSpiralLength;
				pVector[2] = -sin(d);
			}
			else if(t - 2.0 - dSpiralLength >= 0)
			{
				d = (t - 2.0 - dSpiralLength) * dWraps * 2 * M_PI / dSpiralLength; // d = radians
				pVector[0] = cos(d);
				pVector[1] = dHeight * (dSpiralLength - (t - 2.0 - dSpiralLength)) / dSpiralLength;
				pVector[2] = sin(d);
			}
			else
			{
				d = (t - dSpiralLength) / 2.0; // 2 = diameter
				pVector[0] = 2.0 * d - 1.0;
				pVector[1] = dHeight;
				pVector[2] = 0;
			}
		}
	}

	// Print the data
	data.print(cout);
}

void ImageTranslatedOverNoise(GArgReader& args)
{
	const char* szFilenameIn = args.pop_string();

	// Parse options
	unsigned int nSeed = getpid() * (unsigned int)time(NULL);
	bool reduced = false;
	while(args.next_is_flag())
	{
		if(args.if_pop("-seed"))
			nSeed = args.pop_uint();
		else if(args.if_pop("-reduced"))
			reduced = true;
		else
			ThrowError("Invalid option: ", args.peek());
	}

	// Load the image
	GImage imageFace;
	imageFace.loadPng(szFilenameIn);

	// Generate the noise image
	GRand prng(nSeed);
	GImage imageNoise;
	imageNoise.setSize(imageFace.width() * 2, imageFace.height() * 2);
	for(int y = 0; y < (int)imageNoise.height(); y++)
	{
		for(int x = 0; x < (int)imageNoise.width(); x++)
			imageNoise.setPixel(x, y, gARGB(0xff, (int)prng.next(256), (int)prng.next(256), (int)prng.next(256)));
	}

	// Make the relation
	GImage imageAll;
	imageAll.setSize(imageNoise.width(), imageNoise.height());
	GArffRelation rel;
	if(reduced)
	{
		rel.addAttribute("x", 0, NULL);
		rel.addAttribute("y", 0, NULL);
	}
	else
	{
		for(int y = 0; y < (int)imageAll.height(); y++)
		{
			for(int x = 0; x < (int)imageAll.width(); x++)
				rel.addAttribute("pix", 0, NULL);
		}
	}

	// Generate the data
	GData data(reduced ? 2 : imageAll.width() * imageAll.height());
	GRect r(0, 0, imageFace.width(), imageFace.height());
	for(int y = 0; y <= (int)imageFace.height(); y++)
	{
		for(int x = 0; x <= (int)imageFace.width(); x++)
		{
			if(reduced)
			{
				double* pVec = data.newRow();
				pVec[0] = x;
				pVec[1] = y;
			}
			else
			{
				imageAll.copy(&imageNoise);
				imageAll.blit(x, y, &imageFace, &r);
				double* pVec = data.newRow();
				for(int yy = 0; yy < (int)imageAll.height(); yy++)
				{
					for(int xx = 0; xx < (int)imageAll.width(); xx++)
						pVec[imageAll.width() * yy + xx] = (double)gGray(imageAll.pixel(xx, yy)) / MAX_GRAY_VALUE;
				}
			}
		}
	}

	// Print the data
	data.print(cout);
}

void SelfIntersectingRibbon(GArgReader& args)
{
	unsigned int points = args.pop_uint();

	// Parse options
	unsigned int nSeed = getpid() * (unsigned int)time(NULL);
	while(args.next_is_flag())
	{
		if(args.if_pop("-seed"))
			nSeed = args.pop_uint();
		else
			ThrowError("Invalid option: ", args.peek());
	}

	GRand prng(nSeed);
	GData data(3);
	for(unsigned int i = 0; i < points; i++)
	{
		double t = M_PI / 4 + (1.5 * M_PI) * (double)i / points;
		double* pPat = data.newRow();
		pPat[0] = sin(t * 2);
		pPat[1] = -2.0 * cos(t);
		pPat[2] = 2.0 * prng.uniform();
	}

	data.print(cout);
}

void WindowedImageData(GArgReader& args)
{
	const char* szFilenameIn = args.pop_string();

	// Load the image
	GImage imageSource;
	imageSource.loadPng(szFilenameIn);

	// Parse options
	bool reduced = false;
	int hstep = 1;
	int vstep = 1;
	int windowWidth = imageSource.width() / 2;
	int windowHeight = imageSource.height() / 2;
	int hole = 0;
	while(args.next_is_flag())
	{
		if(args.if_pop("-reduced"))
			reduced = true;
		else if(args.if_pop("-stepsizes"))
		{
			hstep = args.pop_uint();
			vstep = args.pop_uint();
		}
		else if(args.if_pop("-windowsize"))
		{
			windowWidth = args.pop_uint();
			windowHeight = args.pop_uint();
		}
		else if(args.if_pop("-hole"))
			hole = args.pop_uint();
		else
			ThrowError("Invalid option: ", args.peek());
	}

	// Make the relation
	GImage imageWindow;
	imageWindow.setSize(windowWidth, windowHeight);
	GArffRelation rel;
	if(reduced)
	{
		rel.addAttribute("x", 0, NULL);
		rel.addAttribute("y", 0, NULL);
	}
	else
	{
		for(int y = 0; y < (int)imageWindow.height(); y++)
		{
			for(int x = 0; x < (int)imageWindow.width(); x++)
			{
				rel.addAttribute("r", 0, NULL);
				rel.addAttribute("g", 0, NULL);
				rel.addAttribute("b", 0, NULL);
			}
		}
	}

	// Generate the data
	int centerx = ((int)imageSource.width() - (int)imageWindow.width()) / 2;
	int centery = ((int)imageSource.height() - (int)imageWindow.height()) / 2;
	GData data(reduced ? 2 : imageWindow.width() * imageWindow.height() * 3);
	for(int y = 0; y + (int)imageWindow.height() <= (int)imageSource.height(); y += vstep)
	{
		for(int x = 0; x + (int)imageWindow.width() <= (int)imageSource.width(); x += hstep)
		{
			if(hole > 0 && x > centerx - hole && x < centerx + hole && y > centery - hole && y < centery + hole)
				continue;
			if(reduced)
			{
				double* pVec = data.newRow();
				pVec[0] = x;
				pVec[1] = y;
			}
			else
			{
				GRect r(x, y, imageWindow.width(), imageWindow.height());
				imageWindow.blit(0, 0, &imageSource, &r);
				double* pVec = data.newRow();
				for(int yy = 0; yy < (int)imageWindow.height(); yy++)
				{
					for(int xx = 0; xx < (int)imageWindow.width(); xx++)
					{
						unsigned int pix = imageWindow.pixel(xx, yy);
						*(pVec++) = (double)gRed(pix) / 256;
						*(pVec++) = (double)gGreen(pix) / 256;
						*(pVec++) = (double)gBlue(pix) / 256;
					}
				}
			}
		}
	}

	// Print the data
	data.print(cout);
}

void addCraneToScene(GRayTraceScene& scene, double craneYaw, double ballHeight, double ballRadius, double x, double y, double z)
{
	double cranePitch = 1.0;
	double craneLength = 3.0;
	double craneRadius = 0.18;
	double cableRadius = 0.03;

	// Add the crane
	GRayTracePhysicalMaterial* pCraneMaterial = new GRayTracePhysicalMaterial();
	pCraneMaterial->setColor(GRayTraceMaterial::Diffuse, 0.9, 0.7, 0.3);
	pCraneMaterial->setColor(GRayTraceMaterial::Reflective, 0.45, 0.35, 0.15);
	scene.addMaterial(pCraneMaterial);
	G3DVector craneBottom, craneTop;
	//craneBottom.set(0, 0, 0);
	craneBottom.set(x - craneLength * cos(cranePitch) * sin(craneYaw), y - craneLength * sin(cranePitch), z - craneLength * cos(cranePitch) * (-cos(craneYaw)));
	craneTop.set(x + craneLength * cos(cranePitch) * sin(craneYaw), y + craneLength * sin(cranePitch), z + craneLength * cos(cranePitch) * (-cos(craneYaw)));
	GRayTraceTriMesh* pCrane = GRayTraceTriMesh::makeCylinder(pCraneMaterial, &craneBottom, &craneTop, craneRadius, 12, false/*end caps*/);
	pCrane->computePhongNormals();
	scene.addMesh(pCrane);

	// Add the cable
	GRayTracePhysicalMaterial* pCableMaterial = new GRayTracePhysicalMaterial();
	pCableMaterial->setColor(GRayTraceMaterial::Diffuse, 0.6, 0.6, 0.6);
	scene.addMaterial(pCableMaterial);
	G3DVector ballCenter;
	ballCenter.copy(&craneTop);
	ballCenter.m_vals[1] = y + ballHeight;
	GRayTraceTriMesh* pCable = GRayTraceTriMesh::makeCylinder(pCableMaterial, &ballCenter, &craneTop, cableRadius, 12, false/*end caps*/);
	pCable->computePhongNormals();
	scene.addMesh(pCable);

	// Add the ball
	GRayTracePhysicalMaterial* pBallMaterial = new GRayTracePhysicalMaterial();
	pBallMaterial->setColor(GRayTraceMaterial::Diffuse, 0.1, 0.1, 0.2);
	pBallMaterial->setColor(GRayTraceMaterial::Reflective, 0.2, 0.2, 0.4);
	scene.addMaterial(pBallMaterial);
	GRayTraceSphere* pBall = new GRayTraceSphere(pBallMaterial, ballCenter.m_vals[0], ballCenter.m_vals[1], ballCenter.m_vals[2], ballRadius);
	scene.addObject(pBall);
}

GImage* makeCraneImage(double craneYaw, double ballHeight, int wid, int hgt, double ballRadius, bool front)
{
	if(front)
		craneYaw = M_PI - craneYaw;

	// Make a scene
	GRand prng(0);
	GRayTraceScene scene(&prng);
	scene.setAmbientLight(0.3, 0.3, 0.3);
	scene.addLight(new GRayTraceDirectionalLight(-1.0, 2.0, 1.0, 0.1/*r*/, 0.2/*g*/, 0.3/*b*/, 0.0/*jitter*/));
	scene.addLight(new GRayTraceDirectionalLight(1.0, 0.0, 0.3, 0.10/*r*/, 0.15/*g*/, 0.05/*b*/, 0.0/*jitter*/));
	scene.addLight(new GRayTraceDirectionalLight(0.0, 0.0, 1.0, 0.15/*r*/, 0.05/*g*/, 0.10/*b*/, 0.0/*jitter*/));
	scene.setBackgroundColor(1/*a*/, 1.0/*r*/, 1.0/*g*/, 1.0/*b*/);
	GRayTraceCamera* pCamera = scene.camera();
	pCamera->setViewAngle(1.0);
	pCamera->setImageSize(wid, hgt);
	if(front)
		pCamera->lookFromPoint()->set(0.0, 0.7, 3.5);
	else
		pCamera->lookFromPoint()->set(0.0, 0.7, 1.0);
	G3DVector cameraDirection(0.0, 0.0, -1.0);
	pCamera->setDirection(&cameraDirection, 0.0);
	//scene.setRenderMode(GRayTraceScene::FAST_RAY_TRACE);
	scene.setRenderMode(GRayTraceScene::QUALITY_RAY_TRACE);

	addCraneToScene(scene, craneYaw, ballHeight, ballRadius, 0.0, 0.0, 0.0);

	// Render the image
	scene.render();
	return scene.releaseImage();
}

GImage* makeThreeCraneImage(double* craneYaw, double* ballHeight, int wid, int hgt, double ballRadius)
{
	craneYaw[0] = M_PI - craneYaw[0];
	craneYaw[1] = M_PI - craneYaw[1];
	craneYaw[2] = M_PI - craneYaw[2];

	// Make a scene
	GRand prng(0);
	GRayTraceScene scene(&prng);
	scene.setAmbientLight(0.3, 0.3, 0.3);
	scene.addLight(new GRayTraceDirectionalLight(-1.0, 2.0, 1.0, 0.1/*r*/, 0.2/*g*/, 0.3/*b*/, 0.0/*jitter*/));
	scene.addLight(new GRayTraceDirectionalLight(1.0, 0.0, 0.3, 0.10/*r*/, 0.15/*g*/, 0.05/*b*/, 0.0/*jitter*/));
	scene.addLight(new GRayTraceDirectionalLight(0.0, 0.0, 1.0, 0.15/*r*/, 0.05/*g*/, 0.10/*b*/, 0.0/*jitter*/));
	scene.setBackgroundColor(1/*a*/, 1.0/*r*/, 1.0/*g*/, 1.0/*b*/);
	GRayTraceCamera* pCamera = scene.camera();
	pCamera->setViewAngle(1.0);
	pCamera->setImageSize(wid, hgt);
	pCamera->lookFromPoint()->set(0.0, 0.7/*0.84*/, 4.2);
	G3DVector cameraDirection(0.0, 0.0, -1.0);
	pCamera->setDirection(&cameraDirection, 0.0);
	//scene.setRenderMode(GRayTraceScene::FAST_RAY_TRACE);
	scene.setRenderMode(GRayTraceScene::QUALITY_RAY_TRACE);

	addCraneToScene(scene, craneYaw[0], ballHeight[0], ballRadius, -1.0, 0.0, 0.0);
	addCraneToScene(scene, craneYaw[1], ballHeight[1], ballRadius, 0.0, 0.0, 0.0);
	addCraneToScene(scene, craneYaw[2], ballHeight[2], ballRadius, 1.0, 0.0, 0.0);

	// Render the image
	scene.render();
	return scene.releaseImage();
}

void CraneDataset(GArgReader& args)
{
	// Parse options
	int wid = 64;
	int hgt = 48;
	int horizFrames = 25;
	int vertFrames = 21;
	string imageFile = "";
	double ballRadius = 0.3;
	bool front = false;
	double blur = 0.0;
	bool gray = false;
	while(args.next_is_flag())
	{
		if(args.if_pop("-front"))
			front = true;
		else if(args.if_pop("-saveimage"))
			imageFile = args.pop_string();
		else if(args.if_pop("-ballradius"))
			ballRadius = args.pop_double();
		else if(args.if_pop("-frames"))
		{
			horizFrames = args.pop_uint();
			vertFrames = args.pop_uint();
		}
		else if(args.if_pop("-size"))
		{
			wid = args.pop_uint();
			hgt = args.pop_uint();
		}
		else if(args.if_pop("-blur"))
			blur = args.pop_double();
		else if(args.if_pop("-gray"))
			gray = true;
		else
			ThrowError("Invalid option: ", args.peek());
	}

	GData data(wid * hgt * (gray ? 1 : 3));
	GImage grid;
	grid.setSize(horizFrames * wid, vertFrames * hgt);
	double dx = (0.75 - (-0.75)) / (horizFrames - 1);
	double dy = (1.25 - 0.0) / (vertFrames - 1);

	for(int vert = 0; vert < vertFrames; vert++)
	{
		double ballHeight = (vertFrames - 1 - vert) * dy;
		for(int horiz = 0; horiz < horizFrames; horiz++)
		{
			double craneYaw = horiz * dx - 0.75;

			GImage* pImage = makeCraneImage(craneYaw, ballHeight, wid, hgt, ballRadius, front);
			if(blur > 0)
				pImage->blur(blur);
			Holder<GImage> hImage(pImage);

			GRect r(0, 0, wid, hgt);
			grid.blit(wid * horiz, hgt * vert, pImage, &r);
			grid.box(wid * horiz, hgt * vert, wid * horiz + wid - 1, hgt * vert + hgt - 1, 0xffb0b0b0);

			double* pRow = data.newRow();
			unsigned int* pPixels = pImage->pixels();
			for(int yy = 0; yy < hgt; yy++)
			{
				for(int xx = 0; xx < wid; xx++)
				{
					if(gray)
						*(pRow++) = gGray(*pPixels);
					else
					{
						*(pRow++) = gRed(*pPixels);
						*(pRow++) = gGreen(*pPixels);
						*(pRow++) = gBlue(*pPixels);
					}
					pPixels++;
				}
			}
		}
	}

	if(imageFile.length() > 0)
		grid.savePng(imageFile.c_str());
	data.print(cout);
}

void cranePath(GArgReader& args)
{
	GData* pActions = GData::loadArff(args.pop_string());
	Holder<GData> hActions(pActions);

	// Parse options
	unsigned int nSeed = getpid() * (unsigned int)time(NULL);
	int wid = 64;
	int hgt = 48;
	int horizFrames = 25;
	int vertFrames = 21;
	double state[2];
	double ballRadius = 0.3;
	bool front = false;
	double blur = 0.0;
	bool gray = false;
	string stateFile = "";
	double transitionNoiseDev = 0.0;
	double observationNoiseDev = 0.0;
	while(args.next_is_flag())
	{
		if(args.if_pop("-front"))
			front = true;
		else if(args.if_pop("-seed"))
			nSeed = args.pop_uint();
		else if(args.if_pop("-ballradius"))
			ballRadius = args.pop_double();
		else if(args.if_pop("-size"))
		{
			wid = args.pop_uint();
			hgt = args.pop_uint();
		}
		else if(args.if_pop("-blur"))
			blur = args.pop_double();
		else if(args.if_pop("-gray"))
			gray = true;
		else if(args.if_pop("-state"))
			stateFile = args.pop_string();
		else if(args.if_pop("-noise"))
		{
			transitionNoiseDev = args.pop_double();
			observationNoiseDev = args.pop_double();
		}
		else
			ThrowError("Invalid option: ", args.peek());
	}

	// Make the observation data
	GRand prng(nSeed);
	GData data(wid * hgt * (gray ? 1 : 3));
	GData stateData(2);
	double dx = (0.75 - (-0.75)) / (horizFrames - 1);
	double dy = (1.25 - 0.0) / (vertFrames - 1);
	double posx = 0.0;
	double posy = 1.25 / 2.0;
	for(size_t i = 0; i < pActions->rows(); i++)
	{
		double ballHeight = 1.25 - posy;
		double craneYaw = posx;
		state[0] = (craneYaw + 0.75) * (horizFrames - 1) / 1.5;
		state[1] = ballHeight * (vertFrames - 1) / 1.25;
		GVec::copy(stateData.newRow(), state, 2);

		// Generate the image vector
		GImage* pImage = makeCraneImage(craneYaw, ballHeight, wid, hgt, ballRadius, front);
		if(blur > 0)
			pImage->blur(blur);
		Holder<GImage> hImage(pImage);
		double* pRow = data.newRow();
		unsigned int* pPixels = pImage->pixels();
		for(int yy = 0; yy < hgt; yy++)
		{
			for(int xx = 0; xx < wid; xx++)
			{
				if(gray)
					*(pRow++) = ClipChan((int)floor((double)gGray(*pPixels) + 0.5 + prng.normal() * 255.0 * observationNoiseDev));
				else
				{
					*(pRow++) = ClipChan((int)floor((double)gRed(*pPixels) + 0.5 + prng.normal() * 255.0 * observationNoiseDev));
					*(pRow++) = ClipChan((int)floor((double)gGreen(*pPixels) + 0.5 + prng.normal() * 255.0 * observationNoiseDev));
					*(pRow++) = ClipChan((int)floor((double)gBlue(*pPixels) + 0.5 + prng.normal() * 255.0 * observationNoiseDev));
				}
				pPixels++;
			}
		}

		// Move
		int a = (int)pActions->row(i)[0];
		double deltax = 0;
		double deltay = 0;
		switch(a)
		{
			case 0: deltax = -1.0; deltay = 0.0; break;
			case 1: deltax = 1.0; deltay = 0.0; break;
			case 2: deltax = 0.0; deltay = -1.0; break;
			case 3: deltax = 0.0; deltay = 1.0; break;
		}
		deltax += prng.normal() * transitionNoiseDev;
		deltay += prng.normal() * transitionNoiseDev;
		deltax *= dx;
		deltay *= dy;
		posx = MIN(0.75, MAX(-0.75, posx + deltax));
		posy = MIN(1.25, MAX(0.0, posy + deltay));
	}

	if(stateFile.length() > 0)
		stateData.saveArff(stateFile.c_str());
	data.print(cout);
}

void threeCranePath(GArgReader& args)
{
	GData* pActions = GData::loadArff(args.pop_string());
	Holder<GData> hActions(pActions);

	// Parse options
	unsigned int nSeed = getpid() * (unsigned int)time(NULL);
	int wid = 64;
	int hgt = 48;
	int horizFrames = 25;
	int vertFrames = 21;
	double state[4];
	double ballRadius = 0.3;
	double blur = 0.0;
	bool gray = false;
	string stateFile = "";
	double transitionNoiseDev = 0.0;
	double observationNoiseDev = 0.0;
	while(args.next_is_flag())
	{
		if(args.if_pop("-seed"))
			nSeed = args.pop_uint();
		else if(args.if_pop("-ballradius"))
			ballRadius = args.pop_double();
		else if(args.if_pop("-size"))
		{
			wid = args.pop_uint();
			hgt = args.pop_uint();
		}
		else if(args.if_pop("-blur"))
			blur = args.pop_double();
		else if(args.if_pop("-gray"))
			gray = true;
		else if(args.if_pop("-state"))
			stateFile = args.pop_string();
		else if(args.if_pop("-noise"))
		{
			transitionNoiseDev = args.pop_double();
			observationNoiseDev = args.pop_double();
		}
		else
			ThrowError("Invalid option: ", args.peek());
	}

	// Make the observation data
	GRand prng(nSeed);
	GData data(wid * hgt * (gray ? 1 : 3));
	GData stateData(4);
	double dx = (0.75 - (-0.75)) / (horizFrames - 1);
	double dy = (1.25 - 0.0) / (vertFrames - 1);
	double posx[3]; posx[0] = 0.0; posx[1] = 0.0; posx[2] = 0.0;
	double posy[3]; posy[0] = 1.25 / 2.0; posy[1] = 1.25 / 2.0; posy[2] = 1.25 / 2.0;
	for(size_t i = 0; i < pActions->rows(); i++)
	{
		double ballHeight[3];
		ballHeight[0] = 1.25 - posy[0];
		ballHeight[1] = 1.25 - posy[1];
		ballHeight[2] = 1.25 - posy[2];
		double craneYaw[3];
		craneYaw[0] = posx[0];
		craneYaw[1] = posx[1];
		craneYaw[2] = posx[2];
		state[0] = (craneYaw[0] + 0.75) * (horizFrames - 1) / 1.5;
		state[1] = ballHeight[0] * (vertFrames - 1) / 1.25;
		state[2] = (craneYaw[2] + 0.75) * (horizFrames - 1) / 1.5;
		state[3] = ballHeight[2] * (vertFrames - 1) / 1.25;
		GVec::copy(stateData.newRow(), state, 4);

		// Generate the image vector
		GImage* pImage = makeThreeCraneImage(craneYaw, ballHeight, wid, hgt, ballRadius);
		if(blur > 0)
			pImage->blur(blur);
		Holder<GImage> hImage(pImage);
		double* pRow = data.newRow();
		unsigned int* pPixels = pImage->pixels();
		for(int yy = 0; yy < hgt; yy++)
		{
			for(int xx = 0; xx < wid; xx++)
			{
				if(gray)
					*(pRow++) = ClipChan((int)floor((double)gGray(*pPixels) + 0.5 + prng.normal() * 255.0 * observationNoiseDev));
				else
				{
					*(pRow++) = ClipChan((int)floor((double)gRed(*pPixels) + 0.5 + prng.normal() * 255.0 * observationNoiseDev));
					*(pRow++) = ClipChan((int)floor((double)gGreen(*pPixels) + 0.5 + prng.normal() * 255.0 * observationNoiseDev));
					*(pRow++) = ClipChan((int)floor((double)gBlue(*pPixels) + 0.5 + prng.normal() * 255.0 * observationNoiseDev));
				}
				pPixels++;
			}
		}

		// Move
		int a = (int)pActions->row(i)[0];
		double deltax = 0;
		double deltay = 0;
		switch(a)
		{
			case 0: deltax = -1.0; deltay = 0.0; break;
			case 1: deltax = 1.0; deltay = 0.0; break;
			case 2: deltax = 0.0; deltay = -1.0; break;
			case 3: deltax = 0.0; deltay = 1.0; break;
		}
		deltax += prng.normal() * transitionNoiseDev;
		deltay += prng.normal() * transitionNoiseDev;
		deltax *= dx;
		deltay *= dy;
		posx[0] = MIN(0.75, MAX(-0.75, posx[0] + deltax));
		posy[0] = MIN(1.25, MAX(0.0, posy[0] + deltay));

		a = (int)prng.next(4);
		switch(a)
		{
			case 0: deltax = -1.0; deltay = 0.0; break;
			case 1: deltax = 1.0; deltay = 0.0; break;
			case 2: deltax = 0.0; deltay = -1.0; break;
			case 3: deltax = 0.0; deltay = 1.0; break;
		}
		deltax += prng.normal() * transitionNoiseDev;
		deltay += prng.normal() * transitionNoiseDev;
		deltax *= dx;
		deltay *= dy;
		posx[1] = MIN(0.75, MAX(-0.75, posx[1] + deltax));
		posy[1] = MIN(1.25, MAX(0.0, posy[1] + deltay));

		a = (int)pActions->row(i)[1];
		switch(a)
		{
			case 0: deltax = -1.0; deltay = 0.0; break;
			case 1: deltax = 1.0; deltay = 0.0; break;
			case 2: deltax = 0.0; deltay = -1.0; break;
			case 3: deltax = 0.0; deltay = 1.0; break;
		}
		deltax += prng.normal() * transitionNoiseDev;
		deltay += prng.normal() * transitionNoiseDev;
		deltax *= dx;
		deltay *= dy;
		posx[2] = MIN(0.75, MAX(-0.75, posx[2] + deltax));
		posy[2] = MIN(1.25, MAX(0.0, posy[2] + deltay));
	}

	if(stateFile.length() > 0)
		stateData.saveArff(stateFile.c_str());
	data.print(cout);
}

void ScaleAndRotate(GArgReader& args)
{
	// Load the image
	GImage imageSource;
	imageSource.loadPng(args.pop_string());

	// Parse options
	int rotateFrames = 40;
	int scaleFrames = 15;
	double arc = 2.0 * M_PI;
	string imageFile = "";
	while(args.next_is_flag())
	{
		if(args.if_pop("-saveimage"))
			imageFile = args.pop_string();
		else if(args.if_pop("-frames"))
		{
			rotateFrames = args.pop_uint();
			scaleFrames = args.pop_uint();
		}
		else if(args.if_pop("-arc"))
			arc = args.pop_double();
		else
			ThrowError("Invalid option: ", args.peek());
	}

	// Make the frames
	double finalSize = 0.4;
	GImage grid;
	grid.setSize(rotateFrames * imageSource.width(), scaleFrames * imageSource.height());
	GImage rotated;
	GImage frame;
	frame.setSize(imageSource.width(), imageSource.height());
	GData data(imageSource.width() * imageSource.height());
	for(int y = 0; y < scaleFrames; y++)
	{
		for(int x = 0; x < rotateFrames; x++)
		{
			// Make the frame
			rotated.rotate(&imageSource, imageSource.width() / 2, imageSource.height() / 2, (double)x * arc / rotateFrames);
			GDoubleRect rSrc(1.0, 1.0, imageSource.width() - 2, imageSource.height() - 2);
			double s = (double)y / (scaleFrames - 1);
			GDoubleRect rDest(1.0, 1.0, (1.0 - s) * rSrc.w + s * finalSize * rSrc.w, (1.0 - s) * rSrc.h + s * finalSize * rSrc.h);
			rDest.x += (rSrc.w - rDest.w) / 2.0;
			rDest.y += (rSrc.h - rDest.h) / 2.0;
			frame.clear(0xffffffff);
			frame.blitStretchInterpolate(&rDest, &rSrc, &rotated);

			// Convert to a data row
			double* pRow = data.newRow();
			for(int b = 0; b < (int)imageSource.height(); b++)
			{
				for(int a = 0; a < (int)imageSource.width(); a++)
					*(pRow++) = gGray(frame.pixel(a, b));
			}

			// Save in the grid image
			GRect r(0, 0, imageSource.width(), imageSource.height());
			frame.box(0, 0, imageSource.width() - 1, imageSource.height() - 1, 0xffa0a0a0);
			grid.blit(x * imageSource.width(), y * imageSource.height(), &frame, &r);
		}
	}
	if(imageFile.length() > 0)
		grid.savePng(imageFile.c_str());
	data.print(cout);
}

void gridRandomWalk(GArgReader& args)
{
	GData* pData = GData::loadArff(args.pop_string());
	Holder<GData> hData(pData);
	int width = args.pop_uint();
	int height = pData->rows() / width;
	if((height * width) != (int)pData->rows())
		ThrowError("Expected a dataset with a number of rows that is divisible by width");
	size_t samples = args.pop_uint();
	if(pData->rows() < 2)
		ThrowError("Expected at least two states");

	// Parse options
	unsigned int nSeed = getpid() * (unsigned int)time(NULL);
	int x = width / 2;
	int y = height / 2;
	string obsFile = "observations.arff";
	string actionFile = "actions.arff";
	while(args.next_is_flag())
	{
		if(args.if_pop("-seed"))
			nSeed = args.pop_uint();
		else if(args.if_pop("-start"))
		{
			x = args.pop_uint();
			y = args.pop_uint();
		}
		else if(args.if_pop("-obsfile"))
			obsFile = args.pop_string();
		else if(args.if_pop("-actionfile"))
			actionFile = args.pop_string();
		else
			ThrowError("Invalid option: ", args.peek());
	}

	// Generate the dataset
	GRand prng(nSeed);
	GData dataObs(pData->relation());
	sp_relation pControlRel = new GUniformRelation(1, 4);
	GData dataControl(pControlRel);
	GData dataState(2);
	int colCount = pData->cols();
	while(dataObs.rows() < samples)
	{
		double* pRow = dataState.newRow();
		pRow[0] = x;
		pRow[1] = y;
		GVec::copy(dataObs.newRow(), pData->row(width * y + x), colCount);
		int action = 0;
		while(true)
		{
			action = (int)prng.next(4);
			if(action == 0 && x > 0)
			{
				x--;
				break;
			}
			else if(action == 1 && x < width - 1)
			{
				x++;
				break;
			}
			else if(action == 2 && y > 0)
			{
				y--;
				break;
			}
			else if(action == 3 && y < height - 1)
			{
				y++;
				break;
			}
		}
		dataControl.newRow()[0] = action;
	}
	dataControl.saveArff(actionFile.c_str());
	dataObs.saveArff(obsFile.c_str());
	dataState.print(cout);
}

class MyHtmlParser1 : public GHtml
{
protected:
	GVocabulary* m_pVocab;

public:
	MyHtmlParser1(GVocabulary* pVocab, const char* pDoc, int nSize)
	: GHtml(pDoc, nSize), m_pVocab(pVocab)
	{
	}

	virtual ~MyHtmlParser1() {}

	virtual void onTextChunk(const char* pChunk, int chunkSize)
	{
		m_pVocab->addWordsFromTextBlock(pChunk, chunkSize);
	}
};

class MyHtmlParser2 : public GHtml
{
protected:
	GSparseMatrix* m_pSM;
	int m_row;
	GVocabulary* m_pVocab;
	bool m_binary;

public:
	MyHtmlParser2(const char* pDoc, int nSize, GSparseMatrix* pSM, int row, GVocabulary* pVocab, bool binary)
	: GHtml(pDoc, nSize), m_pSM(pSM), m_row(row), m_pVocab(pVocab), m_binary(binary)
	{
	}

	virtual ~MyHtmlParser2() {}

	virtual void onTextChunk(const char* pChunk, int chunkSize)
	{
		GWordIterator it(pChunk, chunkSize);
		const char* pWord;
		int wordLen;
		while(true)
		{
			if(!it.next(&pWord, &wordLen))
				break;
			int col = m_pVocab->wordIndex(pWord, wordLen);
			if(col >= 0)
			{
				if(m_binary)
					m_pSM->set(m_row, col, 1.0);
				else
					m_pSM->set(m_row, col, m_pSM->get(m_row, col) + m_pVocab->weight(col));
			}
		}
	}
};

void addWordsToVocabFromHtmlFile(GVocabulary* pVocab, const char* szFilename)
{
	size_t len;
	char* pFile = GFile::loadFile(szFilename, &len);
	ArrayHolder<char> hFile(pFile);
	pVocab->newDoc();
	MyHtmlParser1 parser(pVocab, pFile, len);
	while(true)
	{
		if(!parser.parseSomeMore())
			break;
	}
}

void makeHtmlFileVector(GSparseMatrix* pSM, int clss, int row, GVocabulary* pVocab, const char* szFilename, bool binary)
{
	size_t len;
	char* pFile = GFile::loadFile(szFilename, &len);
	ArrayHolder<char> hFile(pFile);
	MyHtmlParser2 parser(pFile, len, pSM, row, pVocab, binary);
	while(true)
	{
		if(!parser.parseSomeMore())
			break;
	}
	if(clss > 0)
		pSM->set(row, pSM->cols() - 1, (double)clss);
}

void addWordsToVocabFromTextFile(GVocabulary* pVocab, const char* szFilename)
{
	size_t len;
	char* pFile = GFile::loadFile(szFilename, &len);
	ArrayHolder<char> hFile(pFile);
	pVocab->newDoc();
	pVocab->addWordsFromTextBlock(pFile, len);
}

void makeTextFileVector(GSparseMatrix* pSM, int clss, int row, GVocabulary* pVocab, const char* szFilename, bool binary)
{
	size_t len;
	char* pFile = GFile::loadFile(szFilename, &len);
	ArrayHolder<char> hFile(pFile);
	GWordIterator it(pFile, len);
	const char* pWord;
	int wordLen;
	while(true)
	{
		if(!it.next(&pWord, &wordLen))
			break;
		int col = pVocab->wordIndex(pWord, wordLen);
		if(col >= 0)
		{
			if(binary)
				pSM->set(row, col, 1.0);
			else
				pSM->set(row, col, pSM->get(row, col) + pVocab->weight(col));
		}
	}
	if(clss > 0)
		pSM->set(row, pSM->cols() - 1, (double)clss);
}

void docsToSparseMatrix(GArgReader& args)
{
	// Parse options
	bool useStemmer = true;
	bool binary = false;
	string outFilename = "docs.sparse";
	string vocabFile = "";
	while(args.next_is_flag())
	{
		if(args.if_pop("-nostem"))
			useStemmer = false;
		else if(args.if_pop("-binary"))
			binary = true;
		else if(args.if_pop("-outfile"))
			outFilename = args.pop_string();
		else if(args.if_pop("-vocabfile"))
			vocabFile = args.pop_string();
		else
			ThrowError("Invalid option: ", args.peek());
	}

	// Parse the vocabulary
	GVocabulary vocab(useStemmer);
	vocab.addTypicalStopWords();
	vector<string> folders;
	while(args.size() > 0)
	{
		const char* szFolder = args.pop_string();
		folders.push_back(szFolder);
		char cwd[300];
		if(!getcwd(cwd, 300))
			ThrowError("Failed to get cwd");
		if(chdir(szFolder) != 0)
			ThrowError("Failed to change directory to: ", szFolder, ", from: ", cwd);
		GDirList dl(false/*recurseSubDirs*/, true/*reportFiles*/, false/*reportDirs*/, false/*reportPaths*/);
		while(true)
		{
			const char* filename = dl.GetNext();
			if(!filename)
				break;
			PathData pd;
			GFile::parsePath(filename, &pd);
			if(_stricmp(filename + pd.extStart, ".txt") == 0)
				addWordsToVocabFromTextFile(&vocab, filename);
			else if(_stricmp(filename + pd.extStart, ".html") == 0 || _stricmp(filename + pd.extStart, ".htm") == 0)
				addWordsToVocabFromHtmlFile(&vocab, filename);
			else
				printf("Skipping file: %s. (Only .txt and .html is supported.)\n", filename);
		}
		if(chdir(cwd) != 0)
			ThrowError("failed to change dir");
	}
	printf("-----\n");

	// Make the sparse matrix
	int colCount = vocab.wordCount();
	if(folders.size() > 1)
		colCount++; // add a column for the class label
	GSparseMatrix sm(vocab.docCount(), colCount);
	size_t row = 0;
	for(int clss = 0; clss < (int)folders.size(); clss++)
	{
		const char* szFolder = folders[clss].c_str();
		char cwd[300];
		if(!getcwd(cwd, 300))
			ThrowError("Failed to get cwd");
		if(chdir(szFolder) != 0)
			ThrowError("Failed to change directory to: ", szFolder, ", from: ", cwd);
		GDirList dl(false/*recurseSubDirs*/, true/*reportFiles*/, false/*reportDirs*/, false/*reportPaths*/);
		while(true)
		{
			const char* filename = dl.GetNext();
			if(!filename)
				break;
			PathData pd;
			GFile::parsePath(filename, &pd);
			if(_stricmp(filename + pd.extStart, ".txt") == 0)
			{
				printf("%d) %s\n", (int)row, filename);
				makeTextFileVector(&sm, clss, row++, &vocab, filename, binary);
			}
			else if(_stricmp(filename + pd.extStart, ".html") == 0 || _stricmp(filename + pd.extStart, ".htm") == 0)
			{
				printf("%d) %s\n", (int)row, filename);
				makeHtmlFileVector(&sm, clss, row++, &vocab, filename, binary);
			}
		}
		if(chdir(cwd) != 0)
			ThrowError("Failed to change dir");
	}

	// Save the vocab file
	if(vocabFile.length() > 0)
	{
		FILE* pFile = fopen(vocabFile.c_str(), "w");
		FileHolder hFile(pFile);
		for(int i = 0; i < vocab.wordCount(); i++)
		{
			const char* szWord = vocab.stats(i).m_szWord;
			fprintf(pFile, "%d) %s\n", i, szWord);
		}
	}

	sm.save(outFilename.c_str());
}

void vectorToImage(GImage* pImage, const double* pVec, int wid, int hgt)
{
	pImage->setSize(wid, hgt);
	unsigned int* pPix = pImage->pixels();
	for(int y = 0; y < hgt; y++)
	{
		for(int x = 0; x < wid; x++)
		{
			int r = ClipChan((int)*pVec++);
			int g = ClipChan((int)*pVec++);
			int b = ClipChan((int)*pVec++);
			*pPix = gARGB(0xff, r, g, b);
			pPix++;
		}
	}
}

void vectorToImage(GArgReader& args)
{
	GData* pData = GData::loadArff(args.pop_string());
	Holder<GData> hData(pData);
	size_t r = args.pop_uint();
	if(r >= pData->rows())
		ThrowError("row index out of range");
	int wid = args.pop_uint();
	int channels = 3;
	int hgt = pData->cols() / (wid * channels);
	if((wid * hgt * channels) != pData->cols())
		ThrowError("Invalid dimensions");
	double* pVec = pData->row(r);
	GImage image;
	vectorToImage(&image, pVec, wid, hgt);
	image.savePng("image.png");
}

void dataToFrames(GArgReader& args)
{
	GData* pData = GData::loadArff(args.pop_string());
	Holder<GData> hData(pData);
	int wid = args.pop_uint();
	int channels = 3;
	int hgt = pData->cols() / (wid * channels);
	if((wid * hgt * channels) != pData->cols())
		ThrowError("Invalid dimensions");
	GImage image;
	GImage master;
	master.setSize(wid, hgt * pData->rows());
	for(unsigned int i = 0; i < pData->rows(); i++)
	{
		vectorToImage(&image, pData->row(i), wid, hgt);
		GRect r(0, 0, wid, hgt);
		master.blit(0, i * hgt, &image, &r);
	}
	master.savePng("frames.png");
}

void sceneRobotSimulationGrid(GArgReader& args)
{
	const char* sceneFilename = args.pop_string();
	int strafeLen = args.pop_uint();
	int zoomLen = args.pop_uint();
	int cameraWid = args.pop_uint();
	int cameraHgt = args.pop_uint();
	GImage scene;
	scene.loadPng(sceneFilename);
	if(scene.height() > scene.width() * cameraHgt / cameraWid)
		ThrowError("Expected a panoramic (wide) scene");
	double maxWid = scene.height() * cameraWid / cameraHgt;
	double maxStrafeStride = (scene.width() - maxWid) / strafeLen;

	GImage frame;
	frame.setSize(cameraWid, cameraHgt);

	GImage master;
	master.setSize(cameraWid * (strafeLen + 1), cameraHgt * (zoomLen + 1));
	for(double y = 0; y <= zoomLen; y++)
	{
		for(double x = 0; x <= strafeLen; x++)
		{
			if(y > 0.3 * zoomLen && x > 0.3 * strafeLen && x < 0.7 * strafeLen)
				continue;
			//double xx = x / strafeLen;
			double yy = y / zoomLen;
			double z = yy;
			double h = z * cameraHgt + (1.0 - z) * scene.height();
			double w = h * cameraWid / cameraHgt;
			double strafeStride = z + (1.0 - z) * maxStrafeStride;
			double left = ((double)scene.width() - w) / 2.0 + (x - strafeLen / 2.0) * strafeStride;
			double top = (scene.height() - h) / 2.0;
			GDoubleRect src(left, top, w, h);
			GDoubleRect dest(0, 0, cameraWid, cameraHgt);
			frame.blitStretchInterpolate(&dest, &src, &scene);
			GRect r(0, 0, cameraWid, cameraHgt);
			master.blit((int)(cameraWid * x), (int)(cameraHgt * y), &frame, &r);
		}
	}
	master.savePng("frames.png");
}

void sceneRobotSimulationPath(GArgReader& args)
{
	const char* sceneFilename = args.pop_string();
	unsigned int frames = args.pop_uint();

	// Parse options
	unsigned int seed = getpid() * (unsigned int)time(NULL);
	double transitionNoise = 0.0;
	double observationNoise = 0.0;
	int cameraWid = 64;
	int cameraHgt = 48;
	int strafeLen = 30;
	int zoomLen = 20;
	string stateFilename = "state.arff";
	string actFilename = "act.arff";
	string obsFilename = "obs.arff";
	while(args.next_is_flag())
	{
		if(args.if_pop("-noise"))
		{
			transitionNoise = args.pop_double();
			observationNoise = args.pop_double();
		}
		else if(args.if_pop("-framesize"))
		{
			cameraWid = args.pop_uint();
			cameraHgt = args.pop_uint();
		}
		else if(args.if_pop("-gridsize"))
		{
			strafeLen = args.pop_uint();
			zoomLen = args.pop_uint();
		}
		else if(args.if_pop("-seed"))
			seed = args.pop_uint();
		else if(args.if_pop("-out"))
		{
			stateFilename = args.pop_string();
			actFilename = args.pop_string();
			obsFilename = args.pop_string();
		}
		else
			ThrowError("Invalid option: ", args.peek());
	}

	GImage scene;
	scene.loadPng(sceneFilename);
	if(scene.height() > scene.width() * cameraHgt / cameraWid)
		ThrowError("Expected a panoramic (wide) scene");
	double maxWid = scene.height() * cameraWid / cameraHgt;
	double maxStrafeStride = (scene.width() - maxWid) / strafeLen;
	GRand prng(seed);

	GImage frame;
	frame.setSize(cameraWid, cameraHgt);
	double x = strafeLen / 2.0;
	double y = 0.0;//zoomLen / 2.0;
	sp_relation pRel = new GUniformRelation(1, 4);
	GData state(2);
	GData act(pRel);
	GData obs(cameraWid * cameraHgt * 3);
	for(unsigned int i = 0; i < frames; i++)
	{
		// Make the state
		double* pVec = state.newRow();
		pVec[0] = x;
		pVec[1] = y;

		// Make the frame
		double yy = y / zoomLen;
		double z = yy;
		double h = z * cameraHgt + (1.0 - z) * scene.height();
		double w = h * cameraWid / cameraHgt;
		double strafeStride = z + (1.0 - z) * maxStrafeStride;
		double left = ((double)scene.width() - w) / 2.0 + (x - strafeLen / 2.0) * strafeStride;
		double top = (scene.height() - h) / 2.0;
		GDoubleRect src(left, top, w, h);
		GDoubleRect dest(0, 0, cameraWid, cameraHgt);
		frame.blitStretchInterpolate(&dest, &src, &scene);

		// Convert to an observation vector
		pVec = obs.newRow();
		unsigned int* pix = frame.pixels();
		for(unsigned int yy = 0; yy < frame.height(); yy++)
		{
			for(unsigned int xx = 0; xx < frame.width(); xx++)
			{
				*(pVec++) = ClipChan((int)(observationNoise * 255 * prng.normal()) + gRed(*pix));
				*(pVec++) = ClipChan((int)(observationNoise * 255 * prng.normal()) + gGreen(*pix));
				*(pVec++) = ClipChan((int)(observationNoise * 255 * prng.normal()) + gBlue(*pix));
				pix++;
			}
		}

		// Do a random action
		double oldx = x;
		double oldy = y;
		int a;
		while(true)
		{
			a = (int)prng.next(4);
			if(a == 0)
				x -= 1.0;
			else if(a == 1)
				x += 1.0;
			else if(a == 2)
				y -= 1.0;
			else if(a == 3)
				y += 1.0;
			x += transitionNoise * prng.normal();
			y += transitionNoise * prng.normal();
			if(x < 0 || y < 0 || x > strafeLen || y > zoomLen)
			{
				x = oldx;
				y = oldy;
				continue;
			}
			if(y > 0.3 * zoomLen && x > 0.3 * strafeLen && x < 0.7 * strafeLen)
			{
				x = oldx;
				y = oldy;
				continue;
			}
			break;
		}
		act.newRow()[0] = a;
	}
	state.saveArff(stateFilename.c_str());
	obs.saveArff(obsFilename.c_str());
	act.saveArff(actFilename.c_str());
}

void mechanicalRabbit(GArgReader& args)
{
	const char* sceneFilename = args.pop_string();
	GTwtDoc doc;
	doc.load(args.pop_string());
	GData* pActions = GData::loadArff(args.pop_string());
	Holder<GData> hActions(pActions);

	// Parse the options
	unsigned int seed = getpid() * (unsigned int)time(NULL);
	string contextFilename = "context.arff";
	string stateFilename = "state.arff";
	string plannedObsFilename = "planned_obs.arff";
	string actualObsFilename = "actual_obs.arff";
	string actionsFilename = "actions.arff";
	double transitionNoise = 0.0;
	double observationNoise = 0.0;
	int cameraWid = 64;
	int cameraHgt = 48;
	int strafeLen = 30;
	int zoomLen = 20;
	while(args.next_is_flag())
	{
		if(args.if_pop("-out"))
		{
			contextFilename = args.pop_string();
			stateFilename = args.pop_string();
			plannedObsFilename = args.pop_string();
			actualObsFilename = args.pop_string();
			actionsFilename = args.pop_string();
		}
		else if(args.if_pop("-seed"))
			seed = args.pop_uint();
		else if(args.if_pop("-noise"))
		{
			transitionNoise = args.pop_double();
			observationNoise = args.pop_double();
		}
		else
			ThrowError("Invalid option: ", args.peek());
	}

	// Generate the context sequence
	GRand prng(seed);
	GRecurrentModel rm(doc.root(), &prng);
	GData obsPlanned(rm.obsDims());
	GData context(rm.contextDims());
	GVec::setAll(rm.context(), 0.0, rm.contextDims());
	GVec::copy(context.newRow(), rm.context(), rm.contextDims());
	rm.predict(obsPlanned.newRow());
	for(size_t i = 0; i < pActions->rows(); i++)
	{
		rm.doAction(pActions->row(i));
		GVec::copy(context.newRow(), rm.context(), rm.contextDims());
		rm.predict(obsPlanned.newRow());
	}

	// Reset the model and the system
	GVec::setAll(rm.context(), 0.0, rm.contextDims());
	GImage scene;
	scene.loadPng(sceneFilename);
	if(scene.height() > scene.width() * cameraHgt / cameraWid)
		ThrowError("Expected a panoramic (wide) scene");
	double maxWid = scene.height() * cameraWid / cameraHgt;
	double maxStrafeStride = (scene.width() - maxWid) / strafeLen;
	GImage frame;
	frame.setSize(cameraWid, cameraHgt);
	double x = strafeLen / 2.0;
	double y = 0.0;
	GData obsActual(rm.obsDims());
	GData state(2);
	sp_relation pRel = new GUniformRelation(1, 4);
	GData actionsActual(pRel);

	// Do the mechanical rabbit
	GTEMPBUF(double, tmpContext, rm.contextDims());
	size_t r = 0;
	size_t gap = 1;
	int safety = 0;
	while(r < context.rows() && gap > 0)
	{
		// Store the actual state
		double* pState = state.newRow();
		pState[0] = x;
		pState[1] = y;

		// Make the frame
		double yy = y / zoomLen;
		double z = yy;
		double h = z * cameraHgt + (1.0 - z) * scene.height();
		double w = h * cameraWid / cameraHgt;
		double strafeStride = z + (1.0 - z) * maxStrafeStride;
		double left = ((double)scene.width() - w) / 2.0 + (x - strafeLen / 2.0) * strafeStride;
		double top = (scene.height() - h) / 2.0;
		GDoubleRect src(left, top, w, h);
		GDoubleRect dest(0, 0, cameraWid, cameraHgt);
		frame.blitStretchInterpolate(&dest, &src, &scene);

		// Convert to an observation vector
		double* pVec = obsActual.newRow();
		double* pObsActual = pVec;
		unsigned int* pix = frame.pixels();
		for(unsigned int yy = 0; yy < frame.height(); yy++)
		{
			for(unsigned int xx = 0; xx < frame.width(); xx++)
			{
				*(pVec++) = ClipChan((int)(observationNoise * 255 * prng.normal()) + gRed(*pix));
				*(pVec++) = ClipChan((int)(observationNoise * 255 * prng.normal()) + gGreen(*pix));
				*(pVec++) = ClipChan((int)(observationNoise * 255 * prng.normal()) + gBlue(*pix));
				pix++;
			}
		}

		// Advance the mechanical rabbit
		rm.calibrate(pObsActual);
		if(r + 1 == context.rows())
			gap--;
		safety++;
		while(safety > 3 || (r + 1 < context.rows() && sqrt(GVec::squaredDistance(rm.context(), context.row(r), rm.contextDims())) < gap))
		{
			r++;
			safety = 0;
		}

		// Pick the best action
		int bestAction = -1;
		double bestSquaredDist = 1e308;
		double dAction;
		for(int i = 0; i < 4; i++)
		{
			GVec::copy(tmpContext, rm.context(), rm.contextDims());
			dAction = (double)i;
			rm.doAction(&dAction);
			double squaredDist = GVec::squaredDistance(context.row(r), rm.context(), rm.contextDims());
			if(bestAction < 0 || squaredDist < bestSquaredDist)
			{
				bestAction = i;
				bestSquaredDist = squaredDist;
			}
			GVec::copy(rm.context(), tmpContext, rm.contextDims());
		}

		// Apply the action to the model
		dAction = (double)bestAction;
		rm.doAction(&dAction);

		// Apply the action to the system
		actionsActual.newRow()[0] = bestAction;
		if(bestAction == 0)
			x -= 1.0;
		else if(bestAction == 1)
			x += 1.0;
		else if(bestAction == 2)
			y -= 1.0;
		else if(bestAction == 3)
			y += 1.0;
		x += transitionNoise * prng.normal();
		y += transitionNoise * prng.normal();
	}

	// Save the results
	context.saveArff(contextFilename.c_str());
	obsPlanned.saveArff(plannedObsFilename.c_str());
	state.saveArff(stateFilename.c_str());
	obsActual.saveArff(actualObsFilename.c_str());
	actionsActual.saveArff(actionsFilename.c_str());
}

void ShowBriefUsage(const char* appName)
{
	UsageNode* pUsageTree = makeGenerateUsageTree();
	Holder<UsageNode> hUsageTree(pUsageTree);
	pUsageTree->print(0, 3, 76, false);
	cout << "__________________________________\nTo see the full usage info, enter:\n    " << appName << " usage\n";
	cout.flush();
}

void ShowInvalidSyntaxError(const char* appName)
{
	cout << "\nInvalid syntax.\n\n";
	ShowBriefUsage(appName);
}

void ShowUsage(const char* appName)
{
	UsageNode* pUsageTree = makeGenerateUsageTree();
	Holder<UsageNode> hUsageTree(pUsageTree);
	pUsageTree->print(0, 3, 76, true);
	cout.flush();
}


int DoIt(int argc, char *argv[])
{
	PathData pd;
	GFile::parsePath(argv[0], &pd);
	const char* appName = argv[0] + pd.fileStart;

	GArgReader args(argc, argv);
	args.pop_string();
	int ret = 0;
	if(args.size() >= 1)
	{
		if(args.if_pop("usage")) ShowUsage(appName);
		else if(args.if_pop("crane")) CraneDataset(args);
		else if(args.if_pop("cranepath")) cranePath(args);
		else if(args.if_pop("cube")) cube(args);
		else if(args.if_pop("datatoframes")) dataToFrames(args);
		else if(args.if_pop("docstosparsematrix")) docsToSparseMatrix(args);
		else if(args.if_pop("entwinedspirals")) EntwinedSpirals(args);
		else if(args.if_pop("fishbowl")) fishBowl(args);
		else if(args.if_pop("gridrandomwalk")) gridRandomWalk(args);
		else if(args.if_pop("imagetranslatedovernoise")) ImageTranslatedOverNoise(args);
		else if(args.if_pop("mechanicalrabbit")) mechanicalRabbit(args);
		else if(args.if_pop("noise")) Noise(args);
		else if(args.if_pop("scalerotate")) ScaleAndRotate(args);
		else if(args.if_pop("scenerobotsimulationgrid")) sceneRobotSimulationGrid(args);
		else if(args.if_pop("scenerobotsimulationpath")) sceneRobotSimulationPath(args);
		else if(args.if_pop("scurve")) SCurve(args);
		else if(args.if_pop("selfintersectingribbon")) SelfIntersectingRibbon(args);
		else if(args.if_pop("swissroll")) SwissRoll(args);
		else if(args.if_pop("threecranepath")) threeCranePath(args);
		else if(args.if_pop("vectortoimage")) vectorToImage(args);
		else if(args.if_pop("windowedimage")) WindowedImageData(args);
		else
		{
			ret = 1;
			ShowInvalidSyntaxError(appName);
		}
	}
	else
	{
		ret = 1;
		ShowInvalidSyntaxError(appName);
	}
	return ret;
}

int main(int argc, char *argv[])
{
	int nRet = 0;
	try
	{
		nRet = DoIt(argc, argv);
	}
	catch(const std::exception& e)
	{
		std::cerr << e.what() << "\n";
		nRet = 100;
	}

	return nRet;
}

