// ----------------------------------------------------------------
// The contents of this file are distributed under the CC0 license.
// See http://creativecommons.org/publicdomain/zero/1.0/
// ----------------------------------------------------------------

#include <stdio.h>
#include <stdlib.h>
#include "../GClasses/GApp.h"
#include "../GClasses/GMacros.h"
#include "../GClasses/GData.h"
#include "../GClasses/GDecisionTree.h"
#include "../GClasses/GRand.h"
#include "../GClasses/GFile.h"
#include "../GClasses/GImage.h"
#include "../GClasses/GOptimizer.h"
#include "../GClasses/GHillClimber.h"
#include "../GSup/G3D.h"
#include "../GClasses/GVec.h"
#include "../GClasses/GTwt.h"
#include "../GClasses/GPlot.h"
#include "../GClasses/GFunction.h"
#include "../GClasses/GLearner.h"
#include "../GClasses/GNeighborFinder.h"
#include "../GClasses/GSystemLearner.h"
#include "plotchart.h"
#include "../wizard/usage.h"
#include <time.h>
#include <iostream>
#ifdef WIN32
#	include <direct.h>
#	include <process.h>
#endif
#include <string>
#include <vector>
#include <exception>

using namespace GClasses;
using std::string;
using std::cout;
using std::vector;

void ShowBriefUsage(const char* appName)
{
	UsageNode* pUsageTree = makePlotUsageTree();
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
	UsageNode* pUsageTree = makePlotUsageTree();
	Holder<UsageNode> hUsageTree(pUsageTree);
	pUsageTree->print(0, 3, 76, true);
	cout.flush();
}

GNeighborFinder* instantiateNeighborFinder(GData* pData, GRand* pRand, GArgReader& args)
{
	// Get the algorithm name
	GNeighborFinder* pNF = NULL;
	const char* alg = args.pop_string();

	// Parse the options
	int cutCycleLen = 0;
	while(args.next_is_flag())
	{
		if(args.if_pop("-cyclecut"))
			cutCycleLen = args.pop_uint();
		else
			ThrowError("Invalid neighbor finder option: ", args.peek());
	}

	// Parse required algorithms
	if(_stricmp(alg, "bruteforce") == 0)
	{
		int neighbors = args.pop_uint();
		pNF = new GBruteForceNeighborFinder(pData, 0, neighbors, NULL, true);
	}
	else if(_stricmp(alg, "kdtree") == 0)
	{
		int neighbors = args.pop_uint();
		pNF = new GKdTree(pData, 0, neighbors, NULL, true);
	}
	else if(_stricmp(alg, "manifold") == 0)
	{
		int neighbors = args.pop_uint();
		int tangentSpaceDims = args.pop_uint();
		double alpha = args.pop_double();
		double beta = args.pop_double();
		pNF = new GManifoldNeighborFinder(pData, neighbors, neighbors * 4, tangentSpaceDims, alpha, beta, false, pRand);
	}
	else if(_stricmp(alg, "system") == 0)
	{
		GData* pControlData = GData::loadArff(args.pop_string());
		Holder<GData> hControlData(pControlData);
		if(pControlData->rows() != pData->rows())
			ThrowError("mismatching number of rows");
		int neighbors = args.pop_uint();
		pNF = new GDynamicSystemNeighborFinder(pData, hControlData.release(), true, neighbors, pRand);
	}
	else
		ThrowError("Unrecognized neighbor finding algorithm: ", alg);

	// Apply CycleCut
	if(cutCycleLen > 0)
	{
		GNeighborFinderCacheWrapper* pNF2 = new GNeighborFinderCacheWrapper(pNF, true);
		pNF2->fillCache();
		pNF2->cutShortcuts(cutCycleLen);
		pNF = pNF2;
	}

	return pNF;
}

void PlotBar(GArgReader& args)
{
	// Load the data
	GData* pData = GData::loadArff(args.pop_string());
	Holder<GData> hData(pData);
	GArffRelation* pRel = (GArffRelation*)pData->relation().get();
	if(pRel->size() != 1 || pRel->areContinuous(0, 1))
		ThrowError("Expected exactly one continuous attribute");
	double* values = new double[pData->rows()];
	ArrayHolder<double> hValues(values);
	for(size_t i = 0; i < pData->rows(); i++)
		values[i] = pData->row(i)[0];

	// Parse options
	bool bLog = false;
	string filename = "plot.png";
	while(args.next_is_flag())
	{
		if(args.if_pop("-out"))
			filename = args.pop_string();
		else if(args.if_pop("-log"))
			bLog = true;
		else
			ThrowError("Invalid option: ", args.peek());
	}

	// Make the chart
	GRand prng(0);
	int minIndex = GVec::indexOfMin(values, pData->rows(), &prng);
	int maxIndex = GVec::indexOfMax(values, pData->rows(), &prng);
	double dMin = values[minIndex];
	double dMax = MAX(1e-12, values[maxIndex]);
	if(dMin > 0 && !bLog && (dMax - dMin) / dMax > 0.05)
		dMin = 0;
	GImage image;
	image.setSize(800, 800);
	image.clear(0xffffffff);
	double xmin = -0.5;
	double ymin = bLog ? log(dMin * 0.7) : dMin - 0.1 * (dMax - dMin);
	double xmax = pData->rows();
	double ymax = bLog ? log(dMax + 0.5 * (dMax - dMin)) : dMax + 0.1 * (dMax - dMin);
	GPlotWindow pw(&image, xmin, ymin, xmax, ymax);
	pw.gridLines(0, (bLog ? -1 : 30), 0xff808080);
	for(size_t i = 0; i < pData->rows(); i++)
	{
		int x1, y1, x2, y2;
			pw.windowToView(i, 0, &x1, &y1);
		if(bLog)
		{
			pw.windowToView(0.5 + i, log(values[i]), &x2, &y2);
			image.boxFill(x1, y2, x2 - x1, MAX(0, (int)image.height() - y2), gAHSV(0xff, (float)i / pData->rows(), 1.0f, 0.5f));
		}
		else
		{
			pw.windowToView(0.5 + i, values[i], &x2, &y2);
			if(y2 < y1)
				std::swap(y1, y2);
			image.boxFill(x1, y1, x2 - x1, y2 - y1, gAHSV(0xff, (float)i / pData->rows(), 1.0f, 0.5f));
		}
	}
	GImage* pLabeledImage = pw.labelAxes(0, (bLog ? -1 : 30), 5/*precision*/, 1/*size*/, 0xff000000/*color*/, 0.0/*angle*/);
	Holder<GImage> hLabeledImage(pLabeledImage);
	pLabeledImage->savePng(filename.c_str());
	cout << "Chart saved to " << filename.c_str() << ".\n";
}

class BigOCritic : public GTargetFunction
{
protected:
	GData* m_pData;
	int m_attr;

public:
	BigOCritic(GData* pData, int attr)
	: GTargetFunction(3), m_pData(pData), m_attr(attr)
	{
	}

	virtual ~BigOCritic()
	{
	}

	virtual bool isStable() { return true; }
	virtual bool isConstrained() { return false; }

	virtual void initVector(double* pVector)
	{
		pVector[0] = 1.0;
		pVector[1] = 1.0;
		pVector[2] = 0.0;
	}

	virtual double computeError(const double* pVector)
	{
		double err = 0;
		for(size_t i = 0; i < m_pData->rows(); i++)
		{
			double* pPat = m_pData->row(i);
			if(pPat[0] != UNKNOWN_REAL_VALUE && pPat[m_attr] != UNKNOWN_REAL_VALUE)
			{
				double d = pVector[0] * (pow(pPat[0], pVector[1]) + pVector[2]) - pPat[m_attr];
				err += (d * d);
			}
		}
		return err;
	}
};

void EstimateBigO(GArgReader& args)
{
	// Load the data
	GData* pData = GData::loadArff(args.pop_string());
	Holder<GData> hData(pData);
	GArffRelation* pRel = (GArffRelation*)pData->relation().get();
	if(pRel->size() < 2)
		ThrowError("Expected at least two attributes");
	if(!pRel->areContinuous(0, pRel->size()))
		ThrowError("Expected all continuous attributes");

	// Regress t=an^b+c for each algorithm
	cout.precision(8);
	for(int i = 1; i < pRel->size(); i++)
	{
		BigOCritic critic(pData, i);
		GMomentumGreedySearch search(&critic);
		search.searchUntil(500, 50, 0.0001);
		double* pVec = search.currentVector();
		cout << pRel->attrName(i) << ": t=" << pVec[0] << " * (n^" << pVec[1] << " + " << pVec[2] << ")\n";
	}
}

void PlotEquation(GArgReader& args)
{
	// Parse options
	string filename = "plot.png";
	int width = 1024;
	int height = 1024;
	double xmin = -10;
	double ymin = -10;
	double xmax = 10;
	double ymax = 10;
	while(args.next_is_flag())
	{
		if(args.if_pop("-out"))
			filename = args.pop_string();
		else if(args.if_pop("-size"))
		{
			width = args.pop_uint();
			height = args.pop_uint();
		}
		else if(args.if_pop("-range"))
		{
			xmin = args.pop_double();
			ymin = args.pop_double();
			xmax = args.pop_double();
			ymax = args.pop_double();
		}
		else
			ThrowError("Invalid option: ", args.peek());
	}

	// Accumulate the expression
	string expr;
	while(args.size() > 0)
		expr += args.pop_string();

	// Parse the expression
	GFunctionParser mfp(expr.c_str());

	// Make the chart
	GImage image;
	image.setSize(width, height);
	image.clear(0xffffffff);
	GPlotWindow pw(&image, xmin, ymin, xmax, ymax);
	pw.gridLines(30, 30, 0xffa0a0a0);

	// Plot all the functions
	char szFuncName[32];
	unsigned int colors[6];
	colors[0] = 0xff000080;
	colors[1] = 0xff800000;
	colors[2] = 0xff008000;
	colors[3] = 0xff808000;
	colors[4] = 0xff800080;
	colors[5] = 0xff008080;
	for(int i = 1; true; i++)
	{
		// Find the function
		sprintf(szFuncName, "f%d", i);
		GFunction* pFunc = mfp.getFunctionNoThrow(szFuncName);
		if(!pFunc)
		{
			if(i == 1)
				ThrowError("There is no function named \"f1\". Nothing to plot.");
			break;
		}
		if(pFunc->m_expectedParams != 1)
			ThrowError("The function ", szFuncName, " takes ", gformat(pFunc->m_expectedParams), " parameters. Expected a function with 1 parameter");

		// Plot it
		unsigned int col = colors[i % 6];
		double dx = pw.pixelWidth();
		vector<double> params;
		double x = xmin;
		params.push_back(x);
		double y = pFunc->call(params);
		while(x <= xmax)
		{
			double xPrev = x;
			double yPrev = y;
			x += dx;
			params[0] = x;
			y = pFunc->call(params);
			if(y > -1e100 && y < 1e100 && yPrev > -1e100 && yPrev < 1e100)
				pw.line(xPrev, yPrev, x, y, col);
		}
	}

	GImage* pLabeledImage = pw.labelAxes(30, 30, 5/*precision*/, 1/*size*/, 0xff000000/*color*/, 0.0/*angle*/);
	Holder<GImage> hLabeledImage(pLabeledImage);
	pLabeledImage->savePng(filename.c_str());
	cout << "Plot saved to " << filename.c_str() << ".\n";
}

void PlotScatter(GArgReader& args)
{
	// Load the data
	GData* pData = GData::loadArff(args.pop_string());
	Holder<GData> hData(pData);

	// Parse options
	GRand prng(0);
	PlotChartMaker pcm(pData->relation(), pData);
	string filename = "plot.png";
	bool horizAxisLabels = true;
	bool vertAxisLabels = true;
	bool showLines = false;
	while(args.next_is_flag())
	{
		if(args.if_pop("-size"))
		{
			int wid = args.pop_uint();
			int hgt = args.pop_uint();
			pcm.SetSize(wid, hgt);
		}
		else if(args.if_pop("-lines"))
			showLines = true;
		else if(args.if_pop("-logx"))
			pcm.SetLogX();
		else if(args.if_pop("-logy"))
			pcm.SetLogY();
		else if(args.if_pop("-nohorizaxislabels"))
			horizAxisLabels = false;
		else if(args.if_pop("-novertaxislabels"))
			vertAxisLabels = false;
		else if(args.if_pop("-pointradius"))
			pcm.SetPointRadius((float)args.pop_double());
		else if(args.if_pop("-textsize"))
			pcm.setTextSize((float)args.pop_double());
		else if(args.if_pop("-linethickness"))
			pcm.SetLineThickness((float)args.pop_double());
		else if(args.if_pop("-maxgridlines"))
		{
			int h = args.pop_uint();
			int v = args.pop_uint();
			pcm.setMaxGridLines(h, v);
		}
		else if(args.if_pop("-mesh"))
			pcm.SetMeshRowSize(args.pop_uint());
		else if(args.if_pop("-aspect"))
			pcm.setAspect();
		else if(args.if_pop("-range"))
		{
			double xmin = args.pop_double();
			double ymin = args.pop_double();
			double xmax = args.pop_double();
			double ymax = args.pop_double();
			pcm.SetCustomRange(xmin, ymin, xmax, ymax);
		}
		else if(args.if_pop("-chartcolors"))
		{
			unsigned int cBackground = hexToRgb(args.pop_string());
			unsigned int cText = hexToRgb(args.pop_string());
			unsigned int cGrid = hexToRgb(args.pop_string());
			pcm.SetChartColors(cBackground, cText, cGrid);
		}
		else if(args.if_pop("-linecolors"))
		{
			unsigned int c1 = hexToRgb(args.pop_string());
			unsigned int c2 = hexToRgb(args.pop_string());
			unsigned int c3 = hexToRgb(args.pop_string());
			unsigned int c4 = hexToRgb(args.pop_string());
			pcm.SetPlotColors(c1, c2, c3, c4);
		}
		else if(args.if_pop("-spectrum"))
			pcm.UseSpectrumColors();
		else if(args.if_pop("-specmod"))
			pcm.UseSpectrumColors(args.pop_uint());
		else if(args.if_pop("-out"))
			filename = args.pop_string();
		else if(args.if_pop("-neighbors"))
			pcm.showNeighbors(instantiateNeighborFinder(pData, &prng, args));
		else
			ThrowError("Invalid option: ", args.peek());
	}
	pcm.ShowAxisLabels(horizAxisLabels, vertAxisLabels);
	if(!showLines)
		pcm.noLines();

	// Make the chart
	GImage* pImage = pcm.MakeChart();
	Holder<GImage> hImage(pImage);
	pImage->savePng(filename.c_str());
	cout << "Plot saved to " << filename.c_str() << ".\n";
}

void makeHistogram(GArgReader& args)
{
	// Load the data
	GData* pData = GData::loadArff(args.pop_string());
	Holder<GData> hData(pData);

	// Parse options
	int wid = 800;
	int hgt = 800;
	int attr = 0;
	string filename = "plot.png";
	double xmin = UNKNOWN_REAL_VALUE;
	double xmax = UNKNOWN_REAL_VALUE;
	double ymax = UNKNOWN_REAL_VALUE;
	while(args.next_is_flag())
	{
		if(args.if_pop("-attr"))
			attr = args.pop_uint();
		else if(args.if_pop("-size"))
		{
			wid = args.pop_uint();
			hgt = args.pop_uint();
		}
		else if(args.if_pop("-out"))
			filename = args.pop_string();
		else if(args.if_pop("-range"))
		{
			xmin = args.pop_double();
			xmax = args.pop_double();
			ymax = args.pop_double();
		}
		else
			ThrowError("Invalid option: ", args.peek());
	}
	if(attr < 0 || attr >= pData->relation()->size())
		ThrowError("attr out of range");

	// Make the histogram
	GImage image;
	image.setSize(wid, hgt);
	image.clear(0xffffffff);
	if(pData->relation()->valueCount(attr) == 0)
	{
		// Make the histogram
		double mean = pData->mean(attr);
		double median = pData->median(attr);
		double dev = sqrt(pData->variance(attr, mean));
		double min = (xmin == UNKNOWN_REAL_VALUE ? median - 4 * dev : xmin);
		double range = (xmax == UNKNOWN_REAL_VALUE ? 8 * dev : xmax - xmin);
		double height = (ymax == UNKNOWN_REAL_VALUE ? 0.6 : ymax);

		size_t buckets = MIN((size_t)image.width(), (size_t)floor(sqrt((double)pData->rows())));
		size_t* hist = new size_t[buckets];
		ArrayHolder<size_t> hHist(hist);
		for(size_t i = 0; i < buckets; i++)
			hist[i] = 0;
		for(size_t i = 0; i < pData->rows(); i++)
		{
			double d = pData->row(i)[attr];
			size_t index = (size_t)floor((d - min) / range * buckets);
			if(index < buckets)
				hist[index]++;
		}

		// Plot it
		size_t max = 0;
		for(size_t i = 1; i < buckets; i++)
		{
			if(hist[i] > hist[max])
				max = i;
		}
		for(int i = 0; i < (int)image.width(); i++)
		{
			double d = (double)i * buckets / wid;
			int b1 = MAX(0, (int)floor(d));
			int b2 = MIN((int)buckets - 1, b1 + 1);
			d -= floor(d);
			d = GMath::softStep(d, 2.0);
			int h = (int)(((1.0 - d) * hist[b1] + d * hist[b2]) * image.height() / (height * pData->rows() * range / buckets));
			image.line(i, image.height(), i, image.height() - h, 0xff00a080);
		}

		// Draw the grid
		GPlotWindow pw(&image, min, 0, min + range, height);
		pw.gridLines(40, 40, 0xff808080);

		// Draw the labels
		GImage* pLabeledImage = pw.labelAxes(30, 30, 5/*precision*/, 1/*size*/, 0xff000000/*color*/, 0.0/*angle*/);
		Holder<GImage> hLabeledImage(pLabeledImage);

		// Save the image
		pLabeledImage->savePng(filename.c_str());
		cout << "Histogram saved to " << filename.c_str() << ".\n";
	}
	else
	{
		int buckets = pData->relation()->valueCount(attr);
		GTEMPBUF(double, hist, buckets);
		GVec::setAll(hist, 0.0, buckets);
		for(size_t i = 0; i < pData->rows(); i++)
		{
			int b = (int)pData->row(i)[attr];
			if(b >= 0 && b < buckets)
				hist[b]++;
		}

		// Plot it
		int max = 0;
		for(int i = 1; i < buckets; i++)
		{
			if(hist[i] > hist[max])
				max = i;
		}
		for(int i = 0; i < (int)image.width(); i++)
		{
			int b = i * buckets / image.width();
			int h = (int)(hist[b] * image.height() / hist[max]);
			image.line(i, image.height(), i, image.height() - h, (((b & 1) == 0) ? 0xff400000 : 0xff008040));
		}
		image.savePng(filename.c_str());
		cout << "Histogram saved to " << filename.c_str() << ".\n";
	}
}

void MakeAttributeSummaryGraph(GRelation* pRelation, GData* pData, GImage* pImage, int attr)
{
	if(pRelation->valueCount(attr) == 0)
	{
		int buckets = MAX(2, (int)floor(sqrt((double)pData->rows()) + 0.5));
		GTEMPBUF(double, hist, buckets);
		GVec::setAll(hist, 0.0, buckets);
		double min, range;
		pData->minAndRange(attr, &min, &range);
		for(size_t i = 0; i < pData->rows(); i++)
		{
			double val = pData->row(i)[attr];
			int b = (int)floor((val - min) * buckets / range);
			if(b >= 0 && b < buckets)
				hist[b]++;
		}

		// Plot it
		pImage->clear(0xffffffff);
		int max = 0;
		for(int i = 1; i < buckets; i++)
		{
			if(hist[i] > hist[max])
				max = i;
		}
		for(int i = 0; i < (int)pImage->width(); i++)
		{
			double d = (double)i * buckets / pImage->width();
			int b1 = MAX(0, (int)floor(d));
			int b2 = MIN(buckets - 1, b1 + 1);
			d -= floor(d);
			d = GMath::softStep(d, 2.0);
			int h = (int)(((1.0 - d) * hist[b1] + d * hist[b2]) * pImage->height() / hist[max]);
			pImage->line(i, pImage->height(), i, pImage->height() - h, 0xff000080);
		}
	}
	else
	{
		int buckets = pRelation->valueCount(attr);
		GTEMPBUF(double, hist, buckets);
		GVec::setAll(hist, 0.0, buckets);
		for(size_t i = 0; i < pData->rows(); i++)
		{
			int b = (int)pData->row(i)[attr];
			if(b >= 0 && b < buckets)
				hist[b]++;
		}

		// Plot it
		pImage->clear(0xffffffff);
		int max = 0;
		for(int i = 1; i < buckets; i++)
		{
			if(hist[i] > hist[max])
				max = i;
		}
		for(int i = 0; i < (int)pImage->width(); i++)
		{
			int b = i * buckets / pImage->width();
			int h = (int)(hist[b] * pImage->height() / hist[max]);
			pImage->line(i, pImage->height(), i, pImage->height() - h, (((b & 1) == 0) ? 0xff400000 : 0xff008040));
		}
	}
}

void MakeCorrelationGraph(GRelation* pRelation, GData* pData, GImage* pImage, int attrx, int attry, double jitter, GRand* pRand)
{
	pImage->clear(0xffffffff);
	double xmin, ymin, xmax, ymax;
	bool bothNominal = true;
	if(pRelation->valueCount(attrx) == 0) //Continuous x attribute
	{
		pData->minAndRange(attrx, &xmin, &xmax);
		xmax += xmin;
		bothNominal = false;
	}
	else //Discrete x attribute
	{
		xmin = -0.5;
		xmax = pRelation->valueCount(attrx) - 0.5;
	}
	if(pRelation->valueCount(attry) == 0) //Continuous y attribute
	{
		pData->minAndRange(attry, &ymin, &ymax);
		ymax += ymin;
		bothNominal = false;
	}
	else //Discrete y atrribute
	{
		ymin = -0.5;
		ymax = pRelation->valueCount(attry) - 0.5;
	}
	if(bothNominal)
	{
		GPlotWindow pw(pImage, 0.0, 0.0, 1.0, 1.0);
		double left = 0.0;
		double right = 0.0;
		size_t tot = pData->rows();
		for(int i = 0; i < pRelation->valueCount(attrx); i++)
		{
			GData tmp(pData->relation());
			pData->splitByDiscreteValue(&tmp, attrx, i);
			right += (double)tmp.rows() / tot;
			double bot = 0.0;
			double top = 0.0;
			for(int j = 0; j < pRelation->valueCount(attry); j++)
			{
				top += (double)tmp.countValue(attry, j) / tmp.rows();
				int l, b, r, t;
				pw.windowToView(left, bot, &l, &b);
				pw.windowToView(right, top, &r, &t);
				pImage->boxFill(l, t, r - l, b - t, gAHSV(0xff, 0.9f * j / pRelation->valueCount(attry), 0.6f + ((i & 1) ? 0.0f : 0.4f), 0.4f + ((i & 1) ? 0.4f : 0.0f)));
				bot = top;
			}
			pData->mergeVert(&tmp);
			left = right;
		}
	}
	else
	{
		GPlotWindow pw(pImage, xmin, ymin, xmax, ymax);
		int step = MAX(1, (int)(pData->rows() / (pImage->width() * pImage->height())));
		for(size_t i = 0; i < pData->rows(); i += step)
		{
			double* pPat = pData->row(i);
			pw.point(pPat[attrx] + pRand->normal() * jitter * (xmax - xmin), pPat[attry] + pRand->normal() * jitter * (ymax - ymin), 0xff000080);
		}
	}
}

void MakeCorrelationLabel(GArffRelation* pRelation, GData* pData, GImage* pImage, int attr, unsigned int bgCol)
{
	pImage->clear(bgCol);
	if(pRelation->valueCount(attr) == 0)
	{
		pImage->text(pRelation->attrName(attr), 0, 0, 1.0f, 0xff400000);
		double min, max;
		pData->minAndRange(attr, &min, &max);

		for(int i = 0; i < 2; i++)
		{
			GImage image2;
			image2.setSize(pImage->width() - 16, 16);
			image2.clear(0);
			int xx = 0;
			char szValue[64];
			sprintf(szValue, "%.4lg", (i == 0 ? min : max));
			int eatspace = pImage->width() - 16 - GImage::measureTextWidth(szValue, 1.0f);
			xx += eatspace;
			image2.text(szValue, xx, 0, 1.0f, 0xff400000);
			GImage image3;
			image3.rotateClockwise90(&image2);
			GRect r3(0, 0, image3.width(), image3.height());
			pImage->blitAlpha((i == 0 ? 0 : pImage->width() - 16), 16, &image3, &r3);
		}
	}
	else
	{
		pImage->text(pRelation->attrName(attr), 0, 0, 1.0f, 0xff000040);
		GImage image2;
		image2.setSize(pImage->width() - 16, 16);
		image2.clear(0);
		GRect r2(0, 0, pImage->width() - 16, 16);

		int valueCount = pRelation->valueCount(attr);
		for(int i = 0; i < valueCount; i++)
		{
			GImage image2;
			image2.setSize(pImage->width() - 16, 16);
			image2.clear(0);
			int xx = 0;
			string sValue;
			pRelation->attrValue(&sValue, attr, i);
			int eatspace = pImage->width() - 16 - GImage::measureTextWidth(sValue.c_str(), 1.0f);
			xx += eatspace;
			image2.text(sValue.c_str(), xx, 0, 1.0f, 0xff000040);
			GImage image3;
			image3.rotateClockwise90(&image2);
			GRect r3(0, 0, image3.width(), image3.height());
			int span = pImage->width() / valueCount;
			int start = MAX(0, (span - 16) / 2);
			pImage->blitAlpha(start + span * i, 16, &image3, &r3);
		}
	}


}

void PlotCorrelations(GArgReader& args)
{
	// Load the data
	GData* pData = GData::loadArff(args.pop_string());
	Holder<GData> hData(pData);
	GArffRelation* pRel = (GArffRelation*)pData->relation().get();

	// Parse options
	string filename = "plot.png";
	int cellsize = 120;
	int bordersize = 4;
	unsigned int bgCol = 0xffd0d0e0;
	int maxAttrs = 30;
	double jitter = 0.03;
	while(args.next_is_flag())
	{
		if(args.if_pop("-out"))
			filename = args.pop_string();
		else if(args.if_pop("-cellsize"))
			cellsize = args.pop_uint();
		else if(args.if_pop("-jitter"))
			jitter = args.pop_double();
		else if(args.if_pop("-maxattrs"))
			maxAttrs = args.pop_uint();
		else
			ThrowError("Invalid option: ", args.peek());
	}

	// Make the chart
	GImage imageBig;
	int wid = (MIN(maxAttrs, pRel->size()) + 1) * (cellsize + bordersize);
	imageBig.setSize(wid, wid);
	imageBig.clear(bgCol);
	GRand prng(getpid() * (unsigned int)time(NULL));
	GImage imageCell;
	GImage imageCell2;
	imageCell.setSize(cellsize, cellsize);
	for(int i = 0; i < pRel->size() && i < maxAttrs; i++)
	{
		MakeCorrelationLabel(pRel, pData, &imageCell, i, bgCol);
		GRect r(0, 0, cellsize, cellsize);
		imageBig.blit((i + 1) * (cellsize + bordersize), 0, &imageCell, &r);
		imageCell2.rotateCounterClockwise90(&imageCell);
		imageBig.blit(0, (i + 1) * (cellsize + bordersize), &imageCell2, &r);
	}
	for(int y = 0; y < pRel->size() && y < maxAttrs; y++)
	{
		for(int x = 0; x < pRel->size() && x < maxAttrs; x++)
		{
			if(x == y)
				MakeAttributeSummaryGraph(pRel, pData, &imageCell, x);
			else
				MakeCorrelationGraph(pRel, pData, &imageCell, x, y, jitter, &prng);
			GRect r(0, 0, cellsize, cellsize);
			imageBig.blit((x + 1) * (cellsize + bordersize), (y + 1) * (cellsize + bordersize), &imageCell, &r);
		}
	}
	imageBig.savePng(filename.c_str());
	cout << "Output saved to " << filename.c_str() << ".\n";
}

class Compare3DPointsByDistanceFromCameraFunctor
{
protected:
	GCamera* m_pCamera;

public:
	Compare3DPointsByDistanceFromCameraFunctor(GCamera* pCamera)
	: m_pCamera(pCamera)
	{
	}

	// returns false if pA is closer than pB
	bool operator() (const double* pA, const double* pB) const
	{
		G3DVector a, b, c, d;
		a.m_vals[0] = pA[0];
		a.m_vals[1] = pA[1];
		a.m_vals[2] = pA[2];
		b.m_vals[0] = pB[0];
		b.m_vals[1] = pB[1];
		b.m_vals[2] = pB[2];
		m_pCamera->project(&a, &c);
		m_pCamera->project(&b, &d);
		return (c.m_vals[2] > d.m_vals[2]);
	}
};

void toImageCoords(GImage* pImage, GCamera* pCamera, G3DVector* pIn, G3DVector* pOut)
{
	pCamera->project(pIn, pOut);

	// Flip the Y value, because positive is down in image coordinates
	pOut->m_vals[1] = pImage->height() - 1 - pOut->m_vals[1];
}

void Plot3d(GImage* pImage, GData* pData, unsigned int bgCol, float pointRadius, double cameraDist, G3DVector* pCameraDirection, bool box, bool labels)
{
	GCamera camera(pImage->width(), pImage->height());
	camera.setViewAngle(M_PI / 3);
	G3DVector mean;
	mean.m_vals[0] = pData->mean(0);
	mean.m_vals[1] = pData->mean(1);
	mean.m_vals[2] = pData->mean(2);
	G3DVector min, max, range;
	pData->minAndRangeUnbiased(0, &min.m_vals[0], &range.m_vals[0]);
	pData->minAndRangeUnbiased(1, &min.m_vals[1], &range.m_vals[1]);
	pData->minAndRangeUnbiased(2, &min.m_vals[2], &range.m_vals[2]);
	max.copy(&range);
	max.add(&min);
	G3DReal dist = sqrt(min.squaredDist(&max)) * cameraDist;
	G3DVector* pCameraPos = camera.lookFromPoint();
	pCameraPos->copy(pCameraDirection);
	pCameraPos->multiply(-1);
	pCameraPos->normalize();
	pCameraPos->multiply(dist);
	pCameraPos->add(&mean);
	camera.setDirection(pCameraDirection, 0.0);

	G3DVector point, coords, point2, coords2;
	pImage->clear(bgCol);

	// Draw box
	if(box)
	{
		min.subtract(&mean);
		min.multiply(1.1);
		min.add(&mean);
		max.subtract(&mean);
		max.multiply(1.1);
		max.add(&mean);
		range.multiply(1.1);
		int x, y, z;
		for(z = 0; z < 2; z++)
		{
			for(y = 0; y < 2; y++)
			{
				for(x = 0; x < 2; x++)
				{
					if(x == 0)
					{
						point.set(min.m_vals[0], min.m_vals[1] + y * range.m_vals[1], min.m_vals[2] + z * range.m_vals[2]);
						point2.set(max.m_vals[0], min.m_vals[1] + y * range.m_vals[1], min.m_vals[2] + z * range.m_vals[2]);
						toImageCoords(pImage, &camera, &point, &coords);
						toImageCoords(pImage, &camera, &point2, &coords2);
						pImage->line((int)coords.m_vals[0], (int)coords.m_vals[1], (int)coords2.m_vals[0], (int)coords2.m_vals[1], 0xff808080);
					}
					if(y == 0)
					{
						point.set(min.m_vals[0] + x * range.m_vals[0], min.m_vals[1], min.m_vals[2] + z * range.m_vals[2]);
						point2.set(min.m_vals[0] + x * range.m_vals[0], max.m_vals[1], min.m_vals[2] + z * range.m_vals[2]);
						toImageCoords(pImage, &camera, &point, &coords);
						toImageCoords(pImage, &camera, &point2, &coords2);
						pImage->line((int)coords.m_vals[0], (int)coords.m_vals[1], (int)coords2.m_vals[0], (int)coords2.m_vals[1], 0xff808080);
					}
					if(z == 0)
					{
						point.set(min.m_vals[0] + x * range.m_vals[0], min.m_vals[1] + y * range.m_vals[1], min.m_vals[2]);
						point2.set(min.m_vals[0] + x * range.m_vals[0], min.m_vals[1] + y * range.m_vals[1], max.m_vals[2]);
						toImageCoords(pImage, &camera, &point, &coords);
						toImageCoords(pImage, &camera, &point2, &coords2);
						pImage->line((int)coords.m_vals[0], (int)coords.m_vals[1], (int)coords2.m_vals[0], (int)coords2.m_vals[1], 0xff808080);
					}
				}
			}
		}

		// Draw axis labels
		if(labels)
		{
			{
				char tmp[32];
				GPlotLabelSpacer pls(min.m_vals[0], max.m_vals[0], 10);
				for(int i = 0; i < pls.count(); i++)
				{
					point.set(pls.label(i), min.m_vals[1], min.m_vals[2]);
					toImageCoords(pImage, &camera, &point, &coords);
					pImage->dot((float)coords.m_vals[0], (float)coords.m_vals[1], 3.0, 0xff404040, bgCol);
					sprintf(tmp, "%.5lg", pls.label(i));
					pImage->text(tmp, (int)coords.m_vals[0] + 4, (int)coords.m_vals[1] - 4, 1.0f, 0xff404040);
				}
			}
			{
				char tmp[32];
				GPlotLabelSpacer pls(min.m_vals[1], max.m_vals[1], 10);
				for(int i = 0; i < pls.count(); i++)
				{
					point.set(min.m_vals[0], pls.label(i), min.m_vals[2]);
					toImageCoords(pImage, &camera, &point, &coords);
					pImage->dot((float)coords.m_vals[0], (float)coords.m_vals[1], 3.0, 0xff404040, bgCol);
					sprintf(tmp, "%.5lg", pls.label(i));
					pImage->text(tmp, (int)coords.m_vals[0] + 4, (int)coords.m_vals[1] - 4, 1.0f, 0xff404040);
				}
			}
			{
				char tmp[32];
				GPlotLabelSpacer pls(min.m_vals[2], max.m_vals[2], 10);
				for(int i = 0; i < pls.count(); i++)
				{
					point.set(min.m_vals[0], min.m_vals[1], pls.label(i));
					toImageCoords(pImage, &camera, &point, &coords);
					pImage->dot((float)coords.m_vals[0], (float)coords.m_vals[1], 3.0, 0xff404040, bgCol);
					sprintf(tmp, "%.5lg", pls.label(i));
					pImage->text(tmp, (int)coords.m_vals[0] + 4, (int)coords.m_vals[1] - 4, 1.0f, 0xff404040);
				}
			}
		}
	}

	// Plot the points
	Compare3DPointsByDistanceFromCameraFunctor comparator(&camera);
	GData copy(4);
	copy.newRows(pData->rows());
	copy.copyColumns(0, pData, 0, 3);
	for(size_t i = 0; i < copy.rows(); i++)
		copy.row(i)[3] = i;
	copy.sort(comparator);
	for(size_t i = 0; i < copy.rows(); i++)
	{
		double* pVec = copy.row(i);
		point.set(pVec[0], pVec[1], pVec[2]);
		toImageCoords(pImage, &camera, &point, &coords);
		float radius = pointRadius / (float)coords.m_vals[2];
		pImage->dot((float)coords.m_vals[0], (float)coords.m_vals[1], radius, gAHSV(0xff, 0.8f * (float)pVec[3] / copy.rows(), 1.0f, 0.5f), bgCol);
	}
}

void Plot3dMulti(GArgReader& args)
{
	// Load
	GData* pData = GData::loadArff(args.pop_string());
	Holder<GData> hData(pData);
	GArffRelation* pRel = (GArffRelation*)pData->relation().get();

	// Parse options
	unsigned int nSeed = getpid() * (unsigned int)time(NULL);
	int horizPlots = 1;
	int vertPlots = 1;
	int wid = 1000;
	int hgt = 1000;
	string filename = "plot.png";
	float pointRadius = 40.0f;
	double cameraDistance = 1.5;
	bool box = true;
	bool labels = true;
	unsigned int cBackground = 0xffffffff;
	G3DVector cameraDirection;
	cameraDirection.set(0.6, -0.3, -0.8);
	bool blast = false;
	while(args.next_is_flag())
	{
		if(args.if_pop("-blast"))
			blast = true;
		else if(args.if_pop("-seed"))
			nSeed = args.pop_uint();
		else if(args.if_pop("-out"))
			filename = args.pop_string();
		else if(args.if_pop("-size"))
		{
			wid = args.pop_uint();
			hgt = args.pop_uint();
		}
		else if(args.if_pop("-pointradius"))
			pointRadius = (float)args.pop_double();
		else if(args.if_pop("-bgcolor"))
			cBackground = hexToRgb(args.pop_string());
		else if(args.if_pop("-cameradistance"))
			cameraDistance = args.pop_double();
		else if(args.if_pop("-cameradirection"))
		{
			cameraDirection.m_vals[0] = args.pop_double();
			cameraDirection.m_vals[1] = args.pop_double();
			cameraDirection.m_vals[2] = args.pop_double();
		}
		else if(args.if_pop("-nobox"))
			box = false;
		else if(args.if_pop("-nolabels"))
			labels = false;
		else
			ThrowError("Invalid option: ", args.peek());
	}
	if(blast)
	{
		pointRadius /= 5;
		wid /= 5;
		hgt /= 5;
		horizPlots *= 5;
		vertPlots *= 5;
	}

	// Check values
	if(pRel->size() != 3)
		ThrowError("Sorry, only data with 3 dims is currently supported");
	if(!pRel->areContinuous(0,3))
		ThrowError("Sorry, only continuous attributes are currently supported");

	// Make plots
	GRand prng(nSeed);
	GImage masterImage;
	masterImage.setSize(horizPlots * wid, vertPlots * hgt);
	GImage tmpImage;
	tmpImage.setSize(wid, hgt);
	for(int y = 0; y < vertPlots; y++)
	{
		for(int x = 0; x < horizPlots; x++)
		{
			if(blast)
			{
				cameraDirection.m_vals[0] = prng.normal();
				cameraDirection.m_vals[1] = prng.normal();
				cameraDirection.m_vals[2] = prng.normal();
				cameraDirection.normalize();
				cout << "row " << y << ", col " << x << ", cameradirection " << cameraDirection.m_vals[0] << " " << cameraDirection.m_vals[1] << " " << cameraDirection.m_vals[2] << "\n";
			}
			Plot3d(&tmpImage, pData, cBackground, pointRadius, cameraDistance, &cameraDirection, box, labels);
			GRect r(0, 0, wid, hgt);
			masterImage.blit(wid * x, hgt * y, &tmpImage, &r);
		}
	}
	masterImage.savePng(filename.c_str());
	cout << "Plot saved to " << filename.c_str() << ".\n";
}

void PrintStats(GArgReader& args)
{
	// Load
	const char* szFilename = args.pop_string();
	GData* pData = GData::loadArff(szFilename);
	Holder<GData> hData(pData);
	GArffRelation* pRel = (GArffRelation*)pData->relation().get();

	// Print some quick stats
	cout.precision(5);
	cout << "Filename: " << szFilename << "\n";
	cout << "Patterns: " << pData->rows() << "\n";
	int continuousAttrs = 0;
	for(int i = 0; i < pRel->size(); i++)
	{
		if(pRel->valueCount(i) == 0)
			continuousAttrs++;
	}
	cout << "Attributes: " << pRel->size() << " (Continuous:" << continuousAttrs << ", Nominal:" << pRel->size() - continuousAttrs << ")\n";
	int stepSize = pRel->size() / 20;
	if(stepSize < 4)
		stepSize = 1;
	else
		cout << "Displaying every " << stepSize << "th attribute...\n";
	for(int i = 0; i < pRel->size(); i += stepSize)
	{
		cout << "  " << i << ") " << pRel->attrName(i) << ", ";
		if(pRel->valueCount(i) == 0)
		{
			cout << "Type: Continuous, ";
			double d1, d2, d3;
			d1 = pData->mean(i);
			d2 = pData->variance(i, d1);
			d3 = pData->median(i);
			cout << "Mean:" << d1 << ", Dev:" << sqrt(d2) << ", Median:" << d3 << ", ";
			pData->minAndRange(i, &d1, &d2);
			cout << "Min:" << d1 << ", Max:" << d1 + d2 << ", ";
			cout << "Missing:" << pData->countValue(i, UNKNOWN_REAL_VALUE) << "\n";
		}
		else
		{
			cout << "Type: Nominal, ";
			cout << "Values:" << pRel->valueCount(i) << ", ";
			int nMostCommonVal = (int)pData->baselineValue(i);
			int mostCommonOccurrences = pData->countValue(i, nMostCommonVal);
			string s;
			pRel->attrValue(&s, i, nMostCommonVal);
			cout << "Most Common:" << s << " (" << ((double)mostCommonOccurrences * 100.0 / pData->rows()) << "%), ";
			cout << "Entropy: " << pData->entropy(i) << ", ";
			cout << "Missing:" << pData->countValue(i, UNKNOWN_DISCRETE_VALUE) << "\n";
			if(pRel->valueCount(i) < 9)
			{
				for(int j = 0; j < pRel->valueCount(i); j++)
				{
					s.clear();
					pRel->attrValue(&s, i, j);
					cout << "     " << ((double)pData->countValue(i, j) * 100.0 / pData->rows()) << "% " << s << "\n";
				}
			}
		}
	}
}

void printDecisionTree(GArgReader& args)
{
	// Load the model
	GTwtDoc doc;
	if(args.size() < 1)
		ThrowError("Model not specified.");
	doc.load(args.pop_string());
	GLearnerLoader ll(true);
	if(_stricmp(doc.root()->field("class")->asString(), "GDecisionTree") != 0)
		ThrowError("That model is not a decision tree");
	GRand prng(0);
	GSupervisedLearner* pModeler = ll.loadModeler(doc.root(), &prng);
	Holder<GSupervisedLearner> hModeler(pModeler);

	GData* pData = NULL;
	if(args.size() > 0)
		pData = GData::loadArff(args.pop_string());
	Holder<GData> hData(pData);

	// Print the decision tree
	((GDecisionTree*)pModeler)->print(cout, pData ? (GArffRelation*)pData->relation().get() : NULL);
}

void model(GArgReader& args)
{
	// Load the model
	GTwtDoc doc;
	if(args.size() < 1)
		ThrowError("Model not specified");
	doc.load(args.pop_string());
	GLearnerLoader ll(true);
	GRand prng(0);
	GSupervisedLearner* pModeler = ll.loadModeler(doc.root(), &prng);
	Holder<GSupervisedLearner> hModeler(pModeler);

	// Load the data
	if(args.size() < 1)
		ThrowError("Expected the filename of a dataset");
	GData* pData = GData::loadArff(args.pop_string());
	Holder<GData> hData(pData);
	if(pData->cols() != pModeler->featureDims() + pModeler->labelDims())
		ThrowError("Model was trained with a different number of attributes than in this data");

	// Get other parameters
	unsigned int attrx = args.pop_uint();
	if(pData->relation()->valueCount(attrx) != 0)
		ThrowError("Sorry, currently only continuous attributes can be plotted");
	unsigned int attry = args.pop_uint();
	if(pData->relation()->valueCount(attry) != 0)
		ThrowError("Sorry, currently only continuous attributes can be plotted");
	int featureDims = pModeler->featureDims();
	if(attrx >= (unsigned int)featureDims || attry >= (unsigned int)featureDims)
		ThrowError("feature attribute out of range");

	// Parse options
	int width = 400;
	int height = 400;
	int labelDim = 0;
	float dotRadius = 3.0f;
	string filename = "plot.png";

	// Compute label range
	double labelMin = 0.0;
	double labelRange = (double)pData->relation()->valueCount(featureDims + labelDim);
	if(labelRange == 0.0)
		pData->minAndRangeUnbiased(featureDims + labelDim, &labelMin, &labelRange);

	// Plot the data
	double xmin, xrange, ymin, yrange;
	pData->minAndRangeUnbiased(attrx, &xmin, &xrange);
	pData->minAndRangeUnbiased(attry, &ymin, &yrange);
	GImage image;
	image.setSize(width, height);
	GPlotWindow pw(&image, xmin, ymin, xmin + xrange, ymin + yrange);
	GTEMPBUF(double, features, pData->cols());
	double* labels = features + featureDims;
	unsigned int* pPix = image.pixels();
	double xx, yy;
	for(int y = 0; y < height; y++)
	{
		cout << ((float)y * 100.0f / height) << "%       \r";
		cout.flush();
		for(int x = 0; x < width; x++)
		{
			pw.viewToWindow(x, y, &xx, &yy);
			size_t r = 0;
			size_t g = 0;
			size_t b = 0;
			for(size_t i = 0; i < pData->rows(); i++)
			{
				GVec::copy(features, pData->row(i), featureDims);
				features[attrx] = xx;
				features[attry] = yy;
				pModeler->predict(features, labels);
				unsigned int hue = gAHSV(0xff, MAX(0.0f, MIN(1.0f, (float)((labels[labelDim] - labelMin) / labelRange))), 1.0f, 0.5f);
				r += gRed(hue);
				g += gGreen(hue);
				b += gBlue(hue);
			}
			r /= pData->rows();
			g /= pData->rows();
			b /= pData->rows();
			*pPix = gARGB(0xff, ClipChan(r), ClipChan(g), ClipChan(b));
			pPix++;
		}
	}
	cout << "                \n";
	cout.flush();

	// Plot the data
	for(size_t i = 0; i < pData->rows(); i++)
	{
		double* pRow = pData->row(i);
		pw.dot(pRow[attrx], pRow[attry], dotRadius, gAHSV(0xff, MAX(0.0f, MIN(1.0f, (float)((pRow[featureDims + labelDim] - labelMin) / labelRange))), 1.0, 1.0), 0xff000000);
	}

	image.savePng(filename.c_str());
	cout << "Output saved to " << filename.c_str() << ".\n";
}

void rowToImage(GArgReader& args)
{
	GData* pData = GData::loadArff(args.pop_string());
	Holder<GData> hData(pData);
	unsigned int r = args.pop_uint();
	if(r > pData->rows())
		ThrowError("row index out of range");
	unsigned int width = args.pop_uint();

	string filename = "plot.png";
	int channels = 3;
	double range = 255.0;

	int cols = pData->cols();
	if((cols % (channels * width)) != 0)
		ThrowError("The row has ", gformat(cols), " dims, which is not a multiple of ", gformat(channels), " channels times ", gformat(width), " pixels wide");
	double* pRow = pData->row(r);
	unsigned int height = cols / (channels * width);
	GImage image;
	GVec::toImage(pRow, &image, width, height, channels, range);
	image.savePng(filename.c_str());
	cout << "Image saved to " << filename.c_str() << ".\n";
}

void systemFrames(GArgReader& args)
{
	GTwtDoc doc;
	doc.load(args.pop_string());
	GData* pActions = GData::loadArff(args.pop_string());
	Holder<GData> hActions(pActions);
	GData* pObs = NULL;
	Holder<GData> hObs(NULL);

	// Parse options
	unsigned int seed = getpid() * (unsigned int)time(NULL);
	bool calibrate = false;
	int frameWidth = 256;
	int stepsPerFrame = 1;
	double scalePredictions = 1.0;
	string outFilename = "frames.png";
	while(args.next_is_flag())
	{
		if(args.if_pop("-seed"))
			seed = args.pop_uint();
		else if(args.if_pop("-calibrate"))
			calibrate = true;
		else if(args.if_pop("-framewidth"))
			frameWidth = args.pop_uint();
		else if(args.if_pop("-stepsperframe"))
			stepsPerFrame = args.pop_uint();
		else if(args.if_pop("-scalepredictions"))
			scalePredictions = args.pop_double();
		else if(args.if_pop("-out"))
			outFilename = args.pop_string();
		else if(args.if_pop("-observations"))
		{
			pObs = GData::loadArff(args.pop_string());
			hObs.reset(pObs);
		}
		else
			ThrowError("Invalid option: ", args.peek());
	}

	// Instantiate the model
	GRand prng(seed);
	GRecurrentModel rm(doc.root(), &prng);
	GImage* pImage = rm.frames(pActions, pObs, calibrate, frameWidth, stepsPerFrame, scalePredictions);
	Holder<GImage> hImage(pImage);
	pImage->savePng(outFilename.c_str());
	cout << "Frames saved to " << outFilename.c_str() << ".\n";
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
		else if(args.if_pop("3d")) Plot3dMulti(args);
		else if(args.if_pop("bar")) PlotBar(args);
		else if(args.if_pop("bigo")) EstimateBigO(args);
		else if(args.if_pop("equation")) PlotEquation(args);
		else if(args.if_pop("histogram")) makeHistogram(args);
		else if(args.if_pop("model")) model(args);
		else if(args.if_pop("overview")) PlotCorrelations(args);
		else if(args.if_pop("rowtoimage")) rowToImage(args);
		else if(args.if_pop("printdecisiontree")) printDecisionTree(args);
		else if(args.if_pop("scatter")) PlotScatter(args);
		else if(args.if_pop("stats")) PrintStats(args);
		else if(args.if_pop("systemframes")) systemFrames(args);
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

