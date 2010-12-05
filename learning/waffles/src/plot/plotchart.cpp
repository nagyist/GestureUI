// --------------------------------------------------------
// This demo file is dedicated to the Public Domain. See:
// http://creativecommons.org/licenses/publicdomain
// --------------------------------------------------------

#include "plotchart.h"
#ifdef WIN32
#else // WIN32
#include <unistd.h>
#endif // !WIN32
#include "../GClasses/GData.h"
#include "../GClasses/GBits.h"
#include "../GClasses/GFile.h"
#include "../GClasses/GTime.h"
#include "../GClasses/GMacros.h"
#include "../GClasses/GMath.h"
#include "../GClasses/GHashTable.h"
#include "../GClasses/GThread.h"
#include "../GClasses/GHillClimber.h"
#include "../GClasses/GRect.h"
#include "../GClasses/GPlot.h"
#include "../GClasses/GNeighborFinder.h"
#include <math.h>

using namespace GClasses;

#define ALMOST1 .999999999

class DoubleRect
{
public:
	double x, y, w, h;

	DoubleRect()
	{
	}

	DoubleRect(double dX, double dY, double dW, double dH)
	{
		x = dX;
		y = dY;
		w = dW;
		h = dH;
	}
};



PlotChartMaker::PlotChartMaker(sp_relation& pRelation, GData* pData)
{
	m_pRelation = pRelation;
	m_nOutputCount = pRelation->size() - 1;
	m_pData = pData;
	m_bCustomRange = false;
	m_logx = false;
	m_logy = false;
	m_showHorizAxisLabels = true;
	m_showVertAxisLabels = true;
	m_showLines = true;
	m_fLineThickness = 3.0;
	m_meshSize = 0;
	m_textSize = 2;
	m_fPointRadius = 7.0;
	m_nChartWidth = 1024;
	m_nChartHeight = 1024;
	m_nVertLabelPrecision = 4;
	m_nHorizLabelPrecision = 4;
	m_cBackground = 0xffffffff;
	m_cText = 0xff000000;
	m_cGrid = 0xff808080;
	m_cPlot[0] = 0xff0000a0;
	m_cPlot[1] = 0xffa00000;
	m_cPlot[2] = 0xff008000;
	m_cPlot[3] = 0xff504010;
	m_useSpectrumColors = false;
	m_spectrumModulus = 0;
	m_aspect = false;
	m_pNeighborFinder = NULL;
	m_maxGridLinesH = 30;
	m_maxGridLinesV = 30;
}

PlotChartMaker::~PlotChartMaker()
{
	delete(m_pNeighborFinder);
}

void PlotChartMaker::ComputeChartRange()
{
	m_pData->minAndRange(0, &m_dXMin, &m_dXMax);
	m_dXMax += m_dXMin;
	m_pData->minAndRange(1, &m_dYMin, &m_dYMax);
	m_dYMax += m_dYMin;
	for(int i = 2; i < m_pRelation->size(); i++)
	{
		double m, r;
		m_pData->minAndRange(i, &m, &r);
		r += m;
		m_dYMin = MIN(m_dYMin, m);
		m_dYMax = MAX(m_dYMax, r);
	}

	// Add some border space
	double border = 0.3;
	if(m_logx)
	{
		m_dXMin = MAX(1e-12, m_dXMin);
		m_dXMax = MAX(1e-10, m_dXMax);
		double d = pow(m_dXMax / m_dXMin, border);
		m_dXMax *= d;
		m_dXMin /= d;
	}
	else
	{
		double d = border * (m_dXMax - m_dXMin);
		m_dXMax += d;
		m_dXMin -= d;
	}
	if(m_logy)
	{
		m_dYMin = MAX(1e-12, m_dYMin);
		m_dYMax = MAX(1e-10, m_dYMax);
		double d = pow(m_dYMax / m_dYMin, border);
		m_dYMax *= d;
		m_dYMin /= d;
	}
	else
	{
		double d = border * (m_dYMax - m_dYMin);
		m_dYMax += d;
		m_dYMin -= d;
	}
}

unsigned int PlotChartMaker::GetLineColor(int line, size_t index)
{
	if(m_useSpectrumColors)
	{
		if(m_spectrumModulus > 0)
			return gAHSV(0xff, (float)(index % m_spectrumModulus) * 0.8f / m_spectrumModulus, 1.0f, 0.5f);
		else
			return gAHSV(0xff, (float)index * 0.8f / m_pData->rows(), 1.0f, 0.5f);
	}
	else
	{
		if(m_pRelation->size() > 5)
			return gAHSV(0xff, (float)line / (m_pRelation->size() - 1), 1.0f, 0.5f);
		else
			return m_cPlot[line];
	}
}

double PlotChartMaker::xval(double x)
{
	if(m_logx)
		return log(x);
	else
		return x;
}

double PlotChartMaker::yval(double y)
{
	if(m_logy)
		return log(y);
	else
		return y;
}

GImage* PlotChartMaker::MakeChart()
{
	// Compute the size of the image
	if(m_nOutputCount <= 0)
		throw "There are no output values to chart";

	GImage* pImage = new GImage();
	Holder<GImage> hImage(pImage);
	pImage->setSize(m_nChartWidth, m_nChartHeight);
	pImage->clear(m_cBackground);

	if(!m_bCustomRange)
		ComputeChartRange();
	if(m_aspect)
	{
		if((m_dXMax - m_dXMin) / (m_dYMax - m_dYMin) < (double)m_nChartWidth / m_nChartHeight)
		{
			double d = 0.5 * ((m_dYMax - m_dYMin) * ((double)m_nChartWidth / m_nChartHeight) - (m_dXMax - m_dXMin));
			m_dXMin -= d;
			m_dXMax += d;
		}
		else
		{
			double d = 0.5 * ((m_dXMax - m_dXMin) * ((double)m_nChartHeight / m_nChartWidth) - (m_dYMax - m_dYMin));
			m_dYMin -= d;
			m_dYMax += d;
		}
	}
	if(m_logx)
	{
		m_dXMin = MAX(1e-24, m_dXMin);
		m_dXMax = MAX(m_dXMin + 1e-12, m_dXMax);
		m_dXMin = log(m_dXMin);
		m_dXMax = log(m_dXMax);
	}
	if(m_logy)
	{
		m_dYMin = MAX(1e-24, m_dYMin);
		m_dYMax = MAX(m_dYMin + 1e-12, m_dYMax);
		m_dYMin = log(m_dYMin);
		m_dYMax = log(m_dYMax);
	}
	GPlotWindow pw(pImage, m_dXMin, m_dYMin, m_dXMax, m_dYMax);

	// Plot the grid lines
	pw.gridLines((m_showHorizAxisLabels ? (m_logx ? -1 : m_maxGridLinesH) : 0), (m_showVertAxisLabels ? (m_logy ? -1 : m_maxGridLinesV) : 0), m_cGrid);

	// Plot the lines
	if(m_showLines && m_meshSize == 0)
	{
		for(int attr = 1; attr < m_pRelation->size(); attr++)
		{
			double prevx = UNKNOWN_REAL_VALUE;
			double prevy = UNKNOWN_REAL_VALUE;
			for(size_t pat = 0; pat < m_pData->rows(); pat++)
			{
				double* pPat = m_pData->row(pat);
				if(pPat[0] != UNKNOWN_REAL_VALUE && pPat[attr] != UNKNOWN_REAL_VALUE)
				{
					if(prevx != UNKNOWN_REAL_VALUE && prevy != UNKNOWN_REAL_VALUE)
					{
						pw.fatLine(xval(prevx), yval(prevy), xval(pPat[0]), yval(pPat[attr]), m_fLineThickness, GetLineColor(attr - 1, pat));
					}
					prevx = pPat[0];
					prevy = pPat[attr];
				}
			}
		}
	}

	// Mesh lines (This is an obscure feature that I don't think anyone else will ever use.)
	if(m_meshSize > 0)
	{
		for(int attr = 1; attr < m_pRelation->size(); attr++)
		{
			for(size_t pat = 0; pat < m_pData->rows(); pat++)
			{
				int x = pat % m_meshSize;
				int y = pat / m_meshSize;
				double* pPat = m_pData->row(pat);
				if(pPat[0] != UNKNOWN_REAL_VALUE && pPat[attr] != UNKNOWN_REAL_VALUE)
				{
					if(x > 0)
					{
						double* pPatNeighbor = m_pData->row(pat - 1);
						if(pPatNeighbor[0] != UNKNOWN_REAL_VALUE && pPatNeighbor[attr] != UNKNOWN_REAL_VALUE)
							pw.fatLine(xval(pPat[0]), yval(pPat[attr]), xval(pPatNeighbor[0]), yval(pPatNeighbor[attr]), m_fLineThickness, GetLineColor(attr - 1, pat));
					}
					if(x < m_meshSize - 1)
					{
						double* pPatNeighbor = m_pData->row(pat + 1);
						if(pPatNeighbor[0] != UNKNOWN_REAL_VALUE && pPatNeighbor[attr] != UNKNOWN_REAL_VALUE)
							pw.fatLine(xval(pPat[0]), yval(pPat[attr]), xval(pPatNeighbor[0]), yval(pPatNeighbor[attr]), m_fLineThickness, GetLineColor(attr - 1, pat));
					}
					if(y > 0)
					{
						double* pPatNeighbor = m_pData->row(pat - m_meshSize);
						if(pPatNeighbor[0] != UNKNOWN_REAL_VALUE && pPatNeighbor[attr] != UNKNOWN_REAL_VALUE)
							pw.fatLine(xval(pPat[0]), yval(pPat[attr]), xval(pPatNeighbor[0]), yval(pPatNeighbor[attr]), m_fLineThickness, GetLineColor(attr - 1, pat));
					}
					if(y < m_meshSize - 1)
					{
						double* pPatNeighbor = m_pData->row(pat + m_meshSize);
						if(pPatNeighbor[0] != UNKNOWN_REAL_VALUE && pPatNeighbor[attr] != UNKNOWN_REAL_VALUE)
							pw.fatLine(xval(pPat[0]), yval(pPat[attr]), xval(pPatNeighbor[0]), yval(pPatNeighbor[attr]), m_fLineThickness, GetLineColor(attr - 1, pat));
					}
				}
			}
		}
	}

	// Draw neighbor connections
	if(m_pNeighborFinder)
	{
		GTEMPBUF(size_t, hood, m_pNeighborFinder->neighborCount());
		for(int attr = 1; attr < m_pRelation->size(); attr++)
		{
			for(size_t pat = 0; pat < m_pData->rows(); pat++)
			{
				double* pPat = m_pData->row(pat);
				if(pPat[0] != UNKNOWN_REAL_VALUE && pPat[attr] != UNKNOWN_REAL_VALUE)
				{
					m_pNeighborFinder->neighbors(hood, pat);
					for(int j = 0; j < m_pNeighborFinder->neighborCount(); j++)
					{
						if(hood[j] >= m_pData->rows())
							continue;
						double* pOther = m_pData->row(hood[j]);
						if(pOther[0] != UNKNOWN_REAL_VALUE && pOther[attr] != UNKNOWN_REAL_VALUE)
							pw.line(xval(pPat[0]), yval(pPat[attr]), xval(pOther[0]), yval(pOther[attr]), GetLineColor(attr - 1, pat));
					}
				}
			}
		}
	}

	// Plot the dots
	for(int attr = 1; attr < m_pRelation->size(); attr++)
	{
		for(size_t pat = 0; pat < m_pData->rows(); pat++)
		{
			double* pPat = m_pData->row(pat);
			if(pPat[0] != UNKNOWN_REAL_VALUE && pPat[attr] != UNKNOWN_REAL_VALUE)
				pw.dot(xval(pPat[0]), yval(pPat[attr]), m_fPointRadius, GetLineColor(attr - 1, pat), m_cBackground);
		}
	}

	// Label the axes
	return pw.labelAxes((m_showHorizAxisLabels ? (m_logx ? -1 : m_maxGridLinesH) : 0), (m_showVertAxisLabels ? (m_logy ? -1 : m_maxGridLinesV) : 0), m_nHorizLabelPrecision, m_textSize, 0xff000000, 45.0 * (M_PI / 180));
}

