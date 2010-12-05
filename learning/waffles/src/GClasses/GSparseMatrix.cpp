/*
	Copyright (C) 2006, Mike Gashler

	This library is free software; you can redistribute it and/or
	modify it under the terms of the GNU Lesser General Public
	License as published by the Free Software Foundation; either
	version 2.1 of the License, or (at your option) any later version.

	see http://www.gnu.org/copyleft/lesser.html
*/

#include "GMacros.h"
#include "GSparseMatrix.h"
#include <math.h>
#include "GData.h"
#include "GVec.h"
#include "GFile.h"
#include "GRand.h"
#include <fstream>

namespace GClasses {

GSparseMatrix::GSparseMatrix(unsigned int rows, unsigned int cols)
: m_rows(rows), m_cols(cols)
{
	m_pRows = new Map[rows];
	m_pCols = new Map[cols];
}

GSparseMatrix::~GSparseMatrix()
{
	delete[] m_pRows;
	delete[] m_pCols;
}

// static
GSparseMatrix* GSparseMatrix::parse(const char* pFile, size_t len)
{
	// Parse the row count
	if(*pFile < '0' || *pFile > '9')
		ThrowError("expected a row count on line 1");
	unsigned int rows = atoi(pFile);
	while(*pFile >= '0' && *pFile <= '9')
		pFile++;
	if(*pFile != ',')
		ThrowError("expected a comma on line 1");
	pFile++;

	// Parse the column count
	if(*pFile < '0' || *pFile > '9')
		ThrowError("expected a column count on line 1");
	unsigned int cols = atoi(pFile);
	while(*pFile >= '0' && *pFile <= '9')
		pFile++;
	while(*pFile != '\n' && *pFile != '\0')
		pFile++;
	if(*pFile != '\n' && rows > 0)
		ThrowError("expected a newline character on line 1");
	pFile++;
	while(*pFile <= ' ' && *pFile != '\0' && *pFile != '\n')
		pFile++;

	// Parse the data
	GSparseMatrix* pSM = new GSparseMatrix(rows, cols);
	Holder<GSparseMatrix> hSM(pSM);
	for(unsigned int r = 0; r < rows; r++)
	{
		while(true)
		{
			if(*pFile < '0' || *pFile > '9')
				break;
			unsigned int c = atoi(pFile);
			while(*pFile >= '0' && *pFile <= '9')
				pFile++;
			if(*pFile != '=')
				ThrowError("expected an equal sign on line ", gformat(r + 2));
			pFile++;
			double v = atof(pFile);
			pSM->set(r, c, v);
			while(*pFile > ' ' && *pFile != ',')
				pFile++;
			if(*pFile == ',')
				pFile++;
		}
		if(*pFile != '\r' && *pFile != '\n')
			ThrowError("unexpected character on line ", gformat(r + 2));
		while(*pFile != '\n' && *pFile != '\0')
			pFile++;
		if(*pFile != '\n' && r < rows - 1)
			ThrowError("expected a newline character on line ", gformat(r + 2));
		pFile++;
		while(*pFile <= ' ' && *pFile != '\0' && *pFile != '\n')
			pFile++;
	}
	return hSM.release();
}

// static
GSparseMatrix* GSparseMatrix::load(const char* filename)
{
	size_t len;
	char* szFile = GFile::loadFile(filename, &len);
	ArrayHolder<char> hFile(szFile);
	return parse(szFile, len);
}

void GSparseMatrix::save(const char* filename)
{
	std::ofstream os;
	os.exceptions(std::ios::failbit|std::ios::badbit);
	try
	{
		os.open(filename, std::ios::binary);
	}
	catch(const std::exception&)
	{
		ThrowError("Error creating file: ", filename);
	}
	print(os);
}

void GSparseMatrix::print(std::ostream& stream)
{
	stream.precision(7);
	stream << m_rows << " " << m_cols << "\n";
	for(unsigned int r = 0; r < m_rows; r++)
	{
		Iter it = rowBegin(r);
		Iter end = rowEnd(r);
		if(it != end)
		{
			stream << it->first << "=" << it->second;
			it++;
		}
		while(it != end)
		{
			stream << "," << it->first << "=" << it->second;
			it++;
		}
		stream << "\n";
	}
}

void GSparseMatrix::fullRow(double* pOutFullRow, unsigned int row)
{
	GVec::setAll(pOutFullRow, 0.0, m_cols);
	Iter end = rowEnd(row);
	for(Iter it = rowBegin(row); it != end; it++)
		pOutFullRow[it->first] = it->second;
}

double GSparseMatrix::get(unsigned int row, unsigned int col)
{
	GAssert(row < m_rows && col < m_cols); // out of range
	if(m_rows >= m_cols)
	{
		It it = m_pRows[row].find(col);
		if(it == m_pRows[row].end())
			return 0;
		return it->second;
	}
	else
	{
		It it = m_pCols[col].find(row);
		if(it == m_pCols[col].end())
			return 0;
		return it->second;
	}
}

void GSparseMatrix::set(unsigned int row, unsigned int col, double val)
{
	GAssert(row < m_rows && col < m_cols); // out of range
	if(val == 0)
	{
		m_pRows[row].erase(col);
		m_pCols[col].erase(row);
	}
	else
	{
		m_pRows[row][col] = val;
		m_pCols[col][row] = val;
	}
}

void GSparseMatrix::multiply(double scalar)
{
	for(unsigned int r = 0; r < m_rows; r++)
	{
		It end = m_pRows[r].end();
		for(It it = m_pRows[r].begin(); it != end; it++)
			it->second *= scalar;
	}
	for(unsigned int c = 0; c < m_cols; c++)
	{
		It end = m_pCols[c].end();
		for(It it = m_pCols[c].begin(); it != end; it++)
			it->second *= scalar;
	}
}

void GSparseMatrix::copyFrom(GSparseMatrix* that)
{
	unsigned int rows = MIN(m_rows, that->rows());
	for(unsigned int r = 0; r < rows; r++)
	{
		Iter end = that->rowEnd(r);
		unsigned int pos = 0;
		for(Iter it = that->rowBegin(r); it != end && pos < m_cols; it++)
		{
			set(r, it->first, it->second);
			pos++;
		}
	}
}

void GSparseMatrix::copyFrom(GData* that)
{
	unsigned int rows = MIN(m_rows, (unsigned int)that->rows());
	unsigned int cols = MIN(m_cols, (unsigned int)that->cols());
	for(unsigned int r = 0; r < rows; r++)
	{
		double* pRow = that->row(r);
		for(unsigned int c = 0; c < cols; c++)
		{
			set(r, c, *pRow);
			pRow++;
		}
	}
}

GData* GSparseMatrix::toFullMatrix()
{
	GData* pData = new GData(m_cols);
	pData->newRows(m_rows);
	for(unsigned int r = 0; r < m_rows; r++)
	{
		double* pRow = pData->row(r);
		GVec::setAll(pRow, 0.0, m_cols);
		Iter end = rowEnd(r);
		for(Iter it = rowBegin(r); it != end; it++)
			pRow[it->first] = it->second;
	}
	return pData;
}

void GSparseMatrix::swapColumns(unsigned int a, unsigned int b)
{
	for(unsigned int r = 0; r < m_rows; r++)
	{
		double aa = get(r, a);
		double bb = get(r, b);
		set(r, a, bb);
		set(r, b, aa);
	}
}

void GSparseMatrix::swapRows(unsigned int a, unsigned int b)
{
	for(unsigned int c = 0; c < m_cols; c++)
	{
		double aa = get(a, c);
		double bb = get(b, c);
		set(a, c, bb);
		set(b, c, aa);
	}
}

void GSparseMatrix::shuffle(GRand* pRand)
{
	for(unsigned int n = m_rows; n > 0; n--)
		swapRows((unsigned int)pRand->next(n), n - 1);
}

GSparseMatrix* GSparseMatrix::subMatrix(int row, int col, int height, int width)
{
	if(row < 0 || col < 0 || row + height >= (int)m_rows || col + width >= (int)m_cols || height < 0 || width < 0)
		ThrowError("out of range");
	GSparseMatrix* pSub = new GSparseMatrix(height, width);
	for(int y = 0; y < height; y++)
	{
		for(int x = 0; x < width; x++)
			pSub->set(y, x, get(row + y, col + x));
	}
	return pSub;
}

void GSparseMatrix::transpose()
{
	std::swap(m_pRows, m_pCols);
	std::swap(m_rows, m_cols);
}

double GSparseMatrix_pythag(double a, double b)
{
	double at = ABS(a);
	double bt = ABS(b);
	if(at > bt)
	{
		double ct = bt / at;
		return at * sqrt(1.0 + ct * ct);
	}
	else if(bt > 0.0)
	{
		double ct = at / bt;
		return bt * sqrt(1.0 + ct * ct);
	}
	else
		return 0.0;
}

double GSparseMatrix_takeSign(double a, double b)
{
	return (b >= 0.0 ? ABS(a) : -ABS(a));
}

void GSparseMatrix::singularValueDecomposition(GSparseMatrix** ppU, double** ppDiag, GSparseMatrix** ppV, int maxIters)
{
	if(rows() >= cols())
		singularValueDecompositionHelper(ppU, ppDiag, ppV, maxIters);
	else
	{
		transpose();
		singularValueDecompositionHelper(ppV, ppDiag, ppU, maxIters);
		transpose();
		(*ppV)->transpose();
		(*ppU)->transpose();
	}
}

void GSparseMatrix::singularValueDecompositionHelper(GSparseMatrix** ppU, double** ppDiag, GSparseMatrix** ppV, int maxIters)
{
	int m = rows();
	int n = cols();
	if(m < n)
		ThrowError("Expected at least as many rows as columns");
	int i, j, k;
	int l = 0;
	int q, iter;
	double c, f, h, s, x, y, z;
	double norm = 0.0;
	double g = 0.0;
	double scale = 0.0;
	GSparseMatrix* pU = new GSparseMatrix(m, m);
	Holder<GSparseMatrix> hU(pU);
	pU->copyFrom(this);
	double* pSigma = new double[n];
	ArrayHolder<double> hSigma(pSigma);
	GSparseMatrix* pV = new GSparseMatrix(n, n);
	Holder<GSparseMatrix> hV(pV);
	GTEMPBUF(double, temp, n + m);
	double* temp2 = temp + n;

	// Householder reduction to bidiagonal form
	for(int i = 0; i < n; i++)
	{
		// Left-hand reduction
		temp[i] = scale * g;
		l = i + 1;
		g = 0.0;
		s = 0.0;
		scale = 0.0;
		if(i < m)
		{
			Iter kend = pU->colEnd(i);
			for(Iter kk = pU->colBegin(i); kk != kend; kk++)
			{
				if(kk->first >= (unsigned int)i)
					scale += ABS(kk->second);
			}
			if(scale != 0.0)
			{
				for(Iter kk = pU->colBegin(i); kk != kend; kk++)
				{
					if(kk->first >= (unsigned int)i)
					{
						double t = kk->second / scale;
						pU->set(kk->first, i, t);
						s += (t * t);
					}
				}
				f = pU->get(i, i);
				g = -GSparseMatrix_takeSign(sqrt(s), f);
				h = f * g - s;
				pU->set(i, i, f - g);
				if(i != n - 1)
				{
					for(j = l; j < n; j++)
					{
						s = 0.0;
						for(Iter kk = pU->colBegin(i); kk != kend; kk++)
						{
							if(kk->first >= (unsigned int)i)
								s += kk->second * pU->get(kk->first, j);
						}
						f = s / h;
						for(Iter kk = pU->colBegin(i); kk != kend; kk++)
						{
							if(kk->first >= (unsigned int)i)
								pU->set(kk->first, j, pU->get(kk->first, j) + f * kk->second);
						}
					}
				}
				for(Iter kk = pU->colBegin(i); kk != kend; kk++)
				{
					if(kk->first >= (unsigned int)i)
						pU->set(kk->first, i, pU->get(kk->first, i) * scale);
				}
			}
		}
		pSigma[i] = scale * g;

		// Right-hand reduction
		g = 0.0;
		s = 0.0;
		scale = 0.0;
		if(i < m && i != n - 1) 
		{
			Iter kend = pU->rowEnd(i);
			for(Iter kk = pU->rowBegin(i); kk != kend; kk++)
			{
				if(kk->first >= (unsigned int)n)
					break;
				if(kk->first >= (unsigned int)l)
					scale += ABS(kk->second);
			}
			if(scale != 0.0) 
			{
				for(Iter kk = pU->rowBegin(i); kk != kend; kk++)
				{
					if(kk->first >= (unsigned int)n)
						break;
					if(kk->first >= (unsigned int)l)
					{
						double t = kk->second / scale;
						pU->set(i, kk->first, t);
						s += (t * t);
					}
				}
				f = pU->get(i, l);
				g = -GSparseMatrix_takeSign(sqrt(s), f);
				h = f * g - s;
				pU->set(i, l, f - g);
				for(k = l; k < n; k++)
					temp[k] = pU->get(i, k) / h;
				if(i != m - 1) 
				{
					for(j = l; j < m; j++) 
					{
						s = 0.0;
						for(Iter kk = pU->rowBegin(i); kk != kend; kk++)
						{
							if(kk->first >= (unsigned int)n)
								break;
							if(kk->first >= (unsigned int)l)
								s += pU->get(j, kk->first) * kk->second;
						}
						Iter kend2 = pU->rowEnd(j);
						for(Iter kk = pU->rowBegin(j); kk != kend2; kk++)
						{
							if(kk->first >= (unsigned int)n)
								break;
							if(kk->first >= (unsigned int)l)
								pU->set(j, kk->first, pU->get(j, kk->first) + s * temp[kk->first]);
						}
					}
				}
				for(Iter kk = pU->rowBegin(i); kk != kend; kk++)
				{
					if(kk->first >= (unsigned int)n)
						break;
					if(kk->first >= (unsigned int)l)
						pU->set(i, kk->first, kk->second * scale);
				}
			}
		}
		norm = MAX(norm, ABS(pSigma[i]) + ABS(temp[i]));
	}

	// Accumulate right-hand transform
	for(int i = n - 1; i >= 0; i--)
	{
		if(i < n - 1)
		{
			if(g != 0.0)
			{
				Iter jend = pU->rowEnd(i);
				for(Iter jj = pU->rowBegin(i); jj != jend; jj++)
				{
					if(jj->first >= (unsigned int)n)
						break;
					if(jj->first >= (unsigned int)l)
						pV->set(i, jj->first, (jj->second / pU->get(i, l)) / g); // (double-division to avoid underflow)
				}
				for(j = l; j < n; j++)
				{
					s = 0.0;
					Iter kend = pU->rowEnd(i);
					for(Iter kk = pU->rowBegin(i); kk != kend; kk++)
					{
						if(kk->first >= (unsigned int)n)
							break;
						if(kk->first >= (unsigned int)l)
							s += kk->second * pV->get(j, kk->first);
					}
					kend = pV->rowEnd(i);
					for(Iter kk = pV->rowBegin(i); kk != kend; kk++)
					{
						if(kk->first >= (unsigned int)n)
							break;
						if(kk->first >= (unsigned int)l)
							pV->set(j, kk->first, pV->get(j, kk->first) + s * kk->second);
					}
				}
			}
			for(j = l; j < n; j++)
			{
				pV->set(i, j, 0.0);
				pV->set(j, i, 0.0);
			}
		}
		pV->set(i, i, 1.0);
		g = temp[i];
		l = i;
	}

	// Accumulate left-hand transform
	for(i = n - 1; i >= 0; i--)
	{
		l = i + 1;
		g = pSigma[i];
		if(i < n - 1)
		{
			for(j = l; j < n; j++)
				pU->set(i, j, 0.0);
		}
		if(g != 0.0)
		{
			g = 1.0 / g;
			if(i != n - 1)
			{
				for(j = l; j < n; j++)
				{
					s = 0.0;
					Iter kend = pU->colEnd(i);
					for(Iter kk = pU->colBegin(i); kk != kend; kk++)
					{
						if(kk->first >= (unsigned int)l)
							s += kk->second * pU->get(kk->first, j);
					}
					f = (s / pU->get(i, i)) * g;
					if(f != 0.0)
					{
						for(Iter kk = pU->colBegin(i); kk != kend; kk++)
						{
							if(kk->first >= (unsigned int)i)
								pU->set(kk->first, j, pU->get(kk->first, j) + f * kk->second);
						}
					}
				}
			}
			for(j = i; j < m; j++)
				pU->set(j, i, pU->get(j, i) * g);
		} 
		else 
		{
			for(j = i; j < m; j++)
				pU->set(j, i, 0.0);
		}
		pU->set(i, i, pU->get(i, i) + 1.0);
	}

	// Diagonalize the bidiagonal matrix
	for(k = n - 1; k >= 0; k--) // For each singular value
	{
		for(iter = 1; iter <= maxIters; iter++)
		{
			// Test for splitting
			bool flag = true;
			for(l = k; l >= 0; l--)
			{
				q = l - 1;
				if(ABS(temp[l]) + norm == norm)
				{
					flag = false;
					break;
				}
				if(ABS(pSigma[q]) + norm == norm)
					break;
			}

			if(flag)
			{
				c = 0.0;
				s = 1.0;
				for(i = l; i <= k; i++)
				{
					f = s * temp[i];
					temp[i] *= c;
					if(ABS(f) + norm == norm)
						break;
					g = pSigma[i];
					h = GSparseMatrix_pythag(f, g);
					pSigma[i] = h;
					h = 1.0 / h;
					c = g * h;
					s = -f * h;
					Iter jendi = pU->colEnd(i);
					Iter jendq = pU->colEnd(q);
					Iter jji = pU->colBegin(i);
					Iter jjq = pU->colBegin(q);
					int tpos;
					for(tpos = 0; jji != jendi || jjq != jendq; tpos++)
					{
						if(jjq == jendq || (jji != jendi && jji->first < jjq->first))
						{
							temp2[tpos] = jji->first;
							jji++;
						}
						else
						{
							temp2[tpos] = jjq->first;
							if(jji != jendi && jjq->first == jji->first)
								jji++;
							jjq++;
						}
					}
					for(int tpos2 = 0; tpos2 < tpos; tpos2++)
					{
						y = pU->get((unsigned int)temp2[tpos2], q);
						z = pU->get((unsigned int)temp2[tpos2], i);
						pU->set((unsigned int)temp2[tpos2], q, y * c + z * s);
						pU->set((unsigned int)temp2[tpos2], i, z * c - y * s);
					}
				}
			}

			z = pSigma[k];
			if(l == k)
			{
				// Detect convergence
				if(z < 0.0)
				{
					// Singular value should be positive
					pSigma[k] = -z;
					for(j = 0; j < n; j++)
						pV->set(k, j, pV->get(k, j) * -1.0);
				}
				break;
			}
			if(iter >= maxIters)
				ThrowError("failed to converge");

			// Shift from bottom 2x2 minor
			x = pSigma[l];
			q = k - 1;
			y = pSigma[q];
			g = temp[q];
			h = temp[k];
			f = ((y - z) * (y + z) + (g - h) * (g + h)) / (2.0 * h * y);
			g = GSparseMatrix_pythag(f, 1.0);
			f = ((x - z) * (x + z) + h * ((y / (f + GSparseMatrix_takeSign(g, f))) - h)) / x;

			// QR transform
			c = 1.0;
			s = 1.0;
			for(j = l; j <= q; j++)
			{
				i = j + 1;
				g = temp[i];
				y = pSigma[i];
				h = s * g;
				g = c * g;
				z = GSparseMatrix_pythag(f, h);
				temp[j] = z;
				c = f / z;
				s = h / z;
				f = x * c + g * s;
				g = g * c - x * s;
				h = y * s;
				y = y * c;
				Iter pendi = pV->rowEnd(i);
				Iter pendj = pV->rowEnd(j);
				Iter ppi = pV->rowBegin(i);
				Iter ppj = pV->rowBegin(j);
				int tpos;
				for(tpos = 0; ppi != pendi || ppj != pendj; tpos++)
				{
					if(ppj == pendj || (ppi != pendi && ppi->first < ppj->first))
					{
						temp2[tpos] = ppi->first;
						ppi++;
					}
					else
					{
						temp2[tpos] = ppj->first;
						if(ppi != pendi && ppj->first == ppi->first)
							ppi++;
						ppj++;
					}
				}
				for(int tpos2 = 0; tpos2 < tpos; tpos2++)
				{
					x = pV->get(j, (unsigned int)temp2[tpos2]);
					z = pV->get(i, (unsigned int)temp2[tpos2]);
					pV->set(j, (unsigned int)temp2[tpos2], x * c + z * s);
					pV->set(i, (unsigned int)temp2[tpos2], z * c - x * s);
				}
				z = GSparseMatrix_pythag(f, h);
				pSigma[j] = z;
				if(z != 0.0)
				{
					z = 1.0 / z;
					c = f * z;
					s = h * z;
				}
				f = c * g + s * y;
				x = c * y - s * g;
				pendi = pU->colEnd(i);
				pendj = pU->colEnd(j);
				ppi = pU->colBegin(i);
				ppj = pU->colBegin(j);
				for(tpos = 0; ppi != pendi || ppj != pendj; tpos++)
				{
					if(ppj == pendj || (ppi != pendi && ppi->first < ppj->first))
					{
						temp2[tpos] = ppi->first;
						ppi++;
					}
					else
					{
						temp2[tpos] = ppj->first;
						if(ppi != pendi && ppj->first == ppi->first)
							ppi++;
						ppj++;
					}
				}
				for(int tpos2 = 0; tpos2 < tpos; tpos2++)
				{
					y = pU->get((unsigned int)temp2[tpos2], j);
					z = pU->get((unsigned int)temp2[tpos2], i);
					pU->set((unsigned int)temp2[tpos2], j, y * c + z * s);
					pU->set((unsigned int)temp2[tpos2], i, z * c - y * s);
				}
			}
			temp[l] = 0.0;
			temp[k] = f;
			pSigma[k] = x;
		}
	}

	// Sort the singular values from largest to smallest
	for(i = 1; i < n; i++)
	{
		for(j = i; j > 0; j--)
		{
			if(pSigma[j - 1] >= pSigma[j])
				break;
			pU->swapColumns(j - 1, j);
			pV->swapRows(j - 1, j);
			std::swap(pSigma[j - 1], pSigma[j]);
		}
	}

	// Return results
	*ppU = hU.release();
	*ppDiag = hSigma.release();
	*ppV = hV.release();
}

#ifndef NO_TEST_CODE
// static
void GSparseMatrix::test()
{
	// Make the data
	GSparseMatrix sm(4, 4);
	sm.set(0, 0, 2.0); sm.set(0, 2, 3.0);
	sm.set(1, 0, 1.0); sm.set(2, 3, -2.0);
	sm.set(2, 2, 5.0);
	sm.set(3, 1, -3.0); sm.set(3, 3, -1.0);
	GData* fm = sm.toFullMatrix();
	Holder<GData> hFM(fm);

	// Do it with the full matrix
	GData* pU;
	double* pDiag;
	GData* pV;
	fm->singularValueDecomposition(&pU, &pDiag, &pV);
	Holder<GData> hU(pU);
	ArrayHolder<double> hDiag(pDiag);
	Holder<GData> hV(pV);

	// Do it with the sparse matrix
	GSparseMatrix* pSU;
	double* pSDiag;
	GSparseMatrix* pSV;
	sm.singularValueDecomposition(&pSU, &pSDiag, &pSV);
	Holder<GSparseMatrix> hSU(pSU);
	ArrayHolder<double> hSDiag(pSDiag);
	Holder<GSparseMatrix> hSV(pSV);

	// Check the results
	GData* pV2 = pSV->toFullMatrix();
	Holder<GData> hV2(pV2);
	double err = pV2->sumSquaredDifference(*pV, false);
	if(err > 1e-6)
		ThrowError("Failed");
}
#endif

} // namespace GClasses

