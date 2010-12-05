/*
	Copyright (C) 2006, Mike Gashler

	This library is free software; you can redistribute it and/or
	modify it under the terms of the GNU Lesser General Public
	License as published by the Free Software Foundation; either
	version 2.1 of the License, or (at your option) any later version.

	see http://www.gnu.org/copyleft/lesser.html
*/

#include "GData.h"
#include "GMacros.h"
#include "GMath.h"
#include "GDistribution.h"
#include "GVec.h"
#include "GFile.h"
#include "GHeap.h"
#include "GTwt.h"
#include "GHashTable.h"
#include <math.h>
#include "GBits.h"
#include "GLearner.h"
#include "GRand.h"
#include <algorithm>
#include <fstream>
#include <sstream>

using namespace GClasses;
using std::vector;
using std::string;
using std::ostream;

// static
smart_ptr<GRelation> GRelation::fromTwt(GTwtNode* pNode)
{
	if(pNode->fieldIfExists("attrs"))
	{
		sp_relation sp;
		sp = new GUniformRelation(pNode);
		return sp;
	}
	else
	{
		sp_relation sp;
		sp = new GMixedRelation(pNode);
		return sp;
	}
}

// virtual
void GRelation::print(ostream& stream, GData* pData, int precision)
{
	stream.precision(precision);

	// Write the relation title
	stream << "@RELATION Untitled\n\n";

	// Write the attributes
	for(int i = 0; i < size(); i++)
	{
		stream << "@ATTRIBUTE Attr" << i << "\t";
		if(valueCount(i) == 0)
			stream << "real";
		else
		{
			stream << "{";
			for(int j = 0; j < valueCount(i); j++)
			{
				if(j > 0)
					stream << ",";
				printValue(stream, j, i);
			}
			stream << "}";
		}
		stream << "\n";
	}

	// Write the data
	stream << "\n@DATA\n";
	if(!pData)
		return;
	for(size_t i = 0; i < pData->rows(); i++)
		printRow(stream, pData->row(i), ",");
}

// virtual
void GRelation::attrValue(string* pOutString, int attr, double value)
{
	std::ostringstream os;
	printValue(os, value, attr);
	*pOutString += os.str();
}

// virtual
void GRelation::printValue(ostream& stream, double value, int column)
{
	int valCount = valueCount(column);
	if(valCount == 0)
	{
		if(value == UNKNOWN_REAL_VALUE)
			stream << "?";
		else
			stream << value;
	}
	else
	{
		int val = (int)value;
		if(val < 0)
			stream << "?";
		else if(val >= valCount)
			ThrowError("value out of range");
		else if(val < 26)
		{
			char tmp[2];
			tmp[0] = 'a' + val;
			tmp[1] = '\0';
			stream << tmp;
		}
		else
			stream << "_" << val;
	}
}

void GRelation::printRow(ostream& stream, double* pRow, const char* separator)
{
	for(int j = 0; j < size(); j++)
	{
		if(j > 0)
			stream << separator;
		printValue(stream, pRow[j], j);
	}
	stream << "\n";
}

int GRelation::countRealSpaceDims(int nFirstAttr, int nAttrCount)
{
	int nDims = 0;
	int i, nValues;
	for(i = 0; i < nAttrCount; i++)
	{
		nValues = valueCount(nFirstAttr + i);
		if(nValues == 0)
			nDims += 2;
		else
			nDims += nValues;
	}
	return nDims;
}

void GRelation::toRealSpace(const double* pIn, double* pOut, int nFirstAttr, int nAttrCount)
{
	int nDims = 0;
	int i, j, k, nValues;
	for(i = 0; i < nAttrCount; i++)
	{
		nValues = valueCount(nFirstAttr + i);
		if(nValues == 0)
		{
			pOut[nDims++] = pIn[i];
			if(pIn[i] == UNKNOWN_REAL_VALUE)
				pOut[nDims++] = UNKNOWN_REAL_VALUE;
			else
				pOut[nDims++] = pIn[i] * pIn[i];
		}
		else
		{
			k = nDims;
			for(j = 0; j < nValues; j++)
				pOut[nDims++] = 0;
			if(pIn[i] >= 0) // For unknown discrete values, set to all zeros.
			{
				GAssert(pIn[i] >= 0 && pIn[i] < nValues);
				pOut[k + (int)pIn[i]] = 1;
			}
		}
	}
}

void GRelation::fromRealSpace(const double* pIn, double* pOut, int nFirstAttr, int nAttrCount, GRand* pRand)
{
	int nDims = 0;
	int i, nValues;
	for(i = 0; i < nAttrCount; i++)
	{
		nValues = valueCount(nFirstAttr + i);
		if(nValues == 0)
		{
			pOut[i] = pIn[nDims++];
			nDims++;
		}
		else
		{
			pOut[i] = (double)GVec::indexOfMax(&pIn[nDims], nValues, pRand);
			nDims += nValues;
		}
	}
}

void GRelation::fromRealSpace(const double* pIn, GPrediction* pOut, int nFirstAttr, int nAttrCount)
{
	int nDims = 0;
	int i, nValues;
	for(i = 0; i < nAttrCount; i++)
	{
		nValues = valueCount(nFirstAttr + i);
		if(nValues == 0)
		{
			GNormalDistribution* pNorm = pOut[i].makeNormal();
			double mean = pIn[nDims++];
			double variance = pIn[nDims++] - (mean * mean);
			pNorm->setMeanAndVariance(mean, variance);
		}
		else
		{
			GCategoricalDistribution* pCat = pOut[i].makeCategorical();
			pCat->setValues(nValues, &pIn[nDims]);
			nDims += nValues;
		}
	}
}

void GRelation::save(GData* pData, const char* szFilename, int precision)
{
	std::ofstream stream;
	stream.exceptions(std::ios::failbit|std::ios::badbit);
	try
	{
		stream.open(szFilename, std::ios::binary);
	}
	catch(const std::exception&)
	{
		ThrowError("Error creating file: ", szFilename);
	}
	print(stream, pData, precision);
}










GUniformRelation::GUniformRelation(GTwtNode* pNode)
{
	m_attrCount = (int)pNode->field("attrs")->asInt();
	m_valueCount = (int)pNode->field("vals")->asInt();
}

GTwtNode* GUniformRelation::toTwt(GTwtDoc* pDoc)
{
	GTwtNode* pNode = pDoc->newObj();
	pNode->addField(pDoc, "attrs", pDoc->newInt(m_attrCount));
	pNode->addField(pDoc, "vals", pDoc->newInt(m_valueCount));
	return pNode;
}

// virtual
void GUniformRelation::deleteAttribute(int index)
{
	if(m_attrCount <= 0)
		ThrowError("Index out of range");
	m_attrCount--;
}




GMixedRelation::GMixedRelation()
{
}

GMixedRelation::GMixedRelation(GTwtNode* pNode)
{
	m_valueCounts.clear();
	GTwtNode* pValueCounts = pNode->field("valueCounts");
	size_t itemCount = pValueCounts->itemCount();
	for(size_t i = 0; i < itemCount; i++)
		m_valueCounts.push_back((int)pValueCounts->item(i)->asInt());
}

GMixedRelation::GMixedRelation(GRelation* pCopyMe)
{
	copy(pCopyMe);
}

GMixedRelation::GMixedRelation(GRelation* pCopyMe, int firstAttr, int attrCount)
{
	addAttrs(pCopyMe, firstAttr, attrCount);
}

// virtual
GMixedRelation::~GMixedRelation()
{
}

GTwtNode* GMixedRelation::toTwt(GTwtDoc* pDoc)
{
	GTwtNode* pNode = pDoc->newObj();
 	GTwtNode* pValueCounts = pNode->addField(pDoc, "valueCounts", pDoc->newList((int)m_valueCounts.size()));
	size_t i;
	for(i = 0; i < m_valueCounts.size(); i++)
		pValueCounts->setItem((int)i, pDoc->newInt(m_valueCounts[i]));
	return pNode;
}

// virtual
GRelation* GMixedRelation::clone()
{
	GMixedRelation* pNewRelation = new GMixedRelation();
	pNewRelation->addAttrs(this, 0, size());
	return pNewRelation;
}

void GMixedRelation::addAttrs(GRelation* pCopyMe, int firstAttr, int attrCount)
{
	if(attrCount < 0)
		attrCount = pCopyMe->size() - firstAttr;
	int i;
	for(i = 0; i < attrCount; i++)
		copyAttr(pCopyMe, firstAttr + i);
}

void GMixedRelation::addAttrs(int attrCount, int valueCount)
{
	for(int i = 0; i < attrCount; i++)
		addAttr(valueCount);
}

void GMixedRelation::copy(GRelation* pCopyMe)
{
	flush();
	if(pCopyMe)
		addAttrs(pCopyMe);
}

// virtual
void GMixedRelation::flush()
{
	m_valueCounts.clear();
}

void GMixedRelation::addAttr(int nValues)
{
	if(type() == ARFF)
		((GArffRelation*)this)->addAttribute("No_name", nValues, NULL);
	else
		m_valueCounts.push_back(nValues);
}

// virtual
void GMixedRelation::copyAttr(GRelation* pThat, int nAttr)
{
	addAttr(pThat->valueCount(nAttr));
}

// virtual
bool GMixedRelation::areContinuous(int first, int count)
{
	int c = first;
	for(int i = 0; i < count; i++)
	{
		if(valueCount(c) != 0)
			return false;
		c++;
	}
	return true;
}

// virtual
bool GMixedRelation::areNominal(int first, int count)
{
	int c = first;
	for(int i = 0; i < count; i++)
	{
		if(valueCount(c) == 0)
			return false;
		c++;
	}
	return true;
}

// virtual
void GMixedRelation::swapAttributes(int nAttr1, int nAttr2)
{
	std::swap(m_valueCounts[nAttr1], m_valueCounts[nAttr2]);
}

// virtual
void GMixedRelation::deleteAttribute(int nAttr)
{
	m_valueCounts.erase(m_valueCounts.begin() + nAttr);
}

// virtual
void GMixedRelation::setAttrValueCount(int nAttr, int nValues)
{
	m_valueCounts[nAttr] = nValues;
}

// ------------------------------------------------------------------

GArffRelation::GArffRelation()
{
}

GArffRelation::~GArffRelation()
{
}

// virtual
void GArffRelation::flush()
{
	GAssert(m_attrs.size() == m_valueCounts.size());
	m_attrs.clear();
	GMixedRelation::flush();
}

// virtual
GRelation* GArffRelation::clone()
{
	GArffRelation* pNewRelation = new GArffRelation();
	pNewRelation->addAttrs(this);
	pNewRelation->setName(name());
	return pNewRelation;
}

void GArffRelation::addAttributeInternal(const char* pName, int nameLen, int valueCount)
{
	GAssert(m_attrs.size() == m_valueCounts.size());
	size_t index = m_valueCounts.size();
	m_valueCounts.push_back(valueCount);
	m_attrs.resize(index + 1);
	m_attrs[index].m_name.append(pName, nameLen);
}

void GArffRelation::addAttributeInternal(const char* pName, int nameLen, const char* pValues, int valuesLen)
{
	// Count the commas and add the attribute
	int valueCount = 0;
	for(int n = 0; n < valuesLen; n++)
	{
		if(pValues[n] == ',')
			valueCount++;
	}
	valueCount++;
	addAttributeInternal(pName, nameLen, valueCount);

	// Parse the values
	int index = (int)m_valueCounts.size() - 1;
	int start = 0;
	int count = 0;
	for(int n = 0; n <= valuesLen; n++)
	{
		if(pValues[n] == ',' || pValues[n] == '}')
		{
			int nStart = start;
			int nEnd = n;
			while(nStart < nEnd && pValues[nStart] <= ' ')
				nStart++;
			while(nStart < nEnd && pValues[nEnd - 1] <= ' ')
				nEnd--;
			string s = "";
			s.append(pValues + nStart, nEnd - nStart);
			m_attrs[index].m_values.push_back(s);
			count++;
			start = n + 1;
		}
	}
}

void GArffRelation::addAttribute(const char* szName, int nValues, vector<const char*>* pValues)
{
	GAssert(m_attrs.size() == m_valueCounts.size());
	int index = (int)m_valueCounts.size();
	m_valueCounts.push_back(nValues);
	m_attrs.resize(index + 1);
	if(szName)
		m_attrs[index].m_name = szName;
	if(pValues)
	{
		if(nValues != (int)pValues->size())
			ThrowError("mismatching value counts");
		int i;
		for(i = 0; i < nValues; i++)
			m_attrs[index].m_values.push_back((*pValues)[i]);
	}
}

// virtual
void GArffRelation::copyAttr(GRelation* pThat, int nAttr)
{
	if(pThat->type() == ARFF)
	{
		int index = (int)m_valueCounts.size();
		GArffRelation* pOther = (GArffRelation*)pThat;
		addAttributeInternal(pOther->m_attrs[nAttr].m_name.c_str(), (int)pOther->m_attrs[nAttr].m_name.length(), pOther->m_valueCounts[nAttr]);
		for(size_t i = 0; i < pOther->m_attrs[nAttr].m_values.size(); i++)
			m_attrs[index].m_values.push_back(pOther->m_attrs[nAttr].m_values[i]);
	}
	else
		addAttribute(NULL, pThat->valueCount(nAttr), NULL);
}

void GArffRelation::setName(const char* szName)
{
	m_name = szName;
}

void GArffRelation::parseAttribute(const char* szFile, size_t nLen, int nLine)
{
	// Eat whitespace
	while(nLen > 0 && *szFile <= ' ')
	{
		if(*szFile == '\n')
			ThrowError("Unexpected end of file at line ", gformat(nLine));
		szFile++;
		nLen--;
	}
	if(nLen == 0)
		ThrowError("Unexpected end of file at line ", gformat(nLine));

	// Parse the name
	int nQuotes = 0;
	if(szFile[0] == '\'' || szFile[0] == '"')
		nQuotes = 1;
	size_t nPos = 1;
	for( ; nPos < nLen && (((unsigned char)szFile[nPos] > (unsigned char)' ' && szFile[nPos] != '{') || nQuotes > 0); nPos++)
	{
		if(szFile[nPos] == '\'' || szFile[nPos] == '"')
			nQuotes--;
	}
	const char* pName = szFile;
	size_t nNameLen = nPos;

	// Eat whitespace
	while(nPos < nLen && szFile[nPos] <= ' ')
	{
		if(szFile[nPos] == '\n')
			ThrowError("Expected a value type at line ", gformat(nLine));
		nPos++;
	}
	if(nPos >= nLen)
		ThrowError("Expected a value type at line ", gformat(nLine));

	// Check for CONTINUOUS
	if(nLen - nPos >= 10 && _strnicmp(&szFile[nPos], "CONTINUOUS", 10) == 0)
	{
		addAttributeInternal(pName, (int)nNameLen, 0);
		return;
	}
	if(nLen - nPos >= 7 && _strnicmp(&szFile[nPos], "INTEGER", 7) == 0)
	{
		addAttributeInternal(pName, (int)nNameLen, 0);
		return;
	}
	if(nLen - nPos >= 7 && _strnicmp(&szFile[nPos], "NUMERIC", 7) == 0)
	{
		addAttributeInternal(pName, (int)nNameLen, 0);
		return;
	}
	if(nLen - nPos >= 4 && _strnicmp(&szFile[nPos], "REAL", 4) == 0)
	{
		addAttributeInternal(pName, (int)nNameLen, 0);
		return;
	}

	// Parse the values
	if(nPos >= nLen || szFile[nPos] != '{')
		ThrowError("Unrecognized value type at line ", gformat(nLine));
	nPos++;

	// Find the end of the values
	size_t n;
	for(n = nPos; n < nLen && szFile[n] != '}' && szFile[n] != '\n'; n++)
	{
	}
	if(n >= nLen || szFile[n] != '}')
		ThrowError("Expected a '}' at line ", gformat(nLine));

	addAttributeInternal(pName, (int)nNameLen, szFile + nPos, (int)(n - nPos));
}

// virtual
void GArffRelation::print(ostream& stream, GData* pData, int precision)
{
	// Write the relation title
	stream.precision(precision);
	stream << "@RELATION ";
	const char* szName = name();
	if(!szName)
		szName = "Untitled";
	stream << szName << "\n\n";

	// Write the attributes
	for(int i = 0; i < size(); i++)
	{
		stream << "@ATTRIBUTE ";
		stream << attrName(i) << "\t";
		if(valueCount(i) == 0)
			stream << "real";
		else
		{
			stream << "{";
			for(int j = 0; j < valueCount(i); j++)
			{
				if(j > 0)
					stream << ",";
				printValue(stream, j, i);
			}
			stream << "}";
		}
		stream << "\n";
	}

	// Write the data
	stream << "\n@DATA\n";
	if(!pData)
		return;
	for(size_t i = 0; i < pData->rows(); i++)
		printRow(stream, pData->row(i), ",");
}

// virtual
void GArffRelation::printValue(ostream& stream, double value, int column)
{
	int valCount = valueCount(column);
	if(valCount == 0)
	{
		if(value == UNKNOWN_REAL_VALUE)
			stream << "?";
		else
			stream << value;
	}
	else
	{
		int val = (int)value;
		if(val < 0)
			stream << "?";
		else if(val >= valCount)
			ThrowError("value out of range");
		else if(m_attrs[column].m_values.size() > 0)
			stream << m_attrs[column].m_values[val];
		else if(val < 26)
		{
			char tmp[2];
			tmp[0] = 'a' + val;
			tmp[1] = '\0';
			stream << tmp;
		}
		else
			stream << "_" << val;
	}
}

int GArffRelation::findEnumeratedValue(int nAttr, const char* szValue)
{
	size_t nValueCount = valueCount(nAttr);
	size_t actualValCount = m_attrs[nAttr].m_values.size();
	if(nValueCount > actualValCount)
		ThrowError("some values have no names");
	size_t i;
	for(i = 0; i < nValueCount; i++)
	{
		if(_stricmp(m_attrs[nAttr].m_values[i].c_str(), szValue) == 0)
			return (int)i;
	}
	return UNKNOWN_DISCRETE_VALUE;
}

inline bool IsRealValue(const char* szValue)
{
	if(*szValue == '-')
		szValue++;
	if(*szValue == '.')
		szValue++;
	if(*szValue >= '0' && *szValue <= '9')
		return true;
	return false;
}

const char* GArffRelation::attrName(int nAttr)
{
	return m_attrs[nAttr].m_name.c_str();
}

void GArffRelation::addAttrValue(int nAttr, const char* szValue)
{
	m_attrs[nAttr].m_values.push_back(szValue);
}

// virtual
void GArffRelation::setAttrValueCount(int nAttr, int nValues)
{
	m_attrs[nAttr].m_values.clear();
	GMixedRelation::setAttrValueCount(nAttr, nValues);
}

// virtual
void GArffRelation::swapAttributes(int nAttr1, int nAttr2)
{
	GMixedRelation::swapAttributes(nAttr1, nAttr2);
	std::swap(m_attrs[nAttr1], m_attrs[nAttr2]);
}

// virtual
void GArffRelation::deleteAttribute(int nAttr)
{
	m_attrs.erase(m_attrs.begin() + nAttr);
	GMixedRelation::deleteAttribute(nAttr);
}

double GArffRelation::parseValue(int attr, const char* val)
{
	int values = valueCount(attr);
	if(values == 0)
	{
		if(strcmp(val, "?") == 0)
			return UNKNOWN_REAL_VALUE;
		else
			return atof(val);
	}
	else
	{
		if(strcmp(val, "?") == 0)
			return UNKNOWN_DISCRETE_VALUE;
		else
		{
			int v = -1;
			for(int j = 0; j < values; j++)
			{
				if(_stricmp(val, m_attrs[attr].m_values[j].c_str()) == 0)
				{
					v = j;
					break;
				}
			}
			if(v < 0)
				v = atoi(val);
			return (double)v;
		}
	}
}

// ------------------------------------------------------------------

GData::GData(sp_relation& pRelation, GHeap* pHeap)
: m_pRelation(pRelation), m_pHeap(pHeap)
{
}

GData::GData(int attrs, GHeap* pHeap)
: m_pHeap(pHeap)
{
	m_pRelation = new GUniformRelation(attrs, 0);
}

GData::GData(GTwtNode* pNode, GHeap* pHeap)
: m_pHeap(pHeap)
{
	m_pRelation = GRelation::fromTwt(pNode->field("rel"));
	GTwtNode* pPats = pNode->field("pats");
	size_t size = pPats->itemCount();
	reserve(size);
	size_t dims = (size_t)m_pRelation->size();
	double* pPat;
	for(size_t i = 0; i < size; i++)
	{
		GTwtNode* pRow = pPats->item(i);
		if(pRow->itemCount() != dims)
			ThrowError("Row ", gformat(i), " has an unexpected number of values");
		pPat = newRow();
		for(size_t j = 0; j < dims; j++)
		{
			*pPat = pRow->item(j)->asDouble();
			pPat++;
		}
	}
}

GData::~GData()
{
	flush();
}

void GData::flush()
{
	if(!m_pHeap)
	{
		for(vector<double*>::iterator it = m_rows.begin(); it != m_rows.end(); it++)
			delete[] (*it);
	}
	m_rows.clear();
}

// static
GData* GData::loadArff(const char* szFilename)
{
	size_t nLen;
	char* szFile = GFile::loadFile(szFilename, &nLen);
	ArrayHolder<char> hFile(szFile);
	return parseArff(szFile, nLen);
}

// static
GData* GData::loadCsv(const char* szFilename, char separator, bool columnNamesInFirstRow, bool tolerant)
{
	size_t nLen;
	char* szFile = GFile::loadFile(szFilename, &nLen);
	ArrayHolder<char> hFile(szFile);
	return parseCsv(szFile, nLen, separator, columnNamesInFirstRow, tolerant);
}

void GData::saveArff(const char* szFilename)
{
	m_pRelation->save(this, szFilename, 14);
}

// static
void GData::parseDataRow(GArffRelation* pRelation, GData* pDataSet, const char* szFile, size_t nLen, int nLine)
{
	char szBuf[512];
	int nAttributeCount = pRelation->size();
	double* pData = new double[nAttributeCount];
	pDataSet->m_rows.push_back(pData);
	int col = 0;
	int n;
	for(n = 0; n < nAttributeCount; n++)
	{
		// Eat whitespace
		while(nLen > 0 && *szFile <= ' ')
		{
			if(*szFile == '\n')
				break;
			szFile++;
			nLen--;
		}

		// Parse the next value
		size_t nPos;
		for(nPos = 0; nPos < nLen; nPos++)
		{
			if(szFile[nPos] == ',' || szFile[nPos] == '\t' || szFile[nPos] == '\n')
				break;
		}
		if(szFile[nPos] == ',' && n >= nAttributeCount - 1)
			ThrowError("Too many values on line ", gformat(nLine));
		if(szFile[nPos] == '\n' && n < nAttributeCount - 1)
			ThrowError("Expected more values on line ", gformat(nLine));
		if(pRelation->valueCount(n) >= 0) // if it's not a comment attribute
		{
			int nEnd;
			for(nEnd = nPos; nEnd > 0 && szFile[nEnd - 1] <= ' '; nEnd--)
			{
			}
			memcpy(szBuf, szFile, nEnd);
			szBuf[nEnd] = '\0';
			if(pRelation->valueCount(n) == 0)
			{
				if(nEnd == 0 || (nEnd == 1 && szBuf[0] == '?'))
					pData[col++] = UNKNOWN_REAL_VALUE;
				else
				{
					// Parse a continuous value
					if(IsRealValue(szBuf))
						pData[col++] = atof(szBuf);
					else
						ThrowError("Expected a numeric value at line ", gformat(nLine), " attribute ", gformat(n + 1));
				}
			}
			else
			{
				if(nEnd == 0 || (nEnd == 1 && szBuf[0] == '?'))
					pData[col++] = UNKNOWN_DISCRETE_VALUE;
				else
				{
					// Parse an enumerated value
					int nVal = pRelation->findEnumeratedValue(n, szBuf);
					if(nVal < 0)
						ThrowError("Unrecognized enumeration value at line ", gformat(nLine), " attribute ", gformat(n + 1));
					pData[col++] = nVal;
				}
			}
		}

		// Advance past the attribute
		if(nPos < nLen)
			nPos++;
		while(nPos > 0)
		{
			szFile++;
			nPos--;
			nLen--;
		}
	}
	GAssert(col == nAttributeCount); // something got off count
}

// static
GData* GData::parseArff(const char* szFile, size_t nLen)
{
	// Parse the relation name
	size_t nPos = 0;
	int nLine = 1;
	GArffRelation* pRelation = new GArffRelation();
	sp_relation sp_rel;
	sp_rel = pRelation;
	while(true)
	{
		// Skip Whitespace
		while(nPos < nLen && szFile[nPos] <= ' ')
		{
			if(szFile[nPos] == '\n')
				nLine++;
			nPos++;
		}
		if(nPos >= nLen)
			ThrowError("Expected @RELATION at line ", gformat(nLine));

		// Check for comments
		if(szFile[nPos] == '%')
		{
			for(nPos++; nPos < nLen && szFile[nPos] != '\n'; nPos++)
			{
			}
			continue;
		}

		// Parse Relation
		if(nLen - nPos < 9 || _strnicmp(&szFile[nPos], "@RELATION", 9) != 0)
			ThrowError("Expected @RELATION at line ", gformat(nLine));
		nPos += 9;

		// Skip Whitespace
		while(szFile[nPos] <= ' ' && nPos < nLen)
		{
			if(szFile[nPos] == '\n')
				nLine++;
			nPos++;
			break;
		}
		if(nPos >= nLen)
			ThrowError("Expected relation name at line ", gformat(nLine));

		// Parse Name
		size_t nNameStart = nPos;
		while(szFile[nPos] > ' ' && nPos < nLen)
			nPos++;
		string s;
		s.assign(szFile + nNameStart, nPos - nNameStart);
		pRelation->setName(s.c_str());
		break;
	}

	// Parse the attribute section
	while(true)
	{
		// Skip Whitespace
		while(nPos < nLen && szFile[nPos] <= ' ')
		{
			if(szFile[nPos] == '\n')
				nLine++;
			nPos++;
		}
		if(nPos >= nLen)
			ThrowError("Expected @ATTRIBUTE or @DATA at line ", gformat(nLine));

		// Check for comments
		if(szFile[nPos] == '%')
		{
			for(nPos++; szFile[nPos] != '\n' && nPos < nLen; nPos++)
			{
			}
			continue;
		}

		// Check for @DATA
		if(nLen - nPos < 5) // 10 = strlen("@DATA")
			ThrowError("Expected @DATA at line ", gformat(nLine));
		if(_strnicmp(&szFile[nPos], "@DATA", 5) == 0)
		{
			nPos += 5;
			break;
		}

		// Parse @ATTRIBUTE
		if(nLen - nPos < 10) // 10 = strlen("@ATTRIBUTE")
			ThrowError("Expected @ATTRIBUTE at line ", gformat(nLine));
		if(_strnicmp(&szFile[nPos], "@ATTRIBUTE", 10) != 0)
			ThrowError("Expected @ATTRIBUTE or @DATA at line ", gformat(nLine));
		nPos += 10;
		pRelation->parseAttribute(&szFile[nPos], nLen - nPos, nLine);

		// Move to next line
		for(nPos++; szFile[nPos] != '\n' && nPos < nLen; nPos++)
		{
		}
	}

	// Parse the data section
	GData* pData = new GData(sp_rel);
	Holder<GData> hData(pData);
	while(true)
	{
		// Skip Whitespace
		while(nPos < nLen && szFile[nPos] <= ' ')
		{
			if(szFile[nPos] == '\n')
				nLine++;
			nPos++;
		}
		if(nPos >= nLen)
			break;

		// Check for comments
		if(szFile[nPos] == '%')
		{
			for(nPos++; szFile[nPos] != '\n' && nPos < nLen; nPos++)
			{
			}
			continue;
		}

		
		// Parse the data line
		parseDataRow(pRelation, pData, &szFile[nPos], nLen - nPos, nLine);

		// Move to next line
		for(nPos++; szFile[nPos] != '\n' && nPos < nLen; nPos++)
		{
		}
		continue;
	}

	return hData.release();
}

class ImportRow
{
public:
	vector<const char*> m_elements;
};

// static
GData* GData::parseCsv(const char* pFile, size_t len, char separator, bool columnNamesInFirstRow, bool tolerant)
{
	// Extract the elements
	GHeap heap(2048);
	vector<ImportRow> rows;
	int elementCount = -1;
	int nFirstDataLine = 1;
	int nLine = 1;
	size_t nPos = 0;
	while(true)
	{
		// Skip Whitespace
		while(nPos < len && pFile[nPos] <= ' ' && pFile[nPos] != separator)
		{
			if(pFile[nPos] == '\n')
				nLine++;
			nPos++;
		}
		if(nPos >= len)
			break;

		// Count the elements
		if(elementCount < 0)
		{
			if(separator == '\0')
			{
				// Elements are separated by an arbitrary amount of whitespace, element values contain no whitespace, and there are no missing elements
				size_t i = nPos;
				elementCount = 0;
				while(true)
				{
					elementCount++;
					while(i < len && pFile[i] > ' ')
						i++;
					while(i < len && pFile[i] <= ' ' && pFile[i] != '\n')
						i++;
					if(pFile[i] == '\n')
						break;
				}
			}
			else
			{
				// Elements are separated by the specified character
				nFirstDataLine = nLine;
				elementCount = 1;
				for(int i = 0; pFile[nPos + i] != '\n' && pFile[nPos + i] != '\0'; i++)
				{
					if(pFile[nPos + i] == separator)
						elementCount++;
				}
			}
		}

		// Extract the elements from the row
		rows.resize(rows.size() + 1);
		ImportRow& row = rows[rows.size() - 1];
		while(true)
		{
			// Skip Whitespace
			while(nPos < len && pFile[nPos] <= ' ' && pFile[nPos] != separator)
			{
				if(pFile[nPos] == '\n')
					break;
				nPos++;
			}
			if(nPos >= len || pFile[nPos] == '\n')
				break;

			// Extract the element
			int i, l;
			if(separator == '\0')
			{
				for(l = 0; pFile[nPos + l] > ' '; l++)
				{
				}
				for(i = l; pFile[nPos + i] <= ' ' && pFile[nPos + i] != '\n' && pFile[nPos + i] != '\0'; i++)
				{
				}
			}
			else
			{
				for(i = 0; pFile[nPos + i] != separator && pFile[nPos + i] != '\n' && pFile[nPos + i] != '\0'; i++)
				{
				}
				for(l = i; l > 0 && pFile[nPos + l - 1] <= ' '; l--)
				{
				}
			}
			char* el = heap.add(pFile + nPos, l);
			row.m_elements.push_back(el);
			if((int)row.m_elements.size() > elementCount)
				break;
			nPos += i;
			if(separator != '\0' && pFile[nPos] == separator)
				nPos++;
		}
		if(tolerant)
		{
			while((int)row.m_elements.size() < elementCount)
				row.m_elements.push_back("");
		}
		else
		{
			if(row.m_elements.size() != (size_t)elementCount)
				ThrowError("Line ", gformat(nLine), " has a different number of elements than line ", gformat(nFirstDataLine));
		}

		// Move to next line
		for(; nPos < len && pFile[nPos] != '\n'; nPos++)
		{
		}
		continue;
	}

	// Parse it all
	GArffRelation* pRelation = new GArffRelation();
	sp_relation pRel;
	pRel = pRelation;
	GData* pData = new GData(pRel);
	Holder<GData> hData(pData);
	size_t rowCount = rows.size();
	if(columnNamesInFirstRow)
		rowCount--;
	pData->reserve(rowCount);
	for(size_t i = 0; i < rowCount; i++)
		pData->m_rows.push_back(new double[elementCount]);
	for(int attr = 0; attr < elementCount; attr++)
	{
		// Determine if the attribute can be real
		bool real = true;
		for(size_t pat = columnNamesInFirstRow ? 1 : 0; pat < rows.size(); pat++)
		{
			const char* el = rows[pat].m_elements[attr];
			if(el[0] == '\0')
				continue; // unknown value
			if(strcmp(el, "?") == 0)
				continue; // unknown value
			if(GBits::isValidFloat(el, (int)strlen(el)))
				continue;
			real = false;
			break;
		}

		// Make the attribute
		if(real)
		{
			if(columnNamesInFirstRow)
				pRelation->addAttribute(rows[0].m_elements[attr], 0, NULL);
			else
			{
				string attrName = "attr";
				attrName += gformat(attr);
				pRelation->addAttribute(attrName.c_str(), 0, NULL);
			}
			size_t i = 0;
			for(size_t pat = columnNamesInFirstRow ? 1 : 0; pat < rows.size(); pat++)
			{
				const char* el = rows[pat].m_elements[attr];
				double val;
				if(el[0] == '\0')
					val = UNKNOWN_REAL_VALUE;
				else if(strcmp(el, "?") == 0)
					val = UNKNOWN_REAL_VALUE;
				else
					val = atof(el);
				pData->row(i)[attr] = val;
				i++;
			}
		}
		else
		{
			// Make the data
			vector<const char*> values;
			GConstStringHashTable ht(31, true);
			void* pVal;
			uintptr_t n;
			size_t i = 0;
			size_t valueCount = 0;
			for(size_t pat = columnNamesInFirstRow ? 1 : 0; pat < rows.size(); pat++)
			{
				const char* el = rows[pat].m_elements[attr];
				if(el[0] == '\0')
					pData->row(i)[attr] = UNKNOWN_DISCRETE_VALUE;
				else if(strcmp(el, "?") == 0)
					pData->row(i)[attr] = UNKNOWN_DISCRETE_VALUE;
				else
				{
					if(ht.get(el, &pVal))
						n = (uintptr_t)pVal;
					else
					{
						values.push_back(el);
						n = valueCount++;
						ht.add(el, (const void*)n);
					}
					pData->row(i)[attr] = (double)n;
				}
				i++;
			}

			// Make the attribute
			if(columnNamesInFirstRow)
				pRelation->addAttribute(rows[0].m_elements[attr], (int)valueCount, &values);
			else
			{
				string attrName = "attr";
				attrName += gformat(attr);
				pRelation->addAttribute(attrName.c_str(), (int)valueCount, &values);
			}
		}
	}
	return hData.release();
}

GTwtNode* GData::toTwt(GTwtDoc* pDoc)
{
	GTwtNode* pData = pDoc->newObj();
	int attrCount = m_pRelation->size();
	pData->addField(pDoc, "rel", m_pRelation->toTwt(pDoc));
	GTwtNode* pPats = pData->addField(pDoc, "pats", pDoc->newList(rows()));
	GTwtNode* pRow;
	double* pPat;
	for(size_t i = 0; i < rows(); i++)
	{
		pPat = row(i);
		pRow = pPats->setItem(i, pDoc->newList(attrCount));
		for(int j = 0; j < attrCount; j++)
			pRow->setItem(j, pDoc->newDouble(pPat[j]));
	}
	return pData;
}

void GData::col(int index, double* pOutVector)
{
	for(size_t i = 0; i < rows(); i++)
		*(pOutVector++) = row(i)[index];
}

void GData::setCol(int index, const double* pVector)
{
	for(size_t i = 0; i < rows(); i++)
		row(i)[index] = *(pVector++);
}

void GData::add(GData* pThat, bool transpose)
{
	if(transpose)
	{
		size_t c = (size_t)cols();
		if(rows() != (size_t)pThat->cols() || c != pThat->rows())
			ThrowError("expected matrices of same size");
		for(size_t i = 0; i < rows(); i++)
		{
			double* pRow = row(i);
			for(size_t j = 0; j < c; j++)
				*(pRow++) += pThat->row(j)[i];
		}
	}
	else
	{
		int c = cols();
		if(rows() != pThat->rows() || c != pThat->cols())
			ThrowError("expected matrices of same size");
		for(size_t i = 0; i < rows(); i++)
			GVec::add(row(i), pThat->row(i), c);
	}
}

void GData::subtract(GData* pThat, bool transpose)
{
	if(transpose)
	{
		size_t c = (size_t)cols();
		if(rows() != (size_t)pThat->cols() || c != pThat->rows())
			ThrowError("expected matrices of same size");
		for(size_t i = 0; i < rows(); i++)
		{
			double* pRow = row(i);
			for(size_t j = 0; j < c; j++)
				*(pRow++) -= pThat->row(j)[i];
		}
	}
	else
	{
		int c = cols();
		if(rows() != pThat->rows() || c != pThat->cols())
			ThrowError("expected matrices of same size");
		for(size_t i = 0; i < rows(); i++)
			GVec::subtract(row(i), pThat->row(i), c);
	}
}

void GData::multiply(double scalar)
{
	int c = cols();
	for(size_t i = 0; i < rows(); i++)
		GVec::multiply(row(i), scalar, c);
}

void GData::multiply(const double* pVectorIn, double* pVectorOut, bool transpose)
{
	size_t rowCount = rows();
	int colCount = cols();
	if(transpose)
	{
		GVec::setAll(pVectorOut, 0.0, colCount);
		for(size_t i = 0; i < rowCount; i++)
			GVec::addScaled(pVectorOut, *(pVectorIn++), row(i), colCount);
	}
	else
	{
		for(size_t i = 0; i < rowCount; i++)
			*(pVectorOut++) = GVec::dotProduct(row(i), pVectorIn, colCount);
	}
}

// static
GData* GData::multiply(GData& a, GData& b, bool transposeA, bool transposeB)
{
	if(transposeA)
	{
		if(transposeB)
		{
			size_t dims = a.rows();
			if((size_t)b.cols() != dims)
				ThrowError("dimension mismatch");
			size_t w = b.rows();
			size_t h = a.cols();
			GData* pOut = new GData(w);
			pOut->newRows(h);
			for(size_t y = 0; y < h; y++)
			{
				double* pRow = pOut->row(y);
				for(size_t x = 0; x < w; x++)
				{
					double* pB = b[x];
					double sum = 0;
					for(size_t i = 0; i < dims; i++)
						sum += a[i][y] * pB[i];
					*(pRow++) = sum;
				}
			}
			return pOut;
		}
		else
		{
			size_t dims = a.rows();
			if(b.rows() != dims)
				ThrowError("dimension mismatch");
			size_t w = b.cols();
			size_t h = a.cols();
			GData* pOut = new GData(w);
			pOut->newRows(h);
			for(size_t y = 0; y < h; y++)
			{
				double* pRow = pOut->row(y);
				for(size_t x = 0; x < w; x++)
				{
					double sum = 0;
					for(size_t i = 0; i < dims; i++)
						sum += a[i][y] * b[i][x];
					*(pRow++) = sum;
				}
			}
			return pOut;
		}
	}
	else
	{
		if(transposeB)
		{
			size_t dims = (size_t)a.cols();
			if((size_t)b.cols() != dims)
				ThrowError("dimension mismatch");
			size_t w = b.rows();
			size_t h = a.rows();
			GData* pOut = new GData(w);
			pOut->newRows(h);
			for(size_t y = 0; y < h; y++)
			{
				double* pRow = pOut->row(y);
				double* pA = a[y];
				for(size_t x = 0; x < w; x++)
					*(pRow++) = GVec::dotProduct(pA, b[x], dims);
			}
			return pOut;
		}
		else
		{
			size_t dims = (size_t)a.cols();
			if(b.rows() != dims)
				ThrowError("dimension mismatch");
			size_t w = b.cols();
			size_t h = a.rows();
			GData* pOut = new GData(w);
			pOut->newRows(h);
			for(size_t y = 0; y < h; y++)
			{
				double* pRow = pOut->row(y);
				double* pA = a[y];
				for(size_t x = 0; x < w; x++)
				{
					double sum = 0;
					for(size_t i = 0; i < dims; i++)
						sum += pA[i] * b[i][x];
					*(pRow++) = sum;
				}
			}
			return pOut;
		}
	}
}

GData* GData::transpose()
{
	size_t r = rows();
	size_t c = (size_t)cols();
	GData* pTarget = new GData(r);
	pTarget->newRows(c);
	for(size_t i = 0; i < c; i++)
	{
		double* pRow = pTarget->row(i);
		for(size_t j = 0; j < r; j++)
			*(pRow++) = row(j)[i];
	}
	return pTarget;
}

double GData::trace()
{
	size_t min = MIN((size_t)cols(), rows());
	double sum = 0;
	for(size_t n = 0; n < min; n++)
		sum += row(n)[n];
	return sum;
}

int GData::toReducedRowEchelonForm()
{
	int nLead = 0;
	double* pRow;
	size_t rowCount = rows();
	int colCount = cols();
	for(size_t nRow = 0; nRow < rowCount; nRow++)
	{
		// Find the next pivot (swapping rows as necessary)
		size_t i = nRow;
		while(ABS(row(i)[nLead]) < 1e-9)
		{
			if(++i >= rowCount)
			{
				i = nRow;
				if(++nLead >= colCount)
					return nRow;
			}
		}
		if(i > nRow)
			swapRows(i, nRow);

		// Scale the pivot to 1
		pRow = row(nRow);
		GVec::multiply(pRow + nLead, 1.0 / pRow[nLead], colCount - nLead);

		// Elliminate all values above and below the pivot
		for(i = 0; i < rowCount; i++)
		{
			if(i != nRow)
				GVec::addScaled(row(i) + nLead, -row(i)[nLead], pRow + nLead, colCount - nLead);
		}

		nLead++;
	}
	return rowCount;
}

bool GData::gaussianElimination(double* pVector)
{
	if(rows() != (size_t)cols())
		ThrowError("Expected a square matrix");
	double d, dBest;
	double* pRow;
	size_t rowCount = rows();
	int colCount = cols();
	for(size_t nRow = 0; nRow < rowCount; nRow++)
	{
		dBest = 0;
		size_t i;
		for(i = nRow; i < rowCount && ABS(row(i)[nRow]) < 1e-4; i++)
		{
		}
		if(i >= rowCount)
			continue;
		if(i > nRow)
		{
			swapRows(i, nRow);
			d = pVector[i];
			pVector[i] = pVector[nRow];
			pVector[nRow] = d;
		}

		// Scale the pivot to 1
		pRow = row(nRow);
		d = 1.0 / pRow[nRow];
		GVec::multiply(pRow + nRow, d, colCount - nRow);
		pVector[nRow] *= d;

		// Elliminate all values above and below the pivot
		for(i = 0; i < rowCount; i++)
		{
			if(i != nRow)
			{
				d = -row(i)[nRow];
				GVec::addScaled(row(i) + nRow, d, pRow + nRow, colCount - nRow);
				pVector[i] += d * pVector[nRow];
			}
		}
	}

	// Arbitrarily assign null-space values to 1
	for(size_t nRow = 0; nRow < rowCount; nRow++)
	{
		if(row(nRow)[nRow] < 0.5)
		{
			if(ABS(pVector[nRow]) >= 1e-4)
				return false;
			for(size_t i = 0; i < rowCount; i++)
			{
				if(i == nRow)
				{
					pVector[nRow] = 1;
					row(nRow)[nRow] = 1;
				}
				else
				{
					pVector[i] -= row(i)[nRow];
					row(i)[nRow] = 0;
				}
			}
		}
	}
	return true;
}

GData* GData::cholesky()
{
	size_t rowCount = rows();
	size_t colCount = (size_t)cols();
	GData* pOut = new GData(m_pRelation);
	pOut->newRows(rowCount);
	double d;
	for(size_t j = 0; j < rowCount; j++)
	{
		size_t i;
		for(i = 0; i < j; i++)
		{
			d = 0;
			for(size_t k = 0; k < i; k++)
				d += (pOut->row(i)[k] * pOut->row(j)[k]);
			pOut->row(j)[i] = (1.0 / pOut->row(i)[i]) * (row(i)[j] - d);
		}
		d = 0;
		for(size_t k = 0; k < i; k++)
			d += (pOut->row(i)[k] * pOut->row(j)[k]);
		d = row(j)[i] - d;
		if(d < 0)
		{
			if(d > -1e-12)
				d = 0; // it's probably just rounding error
			else
				ThrowError("not positive definite");
		}
		pOut->row(j)[i] = sqrt(d);
		for(i++; i < colCount; i++)
			pOut->row(j)[i] = 0;
	}
	return pOut;
}
/*
void GData::invert()
{
	if(rows() != (size_t)cols())
		ThrowError("only square matrices supported");
	if(rows() == 1)
	{
		row(0)[0] = 1.0 / row(0)[0];
		return;
	}

	// Do LU decomposition (I think this is the Doolittle algorithm)
	int colCount = cols();
	double* pRow = row(0);
	for(int i = 1; i < colCount; i++)
		pRow[i] /= pRow[0];
	for(int i = 1; i < colCount; i++)
	{
		for(int j = i; j < colCount; j++)
		{ // do a column of L
			double sum = 0.0;
			for(int k = 0; k < i; k++)
				sum += row(j)[k] * row(k)[i];
			row(j)[i] -= sum;
		}
		if(i == colCount - 1)
			continue;
		for(int j = i + 1; j < colCount; j++)
		{ // do a row of U
			double sum = 0.0;
			for(int k = 0; k < i; k++)
				sum += row(i)[k] * row(k)[j];
			row(i)[j] = (row(i)[j] - sum) / row(i)[i];
		}
	}

	// Invert L
	for(int i = 0; i < colCount; i++)
	{
		for(int j = i; j < colCount; j++ )
		{
			double x = 1.0;
			if ( i != j )
			{
				x = 0.0;
				for(int k = i; k < j; k++ ) 
					x -= row(j)[k] * row(k)[i];
			}
			row(j)[i] = x / row(j)[j];
		}
	}

	// Invert U
	for(int i = 0; i < colCount; i++)
	{
		for(int j = i; j < colCount; j++ )
		{
			if( i == j )
				continue;
			double sum = 0.0;
			for (int k = i; k < j; k++ )
				sum += row(k)[j] * ((i == k) ? 1.0 : row(i)[k]);
			row(i)[j] = -sum;
		}
	}

	// A^-1 = U^-1 x L^-1
	for(int i = 0; i < colCount; i++ )
	{
		for(int j = 0; j < colCount; j++ )
		{
			double sum = 0.0;
			for(int k = ((i > j) ? i : j); k < colCount; k++)
				sum += ((j == k) ? 1.0 : row(j)[k]) * row(k)[i];
			row(j)[i] = sum;
		}
	}
}
*/
void GData::inPlaceSquareTranspose()
{
	size_t size = rows();
	if(size != (size_t)cols())
		ThrowError("Expected a square matrix");
	for(size_t a = 0; a < size; a++)
	{
		for(size_t b = a + 1; b < size; b++)
			std::swap(row(a)[b], row(b)[a]);
	}
}

double GData_pythag(double a, double b)
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

double GData_takeSign(double a, double b) 
{
	return (b >= 0.0 ? ABS(a) : -ABS(a));
}

void GData::singularValueDecomposition(GData** ppU, double** ppDiag, GData** ppV, bool throwIfNoConverge, int maxIters)
{
	if(rows() >= (size_t)cols())
		singularValueDecompositionHelper(ppU, ppDiag, ppV, throwIfNoConverge, maxIters);
	else
	{
		GData* pTemp = transpose();
		Holder<GData> hTemp(pTemp);
		pTemp->singularValueDecompositionHelper(ppV, ppDiag, ppU, throwIfNoConverge, maxIters);
		(*ppV)->inPlaceSquareTranspose();
		(*ppU)->inPlaceSquareTranspose();
	}
}

double GData_safeDivide(double n, double d)
{
	if(d == 0.0 && n == 0.0)
		return 0.0;
	else
	{
		double t = n / d;
		//GAssert(t > -1e200, "prob");
		return t;
	}
}

void GData::fixNans()
{
	int colCount = cols();
	for(size_t i = 0; i < rows(); i++)
	{
		double* pRow = row(i);
		for(int j = 0; j < colCount; j++)
		{
			if(*pRow >= -1e308 && *pRow < 1e308)
			{
			}
			else
				*pRow = (i == (size_t)j ? 1.0 : 0.0);
			pRow++;
		}
	}
}

void GData::singularValueDecompositionHelper(GData** ppU, double** ppDiag, GData** ppV, bool throwIfNoConverge, int maxIters)
{
	int m = (int)rows();
	int n = cols();
	if(m < n)
		ThrowError("Expected at least as many rows as columns");
	int i, j, k;
	int l = 0;
	int p, q, iter;
	double c, f, h, s, x, y, z;
	double norm = 0.0;
	double g = 0.0;
	double scale = 0.0;
	GData* pU = new GData(m);
	Holder<GData> hU(pU);
	pU->newRows(m);
	pU->setAll(0.0);
	pU->copyColumns(0, this, 0, n);
	double* pSigma = new double[n];
	ArrayHolder<double> hSigma(pSigma);
	GData* pV = new GData(n);
	Holder<GData> hV(pV);
	pV->newRows(n);
	pV->setAll(0.0);
	GTEMPBUF(double, temp, n);

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
			for(k = i; k < m; k++)
				scale += ABS(pU->row(k)[i]);
			if(scale != 0.0)
			{
				for(k = i; k < m; k++)
				{
					pU->row(k)[i] = GData_safeDivide(pU->row(k)[i], scale);
					double t = pU->row(k)[i];
					s += t * t;
				}
				f = pU->row(i)[i];
				g = -GData_takeSign(sqrt(s), f);
				h = f * g - s;
				pU->row(i)[i] = f - g;
				if(i != n - 1)
				{
					for(j = l; j < n; j++)
					{
						s = 0.0;
						for(k = i; k < m; k++) 
							s += pU->row(k)[i] * pU->row(k)[j];
						f = GData_safeDivide(s, h);
						for(k = i; k < m; k++)
							pU->row(k)[j] += f * pU->row(k)[i];
					}
				}
				for(k = i; k < m; k++)
					pU->row(k)[i] *= scale;
			}
		}
		pSigma[i] = scale * g;

		// Right-hand reduction
		g = 0.0;
		s = 0.0;
		scale = 0.0;
		if(i < m && i != n - 1) 
		{
			for(k = l; k < n; k++)
				scale += ABS(pU->row(i)[k]);
			if(scale != 0.0) 
			{
				for(k = l; k < n; k++) 
				{
					pU->row(i)[k] = GData_safeDivide(pU->row(i)[k], scale);
					double t = pU->row(i)[k];
					s += t * t;
				}
				f = pU->row(i)[l];
				g = -GData_takeSign(sqrt(s), f);
				h = f * g - s;
				pU->row(i)[l] = f - g;
				for(k = l; k < n; k++)
					temp[k] = GData_safeDivide(pU->row(i)[k], h);
				if(i != m - 1) 
				{
					for(j = l; j < m; j++) 
					{
						s = 0.0;
						for(k = l; k < n; k++)
							s += pU->row(j)[k] * pU->row(i)[k];
						for(k = l; k < n; k++)
							pU->row(j)[k] += s * temp[k];
					}
				}
				for(k = l; k < n; k++)
					pU->row(i)[k] *= scale;
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
				for(j = l; j < n; j++)
					pV->row(i)[j] = GData_safeDivide(GData_safeDivide(pU->row(i)[j], pU->row(i)[l]), g); // (double-division to avoid underflow)
				for(j = l; j < n; j++)
				{
					s = 0.0;
					for(k = l; k < n; k++)
						s += pU->row(i)[k] * pV->row(j)[k];
					for(k = l; k < n; k++)
						pV->row(j)[k] += s * pV->row(i)[k];
				}
			}
			for(j = l; j < n; j++)
			{
				pV->row(i)[j] = 0.0;
				pV->row(j)[i] = 0.0;
			}
		}
		pV->row(i)[i] = 1.0;
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
				pU->row(i)[j] = 0.0;
		}
		if(g != 0.0)
		{
			g = GData_safeDivide(1.0, g);
			if(i != n - 1)
			{
				for(j = l; j < n; j++)
				{
					s = 0.0;
					for(k = l; k < m; k++)
						s += pU->row(k)[i] * pU->row(k)[j];
					f = GData_safeDivide(s, pU->row(i)[i]) * g;
					for(k = i; k < m; k++)
						pU->row(k)[j] += f * pU->row(k)[i];
				}
			}
			for(j = i; j < m; j++)
				pU->row(j)[i] *= g;
		} 
		else 
		{
			for(j = i; j < m; j++)
				pU->row(j)[i] = 0.0;
		}
		pU->row(i)[i] += 1.0;
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
					h = GData_pythag(f, g);
					pSigma[i] = h;
					h = GData_safeDivide(1.0, h);
					c = g * h;
					s = -f * h;
					for(j = 0; j < m; j++)
					{
						y = pU->row(j)[q];
						z = pU->row(j)[i];
						pU->row(j)[q] = y * c + z * s;
						pU->row(j)[i] = z * c - y * s;
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
						pV->row(k)[j] *= -1.0; 
				}
				break;
			}
			if(throwIfNoConverge && iter >= maxIters)
				ThrowError("failed to converge");

			// Shift from bottom 2x2 minor
			x = pSigma[l];
			q = k - 1;
			y = pSigma[q];
			g = temp[q];
			h = temp[k];
			f = GData_safeDivide(((y - z) * (y + z) + (g - h) * (g + h)), (2.0 * h * y));
			g = GData_pythag(f, 1.0);
			f = GData_safeDivide(((x - z) * (x + z) + h * (GData_safeDivide(y, (f + GData_takeSign(g, f))) - h)), x);

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
				z = GData_pythag(f, h);
				temp[j] = z;
				c = GData_safeDivide(f, z);
				s = GData_safeDivide(h, z);
				f = x * c + g * s;
				g = g * c - x * s;
				h = y * s;
				y = y * c;
				for(p = 0; p < n; p++)
				{
					x = pV->row(j)[p];
					z = pV->row(i)[p];
					pV->row(j)[p] = x * c + z * s;
					pV->row(i)[p] = z * c - x * s;
				}
				z = GData_pythag(f, h);
				pSigma[j] = z;
				if(z != 0.0)
				{
					z = GData_safeDivide(1.0, z);
					c = f * z;
					s = h * z;
				}
				f = c * g + s * y;
				x = c * y - s * g;
				for(p = 0; p < m; p++)
				{
					y = pU->row(p)[j];
					z = pU->row(p)[i];
					pU->row(p)[j] = y * c + z * s;
					pU->row(p)[i] = z * c - y * s;
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
	pU->fixNans();
	pV->fixNans();
	*ppU = hU.release();
	*ppDiag = hSigma.release();
	*ppV = hV.release();
}

GData* GData::pseudoInverse()
{
	GData* pU;
	double* pDiag;
	GData* pV;
	int colCount = cols();
	size_t rowCount = rows();
	if(rowCount < (size_t)colCount)
	{
		GData* pTranspose = transpose();
		Holder<GData> hTranspose(pTranspose);
		pTranspose->singularValueDecompositionHelper(&pU, &pDiag, &pV, false, 80);
	}
	else
		singularValueDecompositionHelper(&pU, &pDiag, &pV, false, 80);
	Holder<GData> hU(pU);
	ArrayHolder<double> hDiag(pDiag);
	Holder<GData> hV(pV);
	GData sigma(rowCount < (size_t)colCount ? rowCount : colCount);
	sigma.newRows(rowCount < (size_t)colCount ? colCount : rowCount);
	sigma.setAll(0.0);
	int m = MIN((int)rowCount, colCount);
	for(int i = 0; i < m; i++)
	{
		if(ABS(pDiag[i]) > 1e-9)
			sigma[i][i] = GData_safeDivide(1.0, pDiag[i]);
		else
			sigma[i][i] = 0.0;
	}
	GData* pT = GData::multiply(*pU, sigma, false, false);
	Holder<GData> hT(pT);
	if(rowCount < (size_t)colCount)
		return GData::multiply(*pT, *pV, false, false);
	else
		return GData::multiply(*pV, *pT, true, true);
}

// static
GData* GData::kabsch(GData* pA, GData* pB)
{
	GData* pCovariance = GData::multiply(*pA, *pB, true, false);
	Holder<GData> hCov(pCovariance);
	GData* pU;
	double* pDiag;
	GData* pV;
	pCovariance->singularValueDecomposition(&pU, &pDiag, &pV);
	Holder<GData> hU(pU);
	delete[] pDiag;
	Holder<GData> hV(pV);
	GData* pK = GData::multiply(*pV, *pU, true, true);
	return pK;
}

// static
GData* GData::align(GData* pA, GData* pB)
{
	size_t columns = pA->cols();
	GTEMPBUF(double, mean, columns);
	pA->centroid(mean);
	GData* pAA = pA->clone();
	Holder<GData> hAA(pAA);
	pAA->centerMeanAtOrigin();
	GData* pBB = pB->clone();
	Holder<GData> hBB(pBB);
	pBB->centerMeanAtOrigin();
	GData* pK = GData::kabsch(pBB, pAA);
	Holder<GData> hK(pK);
	hAA.reset(NULL);
	GData* pAligned = GData::multiply(*pBB, *pK, false, true);
	Holder<GData> hAligned(pAligned);
	hBB.reset(NULL);
	for(vector<double*>::iterator it = pAligned->m_rows.begin(); it != pAligned->m_rows.end(); it++)
		GVec::add(*it, mean, columns);
	return hAligned.release();
}

double GData::determinant()
{
	// Check size
	int n = (int)rows();
	if(n != cols())
		ThrowError("Only square matrices are supported");

	// Convert to a triangular matrix
	double epsilon = 1e-10;
	GData* pC = this->clone();
	Holder<GData> hC(pC);
	GData& C = *pC;
	GTEMPBUF(int, Kp, 2 * n);
	int* Lp = Kp + n;
	int l, ko, lo;
	double po,t0;
	bool nonSingular = true;
	int k = 0;
	while(nonSingular && k < n)
	{
		po = C[k][k];
		lo = k;
		ko = k;
		for(int i = k; i < n; i++)
			for(int j = k; j < n; j++)
				if(ABS(C[i][j]) > ABS(po))
				{
					po = C[i][j];
					lo = i;
					ko = j;
				}
		Lp[k] = lo;
		Kp[k] = ko;
		if(ABS(po) < epsilon)
		{
			nonSingular = false;
			//ThrowError("Failed to compute determinant. Pivot too small.");
		}
		else
		{
			if(lo != k)
			{
				for(int j = k; j < n; j++)
				{
					t0 = C[k][j];
					C[k][j] = C[lo][j];
					C[lo][j] = t0;
				}
			}
			if(ko != k)
			{
				for(int i = 0; i < n; i++)
				{
					t0 = C[i][k];
					C[i][k] = C[i][ko];
					C[i][ko] = t0;
				}
			}
			for(int i = k + 1; i < n; i++)
			{
				C[i][k] /= po;
				for(int j = k + 1; j < n; j++)
					C[i][j] -= C[i][k] * C[k][j];
			}
			k++;
		}
	}
	if(nonSingular && ABS(C[n - 1][n - 1]) < epsilon)
		nonSingular = false;

	// Compute determinant
	if(!nonSingular)
		return 0.0;
	else
	{
		double det = 1.0;
		for(k = 0; k < n; k++)
			det *= C[k][k];
		l = 0;
		for(k = 0; k < n - 1; k++)
		{
			if(Lp[k] != k)
				l++;
			if(Kp[k] != k)
				l++;
		}
		if((l % 2) != 0)
			det = -det;
		return det;
	}
}

void GData::makeIdentity()
{
	size_t rowCount = rows();
	int colCount = cols();
	for(size_t nRow = 0; nRow < rowCount; nRow++)
		GVec::setAll(row(nRow), 0.0, colCount);
	size_t nMin = MIN((size_t)colCount, rowCount);
	for(size_t i = 0; i < nMin; i++)
		row(i)[i] = 1.0;
}

void GData::mirrorTriangle(bool upperToLower)
{
	size_t n = MIN(rows(), (size_t)cols());
	if(upperToLower)
	{
		for(size_t i = 0; i < n; i++)
		{
			for(size_t j = i + 1; j < n; j++)
				row(j)[i] = row(i)[j];
		}
	}
	else
	{
		for(size_t i = 0; i < n; i++)
		{
			for(size_t j = i + 1; j < n; j++)
				row(i)[j] = row(j)[i];
		}
	}
}

double GData::eigenValue(const double* pEigenVector)
{
	// Find the element with the largest magnitude
	int i, nEl = 0;
	int colCount = cols();
	for(i = 1; i < colCount; i++)
	{
		if(ABS(pEigenVector[i]) > ABS(pEigenVector[nEl]))
			nEl = i;
	}
	return GVec::dotProduct(row(nEl), pEigenVector, colCount) / pEigenVector[nEl];
}

void GData::eigenVector(double eigenvalue, double* pOutVector)
{
	GAssert(rows() == (size_t)cols()); // Expected a square matrix
	size_t rowCount = rows();
	for(size_t i = 0; i < rowCount; i++)
		row(i)[i] = row(i)[i] - eigenvalue;
	GVec::setAll(pOutVector, 0.0, rowCount);
	if(!gaussianElimination(pOutVector))
		ThrowError("no solution");
	GVec::normalize(pOutVector, rowCount);
}

GData* GData::eigs(int nCount, double* pEigenVals, GRand* pRand, bool mostSignificant)
{
	int dims = cols();
	if(nCount > dims)
		ThrowError("Can't have more eigenvectors than columns");
	if(rows() != (size_t)dims)
		ThrowError("expected a square matrix");

/*
	// The principle components of the Cholesky (square-root) matrix are the same as
	// the eigenvectors of this matrix.
	GData* pDeviation = cholesky();
	Holder<GData> hDeviation(pDeviation);
	GData* pData = pDeviation->transpose();
	Holder<GData> hData(pData);
	size_t s = pData->rows();
	for(size_t i = 0; i < s; i++)
	{
		double* pRow = pData->newRow();
		GVec::copy(pRow, pData->row(i), dims);
		GVec::multiply(pRow, -1, dims);
	}

	// Extract the principle components
	GData* pOut = new GData(m_pRelation);
	pOut->newRows(nCount);
	for(int i = 0; i < nCount; i++)
	{
		pData->principalComponentAboutOrigin(pOut->row(i), dims, pRand);
		pData->removeComponentAboutOrigin(pOut->row(i), dims);
	}
*/

	// Use the power method to compute the first few eigenvectors. todo: we really should use the Lanczos method instead
	GData* pOut = new GData(m_pRelation);
	pOut->newRows(nCount);
	GData* pA;
	if(mostSignificant)
		pA = clone();
	else
		pA = pseudoInverse();
	Holder<GData> hA(pA);
	GTEMPBUF(double, pTemp, dims);
	for(int i = 0; i < nCount; i++)
	{
		// Use the power method to compute the next eigenvector
		double* pX = pOut->row(i);
		pRand->spherical(pX, dims);
		for(int j = 0; j < 100; j++) // todo: is there a better way to detect convergence?
		{
			pA->multiply(pX, pTemp);
			GVec::copy(pX, pTemp, dims);
			GVec::safeNormalize(pX, dims, pRand);
		}

		// Compute the corresponding eigenvalue
		double lambda = pA->eigenValue(pX);
		if(pEigenVals)
			pEigenVals[i] = lambda;

		// Deflate (subtract out the eigenvector)
		for(int j = 0; j < dims; j++)
		{
			double* pRow = pA->row(j);
			for(int k = 0; k < dims; k++)
			{
				*pRow = *pRow - lambda * pX[j] * pX[k];
				pRow++;
			}
		}
	}

	return pOut;
}
/*
GData* GData::leastSignificantEigenVectors(int nCount, GRand* pRand)
{
	GData* pInv = clone();
	Holder<GData> hInv(pInv);
	pInv->invert();
	GData* pOut = pInv->mostSignificantEigenVectors(nCount, pRand);
	double eigenvalue;
	int i;
	for(i = 0; i < nCount; i++)
	{
		eigenvalue = 1.0 / pInv->eigenValue(pOut->row(i));
		GData* cp = clone();
		Holder<GData> hCp(cp);
		cp->eigenVector(eigenvalue, pOut->row(i));
	}
	return pOut;
}
*/
double* GData::newRow()
{
	int nAttributes = m_pRelation->size();
	double* pNewRow;
	if(m_pHeap)
		pNewRow = (double*)m_pHeap->allocate(sizeof(double) * nAttributes);
	else
		pNewRow = new double[nAttributes];
	m_rows.push_back(pNewRow);
	return pNewRow;
}

void GData::takeRow(double* pRow)
{
	m_rows.push_back(pRow);
}

void GData::newRows(size_t nRows)
{
	reserve(m_rows.size() + nRows);
	for(size_t i = 0; i < nRows; i++)
		newRow();
}

void GData::fromVector(const double* pVec, size_t nRows)
{
	if(rows() < nRows)
		newRows(nRows - rows());
	else
	{
		while(rows() > nRows)
			deleteRow(0);
	}
	int nCols = m_pRelation->size();
	for(size_t r = 0; r < nRows; r++)
	{
		double* pRow = row(r);
		GVec::copy(pRow, pVec, nCols);
		pVec += nCols;
	}
}

void GData::toVector(double* pVec)
{
	int nCols = cols();
	for(size_t i = 0; i < rows(); i++)
	{
		GVec::copy(pVec, row(i), nCols);
		pVec += nCols;
	}
}

void GData::setAll(double val)
{
	int colCount = cols();
	for(size_t i = 0; i < rows(); i++)
		GVec::setAll(row(i), val, colCount);
}

void GData::copy(GData* pThat)
{
	m_pRelation = pThat->m_pRelation;
	flush();
	newRows(pThat->rows());
	copyColumns(0, pThat, 0, m_pRelation->size());
}

GData* GData::clone()
{
	GData* pOther = new GData(relation());
	pOther->newRows(rows());
	pOther->copyColumns(0, this, 0, cols());
	return pOther;
}

void GData::copyRow(const double* pRow)
{
	double* pNewRow = newRow();
	GVec::copy(pNewRow, pRow, m_pRelation->size());
}

void GData::copyColumns(int nDestStartColumn, GData* pSource, int nSourceStartColumn, int nColumnCount)
{
	if(rows() != pSource->rows())
		ThrowError("expected datasets to have the same number of rows");
	size_t count = rows();
	for(size_t i = 0; i < count; i++)
		GVec::copy(row(i) + nDestStartColumn, pSource->row(i) + nSourceStartColumn, nColumnCount);
}

GData* GData::attrSubset(int firstAttr, int attrCount)
{
	if(firstAttr < 0 || attrCount < 0 || firstAttr + attrCount > m_pRelation->size())
		ThrowError("index out of range");
	sp_relation relNew;
	if(relation()->type() == GRelation::UNIFORM)
	{
		GUniformRelation* pNewRelation = new GUniformRelation(attrCount, relation()->valueCount(firstAttr));
		relNew = pNewRelation;
	}
	else
	{
		GMixedRelation* pNewRelation = new GMixedRelation();
		pNewRelation->addAttrs(m_pRelation.get(), firstAttr, attrCount);
		relNew = pNewRelation;
	}
	GData* pNewData = new GData(relNew);
	pNewData->newRows(rows());
	pNewData->copyColumns(0, this, firstAttr, attrCount);
	return pNewData;
}

double* GData::replaceRow(size_t nIndex, double* pRow)
{
	double* pOldRow = m_rows[nIndex];
	m_rows[nIndex] = pRow;
	return pOldRow;
}

void GData::swapRows(size_t a, size_t b)
{
	std::swap(m_rows[a], m_rows[b]);
}

void GData::swapColumns(int nAttr1, int nAttr2)
{
	if(nAttr1 == nAttr2)
		return;
	m_pRelation = m_pRelation->clone();
	m_pRelation->swapAttributes(nAttr1, nAttr2);
	size_t nCount = rows();
	double tmp;
	double* pRow;
	for(size_t i = 0; i < nCount; i++)
	{
		pRow = row(i);
		tmp = pRow[nAttr1];
		pRow[nAttr1] = pRow[nAttr2];
		pRow[nAttr2] = tmp;
	}
}

void GData::deleteColumn(int index)
{
	m_pRelation = m_pRelation->clone();
	m_pRelation->deleteAttribute(index);
	size_t nCount = rows();
	double* pRow;
	int nAttrCountBefore = m_pRelation->size();
	for(size_t i = 0; i < nCount; i++)
	{
		pRow = row(i);
		for(int j = index; j < nAttrCountBefore; j++)
			pRow[j] = pRow[j + 1];
	}
}

double* GData::releaseRow(size_t index)
{
	size_t last = m_rows.size() - 1;
	double* pRow = m_rows[index];
	m_rows[index] = m_rows[last];
	m_rows.pop_back();
	return pRow;
}

void GData::deleteRow(size_t index)
{
	double* pRow = releaseRow(index);
	if(!m_pHeap)
		delete[] pRow;
}

double* GData::releaseRowPreserveOrder(size_t index)
{
	double* pRow = m_rows[index];
	m_rows.erase(m_rows.begin() + index);
	return pRow;
}

void GData::deleteRowPreserveOrder(size_t index)
{
	double* pRow = releaseRowPreserveOrder(index);
	if(!m_pHeap)
		delete[] pRow;
}

void GData::releaseAllRows()
{
	m_rows.clear();
}

// static
GData* GData::mergeHoriz(GData* pSetA, GData* pSetB)
{
	if(pSetA->rows() != pSetB->rows())
		ThrowError("Expected same number of rows");
	GArffRelation* pRel = new GArffRelation();
	sp_relation spRel;
	spRel = pRel;
	GRelation* pRelA = pSetA->relation().get();
	GRelation* pRelB = pSetB->relation().get();
	int nSetADims = pRelA->size();
	int nSetBDims = pRelB->size();
	pRel->addAttrs(pRelA);
	pRel->addAttrs(pRelB);
	GData* pNewSet = new GData(spRel);
	Holder<GData> hNewSet(pNewSet);
	pNewSet->reserve(pSetA->rows());
	double* pNewRow;
	for(size_t i = 0; i < pSetA->rows(); i++)
	{
		pNewRow = pNewSet->newRow();
		GVec::copy(pNewRow, pSetA->row(i), nSetADims);
		GVec::copy(&pNewRow[nSetADims], pSetB->row(i), nSetBDims);
	}
	return hNewSet.release();
}

void GData::shuffle(GRand* pRand)
{
	for(size_t n = m_rows.size(); n > 0; n--)
		std::swap(m_rows[(size_t)pRand->next(n)], m_rows[n - 1]);
}

void GData::shuffle2(GRand* pRand, GData& other)
{
	for(size_t n = m_rows.size(); n > 0; n--)
	{
		size_t r = (size_t)pRand->next(n);
		std::swap(m_rows[r], m_rows[n - 1]);
		std::swap(other.m_rows[r], other.m_rows[n - 1]);
	}
}

void GData::shuffleLikeCards()
{
	for(size_t i = 0; i < rows(); i++)
	{
		size_t n = i;
		while(n & 1)
			n = (n >> 1);
		n = (n >> 1) + rows() / 2;
		std::swap(m_rows[i], m_rows[n]);
	}
}

double GData::entropy(int nColumn)
{
	// Count the number of occurrences of each value
	GAssert(m_pRelation->valueCount(nColumn) > 0); // continuous attributes are not supported
	int nPossibleValues = m_pRelation->valueCount(nColumn);
	GTEMPBUF(int, pnCounts, nPossibleValues);
	int nTotalCount = 0;
	memset(pnCounts, '\0', m_pRelation->valueCount(nColumn) * sizeof(int));
	size_t nRows = rows();
	for(size_t n = 0; n < nRows; n++)
	{
		int nValue = (int)row(n)[nColumn];
		if(nValue < 0)
		{
			GAssert(nValue == UNKNOWN_DISCRETE_VALUE);
			continue;
		}
		GAssert(nValue < nPossibleValues);
		pnCounts[nValue]++;
		nTotalCount++;
	}
	if(nTotalCount == 0)
		return 0;

	// Total up the entropy
	double dEntropy = 0;
	double dRatio;
	for(int n = 0; n < nPossibleValues; n++)
	{
		if(pnCounts[n] > 0)
		{
			dRatio = (double)pnCounts[n] / nTotalCount;
			dEntropy -= dRatio * log(dRatio);
		}
	}
	return M_LOG2E * dEntropy;
}

void GData::splitByPivot(GData* pGreaterOrEqual, int nAttribute, double dPivot)
{
	GAssert(pGreaterOrEqual->m_pHeap == m_pHeap);
	size_t nUnknowns = 0;
	double* pRow;
	size_t n;
	for(n = rows() - 1; n >= nUnknowns && n < rows(); n--)
	{
		pRow = row(n);
		if(pRow[nAttribute] == UNKNOWN_REAL_VALUE)
			std::swap(m_rows[nUnknowns++], m_rows[n++]);
		else if(pRow[nAttribute] >= dPivot)
			pGreaterOrEqual->takeRow(releaseRow(n));
	}

	// Send all the unknowns to the side with more rows
	if(pGreaterOrEqual->rows() > rows() - nUnknowns)
	{
		for(; n < rows(); n--)
			pGreaterOrEqual->takeRow(releaseRow(n));
	}
}

int DoubleRefComparer(void* pThis, void* pA, void* pB)
{
	if(*(double*)pA > *(double*)pB)
		return 1;
	if(*(double*)pA < *(double*)pB)
		return -1;
	return 0;
}

void GData::splitByDiscreteValue(GData* pSingleClass, int nAttr, int nValue)
{
	GAssert(pSingleClass->m_pHeap == m_pHeap);
	bool take = false;
	double* pVec;
	for(size_t i = rows() - 1; i < rows(); i--)
	{
		pVec = row(i);
		if((int)pVec[nAttr] == nValue)
		{
			pSingleClass->takeRow(releaseRow(i));
			take = true;
		}
		else if((int)pVec[nAttr] == UNKNOWN_DISCRETE_VALUE)
		{
			// Do whatever the previous row did
			if(take)
				pSingleClass->takeRow(releaseRow(i));
		}
		else
			take = false;
	}
}

void GData::splitBySize(GData* pOtherData, size_t nOtherRows)
{
	GAssert(pOtherData->m_pHeap == m_pHeap);
	GAssert(nOtherRows >= 0 && nOtherRows <= rows());
	size_t nTargetSize = pOtherData->rows() + nOtherRows;
	while(pOtherData->rows() < nTargetSize)
		pOtherData->takeRow(releaseRow(rows() - 1));
}

void GData::mergeVert(GData* pData)
{
	if(cols() != pData->cols())
		ThrowError("Different number of cols");
	while(pData->rows() > 0)
		takeRow(pData->releaseRow(0));
}

double GData::mean(int nAttribute)
{
	if(nAttribute >= cols() || nAttribute <  0)
		ThrowError("attribute index out of range");
	double sum = 0;
	size_t missing = 0;
	for(vector<double*>::iterator it = m_rows.begin(); it != m_rows.end(); it++)
	{
		if((*it)[nAttribute] == UNKNOWN_REAL_VALUE)
			missing++;
		else
			sum += (*it)[nAttribute];
	}
	size_t count = m_rows.size() - missing;
	if(count > 0)
		return sum / count;
	else
	{
		ThrowError("at least one value is required to compute a mean");
		return 0.0;
	}
}

double GData::median(int nAttribute)
{
	if(nAttribute >= cols() || nAttribute <  0)
		ThrowError("attribute index out of range");
	vector<double> vals;
	vals.reserve(rows());
	for(vector<double*>::iterator it = m_rows.begin(); it != m_rows.end(); it++)
	{
		double d = (*it)[nAttribute];
		if(d != UNKNOWN_REAL_VALUE)
			vals.push_back(d);
	}
	if(vals.size() < 1)
		ThrowError("at least one value is required to compute a median");
	std::sort(vals.begin(), vals.end());
	if(vals.size() & 1)
		return vals[vals.size() / 2];
	else
	{
		size_t half = vals.size() / 2;
		return 0.5 * (vals[half - 1] + vals[half]);
	}
}

void GData::centroid(double* pOutMeans)
{
	int c = cols();
	int n;
	for(n = 0; n < c; n++)
		pOutMeans[n] = mean(n);
}

double GData::variance(int nAttr, double mean)
{
	double d;
	double dSum = 0;
	size_t nMissing = 0;
	for(vector<double*>::iterator it = m_rows.begin(); it != m_rows.end(); it++)
	{
		if((*it)[nAttr] == UNKNOWN_REAL_VALUE)
		{
			nMissing++;
			continue;
		}
		d = (*it)[nAttr] - mean;
		dSum += (d * d);
	}
	size_t nCount = m_rows.size() - nMissing;
	if(nCount > 1)
		return dSum / (nCount - 1);
	else
		return 0; // todo: wouldn't UNKNOWN_REAL_VALUE be better here?
}

void GData::minAndRange(int nAttribute, double* pMin, double* pRange)
{
	double dMin = 1e300;
	double dMax = -1e300;
	for(vector<double*>::iterator it = m_rows.begin(); it != m_rows.end(); it++)
	{
		if((*it)[nAttribute] == UNKNOWN_REAL_VALUE)
			continue;
		if((*it)[nAttribute] < dMin)
			dMin = (*it)[nAttribute];
		if((*it)[nAttribute] > dMax)
			dMax = (*it)[nAttribute];
	}
	if(dMax >= dMin)
	{
		*pMin = dMin;
		*pRange = dMax - dMin;
	}
	else
	{
		*pMin = UNKNOWN_REAL_VALUE;
		*pRange = UNKNOWN_REAL_VALUE;
	}
}

void GData::minAndRangeUnbiased(int nAttribute, double* pMin, double* pRange)
{
	double min, range, d;
	minAndRange(nAttribute, &min, &range);
	d = .5 * (range * (rows() + 1) / (rows() - 1) - range);
	*pMin = (min - d);
	*pRange = (range + d);
}

void GData::normalize(int nAttribute, double dInputMin, double dInputRange, double dOutputMin, double dOutputRange)
{
	GAssert(dInputRange > 0);
	double dScale = dOutputRange / dInputRange;
	for(vector<double*>::iterator it = m_rows.begin(); it != m_rows.end(); it++)
	{
		(*it)[nAttribute] -= dInputMin;
		(*it)[nAttribute] *= dScale;
		(*it)[nAttribute] += dOutputMin;
	}
}

/*static*/ double GData::normalize(double dVal, double dInputMin, double dInputRange, double dOutputMin, double dOutputRange)
{
	GAssert(dInputRange > 0);
	dVal -= dInputMin;
	dVal /= dInputRange;
	dVal *= dOutputRange;
	dVal += dOutputMin;
	return dVal;
}

double GData::baselineValue(int nAttribute)
{
	if(m_pRelation->valueCount(nAttribute) == 0)
		return mean(nAttribute);
	int j, val;
	int nValues = m_pRelation->valueCount(nAttribute);
	GTEMPBUF(int, counts, nValues + 1);
	memset(counts, '\0', sizeof(int) * (nValues + 1));
	for(vector<double*>::iterator it = m_rows.begin(); it != m_rows.end(); it++)
	{
		val = (int)(*it)[nAttribute] + 1;
		counts[val]++;
	}
	val = 1;
	for(j = 2; j <= nValues; j++)
	{
		if(counts[j] > counts[val])
			val = j;
	}
	return (double)(val - 1);
}

void GData::baselineVector(int nOutputCount, double* pOutputs)
{
	int i;
	int nInputCount = m_pRelation->size() - nOutputCount;
	for(i = 0; i < nOutputCount; i++)
		pOutputs[i] = baselineValue(nInputCount + i);
}

bool GData::isAttrHomogenous(int col)
{
	if(m_pRelation->valueCount(col) > 0)
	{
		int d;
		vector<double*>::iterator it = m_rows.begin();
		for( ; it != m_rows.end(); it++)
		{
			d = (int)(*it)[col];
			if(d != UNKNOWN_DISCRETE_VALUE)
			{
				it++;
				break;
			}
		}
		for( ; it != m_rows.end(); it++)
		{
			int t = (int)(*it)[col];
			if(t != d && t != UNKNOWN_DISCRETE_VALUE)
				return false;
		}
	}
	else
	{
		double d;
		vector<double*>::iterator it = m_rows.begin();
		for( ; it != m_rows.end(); it++)
		{
			d = (*it)[col];
			if(d != UNKNOWN_REAL_VALUE)
			{
				it++;
				break;
			}
		}
		for( ; it != m_rows.end(); it++)
		{
			double t = (*it)[col];
			if(t != d && t != UNKNOWN_REAL_VALUE)
				return false;
		}
	}
	return true;
}

bool GData::areLabelsHomogenous(int labelDims)
{
	int featureDims = m_pRelation->size() - labelDims;
	for(int i = 0; i < labelDims; i++)
	{
		if(!isAttrHomogenous(featureDims + i))
			return false;
	}
	return true;
}

void GData::randomlyReplaceMissingValues(int nAttr, GRand* pRand)
{
	// Just return if there are no missing attributes
	double dOk = m_pRelation->valueCount(nAttr) == 0 ? -1e200 : 0;
	double* pRow;
	size_t n, i, nStart;
	size_t nCount = rows();
	for(i = 0; i < nCount; i++)
	{
		pRow = row(i);
		if(pRow[nAttr] < dOk)
			break;
	}
	if(i >= nCount)
		return;

	// Sort (which will put all missing attributes at the front)
	sort(nAttr);

	// Find where the missing values end
	for(nStart = 0; nStart < nCount; nStart++)
	{
		pRow = row(nStart);
		if(pRow[nAttr] >= dOk)
			break;
	}

	// Replace with randomly selected good values
	if(nStart >= nCount)
	{
		for(i = 0; i < nStart; i++)
			row(i)[nAttr] = 0;
	}
	else
	{
		for(i = 0; i < nStart; i++)
		{
			n = (size_t)pRand->next(nCount - nStart) + nStart;
			row(i)[nAttr] = row(n)[nAttr];
		}
	}
}

void GData::principalComponent(double* pOutVector, int dims, const double* pMean, GRand* pRand)
{
	// Initialize the out-vector to a random direction
	pRand->spherical(pOutVector, dims);

	// Iterate
	double* pVector;
	size_t nCount = rows();
	GTEMPBUF(double, pAccumulator, dims);
	double d;
	double mag = 0;
	for(int iters = 0; iters < 200; iters++)
	{
		GVec::setAll(pAccumulator, 0.0, dims);
		for(size_t n = 0; n < nCount; n++)
		{
			pVector = row(n);
			d = GVec::dotProduct(pMean, pVector, pOutVector, dims);
			double* pAcc = pAccumulator;
			const double* pM = pMean;
			for(int j = 0; j < dims; j++)
				*(pAcc++) += d * (*(pVector++) - *(pM++));
		}
		GVec::copy(pOutVector, pAccumulator, dims);
		GVec::safeNormalize(pOutVector, dims, pRand);
		d = GVec::squaredMagnitude(pAccumulator, dims);
		if(iters < 6 || d > mag)
			mag = d;
		else
			break;
	}
}

void GData::principalComponentAboutOrigin(double* pOutVector, int dims, GRand* pRand)
{
	// Initialize the out-vector to a random direction
	pRand->spherical(pOutVector, dims);

	// Iterate
	double* pVector;
	size_t nCount = rows();
	GTEMPBUF(double, pAccumulator, dims);
	double d;
	double mag = 0;
	for(int iters = 0; iters < 200; iters++)
	{
		GVec::setAll(pAccumulator, 0.0, dims);
		for(size_t n = 0; n < nCount; n++)
		{
			pVector = row(n);
			d = GVec::dotProduct(pVector, pOutVector, dims);
			double* pAcc = pAccumulator;
			for(int j = 0; j < dims; j++)
				*(pAcc++) += d * *(pVector++);
		}
		GVec::copy(pOutVector, pAccumulator, dims);
		GVec::safeNormalize(pOutVector, dims, pRand);
		d = GVec::squaredMagnitude(pAccumulator, dims);
		if(iters < 6 || d > mag)
			mag = d;
		else
			break;
	}
}

void GData::principalComponentIgnoreUnknowns(double* pOutVector, int dims, const double* pMean, GRand* pRand)
{
	// Initialize the out-vector to a random direction
	pRand->spherical(pOutVector, dims);

	// Iterate
	double* pVector;
	size_t nCount = rows();
	GTEMPBUF(double, pAccumulator, dims);
	double d;
	double mag = 0;
	for(int iters = 0; iters < 200; iters++)
	{
		GVec::setAll(pAccumulator, 0.0, dims);
		for(size_t n = 0; n < nCount; n++)
		{
			pVector = row(n);
			d = GVec::dotProductIgnoringUnknowns(pMean, pVector, pOutVector, dims);
			double* pAcc = pAccumulator;
			const double* pM = pMean;
			for(int j = 0; j < dims; j++)
			{
				if(*pVector != UNKNOWN_REAL_VALUE)
					(*pAcc) += d * (*pVector - *pM);
				pVector++;
				pAcc++;
				pM++;
			}
		}
		GVec::copy(pOutVector, pAccumulator, dims);
		GVec::safeNormalize(pOutVector, dims, pRand);
		d = GVec::squaredMagnitude(pAccumulator, dims);
		if(iters < 6 || d > mag)
			mag = d;
		else
			break;
	}
}

void GData::weightedPrincipalComponent(double* pOutVector, int dims, const double* pMean, const double* pWeights, GRand* pRand)
{
	// Initialize the out-vector to a random direction
	pRand->spherical(pOutVector, dims);

	// Iterate
	double* pVector;
	size_t nCount = rows();
	GTEMPBUF(double, pAccumulator, dims);
	double d;
	double mag = 0;
	for(int iters = 0; iters < 200; iters++)
	{
		GVec::setAll(pAccumulator, 0.0, dims);
		const double* pW = pWeights;
		for(size_t n = 0; n < nCount; n++)
		{
			pVector = row(n);
			d = GVec::dotProduct(pMean, pVector, pOutVector, dims);
			double* pAcc = pAccumulator;
			const double* pM = pMean;
			for(int j = 0; j < dims; j++)
				*(pAcc++) += (*pW) * d * (*(pVector++) - *(pM++));
			pW++;
		}
		GVec::copy(pOutVector, pAccumulator, dims);
		GVec::safeNormalize(pOutVector, dims, pRand);
		d = GVec::squaredMagnitude(pAccumulator, dims);
		if(iters < 6 || d > mag)
			mag = d;
		else
			break;
	}
}

double GData::eigenValue(const double* pMean, const double* pEigenVector, int dims, GRand* pRand)
{
	// Use the element of the eigenvector with the largest magnitude,
	// because that will give us the least rounding error when we compute the eigenvalue.
	int index = GVec::indexOfMaxMagnitude(pEigenVector, dims, pRand);

	// The eigenvalue is the factor by which the eigenvector is scaled by the covariance matrix,
	// so we compute just the part of the covariance matrix that we need to see how much the
	// max-magnitude element of the eigenvector is scaled.
	double d = 0;
	for(int i = 0; i < dims; i++)
		d += covariance(index, pMean[index], i, pMean[i]) * pEigenVector[i];
	return d / pEigenVector[index];
}

void GData::removeComponent(const double* pMean, const double* pComponent, int dims)
{
	size_t nCount = rows();
	for(size_t i = 0; i < nCount; i++)
	{
		double* pVector = row(i);
		double d = GVec::dotProductIgnoringUnknowns(pMean, pVector, pComponent, dims);
		for(int j = 0; j < dims; j++)
		{
			if(*pVector != UNKNOWN_REAL_VALUE)
				(*pVector) -= d * pComponent[j];
			pVector++;
		}
	}
}

void GData::removeComponentAboutOrigin(const double* pComponent, int dims)
{
	size_t nCount = rows();
	for(size_t i = 0; i < nCount; i++)
	{
		double* pVector = row(i);
		double d = GVec::dotProduct(pVector, pComponent, dims);
		for(int j = 0; j < dims; j++)
		{
			(*pVector) -= d * pComponent[j];
			pVector++;
		}
	}
}

void GData::centerMeanAtOrigin()
{
	int dims = cols();
	GTEMPBUF(double, mean, dims);
	centroid(mean);
	for(vector<double*>::iterator it = m_rows.begin(); it != m_rows.end(); it++)
		GVec::subtract(*it, mean, dims);
}

int GData::countPrincipalComponents(double d, GRand* pRand)
{
	int dims = cols();
	GData tmpData(relation(), heap());
	tmpData.copy(this);
	tmpData.centerMeanAtOrigin();
	GTEMPBUF(double, vec, dims);
	double thresh = d * d * tmpData.sumSquaredDistance(NULL);
	int i;
	for(i = 1; i < dims; i++)
	{
		tmpData.principalComponentAboutOrigin(vec, dims, pRand);
		tmpData.removeComponentAboutOrigin(vec, dims);
		if(tmpData.sumSquaredDistance(NULL) < thresh)
			break;
	}
	return i;
}

double GData::sumSquaredDistance(const double* pPoint)
{
	int dims = relation()->size();
	double err = 0;
	if(pPoint)
	{
		for(size_t i = 0; i < rows(); i++)
			err += GVec::squaredDistance(pPoint, row(i), dims);
	}
	else
	{
		for(size_t i = 0; i < rows(); i++)
			err += GVec::squaredMagnitude(row(i), dims);
	}
	return err;
}

double GData::sumSquaredDifference(GData& that, bool transpose)
{
	if(transpose)
	{
		size_t colCount = (size_t)cols();
		if(rows() != (size_t)that.cols() || colCount != that.rows())
			ThrowError("expected matrices of same size");
		double err = 0;
		for(size_t i = 0; i < rows(); i++)
		{
			double* pRow = row(i);
			for(size_t j = 0; j < colCount; j++)
			{
				double d = *(pRow++) - that[j][i];
				err += (d * d);
			}
		}
		return err;
	}
	else
	{
		if(this->rows() != that.rows() || this->cols() != that.cols())
			ThrowError("different sizes");
		int colCount = cols();
		double d = 0;
		for(size_t i = 0; i < rows(); i++)
			d += GVec::squaredDistance(this->row(i), that[i], colCount);
		return d;
	}
}

double GData::linearCorrelationCoefficient(int attr1, double attr1Origin, int attr2, double attr2Origin)
{
	double sx = 0;
	double sy = 0;
	double sxy = 0;
	double mx, my;
	double* pPat;
	size_t count = rows();
	size_t i;
	for(i = 0; i < count; i++)
	{
		pPat = row(i);
		mx = pPat[attr1] - attr1Origin;
		my = pPat[attr2] - attr2Origin;
		if(pPat[attr1] == UNKNOWN_REAL_VALUE || pPat[attr2] == UNKNOWN_REAL_VALUE)
			continue;
		break;
	}
	if(i >= count)
		return 0;
	double d, x, y;
	int j = 1;
	for(i++; i < count; i++)
	{
		pPat = row(i);
		if(pPat[attr1] == UNKNOWN_REAL_VALUE || pPat[attr2] == UNKNOWN_REAL_VALUE)
			continue;
		x = pPat[attr1] - attr1Origin;
		y = pPat[attr2] - attr2Origin;
		d = (double)j / (j + 1);
		sx += (x - mx) * (x - mx) * d;
		sy += (y - my) * (y - my) * d;
		sxy += (x - mx) * (y -  my) * d;
		mx += (x - mx) / (j + 1);
		my += (y - my) / (j + 1);
	}
	if(sx == 0 || sy == 0 || sxy == 0)
		return 0;
	return (sxy / j) / (sqrt(sx / j) * sqrt(sy / j));
}

double GData::covariance(int nAttr1, double dMean1, int nAttr2, double dMean2)
{
	size_t nRowCount = rows();
	double* pVector;
	double dSum = 0;
	for(size_t i = 0; i < nRowCount; i++)
	{
		pVector = row(i);
		dSum += ((pVector[nAttr1] - dMean1) * (pVector[nAttr2] - dMean2));
	}
	return dSum / (nRowCount - 1);
}

GData* GData::covarianceMatrix()
{
	int colCount = cols();
	GData* pOut = new GData(colCount);
	pOut->newRows(colCount);

	// Compute the deviations
	GTEMPBUF(double, pMeans, colCount);
	int n, i;
	for(i = 0; i < colCount; i++)
		pMeans[i] = mean(i);

	// Compute the covariances for half the matrix
	for(i = 0; i < colCount; i++)
	{
		double* pRow = pOut->row(i);
		pRow += i;
		for(n = i; n < colCount; n++)
			*(pRow++) = covariance(i, pMeans[i], n, pMeans[n]);
	}

	// Fill out the other half of the matrix
	for(i = 1; i < colCount; i++)
	{
		double* pRow = pOut->row(i);
		for(n = 0; n < i; n++)
			*(pRow++) = pOut->row(n)[i];
	}
	return pOut;
}

class Row_Binary_Predicate_Functor
{
protected:
	int m_dim;

public:
	Row_Binary_Predicate_Functor(int dim) : m_dim(dim)
	{
	}

	bool operator() (double* pA, double* pB) const
	{
		return pA[m_dim] < pB[m_dim];
	}
};

void GData::sort(int dim)
{
	Row_Binary_Predicate_Functor comparer(dim);
	std::sort(m_rows.begin(), m_rows.end(), comparer);
}

void GData::reverseRows()
{
	std::reverse(m_rows.begin(), m_rows.end());
}

double GData_PairedTTestHelper(void* pThis, double x)
{
	double v = *(double*)pThis;
	return pow(1.0 + x * x / v, -(v + 1) / 2);
}

void GData::pairedTTest(int* pOutV, double* pOutT, int attr1, int attr2, bool normalize)
{
	double* pPat;
	double a, b, m;
	double asum = 0;
	double asumOfSquares = 0;
	double bsum = 0;
	double bsumOfSquares = 0;
	size_t rowCount = rows();
	for(size_t i = 0; i < rowCount; i++)
	{
		pPat = row(i);
		a = pPat[attr1];
		b = pPat[attr2];
		if(normalize)
		{
			m = (a + b) / 2;
			a /= m;
			b /= m;
		}
		asum += a;
		asumOfSquares += (a * a);
		bsum += b;
		bsumOfSquares += (b * b);
	}
	double amean = asum / rowCount;
	double avariance = (asumOfSquares / rowCount - amean * amean) * rowCount / (rowCount - 1);
	double bmean = bsum / rowCount;
	double bvariance = (bsumOfSquares / rowCount - bmean * bmean) * rowCount / (rowCount - 1);
	double grand = sqrt((avariance + bvariance) / rowCount);
	*pOutV = 2 * (int)rowCount - 2;
	*pOutT = ABS(bmean - amean) / grand;
}

double GData::wilcoxonSignedRanksTest(int attr1, int attr2, double tolerance)
{
	// Make sorted list of differences
	GHeap heap(1024);
	size_t rowCount = rows();
	GData tmp(2, &heap);
	tmp.reserve(rowCount);
	double* pPat;
	double* pPat2;
	for(size_t i = 0; i < rowCount; i++)
	{
		pPat = row(i);
		pPat2 = tmp.newRow();
		pPat2[0] = ABS(pPat[attr2] - pPat[attr1]);
		if(pPat2[0] < tolerance)
			pPat2[1] = 0;
		else if(pPat[attr1] < pPat[attr2])
			pPat2[1] = -1;
		else
			pPat2[1] = 1;
	}
	tmp.sort(0);

	// Convert to ranks
	double prev = UNKNOWN_REAL_VALUE;
	size_t index = 0;
	size_t j;
	double ave;
	for(size_t i = 0; i < rowCount; i++)
	{
		pPat = tmp[i];
		if(ABS(pPat[0] - prev) >= tolerance)
		{
			ave = (double)(index + 1 + i) / 2;
			for(j = index; j < i; j++)
			{
				pPat2 = tmp[j];
				pPat2[0] = ave;
			}
			prev = pPat[0];
			index = i;
		}
	}
	ave = (double)(index + 1 + rowCount) / 2;
	for(j = index; j < rowCount; j++)
	{
		pPat2 = tmp[j];
		pPat2[0] = ave;
	}

	// Sum up the scores
	double a = 0;
	double b = 0;
	for(size_t i = 0; i < rowCount; i++)
	{
		pPat = tmp[i];
		if(pPat[1] > 0)
			a += pPat[0];
		else if(pPat[1] < 0)
			b += pPat[0];
		else
		{
			a += pPat[0] / 2;
			b += pPat[0] / 2;
		}
	}
	return MIN(a, b);
}

size_t GData::countValue(int attribute, double value)
{
	size_t count = 0;
	for(size_t i = 0; i < rows(); i++)
	{
		if(row(i)[attribute] == value)
			count++;
	}
	return count;
}

void GData::ensureDataHasNoMissingReals()
{
	int dims = m_pRelation->size();
	for(size_t i = 0; i < rows(); i++)
	{
		double* pPat = row(i);
		for(int j = 0; j < dims; j++)
		{
			if(m_pRelation->valueCount(j) != 0)
				continue;
			if(pPat[i] == UNKNOWN_REAL_VALUE)
				ThrowError("Missing values in continuous attributes are not supported");
		}
	}
}

void GData::ensureDataHasNoMissingNominals()
{
	int dims = m_pRelation->size();
	for(size_t i = 0; i < rows(); i++)
	{
		double* pPat = row(i);
		for(int j = 0; j < dims; j++)
		{
			if(m_pRelation->valueCount(j) == 0)
				continue;
			if((int)pPat[i] == (int)UNKNOWN_DISCRETE_VALUE)
				ThrowError("Missing values in nominal attributes are not supported");
		}
	}
}

void GData::print(ostream& stream)
{
	m_pRelation->print(stream, this, 14);
}

double GData::measureLabelInfo(int labelDims)
{
	int nInputCount = cols() - labelDims;
	double dInfo = 0;
	int n, nIndex;
	for(n = 0; n < labelDims; n++)
	{
		nIndex = nInputCount + n;
		if(m_pRelation->valueCount(nIndex) == 0)
		{
			double m = mean(nIndex);
			dInfo += variance(nIndex, m);
		}
		else
			dInfo += entropy(nIndex);
	}
	return dInfo;
}

double GData::sumSquaredDiffWithIdentity()
{
	size_t m = MIN(rows(), (size_t)cols());
	double err = 0;
	double d;
	for(size_t i = 0; i < m; i++)
	{
		double* pRow = row(i);
		for(size_t j = 0; j < m; j++)
		{
			d = *(pRow++);
			if(i == j)
				d -= 1;
			err += (d * d);
		}
	}
	return err;
}
/*
bool GData::leastCorrelatedVector(double* pOut, GData* pThat, GRand* pRand)
{
	if(rows() != pThat->rows() || cols() != pThat->cols())
		ThrowError("Expected matrices with the same dimensions");
	GData* pC = GData::multiply(*pThat, *this, false, true);
	Holder<GData> hC(pC);
	GData* pE = GData::multiply(*pC, *pC, true, false);
	Holder<GData> hE(pE);
	double d = pE->sumSquaredDiffWithIdentity();
	if(d < 0.001)
		return false;
	GData* pF = pE->mostSignificantEigenVectors(rows(), pRand);
	Holder<GData> hF(pF);
	GVec::copy(pOut, pF->row(rows() - 1), rows());
	return true;
}
*/
bool GData::leastCorrelatedVector(double* pOut, GData* pThat, GRand* pRand)
{
	if(rows() != pThat->rows() || cols() != pThat->cols())
		ThrowError("Expected matrices with the same dimensions");
	GData* pC = GData::multiply(*pThat, *this, false, true);
	Holder<GData> hC(pC);
	GData* pD = GData::multiply(*pThat, *pC, true, false);
	Holder<GData> hD(pD);
	double d = pD->sumSquaredDifference(*this, true);
	if(d < 1e-9)
		return false;
	pD->subtract(this, true);
	pD->principalComponentAboutOrigin(pOut, rows(), pRand);
	return true;
/*
	GData* pE = GData::multiply(*pD, *pD, true, false);
	Holder<GData> hE(pE);
	GData* pF = pE->mostSignificantEigenVectors(1, pRand);
	Holder<GData> hF(pF);
	GVec::copy(pOut, pF->row(0), rows());
	return true;
*/
}

double GData::dihedralCorrelation(GData* pThat, GRand* pRand)
{
	int colCount = cols();
	if(rows() == 1)
		return ABS(GVec::dotProduct(row(0), pThat->row(0), colCount));
	GTEMPBUF(double, pBuf, rows() + 2 * colCount);
	double* pA = pBuf + rows();
	double* pB = pA + colCount;
	if(!leastCorrelatedVector(pBuf, pThat, pRand))
		return 1.0;
	multiply(pBuf, pA, true);
	if(!pThat->leastCorrelatedVector(pBuf, this, pRand))
		return 1.0;
	pThat->multiply(pBuf, pB, true);
	return ABS(GVec::correlation(pA, pB, colCount));
}

void GData::project(double* pDest, const double* pPoint)
{
	int dims = cols();
	GVec::setAll(pDest, 0.0, dims);
	for(size_t i = 0; i < rows(); i++)
	{
		double* pBasis = row(i);
		GVec::addScaled(pDest, GVec::dotProduct(pPoint, pBasis, dims), pBasis, dims);
	}
}

void GData::project(double* pDest, const double* pPoint, const double* pOrigin)
{
	int dims = cols();
	GVec::copy(pDest, pOrigin, dims);
	for(size_t i = 0; i < rows(); i++)
	{
		double* pBasis = row(i);
		GVec::addScaled(pDest, GVec::dotProduct(pOrigin, pPoint, pBasis, dims), pBasis, dims);
	}
}

#ifndef NO_TEST_CODE
void GData_testMultiply()
{
	GData a(2);
	a.newRows(2);
	GData b(2);
	b.newRows(2);
	a[0][0] = 2; a[0][1] = 17;
	a[1][0] = 7; a[1][1] = 19;
	b[0][0] = 11; b[0][1] = 3;
	b[1][0] = 5; b[1][1] = 13;
	GData* pC;
	pC = GData::multiply(a, b, false, false);
	if(pC->rows() != 2 || pC->cols() != 2)
		ThrowError("wrong size");
	if(pC->row(0)[0] != 107 || pC->row(0)[1] != 227 ||
		pC->row(1)[0] != 172 || pC->row(1)[1] != 268)
		ThrowError("wrong answer");
	delete(pC);
	GData* pA = a.transpose();
	pC = GData::multiply(*pA, b, true, false);
	if(pC->rows() != 2 || pC->cols() != 2)
		ThrowError("wrong size");
	if(pC->row(0)[0] != 107 || pC->row(0)[1] != 227 ||
		pC->row(1)[0] != 172 || pC->row(1)[1] != 268)
		ThrowError("wrong answer");
	delete(pC);
	GData* pB = b.transpose();
	pC = GData::multiply(a, *pB, false, true);
	if(pC->rows() != 2 || pC->cols() != 2)
		ThrowError("wrong size");
	if(pC->row(0)[0] != 107 || pC->row(0)[1] != 227 ||
		pC->row(1)[0] != 172 || pC->row(1)[1] != 268)
		ThrowError("wrong answer");
	delete(pC);
	pC = GData::multiply(*pA, *pB, true, true);
	if(pC->rows() != 2 || pC->cols() != 2)
		ThrowError("wrong size");
	if(pC->row(0)[0] != 107 || pC->row(0)[1] != 227 ||
		pC->row(1)[0] != 172 || pC->row(1)[1] != 268)
		ThrowError("wrong answer");
	delete(pC);
	delete(pA);
	delete(pB);
}

void GData_testCholesky()
{
	GData m1(3);
	m1.newRows(3);
	m1[0][0] = 3;	m1[0][1] = 0;	m1[0][2] = 0;
	m1[1][0] = 1;	m1[1][1] = 4;	m1[1][2] = 0;
	m1[2][0] = 2;	m1[2][1] = 2;	m1[2][2] = 7;
	GData* pM3 = GData::multiply(m1, m1, false, true);
	Holder<GData> hM3(pM3);
	GData* pM4 = pM3->cholesky();
	Holder<GData> hM4(pM4);
	if(m1.sumSquaredDifference(*pM4, false) >= .0001)
		ThrowError("Cholesky decomposition didn't work right");
}

void GData_testInvert()
{
	GData i1(3);
	i1.newRows(3);
	i1[0][0] = 2;	i1[0][1] = -1;	i1[0][2] = 0;
	i1[1][0] = -1;	i1[1][1] = 2;	i1[1][2] = -1;
	i1[2][0] = 0;	i1[2][1] = -1;	i1[2][2] = 2;
//	i1.invert();
	GData* pInv = i1.pseudoInverse();
	Holder<GData> hInv(pInv);
	GData i2(3);
	i2.newRows(3);
	i2[0][0] = .75;	i2[0][1] = .5;	i2[0][2] = .25;
	i2[1][0] = .5;	i2[1][1] = 1;	i2[1][2] = .5;
	i2[2][0] = .25;	i2[2][1] = .5;	i2[2][2] = .75;
	if(pInv->sumSquaredDifference(i2, false) >= .0001)
		ThrowError("Not good enough");
//	i1.invert();
	GData* pInvInv = pInv->pseudoInverse();
	Holder<GData> hInvInv(pInvInv);
	GData* pI3 = GData::multiply(*pInvInv, i2, false, false);
	Holder<GData> hI3(pI3);
	GData i4(3);
	i4.newRows(3);
	i4.makeIdentity();
	if(pI3->sumSquaredDifference(i4, false) >= .0001)
		ThrowError("Not good enough");
}

void GData_testDeterminant()
{
	const double dettest[] =
	{
		1,2,3,4,
		5,6,7,8,
		2,6,4,8,
		3,1,1,2,
	};
	GData d1(4);
	d1.fromVector(dettest, 4);
	double det = d1.determinant();
	if(ABS(det - 72.0) >= .0001)
		ThrowError("wrong");
	const double dettest2[] =
	{
		3,2,
		5,7,
	};
	GData d2(2);
	d2.fromVector(dettest2, 2);
	det = d2.determinant();
	if(ABS(det - 11.0) >= .0001)
		ThrowError("wrong");
	const double dettest3[] =
	{
		1,2,3,
		4,5,6,
		7,8,9,
	};
	GData d3(3);
	d3.fromVector(dettest3, 3);
	det = d3.determinant();
	if(ABS(det - 0.0) >= .0001)
		ThrowError("wrong");
}

void GData_testReducedRowEchelonForm()
{
	const double reducedrowechelonformtest[] =
	{
		1,-1,1,0,2,
		2,-2,0,2,2,
		-1,1,2,-3,1,
		-2,2,1,-3,-1,
	};
	const double reducedrowechelonformanswer[] =
	{
		1,-1,0,1,1,
		0,0,1,-1,1,
		0,0,0,0,0,
		0,0,0,0,0,
	};
	GData r1(5);
	r1.fromVector(reducedrowechelonformtest, 4);
	if(r1.toReducedRowEchelonForm() != 2)
		ThrowError("wrong answer");
	GData r2(5);
	r2.fromVector(reducedrowechelonformanswer, 4);
	if(r1.sumSquaredDifference(r2) > .001)
		ThrowError("wrong answer");
	const double reducedrowechelonformtest2[] =
	{
		-2, -4, 4,
		2, -8, 0,
		8, 4, -12,
	};
	const double reducedrowechelonformanswer2[] =
	{
		1, 0, -4.0/3,
		0, 1, -1.0/3,
		0, 0, 0,
	};
	GData r3(3);
	r3.fromVector(reducedrowechelonformtest2, 3);
	if(r3.toReducedRowEchelonForm() != 2)
		ThrowError("wrong answer");
	GData r4(3);
 	r4.fromVector(reducedrowechelonformanswer2, 3);
	if(r4.sumSquaredDifference(r3) > .001)
		ThrowError("wrong answer");
}

void GData_testPrincipalComponents(GRand& prng)
{
	// Test principal components
	GHeap heap(1000);
	GData data(2, &heap);
	int i;
	for(i = 0; i < 100; i++)
	{
		double* pNewRow = data.newRow();
		pNewRow[0] = prng.uniform();
		pNewRow[1] = 2 * pNewRow[0];
	}
	double mean[2];
	mean[0] = data.mean(0);
	mean[1] = data.mean(1);
	double eig[2];
	data.principalComponent(eig, 2, mean, &prng);
	if(ABS(eig[0] * 2 - eig[1]) > .0001)
		ThrowError("incorrect value");

	// Compute principal components via eigenvectors of covariance matrix, and
	// make sure they're the same
	GMixedRelation rel;
	rel.addAttr(0);
	rel.addAttr(0);
	GData* pM = data.covarianceMatrix();
	Holder<GData> hM(pM);
	double ev;
	GData* pEigenVecs = pM->eigs(1, &ev, &prng, true);
	Holder<GData> hEigenVecs(pEigenVecs);
	if(ABS(pEigenVecs->row(0)[0] * pEigenVecs->row(0)[1] - eig[0] * eig[1]) > .0001)
		ThrowError("answers don't agree");

	// Test most significant eigenvector computation
	GData e1(2);
	e1.newRows(2);
	e1[0][0] = 1;	e1[0][1] = 1;
	e1[1][0] = 1;	e1[1][1] = 4;
	double ev2[2];
	GData* pE2 = e1.eigs(2, ev2, &prng, true);
	Holder<GData> hE2(pE2);
	if(ABS(pE2->row(0)[0] * pE2->row(0)[0] + pE2->row(0)[1] * pE2->row(0)[1] - 1) > .0001)
		ThrowError("answer not normalized");
	if(ABS(pE2->row(0)[0] * pE2->row(0)[1] - .27735) >= .0001)
		ThrowError("wrong answer");
	if(ABS(pE2->row(1)[0] * pE2->row(1)[0] + pE2->row(1)[1] * pE2->row(1)[1] - 1) > .0001)
		ThrowError("answer not normalized");
	if(ABS(pE2->row(1)[0] * pE2->row(1)[1] + .27735) >= .0001)
		ThrowError("wrong answer");

	// Test least significant eigenvector computation and gaussian ellimination
	GData e3(2);
	e3.newRows(2);
	e3[0][0] = 9;	e3[0][1] = 3;
	e3[1][0] = 3;	e3[1][1] = 5;
	GData* pE4 = e3.eigs(2, ev2, &prng, true);
	Holder<GData> hE4(pE4);
	GData* pE5 = e3.eigs(2, ev2, &prng, false);
	Holder<GData> hE5(pE5);
	if(ABS(ABS(pE4->row(0)[0]) - ABS(pE5->row(1)[0])) >= .0001)
		ThrowError("failed");
	if(ABS(ABS(pE4->row(0)[1]) - ABS(pE5->row(1)[1])) >= .0001)
		ThrowError("failed");
	if(ABS(ABS(pE4->row(1)[0]) - ABS(pE5->row(0)[0])) >= .0001)
		ThrowError("failed");
	if(ABS(ABS(pE4->row(1)[1]) - ABS(pE5->row(0)[1])) >= .0001)
		ThrowError("failed");
}

void GData_testDihedralCorrelation(GRand& prng)
{
	// Test dihedral angle computation
	for(int iter = 0; iter < 500; iter++)
	{
		// Make a random set of orthonormal basis vectors
		int dims = 5;
		GData basis(dims);
		basis.newRows(dims);
		for(int i = 0; i < dims; i++)
		{
			prng.spherical(basis[i], dims);
			for(int j = 0; j < i; j++)
			{
				GVec::subtractComponent(basis[i], basis[j], dims);
				GVec::normalize(basis[i], dims);
			}
		}
	
		// Make two planes with a known dihedral angle
		double angle = prng.uniform() * 0.5 * M_PI;
		double angle2 = prng.uniform() * 2.0 * M_PI;
		GData p1(dims);
		p1.newRows(2);
		GVec::setAll(p1[0], 0.0, dims);
		p1[0][0] = cos(angle2);
		p1[0][2] = sin(angle2);
		GVec::setAll(p1[1], 0.0, dims);
		p1[1][0] = -sin(angle2);
		p1[1][2] = cos(angle2);
		GData p2(dims);
		p2.newRows(2);
		GVec::setAll(p2[0], 0.0, dims);
		p2[0][0] = cos(angle);
		p2[0][1] = sin(angle);
		GVec::setAll(p2[1], 0.0, dims);
		p2[1][2] = 1.0;

		// Transform the planes with the basis matrix
		GData p3(dims);
		p3.newRows(2);
		basis.multiply(p1[0], p3[0]);
		basis.multiply(p1[1], p3[1]);
		GData p4(dims);
		p4.newRows(2);
		basis.multiply(p2[0], p4[0]);
		basis.multiply(p2[1], p4[1]);

		// Measure the angle
		double actual = cos(angle);
		double measured = p3.dihedralCorrelation(&p4, &prng);
		if(ABS(measured - actual) > 1e-8)
			ThrowError("failed");
	}

	// Measure the dihedral angle of two 3-hyperplanes in 5-space
	double angle = 0.54321;
	GData bas(5);
	bas.newRows(5);
	bas.makeIdentity();
	bas[2][2] = cos(angle);
	bas[2][4] = sin(angle);
	bas[4][2] = -sin(angle);
	bas[4][4] = cos(angle);
	GData sp1(5);
	sp1.newRows(3);
	sp1.makeIdentity();
	GData* sp3 = GData::multiply(sp1, bas, false, true);
	Holder<GData> hSp3(sp3);
	double cosangle = sp1.dihedralCorrelation(sp3, &prng);
	double measured = acos(cosangle);
	if(ABS(measured - angle) > 1e-8)
		ThrowError("failed");

	// Make sure dihedral angles are computed correctly with parallel planes
	static const double aa[] = {1.0, 0.0, 0.0, 0.0, -1.0, 0.0};
	static const double bb[] = {0.6, 0.8, 0.0, -0.8, 0.6, 0.0};
	GData planeA(3);
	planeA.fromVector(aa, 2);
	GData planeB(3);
	planeB.fromVector(bb, 2);
	cosangle = planeA.dihedralCorrelation(&planeB, &prng);
	if(ABS(cosangle - 1.0) > 1e-8)
		ThrowError("failed");
	cosangle = planeB.dihedralCorrelation(&planeA, &prng);
	if(ABS(cosangle - 1.0) > 1e-8)
		ThrowError("failed");
}

void GData_testSingularValueDecomposition()
{
	GData* pU;
	double* pDiag;
	GData* pV;
	GData M(2);
	M.newRows(2);
	M[0][0] = 4.0; M[0][1] = 3.0;
	M[1][0] = 0.0; M[1][1] = -5.0;
	M.singularValueDecomposition(&pU, &pDiag, &pV);
	Holder<GData> hU(pU);
	ArrayHolder<double> hDiag(pDiag);
	Holder<GData> hV(pV);

	// Test that the diagonal values are correct
	if(ABS(pDiag[0] - sqrt(40.0)) > 1e-8)
		ThrowError("pDiag is not correct");
	if(ABS(pDiag[1] - sqrt(10.0)) > 1e-8)
		ThrowError("pDiag is not correct");

	// Test that U is unitary
	GData* pT1 = GData::multiply(*pU, *pU, false, true);
	Holder<GData> hT1(pT1);
	if(pT1->sumSquaredDiffWithIdentity() > 1e-8)
		ThrowError("U is not unitary");

	// Test that V is unitary
	GData* pT2 = GData::multiply(*pV, *pV, false, true);
	Holder<GData> hT2(pT2);
	if(pT2->sumSquaredDiffWithIdentity() > 1e-8)
		ThrowError("V is not unitary");
}

void GData_testPseudoInverse()
{
	{
		GData M(2);
		M.newRows(2);
		M[0][0] = 1.0; M[0][1] = 1.0;
		M[1][0] = 2.0; M[1][1] = 2.0;
		GData* A = M.pseudoInverse();
		Holder<GData> hA(A);
		GData B(2);
		B.newRows(2);
		B[0][0] = 0.1; B[0][1] = 0.2;
		B[1][0] = 0.1; B[1][1] = 0.2;
		if(A->sumSquaredDifference(B, false) > 1e-8)
			ThrowError("failed");
	}
	{
		GData M(2);
		M.newRows(3);
		M[0][0] = 1.0; M[0][1] = 2.0;
		M[1][0] = 3.0; M[1][1] = 4.0;
		M[2][0] = 5.0; M[2][1] = 6.0;
		GData* A = M.pseudoInverse();
		Holder<GData> hA(A);
		if(A->rows() != 2 || A->cols() != 3)
			ThrowError("wrong size");
		GData B(3);
		B.newRows(2);
		B[0][0] = -16.0/12.0; B[0][1] = -4.0/12.0; B[0][2] = 8.0/12.0;
		B[1][0] = 13.0/12.0; B[1][1] = 4.0/12.0; B[1][2] = -5.0/12.0;
		if(A->sumSquaredDifference(B, false) > 1e-8)
			ThrowError("failed");
	}
	{
		GData M(3);
		M.newRows(2);
		M[0][0] = 1.0; M[0][1] = 3.0; M[0][2] = 5.0;
		M[1][0] = 2.0; M[1][1] = 4.0; M[1][2] = 6.0;
		GData* A = M.pseudoInverse();
		Holder<GData> hA(A);
		if(A->rows() != 3 || A->cols() != 2)
			ThrowError("wrong size");
		GData B(2);
		B.newRows(3);
		B[0][0] = -16.0/12.0; B[0][1] = 13.0/12.0;
		B[1][0] = -4.0/12.0; B[1][1] = 4.0/12.0;
		B[2][0] = 8.0/12.0; B[2][1] = -5.0/12.0;
		if(A->sumSquaredDifference(B, false) > 1e-8)
			ThrowError("failed");
	}
}

void GData_testKabsch()
{
	GRand prng(0);
	GData a(5);
	a.newRows(20);
	for(size_t i = 0; i < 20; i++)
	{
		prng.spherical(a[i], 5);
		GVec::multiply(a[i], prng.uniform() + 0.5, 5);
	}
	GData rot(5);
	static const double rr[] = {
		0.0, 1.0, 0.0, 0.0, 0.0,
		0.0, 0.0, 0.0, 1.0, 0.0,
		0.0, 0.0, 0.8, 0.0, -0.6,
		1.0, 0.0, 0.0, 0.0, 0.0,
		0.0, 0.0, 0.6, 0.0, 0.8
	};
	rot.newRows(5);
	rot.fromVector(rr, 5);
	GData* pB = GData::multiply(a, rot, false, false);
	Holder<GData> hB(pB);
	GData* pK = GData::kabsch(&a, pB);
	Holder<GData> hK(pK);
	if(pK->sumSquaredDifference(rot, true) > 1e-6)
		ThrowError("Failed to recover rotation matrix");
	GData* pC = GData::multiply(*pB, *pK, false, false);
	Holder<GData> hC(pC);
	if(a.sumSquaredDifference(*pC, false) > 1e-6)
		ThrowError("Failed to align data");
}

// static
void GData::test()
{
	GRand prng(0);
	GData_testMultiply();
	GData_testCholesky();
	GData_testInvert();
	GData_testDeterminant();
	GData_testReducedRowEchelonForm();
	GData_testPrincipalComponents(prng);
	GData_testDihedralCorrelation(prng);
	GData_testSingularValueDecomposition();
	GData_testPseudoInverse();
	GData_testKabsch();
}
#endif // !NO_TEST_CODE




GDataArray::GDataArray(sp_relation& pRelation)
: m_pRelation(pRelation)
{
}

GDataArray::GDataArray(int cols)
{
	m_pRelation = new GUniformRelation(cols, 0);
}

GDataArray::~GDataArray()
{
	flush();
}

void GDataArray::flush()
{
	for(vector<GData*>::iterator it = m_sets.begin(); it != m_sets.end(); it++)
		delete(*it);
	m_sets.clear();
}

GData* GDataArray::newSet(int rows)
{
	GData* pData = new GData(m_pRelation);
	m_sets.push_back(pData);
	pData->newRows(rows);
	return pData;
}

void GDataArray::newSets(int count, int rows)
{
	m_sets.reserve(m_sets.size() + count);
	for(int i = 0; i < count; i++)
		newSet(rows);
}

size_t GDataArray::largestSet()
{
	vector<GData*>::iterator it = m_sets.begin();
	size_t biggestRows = (*it)->rows();
	size_t biggestIndex = 0;
	size_t i = 1;
	for(it++; it != m_sets.end(); it++)
	{
		if((*it)->rows() > biggestRows)
		{
			biggestRows = (*it)->rows();
			biggestIndex = i;
		}
		i++;
	}
	return biggestIndex;
}
