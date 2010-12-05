/*
	Copyright (C) 2006, Mike Gashler

	This library is free software; you can redistribute it and/or
	modify it under the terms of the GNU Lesser General Public
	License as published by the Free Software Foundation; either
	version 2.1 of the License, or (at your option) any later version.

	see http://www.gnu.org/copyleft/lesser.html
*/

#ifndef __GDATA_H__
#define __GDATA_H__

#include "GMacros.h"
#include <vector>
#include <string>
#include <algorithm>
#include <iostream>
#include "GHolders.h"

namespace GClasses {

#define UNKNOWN_REAL_VALUE -1e308

// Why do we need a different value for unknown discrete values? Because it's
// common to store discrete values in an integer. Integers can't store -1e308,
// and we can't use -1 for unknown reals b/c it's not an insane value.
#define UNKNOWN_DISCRETE_VALUE -1


class GData;
class GPrediction;
class GRand;
class GHeap;
class GTwtDoc;
class GTwtNode;

/// Holds the metadata for a dataset, including which attributes
/// are continuous or nominal, and how many values each nominal
/// attribute supports.
class GRelation
{
public:
	enum RelationType
	{
		UNIFORM,
		MIXED,
		ARFF,
	};

	GRelation() {}
	virtual ~GRelation() {}

	/// Returns the type of relation
	virtual RelationType type() = 0;

	/// Saves this relation to a text-based format
	virtual GTwtNode* toTwt(GTwtDoc* pDoc) = 0;

	/// Returns the number of attributes (columns)
	virtual int size() = 0;

	/// Returns the number of values in the specified attribute. (Returns 0 for
	/// continuous attributes.)
	virtual int valueCount(int nAttr) = 0;

	/// Returns true of all of the attributes in the specified range are continuous
	virtual bool areContinuous(int first, int count) = 0;

	/// Returns true of all of the attributes in the specified range are nominal
	virtual bool areNominal(int first, int count) = 0;

	/// Makes a copy of this relation
	virtual GRelation* clone() = 0;

	/// Deletes the specified attribute
	virtual void deleteAttribute(int index) = 0;

	/// Swaps two attributes
	virtual void swapAttributes(int nAttr1, int nAttr2) = 0;

	/// Prints as an ARFF file to the specified stream. (pData can be NULL if data is not available)
	virtual void print(std::ostream& stream, GData* pData, int precision);

	/// Appends a textual representation of the specified value of the specified attribut to pOutString
	void attrValue(std::string* pOutString, int attr, double value);

	/// Prints the specified value
	virtual void printValue(std::ostream& stream, double value, int column);

	/// Print a single row in ARFF format
	void printRow(std::ostream& stream, double* pRow, const char* separator);

	/// Counts the size of the corresponding real-space vector
	int countRealSpaceDims(int nFirstAttr, int nAttrCount);

	/// Converts a row (pIn) to a real-space vector (pOut)
	/// (pIn should point to the nFirstAttr'th element, not the first element)
	void toRealSpace(const double* pIn, double* pOut, int nFirstAttr, int nAttrCount);

	/// Converts a real-space vector (pIn) to a row (pOut)
	/// nFirstAttr and nAttrCount refer to the row indexes
	void fromRealSpace(const double* pIn, double* pOut, int nFirstAttr, int nAttrCount, GRand* pRand);

	/// Converts a real-space vector (pIn) to an array of predictions (pOut)
	/// nFirstAttr and nAttrCount refer to the prediction indexes
	void fromRealSpace(const double* pIn, GPrediction* pOut, int nFirstAttr, int nAttrCount);

	/// Loads from a text-based format
	static smart_ptr<GRelation> fromTwt(GTwtNode* pNode);

	/// Saves to a file
	void save(GData* pData, const char* szFilename, int precision);
};

typedef smart_ptr<GRelation> sp_relation;


/// A relation with a minimal memory footprint that assumes
/// all attributes are continuous, or all of them are nominal
/// and have the same number of possible values.
class GUniformRelation : public GRelation
{
protected:
	int m_attrCount;
	int m_valueCount;

public:
	GUniformRelation(int attrCount, int valueCount = 0)
	: m_attrCount(attrCount), m_valueCount(valueCount)
	{
	}

	GUniformRelation(GTwtNode* pNode);

	virtual RelationType type() { return UNIFORM; }
	
	/// Serializes this object
	virtual GTwtNode* toTwt(GTwtDoc* pDoc);
	
	/// Returns the number of attributes (columns)
	virtual int size() { return m_attrCount; }
	
	/// Returns the number of values in each nominal attribute (or 0 if the attributes are continuous)
	virtual int valueCount(int nAttr) { return m_valueCount; }
	
	/// See the comment for GRelation::areContinuous
	virtual bool areContinuous(int first, int count) { return m_valueCount == 0; }
	
	/// See the comment for GRelation::areNominal
	virtual bool areNominal(int first, int count) { return m_valueCount != 0; }

	/// Returns a copy of this object
	virtual GRelation* clone() { return new GUniformRelation(m_attrCount, m_valueCount); }
	
	/// Drop the specified attribute
	virtual void deleteAttribute(int index);
	
	/// Swap two attributes
	virtual void swapAttributes(int nAttr1, int nAttr2) {}
};



class GMixedRelation : public GRelation
{
protected:
	std::vector<int> m_valueCounts;

public:
	/// Makes an empty relation
	GMixedRelation();

	/// Loads from a text file in ".twt" format
	GMixedRelation(GTwtNode* pNode);

	/// Makes a copy of pCopyMe
	GMixedRelation(GRelation* pCopyMe);

	/// Makes a copy of the specified range of pCopyMe
	GMixedRelation(GRelation* pCopyMe, int firstAttr, int attrCount);

	virtual ~GMixedRelation();

	virtual RelationType type() { return MIXED; }

	/// Writes to a text file in ".twt" format
	virtual GTwtNode* toTwt(GTwtDoc* pDoc);

	virtual GRelation* clone();

	/// Deletes all the attributes
	virtual void flush();

	/// If nValues is zero, adds a real attribute. If nValues is > 0, adds
	/// an attribute with "nValues" nominal values
	void addAttr(int nValues);

	/// Adds "attrCount" new attributes, each with "valueCount" values. (Use valueCount=0
	/// for continuous attributes.)
	void addAttrs(int attrCount, int valueCount);

	/// Copies the specified attributes and adds them to this relation.
	/// If attrCount < 0, then it will copy all attributes from firstAttr to the end.
	void addAttrs(GRelation* pCopyMe, int firstAttr = 0, int attrCount = -1);

	/// Flushes this relation and then copies all of the attributes from pCopyMe
	void copy(GRelation* pCopyMe);

	/// Adds a copy of the specified attribute to this relation
	virtual void copyAttr(GRelation* pThat, int nAttr);

	/// Returns the total number of attributes in this relation
	virtual int size()
	{
		return (int)m_valueCounts.size();
	}

	/// Returns the number of nominal values in the specified attribute
	virtual int valueCount(int nAttr)
	{
		return m_valueCounts[nAttr];
	}

	/// Sets the number of values for this attribute
	virtual void setAttrValueCount(int nAttr, int nValues);

	virtual bool areContinuous(int first, int count);

	virtual bool areNominal(int first, int count);

	/// Swaps two columns
	virtual void swapAttributes(int nAttr1, int nAttr2);

	/// Deletes an attribute.
	virtual void deleteAttribute(int nAttr);
};


class GArffAttribute
{
public:
	std::string m_name;
	std::vector<std::string> m_values;
};


/// ARFF = Attribute-Relation File Format. This stores richer information
/// than GRelation. This includes a name, a name for each attribute, and
/// names for each supported nominal value.
class GArffRelation : public GMixedRelation
{
friend class GData;
protected:
	std::string m_name;
	std::vector<GArffAttribute> m_attrs;

public:
	GArffRelation();
	virtual ~GArffRelation();

	virtual RelationType type() { return ARFF; }

	/// Returns a copy of this object
	virtual GRelation* clone();

	/// Deletes all the attributes
	virtual void flush();

	/// Prints as an ARFF file to the specified stream. (pData can be NULL if data is not available)
	virtual void print(std::ostream& stream, GData* pData, int precision);

	/// Prints the specified value
	virtual void printValue(std::ostream& sream, double value, int column);

	/// Adds a new attribute (column) to the relation
	void addAttribute(const char* szName, int nValues, std::vector<const char*>* pValues);

	/// Adds a copy of the specified attribute to this relation
	virtual void copyAttr(GRelation* pThat, int nAttr);

	/// Returns the name of the relation
	const char* name() { return m_name.c_str(); }

	/// Sets the name of this relation
	void setName(const char* szName);

	/// Returns the name of the specified attribute
	const char* attrName(int nAttr);

	/// Adds a new possible value to a nominal attribute
	void addAttrValue(int nAttr, const char* szValue);

	/// Sets the number of values for the specified attribute
	virtual void setAttrValueCount(int nAttr, int nValues);

	/// Swaps two columns
	virtual void swapAttributes(int nAttr1, int nAttr2);

	/// Deletes an attribute
	virtual void deleteAttribute(int nAttr);

	/// Returns the nominal index for the specified attribute with the given value
	int findEnumeratedValue(int nAttr, const char* szValue);

	/// Parses a value
	double parseValue(int attr, const char* val);

protected:
	/// takes ownership of ppValues
	void addAttributeInternal(const char* pName, int nameLen, int valueCount);
	void addAttributeInternal(const char* pName, int nameLen, const char* pValues, int valuesLen);

	void parseAttribute(const char* szFile, size_t nLen, int nLine);
};


/// Represents a matrix or a database table. Elements can be discrete or continuous.
/// References a GRelation object, which stores the meta-information about each column.
class GData
{
protected:
	sp_relation m_pRelation;
	GHeap* m_pHeap;
	std::vector<double*> m_rows;

public:
	/// pRelation is a smart-pointer to a relation, which specifies
	/// the type of each attribute.
	GData(sp_relation& pRelation, GHeap* pHeap = NULL);

	/// attrs specifies the number of attributes. All attributes
	/// will be real-valued.
	GData(int attrs, GHeap* pHeap = NULL);

	/// Loads a dataset from a text-based format
	GData(GTwtNode* pNode, GHeap* pHeap = NULL);

	~GData();

	/// Matrix add. Adds the values in pThat to this. (If transpose
	/// is true, adds the transpose of pThat to this.) Both datasets
	/// must have the same dimensions. Behavior is undefined for nominal columns.
	void add(GData* pThat, bool transpose);

	/// Returns a new dataset that contains a subset of the attributes in this dataset
	GData* attrSubset(int firstAttr, int attrCount);

	/// This computes the square root of this matrix. (If you take the matrix that this
	/// returns and multiply it by its transpose, you should get the original dataset
	/// again.) Behavior is undefined if there are nominal attributes. If this matrix
	/// is not positive definate, it will throw an exception.
	GData* cholesky();

	/// Makes a copy of this dataset
	GData* clone();

	/// Copies the specified column into pOutVector
	void col(int index, double* pOutVector);

	/// Returns the number of columns in the dataset
	int cols() const { return m_pRelation->size(); }

	/// Copies all the data from pThat. (Just references the same relation)
	void copy(GData* pThat);

	/// Copies the specified block of columns from pSource to this dataset. pSource must have
	/// the same number of rows as this dataset.
	void copyColumns(int nDestStartColumn, GData* pSource, int nSourceStartColumn, int nColumnCount);

	/// Adds a copy of the row to the data set
	void copyRow(const double* pRow);

	/// Computes the determinant of this matrix
	double determinant();

	/// Computes the eigenvalue that corresponds to the specified eigenvector
	/// of this matrix
	double eigenValue(const double* pEigenVector);

	/// Computes the eigenvector that corresponds to the specified eigenvalue of
	/// this matrix. Note that this method trashes this matrix, so
	/// make a copy first if you care.
	void eigenVector(double eigenvalue, double* pOutVector);

	/// Computes y in the equation M*y=x (or y=M^(-1)x), where M is this dataset, which
	/// must be a square matrix, and x is pVector as passed in, and y is pVector after
	/// the call. If there are multiple solutions, it finds the one for which all
	/// the variables in the null-space have a value of 1. If there are no solutions,
	/// it returns false. Note that this method trashes this dataset (so make a copy
	/// first if you care).
	bool gaussianElimination(double* pVector);

	/// Returns the heap used to allocate rows for this dataset
	GHeap* heap() { return m_pHeap; }
/*
	/// Inverts the matrix using an LU-decomposition method. Overwrites this matrix.
	void invert();
*/
	/// This computes K=kabsch(A,B), such that K is an n-by-n matrix, where n is pA->cols().
	/// K is the optimal orthonormal rotation matrix to align A and B, such that A(K^T)
	/// minimizes sum-squared error with B, and BK minimizes sum-squared error with A.
	/// (This rotates around the origin, so typically you will want to subtract the
	/// centroid from both pA and pB before calling this.)
	static GData* kabsch(GData* pA, GData* pB);

	/// This uses the Kabsch algorithm to rotate and translate pB in order to minimize
	/// RMS with pA. (pA and pB must have the same number of rows and columns.)
	static GData* align(GData* pA, GData* pB);

	/// Loads an ARFF file and returns the data. This will throw an exception if
	/// there's an error.
	static GData* loadArff(const char* szFilename);

	/// Loads a file in CSV format.
	static GData* loadCsv(const char* szFilename, char separator, bool columnNamesInFirstRow, bool tolerant);

	/// Sets this dataset to an identity matrix. (It doesn't change the
	/// number of columns or rows. It just stomps over existing values.)
	void makeIdentity();

	/// If upperToLower is true, copies the upper triangle of this matrix over the lower triangle
	/// If upperToLower is false, copies the lower triangle of this matrix over the upper triangle
	void mirrorTriangle(bool upperToLower);

	/// Merges two datasets side-by-side. The resulting dataset
	/// will contain the attributes of both datasets. Both pSetA
	/// and pSetB (and the resulting dataset) must have the same
	/// number of rows
	static GData* mergeHoriz(GData* pSetA, GData* pSetB);

	/// Steals all the rows from pData and adds them to this set.
	/// (You still have to delete pData.) Both datasets must have
	/// the same number of columns.
	void mergeVert(GData* pData);

	/// Computes nCount eigenvectors and the corresponding eigenvalues using the power method.
	/// (This method is only accurate if a small number of eigenvalues/vectors are needed.)
	/// If mostSignificant is true, the largest eigenvalues are found. If mostSignificant
	/// is false, the smallest eigenvalues are found.
	GData* eigs(int nCount, double* pEigenVals, GRand* pRand, bool mostSignificant);

	/// Multiplies every element in the dataset by scalar.
	/// Behavior is undefined for nominal columns.
	void multiply(double scalar);

	/// Multiplies this matrix by the column vector pVectorIn to get
	/// pVectorOut. (If transpose is true, then it multiplies the transpose
	/// of this matrix by pVectorIn to get pVectorOut.) pVectorIn should have
	/// the same number of elements as columns (or rows if transpose is true)
	/// and pVectorOut should have the same number of elements as rows (or
	/// cols, if transpose is true.) Note that if transpose is true, it is the
	/// same as if pVectorIn is a row vector and you multiply it by this matrix
	/// to get pVectorOut.
	void multiply(const double* pVectorIn, double* pVectorOut, bool transpose = false);

	/// Matrix multiply. For convenience, you can also specify that neither, one, or both
	/// of the inputs are virtually transposed prior to the multiplication. (If you want
	/// the results to come out transposed, you can use the equality AB=((B^T)(A^T))^T to
	/// figure out how to specify the parameters.)
	static GData* multiply(GData& a, GData& b, bool transposeA, bool transposeB);

	///
	/// Adds a new row to the dataset. (The values in the row are not initialized)
	///
	double* newRow();

	/// Adds "nRows" uninitialized rows of size "nAttributes" to the data set
	void newRows(size_t nRows);

	/// Parses an ARFF file and returns the data. This will throw an exception if
	/// there's an error.
	static GData* parseArff(const char* szFile, size_t nLen);

	/// Imports data from a text file. Determines the meta-data automatically.
	/// Note: This method does not support Mac line-endings. You should first replace all '\r' with '\n' if your data comes from a Mac.
	/// As a special case, if separator is '\0', then it assumes data elements are separated by any number of whitespace characters, that
	/// element values themselves contain no whitespace, and that there are no missing elements. (This is the case when you save a
	/// Matlab matrix to an ascii file.)
	static GData* parseCsv(const char* pFile, size_t len, char separator, bool columnNamesInFirstRow, bool tolerant = false);

	/// Computes the Moore-Penrose pseudoinverse of this matrix (using the SVD method). You
	/// are responsible to delete the matrix this returns.
	GData* pseudoInverse();

	/// Returns a relation object, which holds meta-data about the attributes (columns)
	sp_relation& relation() { return m_pRelation; }

	/// Allocates space for the specified number of patters (to avoid superfluous resizing)
	void reserve(size_t n) { m_rows.reserve(n); }

	/// Returns the number of rows in the dataset
	size_t rows() const { return m_rows.size(); }

	/// Saves the dataset to a file in ARFF format
	void saveArff(const char* szFilename);

	/// Sets the relation for this dataset
	void setRelation(sp_relation& pRelation) { m_pRelation = pRelation; }

	/// Performs SVD on A, where A is this m-by-n matrix.
	/// *ppU will be set to an m-by-m matrix where the columns are the eigenvectors of A(A^T).
	/// *ppDiag will be set to an array of n doubles holding the square roots of the corresponding eigenvalues.
	/// *ppV will be set to an n-by-n matrix where the rows are the eigenvectors of (A^T)A.
	/// You are responsible to delete(*ppU), delete(*ppV), and delete[] *ppDiag.
	void singularValueDecomposition(GData** ppU, double** ppDiag, GData** ppV, bool throwIfNoConverge = false, int maxIters = 80);

	/// Matrix subtract. Subtracts the values in pThat from this. (If transpose
	/// is true, subtracts the transpose of pThat from this.) Both datasets
	/// must have the same dimensions. Behavior is undefined for nominal columns.
	void subtract(GData* pThat, bool transpose);

	/// Returns the sum squared difference between this matrix and an identity matrix
	double sumSquaredDiffWithIdentity();

	/// Adds an already-allocated row to this dataset. The row must
	/// be allocated in the same heap that this dataset uses. (There is no way
	/// for this method to verify that, so be careful.)
	void takeRow(double* pRow);

	/// Converts the matrix to reduced row echelon form
	int toReducedRowEchelonForm();

	/// Copies all the data from this dataset into pVector. pVector must be
	/// big enough to hold rows() x cols() doubles.
	void toVector(double* pVector);

	/// Writes data to a text file in ".twt" format
	GTwtNode* toTwt(GTwtDoc* pDoc);

	/// Returns the sum of the diagonal elements
	double trace();

	/// Returns a dataset that is this dataset transposed. (All columns
	/// in the returned dataset will be continuous.)
	GData* transpose();

	/// Copies the data from pVector over this dataset. nRows specifies the number
	/// of rows of data in pVector.
	void fromVector(const double* pVector, size_t nRows);

	/// Returns a pointer to the specified row
	inline double* row(size_t index) { return m_rows[index]; }

	/// Returns a pointer to the specified row
	double* operator [](size_t index) { return m_rows[index]; }

	/// Sets all elements in this dataset to the specified value
	void setAll(double val);

	/// Copies pVector over the specified column
	void setCol(int index, const double* pVector);

	/// Swaps pRow with the row at nIndex. You're responsible to delete the
	/// row this returns
	double* replaceRow(size_t nIndex, double* pRow);

	/// Swaps the two specified rows
	void swapRows(size_t a, size_t b);

	/// Swaps two columns
	void swapColumns(int nAttr1, int nAttr2);

	/// Deletes a column
	void deleteColumn(int index);

	/// Swaps the specified row with the last row, and then releases it from the dataset.
	/// If this dataset does not have its own heap, then you must delete the row this returns
	double* releaseRow(size_t index);

	/// Swaps the specified row with the last row, and then deletes it.
	void deleteRow(size_t index);

	/// Releases the specified row from the dataset and shifts everything after it up one slot.
	/// If this dataset does not have its own heap, then you must delete the row this returns
	double* releaseRowPreserveOrder(size_t index);

	/// Deletes the specified row and shifts everything after it up one slot
	void deleteRowPreserveOrder(size_t index);

	/// Replaces any occurrences of NAN in the matrix with the corresponding values from
	/// an identity matrix.
	void fixNans();

	/// Deletes all the data
	void flush();

	/// Abandons (leaks) all the rows of data
	void releaseAllRows();

	/// Randomizes the order of the rows
	void shuffle(GRand* pRand);

	/// Shuffles the order of the rows. Also shuffles the rows in "other" in
	/// the same way, such that corresponding rows are preserved.
	void shuffle2(GRand* pRand, GData& other);

	/// This is an inferior way to shuffle the data
	void shuffleLikeCards();

	/// Sorts the data from smallest to largest in the specified dimension
	void sort(int nDimension);

	/// Reverses the row order
	void reverseRows();

	/// Sorts rows according to the specified compare function. (Return true
	/// to indicate thate the first row comes before the second row.)
	template<typename CompareFunc>
	void sort(CompareFunc& pComparator)
	{
		std::sort(m_rows.begin(), m_rows.end(), pComparator);
	}

	/// Splits this set of data into two sets. Values greater-than-or-equal-to
	/// dPivot stay in this data set. Values less than dPivot go into pLessThanPivot
	void splitByPivot(GData* pGreaterOrEqual, int nAttribute, double dPivot);

	/// Moves all rows with the specified value in the specified attribute
	/// into pSingleClass
	void splitByDiscreteValue(GData* pSingleClass, int nAttr, int nValue);

	/// Removes the last nOtherRows rows from this data set and
	/// puts them in pOtherData
	void splitBySize(GData* pOtherData, size_t nOtherRows);

	/// Measures the entropy of the specified attribute
	double entropy(int nColumn);

	/// Finds the min and the range of the values of the specified attribute
	void minAndRange(int nAttribute, double* pMin, double* pRange);

	/// Estimates the actual min and range based on a random sample
	void minAndRangeUnbiased(int nAttribute, double* pMin, double* pRange);

	/// Shifts the data such that the mean occurs at the origin
	void centerMeanAtOrigin();

	/// Computes the arithmetic mean of the values in the specified column
	double mean(int nAttribute);

	/// Computes the median of the values in the specified column
	double median(int nAttribute);

	/// Computes the arithmetic means of all attributes
	void centroid(double* pOutCentroid);

	/// Computes the average variance of a single attribute
	double variance(int nAttr, double mean);

	/// Normalizes the specified attribute values
	void normalize(int nAttribute, double dInputMin, double dInputRange, double dOutputMin, double dOutputRange);

	/// Normalize a value from the input min and range to the output min and range
	static double normalize(double dVal, double dInputMin, double dInputRange, double dOutputMin, double dOutputRange);

	/// Determines the mean (for continuous attributes) or the most common value (for nominal attributes)
	double baselineValue(int nAttribute);

	/// Produce an output row with the most common output (for nominal attributes) and the mean
	/// (for continuous attributes)
	void baselineVector(int nOutputCount, double* pOutputs);

	/// Returns true iff the specified attribute contains homogenous values. (Unknowns are counted as homogenous with anything)
	bool isAttrHomogenous(int col);

	/// Returns true iff each of the last labelDims columns in the data are homogenous
	bool areLabelsHomogenous(int labelDims);

	/// Replaces all missing values by copying a randomly selected non-missing value in the same attribute
	/// (This messes up the order of the rows because it uses sorting to find the missing values, so
	///  you should probably call Shuffle() when you're done calling this method.)
	void randomlyReplaceMissingValues(int nAttr, GRand* pRand);

	/// This is an efficient algorithm for iteratively computing the principal component vector
	/// (the eigenvector of the covariance matrix) of the data. See "EM Algorithms for PCA and SPCA"
	/// by Sam Roweis, 1998 NIPS.
	/// nIterations should be a small constant. 20 seems work well for most applications.
	/// (To compute the next principal component, call RemoveComponent, then call this again.)
	void principalComponent(double* pOutVector, int dims, const double* pMean, GRand* pRand);

	/// Computes the first principal component assuming the mean is already subtracted out of the data
	void principalComponentAboutOrigin(double* pOutVector, int dims, GRand* pRand);

	/// Computes principal components, while ignoring missing values
	void principalComponentIgnoreUnknowns(double* pOutVector, int dims, const double* pMean, GRand* pRand);

	/// Computes the first principal component of the data with each row weighted according to the
	/// vector pWeights. (pWeights must have an element for each row.)
	void weightedPrincipalComponent(double* pOutVector, int dims, const double* pMean, const double* pWeights, GRand* pRand);

	/// After you compute the principal component, you can call this to obtain the eigenvalue that
	/// corresponds to that principal component vector (eigenvector).
	double eigenValue(const double* pMean, const double* pEigenVector, int dims, GRand* pRand);

	/// Removes the component specified by pComponent from the data. (pComponent should already be normalized.)
	/// This might be useful, for example, to remove the first principal component from the data so you can
	/// then proceed to compute the second principal component, and so forth.
	void removeComponent(const double* pMean, const double* pComponent, int dims);

	/// Removes the specified component assuming the mean is zero.
	void removeComponentAboutOrigin(const double* pComponent, int dims);

	/// Computes the minimum number of principal components necessary so that
	/// less than the specified portion of the deviation in the data
	/// is unaccounted for. (For example, if the data projected onto
	/// the first 3 principal components contains 90 percent of the
	/// deviation that the original data contains, then if you pass
	/// the value 0.1 to this method, it will return 3.)
	int countPrincipalComponents(double d, GRand* pRand);

	/// Computes the sum-squared distance between pPoint and all of the points in the dataset.
	/// (If pPoint is NULL, it computes the sum-squared distance with the origin.)
	/// (Note that this is equal to the sum of all the eigenvalues times the number of dimensions,
	/// so you can efficiently compute eigenvalues as the difference in sumSquaredDistance with
	/// the mean after removing the corresponding component, and then dividing by the number of
	/// dimensions. This is more efficient than calling eigenValue.)
	double sumSquaredDistance(const double* pPoint);

	/// Computes the squared distance between this and that. (If transpose is true, computes
	/// the difference between this and the transpose of that.)
	double sumSquaredDifference(GData& that, bool transpose = false);

	/// Computes the linear coefficient between the two specified attributes.
	/// Usually you will want to pass the mean values for attr1Origin and attr2Origin.
	double linearCorrelationCoefficient(int attr1, double attr1Origin, int attr2, double attr2Origin);

	/// Computes the covariance between two attributes
	double covariance(int nAttr1, double dMean1, int nAttr2, double dMean2);

	/// Computes the covariance matrix of the data
	GData* covarianceMatrix();

	/// Performs a paired T-Test with data from the two specified attributes.
	/// pOutV will hold the degrees of freedom. pOutT will hold the T-value.
	/// You can use GMath::tTestAlphaValue to convert these to a P-value.
	void pairedTTest(int* pOutV, double* pOutT, int attr1, int attr2, bool normalize);

	/// Performs the Wilcoxon signed ranks test from the two specified attributes
	/// and returns the T-value. If two values are closer than tolerance, they are
	/// considered to be equal.
	double wilcoxonSignedRanksTest(int attr1, int attr2, double tolerance);

	/// Prints the data to the specified stream
	void print(std::ostream& stream);

	/// Returns the number of ocurrences of the specified value in the specified attribute
	size_t countValue(int attribute, double value);

	/// Throws an exception if this data contains any missing values in a continuous attribute
	void ensureDataHasNoMissingReals();

	/// Throws an exception if this data contains any missing values in a nominal attribute
	void ensureDataHasNoMissingNominals();

	/// Computes the entropy of the labels (or the deviation if continuous)
	double measureLabelInfo(int labelDims);

	/// Computes the vector in this subspace that has the greatest distance from its
	/// projection into pThat subspace. Returns true if the results are computed. Returns
	/// false if the subspaces are so nearly parallel that pOut cannot be computed with
	/// accuracy.
	bool leastCorrelatedVector(double* pOut, GData* pThat, GRand* pRand);

	/// Computes the cosine of the dihedral angle between this subspace and pThat subspace
	double dihedralCorrelation(GData* pThat, GRand* pRand);

	/// Projects pPoint onto this hyperplane (where each row defines
	/// one of the orthonormal basis vectors of this hyperplane)
	/// This computes (A^T)Ap, where A is this matrix, and p is pPoint.
	void project(double* pDest, const double* pPoint);

	/// Projects pPoint onto this hyperplane (where each row defines
	/// one of the orthonormal basis vectors of this hyperplane)
	void project(double* pDest, const double* pPoint, const double* pOrigin);

#ifndef NO_TEST_CODE
	/// Performs unit tests for this class. Throws an exception if there is a failure.
	static void test();
#endif // !NO_TEST_CODE
protected:
	static void parseDataRow(GArffRelation* pRelation, GData* pData, const char* szFile, size_t nLen, int nLine);
	double determinantHelper(int nEndRow, int* pColumnList);
	void inPlaceSquareTranspose();
	void singularValueDecompositionHelper(GData** ppU, double** ppDiag, GData** ppV, bool throwIfNoConverge, int maxIters);
};


/// This is a special holder that guarantees the data set
/// will release all of its data before it is deleted
class GReleaseDataHolder
{
protected:
	GData* m_pData;

public:
	GReleaseDataHolder(GData* pData)
	{
		m_pData = pData;
	}

	~GReleaseDataHolder()
	{
		m_pData->releaseAllRows();
	}
};


class GDataArray
{
protected:
	sp_relation m_pRelation;
	std::vector<GData*> m_sets;

public:
	GDataArray(sp_relation& pRelation);
	GDataArray(int cols);
	~GDataArray();
	std::vector<GData*>& sets() { return m_sets; }

	/// Adds a new dataset to the array and preallocates
	/// the specified number of rows
	GData* newSet(int rows);

	/// Adds count new datasets to the array, and
	/// preallocates the specified number of rows in
	/// each one
	void newSets(int count, int rows);

	/// Deletes all the datasets
	void flush();

	/// Returns the index of the largest data set
	size_t largestSet();
};

} // namespace GClasses

#endif // __GDATA_H__
