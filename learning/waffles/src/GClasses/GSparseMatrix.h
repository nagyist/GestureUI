/*
	Copyright (C) 2006, Mike Gashler

	This library is free software; you can redistribute it and/or
	modify it under the terms of the GNU Lesser General Public
	License as published by the Free Software Foundation; either
	version 2.1 of the License, or (at your option) any later version.

	see http://www.gnu.org/copyleft/lesser.html
*/

#ifndef __GSPARSEMATRIX_H__
#define __GSPARSEMATRIX_H__

#include <map>
#include <iostream>

namespace GClasses {

class GData;
class GRand;

/// This is a helper class used internally by GSparseMatrix.
class GSparseMatrixElement
{
public:
	unsigned int row;
	unsigned int col;
	double val;
};



/// This class simultaneously stores a row-compressed sparse matrix and
/// a column-compressed sparse matrix. It updates both of them
class GSparseMatrix
{
protected:
	typedef std::map<unsigned int,double> Map;
	typedef Map::iterator It;
	unsigned int m_rows, m_cols;
	Map* m_pRows;
	Map* m_pCols;

public:
	GSparseMatrix(unsigned int rows, unsigned int cols);
	~GSparseMatrix();

#ifndef NO_TEST_CODE
	static void test();
#endif

	typedef Map::const_iterator Iter;

	/// Returns the number of rows (as if this matrix were dense)
	unsigned int rows() { return m_rows; }

	/// Returns the number of columns (as if this matrix were dense)
	unsigned int cols() { return m_cols; }

	/// Parses from the specified string
	static GSparseMatrix* parse(const char* pFile, size_t len);

	/// Loads and parses from the specified file
	static GSparseMatrix* load(const char* filename);

	/// Prints in a text format to the specifed stream
	void print(std::ostream& stream);

	/// Saves to the specified file
	void save(const char* filename);

	/// Copies a row into a non-sparse vector
	void fullRow(double* pOutFullRow, unsigned int row);

	/// Returns a const_iterator to the beginning of a row. The iterator
	/// references a pair, such that first is the column, and second is the value.
	Iter rowBegin(unsigned int i) { return m_pRows[i].begin(); }

	/// Returns a const_iterator to the end of a row. The iterator
	/// references a pair, such that first is the column, and second is the value.
	Iter rowEnd(unsigned int i) { return m_pRows[i].end(); }

	/// Returns a const_iterator to the beginning of a column. The iterator
	/// references a pair, such that first is the row, and second is the value.
	Iter colBegin(unsigned int i) { return m_pRows[i].begin(); }

	/// Returns a const_iterator to the end of a column. The iterator
	/// references a pair, such that first is the row, and second is the value.
	Iter colEnd(unsigned int i) { return m_pRows[i].end(); }

	/// Returns the value at the specified position in the matrix. Returns 0
	/// if no element is stored at that position.
	double get(unsigned int row, unsigned int col);

	/// Sets a value at the specified position in the matrix. (Updates
	/// both copies of the matrix. If val is 0, it removes the element
	/// from both copies of the matrix.)
	void set(unsigned int row, unsigned int col, double val);

	/// Copies values from "that" into "this". If the matrices are different
	/// sizes, any non-overlapping elements will be left at zero.
	void copyFrom(GSparseMatrix* that);

	/// Copies values from "that" into "this". If the matrices are different
	/// sizes, any non-overlapping elements will be left at zero.
	void copyFrom(GData* that);

	/// Converts to a full matrix
	GData* toFullMatrix();

	/// Multiplies the matrix by a scalar value
	void multiply(double scalar);

	/// Swaps the two specified columns
	void swapColumns(unsigned int a, unsigned int b);

	/// Swaps the two specified rows
	void swapRows(unsigned int a, unsigned int b);

	/// Shuffles the rows in this matrix
	void shuffle(GRand* pRand);

	/// Returns a sub-matrix of this matrix
	GSparseMatrix* subMatrix(int row, int col, int height, int width);

	/// Transposes this matrix in place. (This operation is very efficient. It only needs to swap two pointers and two values.)
	void transpose();

	/// This method is currently buggy. Don't use it until I remove this comment.
	void singularValueDecomposition(GSparseMatrix** ppU, double** ppDiag, GSparseMatrix** ppV, int maxIters = 30);

protected:
	void singularValueDecompositionHelper(GSparseMatrix** ppU, double** ppDiag, GSparseMatrix** ppV, int maxIters);

};

} // namespace GClasses

#endif // __GSPARSEMATRIX_H__
