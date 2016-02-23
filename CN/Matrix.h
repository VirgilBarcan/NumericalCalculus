#pragma once

#include <string>

class Matrix
{
public:
	/*
	Set the number of lines that the matrix will have
	@param noOfLines - the number of lines
	*/
	void setNoOfLines(int noOfLines) { this->noOfLines = noOfLines; }

	/*
	Get the number of lines that the matrix has
	@return - the number of lines
	*/
	int getNoOfLines() { return this->noOfLines; }

	/*
	Set the number of columns that the matrix will have
	@param noOfColumns - the number of columns
	*/
	void setNoOfColumns(int noOfColumns) { this->noOfColumns = noOfColumns; }

	/*
	Get the number of columns that the matrix has
	@return - the number of columns
	*/
	int getNoOfColumns() { return this->noOfColumns; }

	/*
	Add a new element to the matrix
	This function will be implemented in the derivated classes
	@param line - the line where the element should be placed
	@param column - the column where the element should be placed
	@param value - the value of the element
	*/
	virtual void addElementAt(int line, int column, double value) = 0;

	/*
	Get the value of an element in the matrix
	This function will be implemented in the derivated classes
	@param line - the line where the element should be placed
	@param column - the column where the element should be placed
	@return - the value of the element at the wanted line and column
	*/
	virtual double getElementAt(int line, int column) = 0;

	/*
	Get the transpose of the matrix
	The transpose is the matrix where the element A[i, j] is placed in A[j, i]:
		/ a11 a12 a13 \				 / a11 a21 a31 \
	A = | a21 a22 a23 |		=>	AT = | a12 a22 a32 |
		\ a31 a32 a33 /				 \ a13 a23 a33 /
	@return - the transpose of the matrix
	*/
	virtual Matrix *transpose() = 0;

	/*
	Get the string version of the matrix
	This function will be implemented in the derivated classes
	Needed to print the matrix
	@return - the string version of the matrix
	*/
	virtual std::string toString() = 0;
protected:
	int noOfLines;
	int noOfColumns;
};

