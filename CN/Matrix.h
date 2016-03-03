#pragma once

#include <string>
#include <random>

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
	Get the matrix from the file given as parameter
	@param filePath - the path to the file that contains the matrix
	*/
	virtual void getFromFile(std::string filePath) = 0;

	/*
	Generate the values of the matrix
	The values are generated in the interval [min, max]
	@param min - the minimum value generated
	@param max - the maximum value generated
	*/
	virtual void generateRandomMatrixValues(double min, double max) = 0;

	/*
	Check if the wanted element line and column are within the bounds of the matrix
	@param line - the line to check the size for
	@param column - the column to check the size for
	@return - true if the the line and column are in the bounds of the matrix, false otherwise
	*/
	bool checkBounds(int line, int column) {
		if ((line >= 0 && line < noOfLines) &&
			(column >= 0 && column < noOfColumns))
			return true;
		return false;
	}

	/*
	Check if the given matrices have equal dimensions
	@param matrix1 - the first matrix
	@param matrix2 - the second matrix
	@return - true if the two matrices have equal size, false otherwise
	*/
	bool checkEqualSizes(Matrix *matrix1, Matrix *matrix2) {
		if ((matrix1->getNoOfLines() == matrix2->getNoOfLines()) &&
			(matrix1->getNoOfColumns() == matrix2->getNoOfColumns()))
			return true;
		return false;
	}

	/*
	Check if the given matrices can be multiplied
	In order to multiply two matrices the first one's number of columns has to be equal to the second one's number of lines
	@param matrix1 - the first matrix
	@param matrix2 - the second matrix
	@return - true if the two matrices are multipliable, false otherwise
	*/
	bool checkMultipliable(Matrix *matrix1, Matrix *matrix2) {
		if ((matrix1->getNoOfColumns() == matrix2->getNoOfLines()))
			return true;
		return false;
	}

	void setEpsilon(double epsilon) {
		this->epsilon = epsilon;
	}

	double getEpsilon() {
		return this->epsilon;
	}

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
	Get the transpose of the matrix
	The transpose is the matrix where the element A[i, j] is placed in A[j, i]:
		/ a11 a12 a13 \				 / a11 a21 a31 \
	A = | a21 a22 a23 |		=>	AT = | a12 a22 a32 |
		\ a31 a32 a33 /				 \ a13 a23 a33 /
	@param matrix - the matrix to be transposed
	@return - the transpose of the matrix
	*/
	virtual Matrix *transpose(Matrix *matrix) = 0;

	/*
	Add the current matrix with the matrix given as parameter
	Matrix addition is performed element by element
	@return - the sum of the matrices
	*/
	virtual Matrix *add(Matrix *matrix) = 0;

	/*
	Add the two matrices given as parameters
	Matrix addition is performed element by element
	@param matrix1 - the first matrix
	@param matrix2 - the second matrix
	@return - the sum of the matrices
	*/
	virtual Matrix *add(Matrix *matrix1, Matrix *matrix2) = 0;

	/*
	Multiply the current matrix with the matrix given as parameter
	Matrix multiplication is performed like this:
	Given A and B two matrice, A of size mxn and B of size pxq, where n = p, we define AxB to be:
		the sum of all products formed by taking elements from each line of A and each column of B

		/ a11 a12 a13 \		
	A = | a21 a22 a23 |		
		\ a31 a32 a33 /3x3		

		/ b11 b12 \
	B = | b21 b22 |
		\ b31 b32 /3x2

			  / a11 * b11 + a12 * b21 + a13 * b31  a11 * b12 + a12 * b22 + a13 * b32 \
	M = AxB = | a21 * b11 + a22 * b21 + a23 * b31  a21 * b12 + a22 * b22 + a23 * b32 |
			  \ a31 * b11 + a32 * b21 + a33 * b31  a31 * b12 + a32 * b22 + a33 * b32 /

	@param matrix1 - the first matrix
	@param matrix2 - the second matrix
	@return - the product of the matrices
	*/
	virtual Matrix *multiply(Matrix *matrix) = 0;

	/*
	Multiply the two matrices given as parameters
	Matrix multiplication is performed like this:
	Given A and B two matrice, A of size mxn and B of size pxq, where n = p, we define AxB to be:
	the sum of all products formed by taking elements from each line of A and each column of B

		/ a11 a12 a13 \
	A = | a21 a22 a23 |
		\ a31 a32 a33 /3x3

		/ b11 b12 \
	B = | b21 b22 |
		\ b31 b32 /3x2

			  / a11 * b11 + a12 * b21 + a13 * b31  a11 * b12 + a12 * b22 + a13 * b32 \
	M = AxB = | a21 * b11 + a22 * b21 + a23 * b31  a21 * b12 + a22 * b22 + a23 * b32 |
			  \ a31 * b11 + a32 * b21 + a33 * b31  a31 * b12 + a32 * b22 + a33 * b32 /

	@param matrix1 - the first matrix
	@param matrix2 - the second matrix
	@return - the product of the matrices
	*/
	virtual Matrix *multiply(Matrix *matrix1, Matrix *matrix2) = 0;

	/*
	The QR decomposition algorithm applied for the current matrix (only for square matrices)
	Q is an orthogonal matrix (transpose(Q) * Q = Q * transpose(Q) = In)
	R is a superior triangular matrix
	@param Q - the Q matrix, output parameter
	@param R - the R matrix, output paramater
	*/
	virtual void qrDecomposition(Matrix **Q, Matrix **R) = 0;

	/*
	The QR decomposition algorithm (only for square matrices)
	Q is an orthogonal matrix (transpose(Q) * Q = Q * transpose(Q) = In)
	R is a superior triangular matrix
	@param A - the input matrix
	@param Q - the Q matrix, output parameter
	@param R - the R matrix, output paramater
	*/
	virtual void qrDecomposition(Matrix *A, Matrix **Q, Matrix **R) = 0;

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
	double epsilon;
};

