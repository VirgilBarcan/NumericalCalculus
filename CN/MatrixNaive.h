#pragma once

#include "Matrix.h"

class MatrixNaive : public Matrix
{
public:
	MatrixNaive();
	MatrixNaive(int noOfLines, int noOfColumns);
	~MatrixNaive();

	void addElementAt(int line, int column, double value);
	double getElementAt(int line, int column);

	void getFromFile(std::string filePath);
	void generateRandomMatrixValues(double min, double max);

	static Matrix *identityMatrix(int n);

	Matrix *transpose();
	Matrix *transpose(Matrix *matrix);

	Matrix *add(Matrix *matrix);
	Matrix *add(Matrix *matrix1, Matrix *matrix2);

	Matrix *multiply(Matrix *matrix);
	Matrix *multiply(Matrix *matrix1, Matrix *matrix2);

	void qrDecomposition(Matrix *b, Matrix **Q, Matrix **R);
	void qrDecomposition(Matrix *A, Matrix *b, Matrix **Q, Matrix **R);

	Matrix *clone();
	Matrix *clone(Matrix *M);

	std::string toString();

private:
	double **matrix;

	void instantiateMatrix();
	void deinstantiateMatrix();
};

