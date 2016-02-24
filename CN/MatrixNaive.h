#pragma once

#include "Matrix.h"

class MatrixNaive : Matrix
{
public:
	MatrixNaive(int noOfLines, int noOfColumns);
	~MatrixNaive();

	void addElementAt(int line, int column, double value);
	double getElementAt(int line, int column);

	Matrix *transpose();
	Matrix *transpose(Matrix *matrix);

	Matrix *add(Matrix *matrix);
	Matrix *add(Matrix *matrix1, Matrix *matrix2);

	Matrix *multiply(Matrix *matrix);
	Matrix *multiply(Matrix *matrix1, Matrix *matrix2);

	std::string toString();

private:
	double **matrix;

	void instantiateMatrix();
	void deinstantiateMatrix();
};

