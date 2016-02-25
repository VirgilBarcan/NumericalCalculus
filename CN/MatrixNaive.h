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

