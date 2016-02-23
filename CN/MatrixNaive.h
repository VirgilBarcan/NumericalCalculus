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

	std::string toString();

private:
	double **matrix;

	void instantiateMatrix();
	void deinstantiateMatrix();
	bool checkBounds(int line, int column);
};

