#include "MatrixNaive.h"

MatrixNaive::MatrixNaive(int noOfLines, int noOfColumns)
{
	this->noOfLines = noOfLines;
	this->noOfColumns = noOfColumns;
	instantiateMatrix();
}

MatrixNaive::~MatrixNaive()
{
	deinstantiateMatrix();
}

void MatrixNaive::instantiateMatrix()
{
	this->matrix = new double*[noOfLines];

	for (int i = 0; i < noOfLines; ++i) {
		this->matrix[i] = new double[noOfColumns];
	}
}

void MatrixNaive::deinstantiateMatrix()
{
	for (int i = 0; i < noOfLines; ++i) {
		delete[] matrix[i];
	}

	delete[] matrix;
}

void MatrixNaive::addElementAt(int line, int column, double value)
{
	if (checkBounds(line, column)) {
		matrix[line][column] = value;
	}
	else {
		//maybe throw an exception
	}
}

double MatrixNaive::getElementAt(int line, int column)
{
	if (checkBounds(line, column)) {
		return matrix[line][column];
	}
	else {
		//maybe throw an exception
		return 0.0;
	}
}

Matrix *MatrixNaive::transpose()
{
	return transpose(this);
}

Matrix *MatrixNaive::transpose(Matrix *matrix)
{
	MatrixNaive *T = new MatrixNaive(matrix->getNoOfColumns(), matrix->getNoOfLines());

	for (int line = 0; line < noOfLines; ++line) {
		for (int column = 0; column < noOfColumns; ++column) {
			T->addElementAt(column, line, matrix->getElementAt(line, column));
		}
	}

	return T;
}

Matrix *MatrixNaive::add(Matrix *matrix)
{
	return add(this, matrix);
}

Matrix *MatrixNaive::add(Matrix *matrix1, Matrix *matrix2)
{
	if (checkEqualSizes(matrix1, matrix2)) {
		MatrixNaive *sum = new MatrixNaive(matrix1->getNoOfLines(), matrix1->getNoOfLines());

		for (int line = 0; line < matrix1->getNoOfLines(); ++line) {
			for (int column = 0; column < matrix1->getNoOfColumns(); ++column) {
				sum->addElementAt(line, column, matrix1->getElementAt(line, column) + matrix2->getElementAt(line, column));
			}
		}

		return sum;
	}
	return nullptr;
}

Matrix *MatrixNaive::multiply(Matrix *matrix)
{
	return multiply(this, matrix);
}

Matrix *MatrixNaive::multiply(Matrix *matrix1, Matrix *matrix2)
{
	//TODO: Implement function
	return nullptr;
}

std::string MatrixNaive::toString()
{
	std::string stringVersion;

	for (int line = 0; line < noOfLines; ++line) {
		for (int column = 0; column < noOfColumns; ++column) {
			stringVersion.append(std::to_string(getElementAt(line, column)));
			stringVersion.append(" ");
		}
		stringVersion.append("\n");
	}
	return stringVersion;
}
