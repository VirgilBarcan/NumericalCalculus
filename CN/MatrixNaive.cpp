#include "MatrixNaive.h"

MatrixNaive::MatrixNaive()
{
}

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

	//initialize values
	for (int line = 0; line < noOfLines; ++line)
		for (int column = 0; column < noOfColumns; ++column)
			this->matrix[line][column] = 0;
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

void MatrixNaive::getFromFile(std::string filePath)
{
	//TODO: read from the file the matrix size and values

	//read size from the first 2 lines

	//instantiate the matrix
	//instantiateMatrix();

	//read data from the file and build the matrix
}

void MatrixNaive::generateRandomMatrixValues(double min, double max)
{
	std::uniform_real_distribution<double> uniform_distribution(min, max);
	std::default_random_engine random_engine;

	for (int line = 0; line < this->getNoOfLines(); ++line) {
		for (int column = 0; column < this->getNoOfColumns(); ++column) {
			//generate random value in the interval [min, max]
			double value = uniform_distribution(random_engine);

			this->addElementAt(line, column, value);
		}
	}
}

Matrix* MatrixNaive::identityMatrix(int n) {
	MatrixNaive *identity = new MatrixNaive(n, n);

	for (int i = 0; i < n; ++i) {
		identity->addElementAt(i, i, 1);
	}

	return identity;
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
	if (checkMultipliable(matrix1, matrix2)) {
        MatrixNaive *productMatrix = new MatrixNaive(matrix1->getNoOfLines(), matrix2->getNoOfColumns());
        
		int line1, column1, column2;
		double sum = 0;

        for (line1 = 0; line1 < matrix1->getNoOfLines(); ++line1) {
            for (column2 = 0; column2 < matrix2->getNoOfColumns(); ++column2) {
                for (column1 = 0; column1 < matrix1->getNoOfColumns(); ++column1) {
					double product = matrix1->getElementAt(line1, column1) * matrix2->getElementAt(column1, column2);
					sum += product;
                }
				productMatrix->addElementAt(line1, column2, sum);
				sum = 0.0;
            }
        }

		return productMatrix;
    }
	return nullptr;
}

void MatrixNaive::qrDecomposition(Matrix **Q, Matrix **R)
{
	qrDecomposition(this, Q, R);
}

void MatrixNaive::qrDecomposition(Matrix *A, Matrix **Q, Matrix **R)
{
	//TODO: Implement QR decomposition
	if (A->getNoOfLines() == A->getNoOfColumns()) { //the decomposition can be done
		int n = A->getNoOfLines();
		MatrixNaive *Q_tilda = reinterpret_cast<MatrixNaive*>(MatrixNaive::identityMatrix(n));
		MatrixNaive *u = new MatrixNaive(1, n);

		double sigma;
		double k;
		double beta;

		//our matrix is 0 indexed
		for (int r = 0; r < n - 1; ++r) {
			//build Pr matrix, find beta and u vector

			sigma = 0; //sum with i = r,..., n of a_ir^2
			for (int i = r; i < n; ++i) {
				sigma += pow(A->getElementAt(i, r), 2);
			}

			if (sigma <= getEpsilon())
				break;

			k = sqrt(sigma);

			if (A->getElementAt(r, r) > 0)
				k = -k;

			beta = sigma - k * A->getElementAt(r, r);

			u->addElementAt(1, r, (A->getElementAt(r, r) - k));
			for (int i = r + 1; i < n; ++i) {
				u->addElementAt(1, i, A->getElementAt(i, r));
			}
		}
	}
	else { //the decomposition can not be done
		*Q = nullptr;
		*R = nullptr;
	}
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
