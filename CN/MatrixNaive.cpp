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

void MatrixNaive::setElementAt(int line, int column, double value)
{
	if (checkBounds(line, column)) {
		matrix[line][column] = value;
	}
	else {
		//maybe throw an exception
		printf("setElementAt: (line, column) pair not in bounds: (%d, %d) != (%d, %d)\n", line, column, this->noOfLines, this->noOfColumns);
	}
}

double MatrixNaive::getElementAt(int line, int column)
{
	if (checkBounds(line, column)) {
		return matrix[line][column];
	}
	else {
		//maybe throw an exception
		printf("getElementAt: (line, column) pair not in bounds: (%d, %d) != (%d, %d)\n", line, column, this->noOfLines, this->noOfColumns);
		return 0.0;
	}
}

Matrix *MatrixNaive::getLine(int line) {
	if (line >= 0 && line < getNoOfLines()) {
		MatrixNaive *theLine = new MatrixNaive(1, getNoOfColumns());

		for (int column = 0; column < getNoOfColumns(); ++column) {
			theLine->setElementAt(0, column, getElementAt(line, column));
		}

		return theLine;
	}

	return nullptr;
}

Matrix *MatrixNaive::getColumn(int column) {
	if (column >= 0 && column < getNoOfColumns()) {
		MatrixNaive *theColumn = new MatrixNaive(getNoOfLines(), 1);

		for (int line = 0; line < getNoOfLines(); ++line) {
			theColumn->setElementAt(line, 0, getElementAt(line, column));
		}

		return theColumn;
	}

	return nullptr;
}

void MatrixNaive::getFromFile(std::string filePath)
{
	//TODO: read from the file the matrix size and values

	//read size from the first 2 lines
	//this->noOfLines = lines;
	//this->noOfColumns = columns;

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

			this->setElementAt(line, column, value);
		}
	}
}

Matrix* MatrixNaive::identityMatrix(int n) {
	MatrixNaive *identity = new MatrixNaive(n, n);

	for (int i = 0; i < n; ++i) {
		identity->setElementAt(i, i, 1);
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
			T->setElementAt(column, line, matrix->getElementAt(line, column));
		}
	}

	return T;
}

Matrix *MatrixNaive::inverse()
{
	return inverse(this);
}

Matrix *MatrixNaive::inverse(Matrix *matrix)
{
	MatrixNaive *inverse = new MatrixNaive(matrix->getNoOfLines(), matrix->getNoOfColumns());
	MatrixNaive *b = new MatrixNaive(matrix->getNoOfLines(), 1);
	MatrixNaive *R = new MatrixNaive(matrix->getNoOfLines(), matrix->getNoOfColumns());
	MatrixNaive *B = new MatrixNaive(matrix->getNoOfLines(), matrix->getNoOfColumns());

	bool isNonsingular = matrix->gaussEliminationMethod(b, R, B);

	if (isNonsingular) {
		for (int column = 0; column < matrix->getNoOfColumns(); ++column) {
			Matrix *x = matrix->inverseSubstitutionMethod(reinterpret_cast<Matrix*>(R), reinterpret_cast<Matrix*>(B->getColumn(column)));

			//copy the vector x to the inverse
			for (int line = 0; line < x->getNoOfLines(); ++line) {
				inverse->setElementAt(line, column, x->getElementAt(line, 0));
			}

			delete x;
		}
	}
	else {
		inverse = nullptr;
		printf("matrix is singular!\n");
	}

	delete b;
	delete R;
	delete B;

	return inverse;
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
				sum->setElementAt(line, column,
								  matrix1->getElementAt(line, column) + matrix2->getElementAt(line, column));
			}
		}

		return sum;
	}
	return nullptr;
}

Matrix *MatrixNaive::subtract(Matrix *matrix2) {
	return subtract(this, matrix2);
}

Matrix *MatrixNaive::subtract(Matrix *matrix1, Matrix *matrix2) {
	if (checkEqualSizes(matrix1, matrix2)) {
		MatrixNaive *difference = new MatrixNaive(matrix1->getNoOfLines(), matrix1->getNoOfLines());

		for (int line = 0; line < matrix1->getNoOfLines(); ++line) {
			for (int column = 0; column < matrix1->getNoOfColumns(); ++column) {
				difference->setElementAt(line, column,
								  matrix1->getElementAt(line, column) - matrix2->getElementAt(line, column));
			}
		}

		return difference;
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
				productMatrix->setElementAt(line1, column2, sum);
				sum = 0.0;
            }
        }

		return productMatrix;
    }
	return nullptr;
}

void MatrixNaive::qrDecomposition(Matrix *b, Matrix **Q, Matrix **R)
{
	qrDecomposition(this, b, Q, R);
}

void MatrixNaive::qrDecomposition(Matrix *A, Matrix *b, Matrix **Q, Matrix **R)
{
	Matrix *A_copy = A->clone();
	if (A_copy->getNoOfLines() == A_copy->getNoOfColumns()) { //the decomposition can be done
		int n = A_copy->getNoOfLines();
		MatrixNaive *Q_tilda = reinterpret_cast<MatrixNaive*>(MatrixNaive::identityMatrix(n));
		MatrixNaive *u = new MatrixNaive(1, n);

		double sigma;
		double beta;
		double gamma;
		double k;

		//our matrix is 0 indexed
		for (int r = 0; r < n - 1; ++r) {
			//build Pr matrix, find beta and u vector

			//calculate sigma
			sigma = 0; //sum with i = r,..., n of a_ir^2
			for (int i = r; i < n; ++i) {
				sigma += pow(A_copy->getElementAt(i, r), 2);
			}

			if (sigma <= getEpsilon())
				break; //A is a singular matrix

			k = sqrt(sigma);

			if (A_copy->getElementAt(r, r) > 0)
				k = -k;

			beta = sigma - k * A_copy->getElementAt(r, r);

			u->setElementAt(0, r, (A_copy->getElementAt(r, r) - k));
			for (int i = r + 1; i < n; ++i) {
				u->setElementAt(0, i, A_copy->getElementAt(i, r));
			}

			//transform the j-th column, j = r+1, ..., n
			for (int j = r + 1; j < n; ++j) {
				gamma = 0;
				for (int i = r; i < n; ++i) {
					gamma += u->getElementAt(0, i) * A_copy->getElementAt(i, j);
				}
				gamma /= beta;

				for (int i = r; i < n; ++i) {
					A_copy->setElementAt(i, j, A_copy->getElementAt(i, j) - gamma * u->getElementAt(0, i));
				}
			}

			//transform the r-th column of A
			A_copy->setElementAt(r, r, k);
			for (int i = r + 1; i < n; ++i) {
				A_copy->setElementAt(i, r, 0);
			}

			gamma = 0;
			for (int i = r; i < n; ++i) {
				gamma = u->getElementAt(0, i) * b->getElementAt(i, 0);
			}
			gamma /= beta;

			for (int i = r; i < n; ++i) {
				b->setElementAt(0, i, b->getElementAt(i, 0) - gamma * u->getElementAt(0, i));
			}

			for (int j = 0; j < n; ++j) {
				gamma = 0;
				for (int i = r; i < n; ++i) {
					gamma += u->getElementAt(0, i) * Q_tilda->getElementAt(i, j);
				}
				gamma /= beta;

				for (int i = r; i < n; ++i) {
					Q_tilda->setElementAt(i, j, Q_tilda->getElementAt(i, j) - gamma * u->getElementAt(0, i));
				}
			}
		}

		//copy the results in Q and R
		//Q_tilda is the transpose of Q
		*Q = reinterpret_cast<MatrixNaive*>(Q_tilda->transpose());

		//R can be found in A
		*R = A_copy;
	}
	else { //the decomposition can not be done
		*Q = nullptr;
		*R = nullptr;
	}
}

Matrix *MatrixNaive::clone()
{
	return clone(this);
}

Matrix *MatrixNaive::clone(Matrix *M)
{
	Matrix *C = new MatrixNaive(M->getNoOfLines(), M->getNoOfColumns());

	for (int line = 0; line < M->getNoOfLines(); ++line) {
		for (int column = 0; column < M->getNoOfColumns(); ++column) {
			C->setElementAt(line, column, M->getElementAt(line, column));
		}
	}

	return C;
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

Matrix *MatrixNaive::inverseSubstitutionMethod(Matrix *b) {
	return inverseSubstitutionMethod(this, b);
}

Matrix *MatrixNaive::inverseSubstitutionMethod(Matrix *A, Matrix *b) {
	Matrix *x;

	//check if A is superior triangular
	if (A->isSuperiorTriangular()) {
		x = new MatrixNaive(A->getNoOfLines(), 1);

		for (int i = A->getNoOfLines() - 1; i >= 0; --i) {
			double x_i = b->getElementAt(i, 0);

			for (int j = i + 1; j < A->getNoOfLines(); ++j) {
				x_i -= A->getElementAt(i, j) * x->getElementAt(j, 0);
			}

			x_i /= A->getElementAt(i, i);

			x->setElementAt(i, 0, x_i);
		}

	}
	else {
		x = nullptr;
	}

	return x;
}

bool MatrixNaive::isSuperiorTriangular() {
	return isSuperiorTriangular(this);
}

bool MatrixNaive::isSuperiorTriangular(Matrix *A) {
	bool result = true;

	//check if the matrix is square
	//if the matrix is not square, return false
	if (A->getNoOfLines() == A->getNoOfColumns()) {
		for (int line = 0; line < A->getNoOfLines() && result; ++line) {
			for (int column = 0; column < line; ++column) {
				if (A->getElementAt(line, column) != 0) {
					result = false;
					break;
				}
			}
		}
	}
	else {
		result = false;
	}

	return result;
}

double MatrixNaive::superiorTriangularMatrixDeterminant() {
	return superiorTriangularMatrixDeterminant(this);
}

double MatrixNaive::superiorTriangularMatrixDeterminant(Matrix *A) {
	double determinant = 1;
	if (isSuperiorTriangular(A)) {
		//we can iterate like that because A should be a square matrix
		for (int i = 0; i < A->getNoOfLines(); ++i) {
			determinant *= A->getElementAt(i, i);
		}
	}
	else {
		determinant = 0;
	}

	return determinant;
}

bool MatrixNaive::gaussEliminationMethod(Matrix *b, Matrix *R, Matrix *B) {
	return gaussEliminationMethod(this, b, R, B);
}

bool MatrixNaive::gaussEliminationMethod(Matrix *A, Matrix *b, Matrix *R, Matrix *B) {
	//build the extended matrix A
	MatrixNaive *extendedA = new MatrixNaive(A->getNoOfLines(), 2 * A->getNoOfColumns());

	for (int line = 0; line < A->getNoOfLines(); ++line) {
		for (int column = 0; column < A->getNoOfColumns(); ++column) {
			extendedA->setElementAt(line, column, A->getElementAt(line, column));
		}
	}

	for (int i = 0; i < A->getNoOfLines(); ++i) {
		extendedA->setElementAt(i, i + A->getNoOfLines(), 1);
	}

//	printf("extendedA: \n%s\n", extendedA->toString().c_str());

	//the algorithm implementation
	int l = 0;

	partialPivoting(l, extendedA, b);

//	printf("extendedA: \n%s\n", extendedA->toString().c_str());
	//printf("b: \n%s\n", b->toString().c_str());

	while ((l < A->getNoOfLines() - 1) && (fabs(extendedA->getElementAt(l, l)) > epsilon)) {
		for (int i = l + 1; i < A->getNoOfLines(); ++i) {
			double f = - extendedA->getElementAt(i, l) / extendedA->getElementAt(l, l);

			for (int j = l + 1; j < 2 * A->getNoOfLines(); ++j) {
				extendedA->setElementAt(i, j, extendedA->getElementAt(i, j) + f * extendedA->getElementAt(l, j));
			}
			extendedA->setElementAt(i, l, 0);

//			printf("extendedA': %d \n%s\n", l, extendedA->toString().c_str());
		}
		++l;
		partialPivoting(l, extendedA, b);

//		printf("extendedA'': %d \n%s\n", l, extendedA->toString().c_str());
		//printf("b: \n%s\n", b->toString().c_str());
	}

	if (fabs(extendedA->getElementAt(l, l)) < epsilon) {
		//singular matrix
		R = nullptr;
		B = nullptr;
		return false;
	}
	//nonsingular matrix
	//copy elements from extendedA to R and B
	for (int line = 0; line < extendedA->getNoOfLines(); ++line) {
		for (int column = 0; column < extendedA->getNoOfColumns(); ++column) {
			R->setElementAt(line, column, extendedA->getElementAt(line, column));
		}
	}

	for (int line = 0; line < extendedA->getNoOfLines(); ++line) {
		for (int column = 0; column < extendedA->getNoOfColumns(); ++column) {
			B->setElementAt(line, column, extendedA->getElementAt(line, column + A->getNoOfColumns()));
		}
	}

	return true;
}

void MatrixNaive::partialPivoting(int l, Matrix *A, Matrix *b) {
	int i0 = -1, max = 0;

	//find i0
	for (int i = l; i < A->getNoOfLines(); ++i) {
		if (max < fabs(A->getElementAt(i, l))) {
			max = fabs(A->getElementAt(i, l));
			i0 = i;
		}
	}

	//change lines i0 and l
	if ( i0 != l) {
		double aux;
		for (int i = 0; i < A->getNoOfColumns(); ++i) {
			aux = A->getElementAt(i0, i);
			A->setElementAt(i0, i, A->getElementAt(l, i));
			A->setElementAt(l, i, aux);
		}

		//change components i0 and l of the vector b
		aux = b->getElementAt(i0, 0);
		b->setElementAt(i0, 0, b->getElementAt(l, 0));
		b->setElementAt(l, 0, aux);
	}
}


