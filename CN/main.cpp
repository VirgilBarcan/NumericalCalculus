# include <cstdio>
# include "MatrixNaive.h"

using namespace std;

void testNaiveMatrix() {
	printf("Let's test the MatrixNaive\n");

	int noOfLines = 3;
	int noOfColumns = 3;
	MatrixNaive *naiveMatrix = new MatrixNaive(noOfLines, noOfColumns);

	//add values to the matrix
	double value = 0.0;
	for (int line = 0; line < noOfLines; ++line) {
		for (int column = 0; column < noOfColumns; ++column) {
			naiveMatrix->addElementAt(line, column, value);
			value += 1;
		}
	}
	printf("The matrix is:\n%s\n", naiveMatrix->toString().c_str());

	MatrixNaive *transpose = reinterpret_cast<MatrixNaive*>(naiveMatrix->transpose());
	printf("The transpose is:\n%s\n", transpose->toString().c_str());

	MatrixNaive *product = reinterpret_cast<MatrixNaive*>(naiveMatrix->multiply((Matrix*)transpose));
	printf("The product is:\n%s\n", product->toString().c_str());

	MatrixNaive *randomMatrix = new MatrixNaive(5, 5);
	randomMatrix->generateRandomMatrixValues(1, 100);
	printf("The random matrix is:\n%s\n", randomMatrix->toString().c_str());

	delete naiveMatrix;
	delete transpose;
	delete product;
	delete randomMatrix;
}

MatrixNaive *calculateB(MatrixNaive *s, MatrixNaive *A, int n) {
	MatrixNaive *b = new MatrixNaive(1, n);

	for (int i = 0; i < n; ++i) {
		double value = 0;

		for (int j = 0; j < n; ++j) {
			value += s->getElementAt(0, j) * A->getElementAt(i, j);
		}
		b->addElementAt(0, i, value);
	}

	return b;
}

void ex1() {
	int n = 3;
	
	MatrixNaive *s = new MatrixNaive(1, n);
	MatrixNaive *A = new MatrixNaive(n, n);

	s->generateRandomMatrixValues(1, 10);
	A->generateRandomMatrixValues(1, 10);

	MatrixNaive *b = calculateB(s, A, n);

	printf("s:\n%s\n", s->toString().c_str());
	printf("A:\n%s\n", A->toString().c_str());
	printf("b:\n%s\n", b->toString().c_str());

	delete s;
	delete A;
	delete b;
}

int main() {

	//testNaiveMatrix();

	ex1();

	return 0;
}