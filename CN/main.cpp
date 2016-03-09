# include <cstdio>
# include <ctime>
# include "MatrixNaive.h"
# include "VectorialNorm.h"
# include <armadillo>

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

void ex2() {
	int n = 250;

	MatrixNaive *A = new MatrixNaive(n, n);
	MatrixNaive *b = new MatrixNaive(1, n);
	MatrixNaive *Q = new MatrixNaive(n, n);
	MatrixNaive *R = new MatrixNaive(n, n);

	A->generateRandomMatrixValues(0, 100);
	b->generateRandomMatrixValues(0, 100);

	//check how much it takes to get the QR decomposition
	clock_t begin = clock();
	A->qrDecomposition(b, reinterpret_cast<Matrix**>(&Q), reinterpret_cast<Matrix**>(&R));
	clock_t end = clock();

	printf("A:\n%s\n", A->toString().c_str());
	printf("Q:\n%s\n", Q->toString().c_str());
	printf("R:\n%s\n", R->toString().c_str());

	printf("Q*R:\n%s\n", Q->multiply(R)->toString().c_str());

	double elapsed_time = double(end - begin);
	printf("The QR decomposition took: %f seconds\n", elapsed_time / CLOCKS_PER_SEC);

	arma::mat A2 = arma::randu<arma::mat>(n, n);
	arma::mat Q2, R2;
	
	arma::qr(Q2, R2, A2);

	//cout << "Q2.n_rows: " << Q2.n_rows << endl;
	//cout << "Q2.n_cols: " << Q2.n_cols << endl;

	//cout << "R2.n_rows: " << R2.n_rows << endl;
	//cout << "R2.n_cols: " << R2.n_cols << endl;

	//cout << "A2.n_rows: " << A2.n_rows << endl;
	//cout << "A2.n_cols: " << A2.n_cols << endl;

	A2.print("A2:");
	Q2.print("Q2:");
	R2.print("R2:");


	delete A;
	delete Q;
	delete R;
}

void ex3() {

}

void ex2_2() {
	int n = 3;

	MatrixNaive *A = new MatrixNaive(n, n);
	MatrixNaive *b = new MatrixNaive(1, n);
	MatrixNaive *Q = new MatrixNaive(n, n);
	MatrixNaive *R = new MatrixNaive(n, n);

	A->addElementAt(0, 0, 1); A->addElementAt(0, 1, 3); A->addElementAt(0, 2, 1);
	A->addElementAt(1, 0, 3); A->addElementAt(1, 1, 2); A->addElementAt(1, 2, 3);
	A->addElementAt(2, 0, 2); A->addElementAt(2, 1, 5); A->addElementAt(2, 2, -2);
	printf("A:\n%s\n", A->toString().c_str());

	b->addElementAt(0, 0, 10); b->addElementAt(0, 1, 16); b->addElementAt(0, 2, 6);
	printf("b:\n%s\n", b->toString().c_str());

	printf("EuclideanNorm of b: %f\n", VectorialNorm::EuclideanNorm(b));
	printf("ManhattanNorm of b: %f\n", VectorialNorm::ManhattanNorm(b));
	printf("MaximumNorm of b: %f\n", VectorialNorm::MaximumNorm(b));
	printf("pNorm(1) of b: %f\n", VectorialNorm::pNorm(b, 1));
	printf("pNorm(2) of b: %f\n", VectorialNorm::pNorm(b, 2));
	printf("pNorm(3) of b: %f\n\n", VectorialNorm::pNorm(b, 3));

	A->qrDecomposition(b, reinterpret_cast<Matrix**>(&Q), reinterpret_cast<Matrix**>(&R));

	printf("A:\n%s\n", A->toString().c_str());
	printf("b:\n%s\n", b->toString().c_str());
	printf("Q:\n%s\n", Q->toString().c_str());
	printf("R:\n%s\n", R->toString().c_str());

	delete A;
	delete Q;
	delete R;
}

int main() {

	//testNaiveMatrix();

	//ex1();

	ex2();

	//ex2_2();

	return 0;
}