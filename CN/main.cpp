# include <cstdio>
# include <ctime>
# include "MatrixNaive.h"
# include "VectorialNorm.h"
# include "MatrixNorm.h"
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
			naiveMatrix->setElementAt(line, column, value);
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
	MatrixNaive *b = new MatrixNaive(n, 1);

	for (int i = 0; i < n; ++i) {
		double value = 0;

		for (int j = 0; j < n; ++j) {
			value += s->getElementAt(j, 0) * A->getElementAt(i, j);
		}
		b->setElementAt(i, 0, value);
	}

	return b;
}

MatrixNaive *copyFromArmadilloMatToMatrix(arma::mat A) {
	MatrixNaive *M = new MatrixNaive(A.n_rows, A.n_cols);

	for (int line = 0; line < A.n_rows; ++line) {
		for (int column = 0; column < A.n_cols; ++column) {
			M->setElementAt(line, column, A.at(line, column));
		}
	}

	return M;
}

MatrixNaive *copyFromArmadilloVecToMatrix(arma::vec A) {
	MatrixNaive *M = new MatrixNaive(A.n_rows, 1);

	for (int line = 0; line < A.n_rows; ++line) {
		M->setElementAt(line, 0, A(line));
	}

	return M;
}

arma::vec copyFromMatrixToArmadilloVec(MatrixNaive *A) {
	arma::vec M = arma::ones(A->getNoOfLines());

	for (int line = 0; line < A->getNoOfLines(); ++line) {
		M(line) = A->getElementAt(line, 0);
	}

	return M;
}

void HW2(int n) {
	printf("Homework 2: QR + Inverse substitution method!\n");
	MatrixNaive *A;
	MatrixNaive *Q = new MatrixNaive(n, n);
	MatrixNaive *R = new MatrixNaive(n, n);
	MatrixNaive *s = new MatrixNaive(n, 1);
	MatrixNaive *b;
	MatrixNaive *xHouseholder;

	arma::mat A2 = arma::randu<arma::mat>(n, n);
	arma::mat Q2, R2;
	arma::vec b2;
	arma::vec xQR;

	double elapsed_time_qr_our;
	double elapsed_time_qr_arma;
	double elapsed_time_solve_our;
	double elapsed_time_solve_arma;
	clock_t begin, end;

	//copy the elements from A2 to our A
	A = copyFromArmadilloMatToMatrix(A2);
	//A->setElementAt(0, 0, 1); A->setElementAt(0, 1, 3); A->setElementAt(0, 2, 1);
	//A->setElementAt(1, 0, 3); A->setElementAt(1, 1, 2); A->setElementAt(1, 2, 3);
	//A->setElementAt(2, 0, 2); A->setElementAt(2, 1, 5); A->setElementAt(2, 2, -2);
	//A->generateRandomMatrixValues(0, 100);
	//printf("A:\n%s\n", A->toString().c_str());

	s->generateRandomMatrixValues(0, 100);

	//b->setElementAt(0, 0, 10);
	//b->setElementAt(1, 0, 16);
	//b->setElementAt(2, 0, 6);
	b = calculateB(s, A, n);
	//printf("b:\n%s\n", b->toString().c_str());

	begin = clock();
	//calculate the QR decomposition of A using our implementation
	A->qrDecomposition(b->clone(), reinterpret_cast<Matrix**>(&Q), reinterpret_cast<Matrix**>(&R));
	//printf("Q:\n%s\n", Q->toString().c_str());
	//printf("R:\n%s\n", R->toString().c_str());
	end = clock();
	elapsed_time_qr_our = double(end - begin) / CLOCKS_PER_SEC;

	begin = clock();
	//calculate the solution of Ax = b with our implementation
	xHouseholder = reinterpret_cast<MatrixNaive*>(R->inverseSubstitutionMethod(Q->transpose()->multiply(b)));
	//printf("xHousehoulder:\n%s\n", xHouseholder->toString().c_str());
	end = clock();
	elapsed_time_solve_our = double(end - begin) / CLOCKS_PER_SEC;

	begin = clock();
	//calculate the QR decomposition of A using Armadillo
	arma::qr(Q2, R2, A2);
	end = clock();
	elapsed_time_qr_arma = double(end - begin) / CLOCKS_PER_SEC;

	//calculate the solution of Ax = b with Armadillo
	b2 = copyFromMatrixToArmadilloVec(b);
	begin = clock();
	xQR = arma::solve(A2, b2);
	end = clock();
	elapsed_time_solve_arma = double(end - begin) / CLOCKS_PER_SEC;

	//xQR.print("xQR:");

	//calculate the norms
	//calculate ||Ainit * xHouseholder - binit||
	double first = VectorialNorm::EuclideanNorm(A->multiply(xHouseholder)->subtract(b));
	printf("||Ainit * xHouseholder - binit|| = %f\n", first);

	//TODO: calculate ||Ainit * xQR - binit||
	double second = VectorialNorm::EuclideanNorm(A->multiply(copyFromArmadilloVecToMatrix(xQR))->subtract(b));
	printf("||Ainit * xQR - binit|| = %f\n", second);

	//calculate ||xHouseholder - s|| / ||s||
	double third = VectorialNorm::EuclideanNorm(xHouseholder->subtract(s)) / VectorialNorm::EuclideanNorm(s);
	printf("||xHouseholder - s|| / ||s|| = %f\n", third);

	//TODO: calculate ||xQR - s|| / ||s||
	double fourth = VectorialNorm::EuclideanNorm(copyFromArmadilloVecToMatrix(xQR)->subtract(s)) / VectorialNorm::EuclideanNorm(s);
	printf("||xQR - s|| / ||s|| = %f\n", fourth);

	//print statistics
	printf("\nOur implementation:\n");
	printf("\tQR time: %f seconds\n", elapsed_time_qr_our);
	printf("\tSolve time: %f seconds\n", elapsed_time_solve_our);
	printf("\nArmadillo implementation:\n");
	printf("\tQR time: %f seconds\n", elapsed_time_qr_arma);
	printf("\tSolve time: %f seconds\n", elapsed_time_solve_arma);

	delete A;
	delete b;
	delete s;
	delete Q;
	delete R;
	delete xHouseholder;
}

void HW3(int n) {
	printf("Homework3: Gauss elimination method!\n");
	MatrixNaive *A = new MatrixNaive(n, n);
	MatrixNaive *I = reinterpret_cast<MatrixNaive*>(MatrixNaive::identityMatrix(n));
	MatrixNaive *inverse = new MatrixNaive(n, n);

	//first example
//	A->setElementAt(0, 0, 1); A->setElementAt(0, 1, 0); A->setElementAt(0, 2, 2);
//	A->setElementAt(1, 0, 0); A->setElementAt(1, 1, 1); A->setElementAt(1, 2, 0);
//	A->setElementAt(2, 0, 1); A->setElementAt(2, 1, 1); A->setElementAt(2, 2, 1);

	//second example
//	A->setElementAt(0, 0, 3); A->setElementAt(0, 1, 0); A->setElementAt(0, 2, 1);
//	A->setElementAt(1, 0, 0); A->setElementAt(1, 1, 1); A->setElementAt(1, 2, 1);
//	A->setElementAt(2, 0, 6); A->setElementAt(2, 1, 1); A->setElementAt(2, 2, 4);

	//random matrix
	A->generateRandomMatrixValues(-100, 100);

	printf("A:\n%s\n", A->toString().c_str());

	inverse = reinterpret_cast<MatrixNaive*>(A->inverse());
	printf("A^(-1):\n%s\n", inverse->toString().c_str());

	double norm = MatrixNorm::MaximumColumnSumNorm(A->multiply(inverse)->subtract(I));
	printf("||A*A^(-1) - I|| = %.32f\n", norm);

	delete A;
	delete I;
	delete inverse;
}

int main() {

	//testNaiveMatrix();

	//Homework 2
//	HW2(250);

	//Homework 3
	HW3(50);

	return 0;
}