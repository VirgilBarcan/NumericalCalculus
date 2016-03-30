# include <cstdio>
# include <ctime>
# include "MatrixNaive.h"
# include "MatrixSparse.h"
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

void testSparseMatrix() {
	MatrixSparse *A = new MatrixSparse(3, 3); A->setStoreType(LINE);
	MatrixSparse *B_line = new MatrixSparse(3, 3); B_line->setStoreType(LINE);
	MatrixSparse *B_column = new MatrixSparse(3, 3); B_column->setStoreType(COLUMN);
	MatrixSparse *x = new MatrixSparse(3, 1); x->setStoreType(COLUMN);
	MatrixSparse *S;
	MatrixSparse *P;

	A->setElementAt(0, 0, 1); A->setElementAt(0, 1, 2); A->setElementAt(0, 2, 3);
	A->setElementAt(1, 0, 1); A->setElementAt(1, 1, 1); A->setElementAt(1, 2, 1);
	A->setElementAt(2, 0, 2); A->setElementAt(2, 1, 4); A->setElementAt(2, 2, 3);

	B_line->setElementAt(0, 0, 1); B_line->setElementAt(0, 1, 1); B_line->setElementAt(0, 2, 1);
	B_line->setElementAt(1, 0, 2); B_line->setElementAt(1, 1, 3); B_line->setElementAt(1, 2, 4);
	B_line->setElementAt(2, 0, 1); B_line->setElementAt(2, 1, 2); B_line->setElementAt(2, 2, 1);

	B_column->setElementAt(0, 0, 1); B_column->setElementAt(0, 1, 1); B_column->setElementAt(0, 2, 1);
	B_column->setElementAt(1, 0, 2); B_column->setElementAt(1, 1, 3); B_column->setElementAt(1, 2, 4);
	B_column->setElementAt(2, 0, 1); B_column->setElementAt(2, 1, 2); B_column->setElementAt(2, 2, 1);

	x->setElementAt(0, 0, 1); x->setElementAt(1, 0, 2); x->setElementAt(2, 0, 3);

	printf("x:\n%s\n", x->toString().c_str());

	printf("A:\n%s\n", A->toString().c_str());
	printf("B:\n%s\n", B_line->toString().c_str());

	S = reinterpret_cast<MatrixSparse*>(A->add(B_line));
	printf("A+B:\n%s\n", S->toString().c_str());

	S = reinterpret_cast<MatrixSparse*>(B_column->subtract(B_column));
	printf("B-B:\n%s\n", S->toString().c_str());

	S = reinterpret_cast<MatrixSparse*>(A->add(B_column));
	printf("A+B:\n%s\n", S->toString().c_str());

	P = reinterpret_cast<MatrixSparse*>(A->multiply(B_column));
	printf("A*B:\n%s\n", P->toString().c_str());

	P = reinterpret_cast<MatrixSparse*>(A->multiply(x));
	printf("A*x:\n%s\n", P->toString().c_str());

	printf("A == A? %d\n", A->equals(A));
	printf("A == S? %d\n", A->equals(S));

	delete A;
	delete B_line;
	delete B_column;
	delete x;
	delete S;
	delete P;
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

MatrixSparse *createX(int n) {
	MatrixSparse *x = new MatrixSparse(n, 1);
	x->setStoreType(COLUMN);

	for (int i = 0; i < n; ++i) {
		x->setElementAt(i, 0, i + 1);
	}

	return x;
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
	void *ptr = R->inverseSubstitutionMethod(Q->transpose()->multiply(b));
	if (!ptr)
		return;

	xHouseholder = reinterpret_cast<MatrixNaive*>(ptr);

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
	printf("||Ainit * xHouseholder - binit|| = %.32f\n", first);

	//TODO: calculate ||Ainit * xQR - binit||
	double second = VectorialNorm::EuclideanNorm(A->multiply(copyFromArmadilloVecToMatrix(xQR))->subtract(b));
	printf("||Ainit * xQR - binit|| = %.32f\n", second);

	//calculate ||xHouseholder - s|| / ||s||
	double third = VectorialNorm::EuclideanNorm(xHouseholder->subtract(s)) / VectorialNorm::EuclideanNorm(s);
	printf("||xHouseholder - s|| / ||s|| = %.32f\n", third);

	//TODO: calculate ||xQR - s|| / ||s||
	double fourth = VectorialNorm::EuclideanNorm(copyFromArmadilloVecToMatrix(xQR)->subtract(s)) / VectorialNorm::EuclideanNorm(s);
	printf("||xQR - s|| / ||s|| = %.32f\n", fourth);

	//print statistics
	printf("\nOur implementation:\n");
	printf("\tQR time: %f seconds\n", elapsed_time_qr_our);
	printf("\tSolve time: %f seconds\n", elapsed_time_solve_our);
	printf("\nArmadillo implementation:\n");
	printf("\tQR time: %f seconds\n", elapsed_time_qr_arma);
	printf("\tSolve time: %f seconds\n", elapsed_time_solve_arma);

	if (A) delete A;
	if (b) delete b;
	if (s) delete s;
	if (Q) delete Q;
	if (R) delete R;
	if (xHouseholder) delete xHouseholder;
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

void HW4() {
	double epsilon = 0.0001;
	MatrixSparse *A = new MatrixSparse(); A->setStoreType(LINE);
	MatrixSparse *b_A = new MatrixSparse(); b_A->setStoreType(COLUMN);
	MatrixSparse *B_line = new MatrixSparse(); B_line->setStoreType(LINE);
	MatrixSparse *B_column = new MatrixSparse(); B_column->setStoreType(COLUMN);
	MatrixSparse *b_B = new MatrixSparse(); b_B->setStoreType(COLUMN);
	MatrixSparse *AplusB = new MatrixSparse(); AplusB->setStoreType(LINE);
	MatrixSparse *b_AplusB = new MatrixSparse(); b_AplusB->setStoreType(COLUMN);
	MatrixSparse *AoriB = new MatrixSparse(); AoriB->setStoreType(LINE);
	MatrixSparse *b_AoriB = new MatrixSparse(); b_AoriB->setStoreType(COLUMN);
	MatrixSparse *S;
	MatrixSparse *P;
	MatrixSparse *x;
	MatrixSparse *Ax;
	MatrixSparse *Bx;
	MatrixSparse *AplusBx;
	MatrixSparse *AoriBx;

	A->getFromFile("/home/virgil/Facultate/An3/Sem2/CN/Laborator/4/a_mat.txt");
	b_A->getFromFile("/home/virgil/Facultate/An3/Sem2/CN/Laborator/4/a_vect.txt");

	B_line->getFromFile("/home/virgil/Facultate/An3/Sem2/CN/Laborator/4/b_mat.txt");
	B_column->getFromFile("/home/virgil/Facultate/An3/Sem2/CN/Laborator/4/b_mat.txt");
	b_B->getFromFile("/home/virgil/Facultate/An3/Sem2/CN/Laborator/4/b_vect.txt");

	AplusB->getFromFile("/home/virgil/Facultate/An3/Sem2/CN/Laborator/4/aplusb_mat.txt");
	b_AplusB->getFromFile("/home/virgil/Facultate/An3/Sem2/CN/Laborator/4/aplusb_vect.txt");

	AoriB->getFromFile("/home/virgil/Facultate/An3/Sem2/CN/Laborator/4/aorib_mat.txt");
	b_AoriB->getFromFile("/home/virgil/Facultate/An3/Sem2/CN/Laborator/4/aorib_vect.txt");

//	printf("A[0, 0]: %f\n", A->getElementAt(0, 0));
//	printf("B_c[0, 0]: %f\n", B_column->getElementAt(0, 0));
//	printf("B_c[0, 92]: %f\n", B_column->getElementAt(0, 92));

	//Calculate S = A + B
	S = reinterpret_cast<MatrixSparse*>(A->add(B_line)); S->setEpsilon(epsilon);
//	printf("S[0, 0]: %f\n", S->getElementAt(0, 0));
//	printf("S[0, 186]: %f\n", S->getElementAt(0, 186));
//	printf("S[5, 6]: %f\n", S->getElementAt(5, 6));

	double begin = clock();
	//Calculate P = A * B
	P = reinterpret_cast<MatrixSparse*>(A->multiply(B_column)); P->setEpsilon(epsilon);
	//printf("A*B:\n%s\n", P->toString().c_str());
	double end = clock();
	double elapsed_time = double(end - begin) / CLOCKS_PER_SEC;
	printf("Multiplication took: %f seconds\n", elapsed_time);

	printf("S == AplusB? %s\n", S->equals(AplusB) == true ? "true" : "false");
	printf("P == AoriB? %s\n", P->equals(AoriB)  == true ? "true" : "false");

	//create x
	x = createX(A->getNoOfLines());
	printf("Created x\n");

	//A*x
	Ax = reinterpret_cast<MatrixSparse*>(A->multiply(x)); Ax->setEpsilon(epsilon);
	printf("A*x == b_A? %s\n", Ax->equals(b_A)  == true ? "true" : "false");

	//B*x
	Bx = reinterpret_cast<MatrixSparse*>(B_line->multiply(x)); Bx->setEpsilon(epsilon);
	printf("B*x == b_B? %s\n", Bx->equals(b_B)  == true ? "true" : "false");

	//AplusB*x
	AplusBx = reinterpret_cast<MatrixSparse*>(AplusB->multiply(x)); AplusBx->setEpsilon(epsilon);
	printf("AplusB*x == b_AplusB? %s\n", AplusBx->equals(b_AplusB)  == true ? "true" : "false");

	//AoriB*x
	AoriBx = reinterpret_cast<MatrixSparse*>(AoriB->multiply(x)); AoriBx->setEpsilon(epsilon);
	printf("AoriB*x == b_AoriB? %s\n", AoriBx->equals(b_AoriB)  == true ? "true" : "false");

	delete A;
	delete AplusB;
	delete AoriB;
	delete B_line;
	delete B_column;
	delete S;
	delete P;
	delete x;
	delete Ax;
	delete Bx;
	delete AplusBx;
	delete AoriBx;
}

void HW5(int p) {
	double epsilon = pow(10, -p);
	MatrixSparse *A = new MatrixSparse(); A->setStoreType(LINE); A->setEpsilon(epsilon);
	MatrixSparse *b = new MatrixSparse(); b->setStoreType(COLUMN); A->setEpsilon(epsilon);
	MatrixSparse *x = new MatrixSparse();

	A->getFromFile("/home/virgil/Facultate/An3/Sem2/CN/Laborator/5/m_rar_2016_1_mat.txt");
	b->getFromFile("/home/virgil/Facultate/An3/Sem2/CN/Laborator/5/m_rar_2016_1_vect.txt");

	x = reinterpret_cast<MatrixSparse*>(A->sorMethod(b));
	printf("x=\n%s\n", x->toString().c_str());

	MatrixSparse *Ax = new MatrixSparse();
	Ax = reinterpret_cast<MatrixSparse*>(A->multiply(x));
	//printf("b=\n%s\n", b->toString().c_str());
	//printf("Ax=\n%s\n", Ax->toString().c_str());

	double norm = MatrixNorm::MaximumRowSumNormSparse(reinterpret_cast<MatrixSparse*>(Ax->subtract(b)));
	printf("||Ax - b|| = %.16f\n", norm);

	delete A;
	delete b;
	delete x;
	delete Ax;
}

int main() {

	//testNaiveMatrix();
	//testSparseMatrix();

	//Homework 2
	//HW2(250);

	//Homework 3
	//HW3(50);

	//Homework 4
	//HW4();

	//Homework 5
	HW5(15);

	return 0;
}