# include <cstdio>
# include <ctime>
# include "MatrixNaive.h"
# include "MatrixSparse.h"
# include "VectorialNorm.h"
# include "MatrixNorm.h"
# include "Polynomial.h"

# include "Function.h"
# include "ConstantFunction.h"
# include "PowerFunction.h"
# include "CompositeFunction.h"

//# include <armadillo>

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

/*
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

	// begin = clock();
	//calculate the QR decomposition of A using our implementation
	A->qrDecomposition(b->clone(), reinterpret_cast<Matrix**>(&Q), reinterpret_cast<Matrix**>(&R));
	//printf("Q:\n%s\n", Q->toString().c_str());
	//printf("R:\n%s\n", R->toString().c_str());
	// end = clock();
	elapsed_time_qr_our = double(end - begin) / CLOCKS_PER_SEC;

	// begin = clock();
	//calculate the solution of Ax = b with our implementation
	void *ptr = R->inverseSubstitutionMethod(Q->transpose()->multiply(b));
	if (!ptr)
		return;

	xHouseholder = reinterpret_cast<MatrixNaive*>(ptr);

	//printf("xHousehoulder:\n%s\n", xHouseholder->toString().c_str());
	// end = clock();
	elapsed_time_solve_our = double(end - begin) / CLOCKS_PER_SEC;

	// begin = clock();
	//calculate the QR decomposition of A using Armadillo
	arma::qr(Q2, R2, A2);
	// end = clock();
	elapsed_time_qr_arma = double(end - begin) / CLOCKS_PER_SEC;

	//calculate the solution of Ax = b with Armadillo
	b2 = copyFromMatrixToArmadilloVec(b);
	// begin = clock();
	xQR = arma::solve(A2, b2);
	// end = clock();
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

	//for (int i = 0; i < n; ++i)
	//	A->setElementAt(0, i, 3.0),
	//	A->setElementAt(1, i, 6.0);

//	printf("A:\n%s\n", A->toString().c_str());
//
//	inverse = reinterpret_cast<MatrixNaive*>(A->inverse());
//	printf("A^(-1):\n%s\n", inverse->toString().c_str());
//
//	double norm = MatrixNorm::MaximumColumnSumNorm(A->multiply(inverse)->subtract(I));
//	printf("||A*A^(-1) - I|| = %.32f\n", norm);
//
//	delete A;
//	delete I;
//	delete inverse;
//}

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

	A->getFromFile("C:\\Users\\MoroJr\\Source\\Repos\\NumericalCalculus\\HW4Files/a_mat.txt");
	b_A->getFromFile("C:\\Users\\MoroJr\\Source\\Repos\\NumericalCalculus\\HW4Files/a_vect.txt");

	B_line->getFromFile("C:\\Users\\MoroJr\\Source\\Repos\\NumericalCalculus\\HW4Files/b_mat.txt");
	B_column->getFromFile("C:\\Users\\MoroJr\\Source\\Repos\\NumericalCalculus\\HW4Files/b_mat.txt");
	b_B->getFromFile("C:\\Users\\MoroJr\\Source\\Repos\\NumericalCalculus\\HW4Files/b_vect.txt");

	AplusB->getFromFile("C:\\Users\\MoroJr\\Source\\Repos\\NumericalCalculus\\HW4Files/aplusb_mat.txt");
	b_AplusB->getFromFile("C:\\Users\\MoroJr\\Source\\Repos\\NumericalCalculus\\HW4Files/aplusb_vect.txt");

	AoriB->getFromFile("C:\\Users\\MoroJr\\Source\\Repos\\NumericalCalculus\\HW4Files/aorib_mat.txt");
	b_AoriB->getFromFile("C:\\Users\\MoroJr\\Source\\Repos\\NumericalCalculus\\HW4Files/aorib_vect.txt");

//	printf("A[0, 0]: %f\n", A->getElementAt(0, 0));
//	printf("B_c[0, 0]: %f\n", B_column->getElementAt(0, 0));
//	printf("B_c[0, 92]: %f\n", B_column->getElementAt(0, 92));

	//Calculate S = A + B
	S = reinterpret_cast<MatrixSparse*>(A->add(B_line)); S->setEpsilon(epsilon);
//	printf("S[0, 0]: %f\n", S->getElementAt(0, 0));
//	printf("S[0, 186]: %f\n", S->getElementAt(0, 186));
//	printf("S[5, 6]: %f\n", S->getElementAt(5, 6));

	double // begin = clock();
	//Calculate P = A * B
	P = reinterpret_cast<MatrixSparse*>(A->multiply(B_column)); P->setEpsilon(epsilon);
	//printf("A*B:\n%s\n", P->toString().c_str());
	double // end = clock();
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

	A->getFromFile("/home/virgil/Facultate/An3/Sem2/CN/Laborator/5/m_rar_2016_4_mat.txt");
	b->getFromFile("/home/virgil/Facultate/An3/Sem2/CN/Laborator/5/m_rar_2016_4_vect.txt");

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
*/

/*
//A and B should be vectors with the same size
double scalarProduct(Matrix *A, Matrix *B) {
	double result = 0.0;

	for (int i = 0; i < A->getNoOfLines(); ++i) {
		result += A->getElementAt(i, 0) * B->getElementAt(i, 0);
	}

	return result;
}

bool powerMethod(Matrix *A, double &eigenvalue, Matrix **eigenvector) {
	bool result = true;

	double lambdaK;
	double lambdaKPlusOne;
	int k;
	int k_max = 1000000;

	clock_t begin, end;

	// begin = clock();
	//pick v randomly, but with ||v|| = 1
	MatrixSparse *v = new MatrixSparse(A->getNoOfLines(), 1); v->setStoreType(COLUMN);
	v->generateRandomMatrixValues(0.1, 1);
//	v->setElementAt(0, 0, 1); v->setElementAt(1, 0, 0); v->setElementAt(2, 0, 0);
	// end = clock();
//	printf("generate took: %d\n", (end - begin) / CLOCKS_PER_SEC);

	double norm = VectorialNorm::EuclideanNorm(v);

	for (int line = 0; line < v->getNoOfLines(); ++line) {
		v->setElementAt(line, 0, v->getElementAt(line, 0) / norm);
	}

	//printf("v:\n%s\n", v->toString().c_str());

	// begin = clock();
	MatrixSparse *w = new MatrixSparse(A->getNoOfLines(), 1); w->setStoreType(COLUMN);
	w = reinterpret_cast<MatrixSparse*>(A->multiply(v));
	// end = clock();
//	printf("w = A * v took: %d\n", (end - begin) / CLOCKS_PER_SEC);
	//printf("w:\n%s\n", w->toString().c_str());

	MatrixSparse *lambdaV = new MatrixSparse(A->getNoOfLines(), 1); lambdaV->setStoreType(COLUMN);

	lambdaK = scalarProduct(w, v);
	//printf("lambda = %f\n", lambdaK);

	k = 0;

	do {
		fprintf(stderr, "======= %d %f\n", k, norm);
		//lambdaV = v(k)
		// // begin = clock();
		lambdaV = reinterpret_cast<MatrixSparse*>(v->clone());
		// end = clock();
//		printf("clone took: %d\n", (end - begin) / CLOCKS_PER_SEC);

		//printf("v(k):\n%s\n", lambdaV->toString().c_str());

		// begin = clock();
		//v(k+1) = w/norm(w)
		norm = VectorialNorm::EuclideanNorm(w);
		for (int line = 0; line < v->getNoOfLines(); ++line) {
			v->setElementAt(line, 0, w->getElementAt(line, 0) / norm);
		}
		// end = clock();
//		printf("Euclidean Norm took: %d\n", (end - begin) / CLOCKS_PER_SEC);

		//printf("v(k+1):\n%s\n", v->toString().c_str());

		//w = A * v(k+1)
		// begin = clock();
		w = reinterpret_cast<MatrixSparse*>(A->multiply(v));
		// end = clock();
//		printf("w = A * v took: %d\n", (end - begin) / CLOCKS_PER_SEC);

		//printf("w:\n%s\n", w->toString().c_str());

		//lambda(k+1) = (w, v(k+1))
		lambdaKPlusOne = scalarProduct(w, v);

		//printf("lambda_k = %f\n", lambdaK);
		//printf("lambda_k+1 = %f\n\n\n", lambdaKPlusOne);

		k++;

		//calculate lambda * v(k)
		for (int line = 0; line < lambdaV->getNoOfLines(); ++line) {
			lambdaV->setElementAt(line, 0, lambdaV->getElementAt(line, 0) * lambdaK);
		}

		// begin = clock();
		norm = VectorialNorm::EuclideanNorm(w->subtract(lambdaV));
		// end = clock();
//		printf("Euclidean Norm took: %d\n", (end - begin) / CLOCKS_PER_SEC);

		lambdaK = lambdaKPlusOne;

	} while (norm > A->getNoOfLines() * A->getEpsilon() && k <= k_max);

	printf("norm = %.16f\n", norm);
	printf("k = %d\n", k);
	//printf("v:\n%s\n", v->toString().c_str());

	if (k > k_max) {
		//we don't have an eigenvalue
		result = false;
	}

	if (norm <= A->getNoOfLines() * A->getEpsilon()) {
		//we have an eigenvalue and an eigenvector
		eigenvalue = lambdaKPlusOne;

		for (int line = 0; line < v->getNoOfLines(); ++line) {
			(*eigenvector)->setElementAt(line, 0, v->getElementAt(line, 0));
		}
	}

	delete v;
	delete w;
	delete lambdaV;

	return result;
}

MatrixSparse *generateRandomSymmetric(int noOfLines, int noOfColumns, int min, int max) {
	MatrixSparse *A = new MatrixSparse(noOfLines, noOfColumns);
	A->setStoreType(LINE);
	int maxElements = 10;

	std::uniform_real_distribution<double> uniform_distribution(min, max);
	std::default_random_engine random_engine;

	for (int line = 0; line < noOfLines; ++line) {
		int noOfElements = rand() % maxElements + 1;
		for (int i = 0; i < noOfElements; ++i) {
			//generate random value in the interval [min, max]
			int column = rand() % noOfColumns;
			double value = uniform_distribution(random_engine);

			A->setElementAt(line, column, value); A->setElementAt(column, line, value);
		}
	}

	return A;
}

void HW6(int p) {
	double epsilon = pow(10, -p);
	MatrixSparse *A = new MatrixSparse(); A->setStoreType(LINE); A->setEpsilon(epsilon);
	MatrixSparse *B;

	A->getFromFile("m_rar_sim_2016.txt");	
	printf("Is A symmetric? %s\n", A->isSymmetric() ? "yes" : "no");

	B = generateRandomSymmetric(10, 10, 0, 100); B->setEpsilon(epsilon);
	printf("Is B symmetric? %s\n", B->isSymmetric() ? "yes" : "no");

	double eigenvalue;
	MatrixSparse *eigenvector = new MatrixSparse(A->getNoOfLines(), 1); eigenvector->setStoreType(COLUMN);
	powerMethod(A, eigenvalue, reinterpret_cast<Matrix**>(&eigenvector));

	FILE *f = fopen("eigenvector.txt", "w");
	fprintf(f, "Eigenvalue: %f\n\n", eigenvalue);
	
	printf("Eigenvalue: %f\n\n", eigenvalue);

	fprintf(f, "Eigenvector:\n%s\n", eigenvector->toString().c_str());
	// printf("Eigenvector:\n%s\n", eigenvector->toString().c_str());

	fprintf(f, "A*u:\n%s\n", A->multiply(eigenvector)->toString().c_str());
	//printf("A*u:\n%s\n", A->multiply(eigenvector)->toString().c_str());

	fclose(f);

	delete A;
	delete eigenvector;
}

void HW6_SVD(int p, int n, int _s) {
	double epsilon = pow(10, -8);

	if (p <= n) {
		printf("p should be greater than n\n");
		return;
	}
	if (p * n > 36 * 36) {
		printf("p * n should be less than %d\n", 36 * 36);
		return;
	}

	arma::mat A = arma::randu<arma::mat>(p, n);
	arma::mat U;
	arma::vec s;
	arma::mat V;

	arma::svd(U, s, V, A);

	bool bSMax = false, bSMin = false;
	double smax, smin;

	size_t rang = 0;
	printf("Valorile singulare ale matricii A: ");
	for (int i = 0; i < s.size(); ++i) {
		printf("%.4f ", s.at(i));

		if (fabs(s.at(i) - 1.0f) > epsilon)
		{
			if (!bSMax) {
				bSMax = true;
				smax = s.at(i);
			}

			smin = s.at(i);
			bSMin = true;

			++rang;
		}
	}

	printf("\nRangul matricii A: %d", rang);

	if (bSMin)
		printf("\nNumarul de conditionare al matricii A: %.4f\n", smax / smin);
	else
		printf("\nNumarul de conditionare al matricii A nu se poate afla!");

	arma::mat S(p, n);
	S.zeros();
	for (int i = 0; i < n; ++i)
		S(i, i) = s[i];

	arma::mat prod = U * S * V.t();
	arma::mat diff = A - prod;
	cout << diff << endl;

	MatrixNaive *diffNaive = copyFromArmadilloMatToMatrix(diff);

	double maxNorm = MatrixNorm::MaximumNorm(diffNaive);

	printf("\n||A - USV^T||%c: %.16f\n", char(236), maxNorm);

	
	arma::mat As(p, n);
	As.zeros();
	for (int i = 0; i < _s; ++i) {
		arma::vec ui = U.col(i);
		arma::vec vi = V.col(i);

		As = As + s[i] * ui * vi.t();
	}

	cout << endl << As << endl;
	

	arma::mat diff2 = A - As;
	MatrixNaive *diffNaive2 = copyFromArmadilloMatToMatrix(diff2);
	double maxNorm2 = MatrixNorm::MaximumNorm(diffNaive2);

	printf("||A - As||%c: %.16f\n", char(236), maxNorm2);
}
*/

double epsilon = pow(10, -3);

void printRoots(std::vector<double> &roots) {
	for (auto it = roots.begin(); it != roots.end(); ++it) {
		printf("%.8f ", (*it));
	}
}

bool alreadyFound(std::vector<double> &roots, double root) {
	for (auto it = roots.begin(); it != roots.end(); ++it) {
		if (fabs((*it) - root) <= epsilon)
			return true;
	}
	return false;
}

void MullerMethod(Function *f, double R, int no_of_roots_searched) {
	bool ok = false;
	int k_max = 100;

	std::vector<double> roots;

	srand(time(0));

	do {
		int v = ((R - (-R)) + 1);
		double x0 = -R + rand() % v;
		double x1 = -R + rand() % v;
		double x2 = -R + rand() % v;
		double x3 = x2;

		printf("\nx0 = %.8f\n", x0);
		printf("x1 = %.8f\n", x1);
		printf("x2 = %.8f\n", x2);

		double h0, h1;
		double delta0, delta1;
		double a, b, c;
		double delta_x = 0.0;
		double max;
		double delta;

		int k = 2;

		do {
			printf("\n========= k = %d =========\n", k);

			h0 = x1 - x0;
			h1 = x2 - x1;
			printf("\nh0 = %.8f\n", h0);
			printf("h1 = %.8f\n", h1);

			delta0 = (f->evaluate(x1) - f->evaluate(x0)) / h0;
			delta1 = (f->evaluate(x2) - f->evaluate(x1)) / h1;
			printf("\ndelta0 = %.8f\n", delta0);
			printf("delta1 = %.8f\n", delta1);

			a = (delta1 - delta0) / (h1 + h0);
			b = a * h1 + delta1;
			c = f->evaluate(x2);
			printf("\na = %.8f\n", a);
			printf("b = %.8f\n", b);
			printf("c = %.8f\n", c);

			delta = b * b - 4 * a * c;
			printf("\ndelta = %.8f\n", delta);
			if (delta < 0) {
				ok = false;
				break;
			}

			max = fmax(b + sqrt(delta), b - sqrt(delta));
			printf("\nmax = %.8f\n", max);
			if (fabs(max) < epsilon) {
				ok = false;
				break;
			}

			ok = true;
			delta_x = (2 * c) / max;
			printf("\ndelta_x = %.8f\n", delta_x);

			x3 = x2 - delta_x;

			x0 = x1;
			x1 = x2;
			x2 = x3;
			printf("\nx0 = %.8f\n", x0);
			printf("x1 = %.8f\n", x1);
			printf("x2 = %.8f\n", x2);

			++k;

		} while (fabs(delta_x) >= epsilon && k <= k_max && fabs(delta_x) <= pow(10, 8));

		if (ok && fabs(delta_x) < epsilon) {
			//convergence

			//check that this root has not been found already
			if (!alreadyFound(roots, x3)) {
				roots.push_back(x3);

				if (roots.size() == no_of_roots_searched) {
					ok = true;
				}
				else {
					ok = false;
				}
			}
			else {
				ok = false;
			}
		}
		else {
			//divergence
			ok = false;
		}

	} while (!ok);

	printRoots(roots);
}

void HW7() {
	/*
	//Test Function
	double x = 5;
	Function *c = new ConstantFunction(3);
	printf("c(x) = 3; c(%f) = %f\n", x, c->evaluate(x));

	Function *p = new PowerFunction(3);
	printf("p(x) = x ^ 3; p(%f) = %f\n", x, p->evaluate(x));

	CompositeFunction *f = new CompositeFunction();
	f->addFunction(c);
	f->addFunction(p);
	printf("f(x) = 3 + x ^ 3; f(%f) = %f\n", x, f->evaluate(x));

	Function *cf1 = (*c) + (*c);
	printf("f(x) = c(x) + c(x); f(%f) = %f\n", x, cf1->evaluate(x));

	double val = 3;
	Function *cf2 = (*c) * val;
	printf("f(x) = c(x) * 3; f(%f) = %f\n", x, cf2->evaluate(x));
	
	delete c;
	delete p;
	delete f;
	delete cf1;
	delete cf2;
	*/

	double coefficients[] = { -6, 11, -6, 1 };

	Function *x3 = new PowerFunction(coefficients[3], 3);
	Function *x2 = new PowerFunction(coefficients[2], 2);
	Function *x1 = new PowerFunction(coefficients[1], 1);
	Function *x0 = new ConstantFunction(coefficients[0]);

	Function *f = (*x0) + (*x1);
	f = (*f) + (*x2);
	f = (*f) + (*x3);

	double max = 0;
	for (double i : coefficients) {
		if (fabs(i) > max) {
			max = fabs(i);
		}
	}

	double R = fabs(coefficients[3] + max) / fabs(coefficients[3]);

	MullerMethod(f, R, 3);
}

int main() {

	//testNaiveMatrix();
	//testSparseMatrix();

	//Homework 2
	//HW2(250);

	//Homework 3
	//HW3(10);

	//Homework 4
	//HW4();

	//Homework 5
	//HW5(8);

	/*
	std::vector<std::map<int, double>> *list1 = new std::vector<std::map<int, double>>();
	list1->push_back({ {1, 1}, {2, 2} });
	list1->push_back({ {3, 1}, {4, 2} });

	std::vector<std::map<int, double>> *list2 = new std::vector<std::map<int, double>>(*list1);

	list1->push_back({ {5, 1} });
	list2->push_back({ {6, 1} });

	for (auto it = list1->begin(); it != list1->end(); ++it) {
		for (auto it2 = it->begin(); it2 != it->end(); ++it2) {
			printf("%d %f\n", it2->first, it2->second);
		}
		printf("\n");
	}

	printf("\n\n");

	for (auto it = list2->begin(); it != list2->end(); ++it) {
		for (auto it2 = it->begin(); it2 != it->end(); ++it2) {
			printf("%d %f\n", it2->first, it2->second);
		}
	}
	*/

	//Homework 6
	//HW6(3);

	//int s = 5;
	//HW6_SVD(11, 10, s);

	HW7();

	return 0;
}