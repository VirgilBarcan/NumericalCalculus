# pragma once

# include "Matrix.h"

class MatrixNaive : public Matrix
{
public:
	MatrixNaive();
	MatrixNaive(int noOfLines, int noOfColumns);
	~MatrixNaive();

	void setElementAt(int line, int column, double value) override;
	double getElementAt(int line, int column) override;

	Matrix *getLine(int line) override;
	Matrix *getColumn(int column) override;

	void getFromFile(std::string filePath) override;
	void generateRandomMatrixValues(double min, double max) override;

	static Matrix *identityMatrix(int n);

	Matrix *transpose() override;
	Matrix *transpose(Matrix *matrix) override;

	Matrix *inverse() override;
	Matrix *inverse(Matrix *matrix) override;

	Matrix *add(Matrix *matrix) override;
	Matrix *add(Matrix *matrix1, Matrix *matrix2) override;

	Matrix *subtract(Matrix *matrix2) override;
	Matrix *subtract(Matrix *matrix1, Matrix *matrix2) override;

	Matrix *multiply(Matrix *matrix) override;
	Matrix *multiply(Matrix *matrix1, Matrix *matrix2) override;

	virtual Matrix *strictSuperiorPart() override;
	virtual Matrix *strictSuperiorPart(Matrix *A) override;

	virtual Matrix *diagonal() override;
	virtual Matrix *diagonal(Matrix *A) override;

	virtual Matrix *strictInferiorPart() override;
	virtual Matrix *strictInferiorPart(Matrix *A) override;

	void qrDecomposition(Matrix *b, Matrix **Q, Matrix **R) override;
	void qrDecomposition(Matrix *A, Matrix *b, Matrix **Q, Matrix **R) override;

	Matrix *inverseSubstitutionMethod(Matrix *b) override;
	Matrix *inverseSubstitutionMethod(Matrix *A, Matrix *b) override;

	bool isSuperiorTriangular() override;
	bool isSuperiorTriangular(Matrix *A) override;

	double superiorTriangularMatrixDeterminant() override;
	double superiorTriangularMatrixDeterminant(Matrix *A) override;

	virtual bool isDiagonalZero() override;
	virtual bool isDiagonalZero(Matrix *A) override;

	virtual bool gaussEliminationMethod(Matrix *b, Matrix *R, Matrix *B) override;
	virtual bool gaussEliminationMethod(Matrix *A, Matrix *b, Matrix *R, Matrix *B) override;

	bool isSymmetric() override;
	bool isSymmetric(Matrix *A) override;

	Matrix *clone() override;
	Matrix *clone(Matrix *M) override;

	bool equals(Matrix *M) override;
	bool equals(Matrix *A, Matrix *M) override;

	std::string toString() override;

private:
	double **matrix;

	void instantiateMatrix();
	void deinstantiateMatrix();


	void partialPivoting(int l, Matrix *A, Matrix *b);
};

