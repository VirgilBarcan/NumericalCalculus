//
// Created by virgil on 21.03.2016.
//

#ifndef NUMERICALCALCULUS_MATRIXSPARSE_H
#define NUMERICALCALCULUS_MATRIXSPARSE_H

# include "Matrix.h"
# include <fstream>
# include <vector>
# include <map>

typedef enum StoreType {
    LINE,
    COLUMN
}StoreType;

class MatrixSparse  : public Matrix
{
public:
    MatrixSparse();
    MatrixSparse(int noOfLines, int noOfColumns);
    ~MatrixSparse();

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

    void qrDecomposition(Matrix *b, Matrix **Q, Matrix **R) override;
    void qrDecomposition(Matrix *A, Matrix *b, Matrix **Q, Matrix **R) override;

    Matrix *inverseSubstitutionMethod(Matrix *b) override;
    Matrix *inverseSubstitutionMethod(Matrix *A, Matrix *b) override;

    bool isSuperiorTriangular() override;
    bool isSuperiorTriangular(Matrix *A) override;

    double superiorTriangularMatrixDeterminant() override;
    double superiorTriangularMatrixDeterminant(Matrix *A) override;

    virtual bool gaussEliminationMethod(Matrix *b, Matrix *R, Matrix *B) override;
    virtual bool gaussEliminationMethod(Matrix *A, Matrix *b, Matrix *R, Matrix *B) override;

    Matrix *clone() override;
    Matrix *clone(Matrix *M) override;

    void setStoreType(StoreType storeType);
    StoreType getStoreType();
    std::map<int, double> getListElements(int lineOrColumn);

    bool equals(Matrix *M) override;
    bool equals(Matrix *A, Matrix *M) override;

    std::string toString() override;

private:
    std::vector<std::map<int, double>> *list;
    StoreType storeType;

    void instantiateMatrix();
    void deinstantiateMatrix();

};


#endif //NUMERICALCALCULUS_MATRIXSPARSE_H
