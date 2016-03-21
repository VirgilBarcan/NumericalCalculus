//
// Created by virgil on 21.03.2016.
//

# include "MatrixSparse.h"

MatrixSparse::MatrixSparse()
{
}

MatrixSparse::MatrixSparse(int noOfLines, int noOfColumns)
{
    this->noOfLines = noOfLines;
    this->noOfColumns = noOfColumns;
    instantiateMatrix();
}

MatrixSparse::~MatrixSparse()
{
    deinstantiateMatrix();
}

void MatrixSparse::instantiateMatrix()
{
    this->list = new std::vector<std::map<int, double>>(noOfLines);
}

void MatrixSparse::deinstantiateMatrix()
{
    delete list;
}

void MatrixSparse::setElementAt(int line, int column, double value)
{
    if (checkBounds(line, column)) {
        try {
            double d =(*list).at(line).at(column);
        }
        catch (std::out_of_range exception) {
            (*list)[line][column] = value;
        }
    }
    else {
        //maybe throw an exception
//        printf("setElementAt: (line, column) pair not in bounds: (%d, %d) != (%d, %d)\n", line, column, this->noOfLines, this->noOfColumns);
    }
}

double MatrixSparse::getElementAt(int line, int column)
{
    if (checkBounds(line, column)) {
        double d;
        try {
            d = (*list)[line][column];
        }
        catch(std::out_of_range ex) {
            printf("exception boss!\n");
            d = 0.0;
        }
        return d;
    }
    else {
        //maybe throw an exception
//        printf("getElementAt: (line, column) pair not in bounds: (%d, %d) != (%d, %d)\n", line, column, this->noOfLines, this->noOfColumns);
        return 0.0;
    }
}

Matrix *MatrixSparse::getLine(int line) {
    return nullptr;
}

Matrix *MatrixSparse::getColumn(int column) {
    return nullptr;
}

void MatrixSparse::getFromFile(std::string filePath)
{
    //TODO: read from the file the matrix size and values
    std::ifstream inputFile(filePath);

    //read size from the first 2 lines
    int size;

    inputFile >> size;
    this->noOfLines = size;
    this->noOfColumns = size;

    //instantiate the matrix
    instantiateMatrix();

    //read data from the file and build the matrix
    double value;
    int line, column;
    while (inputFile >> value >> line >> column) {
        printf("l: %d c: %d v: %f\n", line, column, value);
        setElementAt(line, column, value);
    }

    inputFile.close();
}

void MatrixSparse::generateRandomMatrixValues(double min, double max)
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

Matrix* MatrixSparse::identityMatrix(int n) {
    return nullptr;
}

Matrix *MatrixSparse::transpose()
{
    return transpose(this);
}

Matrix *MatrixSparse::transpose(Matrix *matrix)
{
    return nullptr;
}

Matrix *MatrixSparse::inverse()
{
    return inverse(this);
}

Matrix *MatrixSparse::inverse(Matrix *matrix)
{
    return nullptr;
}

Matrix *MatrixSparse::add(Matrix *matrix)
{
    return add(this, matrix);
}

Matrix *MatrixSparse::add(Matrix *matrix1, Matrix *matrix2)
{
    return nullptr;
}

Matrix *MatrixSparse::subtract(Matrix *matrix2) {
    return subtract(this, matrix2);
}

Matrix *MatrixSparse::subtract(Matrix *matrix1, Matrix *matrix2) {
    return nullptr;
}

Matrix *MatrixSparse::multiply(Matrix *matrix)
{
    return multiply(this, matrix);
}

Matrix *MatrixSparse::multiply(Matrix *matrix1, Matrix *matrix2)
{
    return nullptr;
}

void MatrixSparse::qrDecomposition(Matrix *b, Matrix **Q, Matrix **R)
{
    qrDecomposition(this, b, Q, R);
}

void MatrixSparse::qrDecomposition(Matrix *A, Matrix *b, Matrix **Q, Matrix **R)
{
    *Q = nullptr;
    *R = nullptr;
}

Matrix *MatrixSparse::clone()
{
    return clone(this);
}

Matrix *MatrixSparse::clone(Matrix *M)
{
    return nullptr;
}

std::string MatrixSparse::toString()
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

Matrix *MatrixSparse::inverseSubstitutionMethod(Matrix *b) {
    return inverseSubstitutionMethod(this, b);
}

Matrix *MatrixSparse::inverseSubstitutionMethod(Matrix *A, Matrix *b) {
    return nullptr;
}

bool MatrixSparse::isSuperiorTriangular() {
    return isSuperiorTriangular(this);
}

bool MatrixSparse::isSuperiorTriangular(Matrix *A) {
    return nullptr;
}

double MatrixSparse::superiorTriangularMatrixDeterminant() {
    return superiorTriangularMatrixDeterminant(this);
}

double MatrixSparse::superiorTriangularMatrixDeterminant(Matrix *A) {
    return 0;
}

bool MatrixSparse::gaussEliminationMethod(Matrix *b, Matrix *R, Matrix *B) {
    return gaussEliminationMethod(this, b, R, B);
}

bool MatrixSparse::gaussEliminationMethod(Matrix *A, Matrix *b, Matrix *R, Matrix *B) {
    return false;
}


