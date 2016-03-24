//
// Created by virgil on 21.03.2016.
//

# include "MatrixSparse.h"

MatrixSparse::MatrixSparse()
{
    this->storeType = LINE; //the default store type is LINE
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
    if (storeType == LINE) {
        this->list = new std::vector<std::map<int, double>>(noOfLines);
    }
    else {
        this->list = new std::vector<std::map<int, double>>(noOfColumns);
    }
}

void MatrixSparse::deinstantiateMatrix()
{
    delete list;
}

void MatrixSparse::setElementAt(int line, int column, double value)
{
    if (checkBounds(line, column)) {
        if (storeType == LINE) {
            try {
                double d = (*list).at(line).at(column);
                (*list)[line][column] = value;
            }
            catch (std::out_of_range ex) {
//                printf("BOOM! Exception! L\n");
                (*list)[line][column] = value;
            }
        }
        else {
            try {
                double d = (*list).at(column).at(line);
                (*list)[column][line] = value;
            }
            catch (std::out_of_range ex) {
//                printf("BOOM! Exception! C: %d %d\n", column, line);
                (*list)[column][line] = value;
            }
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
        if (storeType == LINE) {
            try {
                d = (*list)[line][column];
            }
            catch(std::out_of_range ex) {
//                printf("exception boss! L\n");
                d = 0.0;
            }
        }
        else {
            try {
                d = (*list)[column][line];
            }
            catch(std::out_of_range ex) {
//                printf("exception boss! C\n");
                d = 0.0;
            }
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

    if (storeType == LINE) {
        this->noOfLines = size;
        this->noOfColumns = size;

        //instantiate the matrix
        instantiateMatrix();

        //read data from the file and build the matrix
        double value;
        int line, column;
        while (inputFile >> value >> line >> column) {
//            printf("l: %d c: %d v: %f\n", line, column, value);
            setElementAt(line, column, value);
        }
    }
    else {
        this->noOfLines = size;
        this->noOfColumns = 1;

        instantiateMatrix();

        //read data from the file and build the matrix
        double value;
        int line = 0;
        while (inputFile >> value) {
//            printf("l: %d c: %d v: %f\n", line, column, value);
            setElementAt(line++, 0, value);
        }
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
    MatrixSparse *I = new MatrixSparse(n, n);

    for (int i = 0; i < n; ++i) {
        I->setElementAt(i, i, 1);
    }

    return I;
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
    if (checkEqualSizes(matrix1, matrix2)) {
        MatrixSparse *sum = new MatrixSparse(matrix1->getNoOfLines(), matrix1->getNoOfLines());

        for (int line = 0; line < matrix1->getNoOfLines(); ++line) {
            //go through the elements of the line
            for (auto p : ((MatrixSparse *) matrix1)->getListElements(line)) { //fancy auto :) (std::pair<int, double>)
                sum->setElementAt(line, p.first,
                                  matrix1->getElementAt(line, p.first) + matrix2->getElementAt(line, p.first));
            }
        }

        for (int line = 0; line < matrix1->getNoOfLines(); ++line) {
            //go through the elements of the line
            for (auto p : ((MatrixSparse *) matrix2)->getListElements(line)) { //auto hits again :) (std::pair<int, double>)
                sum->setElementAt(line, p.first,
                                  matrix1->getElementAt(line, p.first) + matrix2->getElementAt(line, p.first));
            }
        }

        return sum;
    }
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
    if (checkMultipliable(matrix1, matrix2)) {
        MatrixSparse *productMatrix = new MatrixSparse(matrix1->getNoOfLines(), matrix2->getNoOfColumns());
        if (matrix2->getNoOfColumns() == 1) {
            productMatrix->setStoreType(COLUMN);
        }

        int line1, column1, column2;
        double sum = 0;

        if (((MatrixSparse*)matrix1)->getStoreType() == LINE && ((MatrixSparse*)matrix2)->getStoreType() == COLUMN) {
            for (line1 = 0; line1 < matrix1->getNoOfLines(); ++line1) {
//                printf("line1: %d\n", line1);
                for (column2 = 0; column2 < matrix2->getNoOfColumns(); ++column2) {
                    for (auto p : ((MatrixSparse *) matrix2)->getListElements(column2)) {
                        double product = matrix1->getElementAt(line1, p.first) * matrix2->getElementAt(p.first, column2);
                        sum += product;
                    }
                    productMatrix->setElementAt(line1, column2, sum);
                    sum = 0.0;
                }
            }
        }

        return productMatrix;
    }
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


std::map<int, double> MatrixSparse::getListElements(int lineOrColumn) {
    return (*list).at(lineOrColumn);
}

void MatrixSparse::setStoreType(StoreType storeType) {
    this->storeType = storeType;
}

StoreType MatrixSparse::getStoreType() {
    return this->storeType;
}

bool MatrixSparse::equals(Matrix *M) {
    return equals(this, M);
}

bool MatrixSparse::equals(Matrix *A, Matrix *M) {
    //compare the sizes of the matrices
    if (A->getNoOfLines() != M->getNoOfLines() || A->getNoOfLines() != M->getNoOfLines())
        return false;

    if (((MatrixSparse*)A)->getStoreType() == LINE && ((MatrixSparse*)M)->getStoreType() == LINE) {
        //compare each list elements
        for (int line = 0; line < A->getNoOfLines(); ++line) {
            for (auto p : ((MatrixSparse*)A)->getListElements(line)) {
                if (fabs(A->getElementAt(line, p.first) - M->getElementAt(line, p.first)) > epsilon) {
                    printf("%d %d %.32f %.32f %.32f\n", line, p.first,
                           A->getElementAt(line, p.first), M->getElementAt(line, p.first),
                           A->getElementAt(line, p.first) - M->getElementAt(line, p.first));
                    return false;
                }
            }
        }
    }
    if (((MatrixSparse*)A)->getStoreType() == COLUMN && ((MatrixSparse*)M)->getStoreType() == COLUMN) {
        //compare each list elements
        for (int column = 0; column < A->getNoOfColumns(); ++column) {
            for (auto p : ((MatrixSparse*)A)->getListElements(column)) {
                if (fabs(A->getElementAt(p.first, column) - M->getElementAt(p.first, column)) > epsilon) {
                    printf("%d %d %.32f %.32f %.32f\n", column, p.first,
                           A->getElementAt(p.first, column), M->getElementAt(p.first, column),
                           A->getElementAt(p.first, column) - M->getElementAt(p.first, column));
                    return false;
                }
            }
        }
    }

    return true;
}

