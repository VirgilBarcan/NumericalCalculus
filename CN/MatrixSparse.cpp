//
// Created by virgil on 21.03.2016.
//

#include "MatrixSparse.h"
#include "MatrixNorm.h"

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
                //double d = (*list).at(line).at(column);
                //(*list)[line][column] = value;

                (*list).at(line).at(column) = value;
            }
            catch (std::out_of_range ex) {
//                printf("BOOM! Exception! L\n");
                (*list)[line][column] = value;
            }
        }
        else {
            try {
                //double d = (*list).at(column).at(line);
                //(*list)[column][line] = value;

                (*list).at(column).at(line) = value;
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
        double d = 0.0;
        if (storeType == LINE) {
            //try {
            //    d = (*list).at(line).at(column);
            //}
            //catch(std::out_of_range ex) {
//                printf("exception boss! L\n");
            //    d = 0.0;
            //}

            auto it = (*list)[line].find(column);
            if (it != (*list)[line].end()) {
                d = it->second;
            }
        }
        else {
            //try {
            //    d = (*list).at(column).at(line);
            //}
            //catch(std::out_of_range ex) {
//                printf("exception boss! C\n");
            //    d = 0.0;
            //}

            auto it = (*list)[column].find(line);
            if (it != (*list)[column].end()) {
                d = it->second;
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

    //read size from the first line
    int size;

	printf("%d\n", inputFile.is_open());

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
    if (checkEqualSizes(matrix1, matrix2)) {
        if (((MatrixSparse *) matrix1)->getStoreType() == LINE && ((MatrixSparse *) matrix1)->getStoreType() == LINE) {
            MatrixSparse *difference = new MatrixSparse(matrix1->getNoOfLines(), matrix1->getNoOfColumns());
            difference->setStoreType(LINE);

            for (int line = 0; line < matrix1->getNoOfLines(); ++line) {
                //go through the elements of the line
                for (auto p : ((MatrixSparse *) matrix1)->getListElements(line)) { //fancy auto :) (std::pair<int, double>)
                    difference->setElementAt(line, p.first,
                                             matrix1->getElementAt(line, p.first) - matrix2->getElementAt(line, p.first));
                }
            }

            for (int line = 0; line < matrix1->getNoOfLines(); ++line) {
                //go through the elements of the line
                for (auto p : ((MatrixSparse *) matrix2)->getListElements(line)) { //auto hits again :) (std::pair<int, double>)
                    difference->setElementAt(line, p.first,
                                             matrix1->getElementAt(line, p.first) - matrix2->getElementAt(line, p.first));
                }
            }
            return difference;
        }

        if (((MatrixSparse *) matrix1)->getStoreType() == COLUMN && ((MatrixSparse *) matrix1)->getStoreType() == COLUMN) {
            MatrixSparse *difference = new MatrixSparse(matrix1->getNoOfLines(), matrix1->getNoOfColumns());
            difference->setStoreType(COLUMN);

            for (int column = 0; column < matrix1->getNoOfColumns(); ++column) {
                //go through the elements of the column
                for (auto p : ((MatrixSparse *) matrix1)->getListElements(column)) { //fancy auto :) (std::pair<int, double>)
                    difference->setElementAt(p.first, column,
                                             matrix1->getElementAt(p.first, column) - matrix2->getElementAt(p.first, column));
                }
            }

            for (int column = 0; column < matrix1->getNoOfColumns(); ++column) {
                //go through the elements of the line
                for (auto p : ((MatrixSparse *) matrix2)->getListElements(column)) { //auto hits again :) (std::pair<int, double>)
                    difference->setElementAt(p.first, column,
                                             matrix1->getElementAt(p.first, column) - matrix2->getElementAt(p.first, column));
                }
            }
            return difference;
        }
    }
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
                std::map<int, double> aux = ((MatrixSparse *)matrix1)->getListElements(line1);
                for (column2 = 0; column2 < matrix2->getNoOfColumns(); ++column2) {
                    for (auto p : aux) {
                        //printf("%d %d %d\n", line1, column2, p.first);
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

Matrix *MatrixSparse::strictSuperiorPart() {
    return strictSuperiorPart(this);
}

Matrix *MatrixSparse::strictSuperiorPart(Matrix *A) {
    if (!A->isDiagonalZero()) {
        MatrixSparse *x = new MatrixSparse(A->getNoOfLines(), 1);

        for (int line = 0; line < A->getNoOfLines(); ++line) {
            for (int column = line + 1; column < A->getNoOfColumns(); ++column) {
                x->setElementAt(line, column, A->getElementAt(line, column));
            }
        }

        return x;
    }
    else {
        printf("We can't return the superior triangular elements of a non-square matrix!\n");
        return nullptr;
    }
}

Matrix *MatrixSparse::diagonal() {
    return diagonal(this);
}

Matrix *MatrixSparse::diagonal(Matrix *A) {
    if (!A->isDiagonalZero()) {
        MatrixSparse *x = new MatrixSparse(A->getNoOfLines(), 1);

        for (int i = 0; i < A->getNoOfLines(); ++i) {
            x->setElementAt(i, i, A->getElementAt(i, i));
        }

        return x;
    }
    else {
        printf("We can't return the diagonal elements of a non-square matrix!\n");
        return nullptr;
    }
}

Matrix *MatrixSparse::strictInferiorPart() {
    return strictInferiorPart(this);
}

Matrix *MatrixSparse::strictInferiorPart(Matrix *A) {
    if (!A->isDiagonalZero()) {
        MatrixSparse *x = new MatrixSparse(A->getNoOfLines(), 1);

        for (int line = 0; line < A->getNoOfLines(); ++line) {
            for (int column = 0; column < line; ++column) {
                x->setElementAt(line, column, A->getElementAt(line, column));
            }
        }

        return x;
    }
    else {
        printf("We can't return the inferior triangular elements of a non-square matrix!\n");
        return nullptr;
    }
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
    if (((MatrixSparse *) M)->getStoreType() == LINE) {
        MatrixSparse *C = new MatrixSparse(M->getNoOfLines(), M->getNoOfColumns());
        C->setStoreType(LINE);
        C->setList(((MatrixSparse*) M)->getList());
        return C;
    }
    if (((MatrixSparse *) M)->getStoreType() == COLUMN) {
        MatrixSparse *C = new MatrixSparse(M->getNoOfLines(), M->getNoOfColumns());
        C->setStoreType(COLUMN);
        C->setList(((MatrixSparse*) M)->getList());
        return C;
    }
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

void MatrixSparse::setList(std::vector<std::map<int, double>>* list)
{
    this->list = new std::vector<std::map<int, double>>(*list);
}

std::vector<std::map<int, double>> *MatrixSparse::getList()
{
    return this->list;
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

bool MatrixSparse::isDiagonalZero() {
    return isSuperiorTriangular(this);
}

bool MatrixSparse::isDiagonalZero(Matrix *A) {
    //check the matrix is square
    if (A->getNoOfLines() != A->getNoOfColumns())
        return true; //we return like that, even though an exception will be better

    for (int i = 0; i < A->getNoOfLines(); ++i) {
        if (fabs(A->getElementAt(i, i)) < epsilon) {
            return true;
        }
    }

    return false;
}

Matrix *MatrixSparse::sorMethod(Matrix *b) {
    return sorMethod(this, b);
}

Matrix *MatrixSparse::sorMethod(Matrix *A, Matrix *b) {
    if (!A->isDiagonalZero()) {
        MatrixSparse *xc = new MatrixSparse(A->getNoOfLines(), 1); xc->setStoreType(COLUMN);
        MatrixSparse *xp = new MatrixSparse(A->getNoOfLines(), 1); xp->setStoreType(COLUMN);

        std::map<int, double> aux;

//        for (int i = 0; i < A->getNoOfLines(); ++i) {
//            xp->setElementAt(i, 0, i + 1);
//        }

        double delta_x = 0;
        int k = 0;

        do {
            for (int i = 0; i < A->getNoOfLines(); ++i) {
                aux = (reinterpret_cast<MatrixSparse*>(A)->getListElements(i));

                if (0 != aux.size()) {
                    //calculate sum with j = 1, i-1 of a[i][j] * x(k+1)[j]
                    double sum1 = 0;
                    //calculate sum with j = i + 1, n of a[i][j] * x(k)[j]
                    double sum2 = 0;
                    for (auto element : aux) {
                        if (element.first < i) {
                            sum1 += A->getElementAt(i, element.first) * xc->getElementAt(element.first, 0);
                        }

                        if (element.first > i) {
                            sum2 += A->getElementAt(i, element.first) * xp->getElementAt(element.first, 0);
                        }
                    }

                    xc->setElementAt(i, 0,
                                     (-0.2 * xp->getElementAt(i, 0)) +
                                     1.2 * (b->getElementAt(i, 0) - sum1 - sum2) / (A->getElementAt(i, i)));
                }
            }

            delta_x = MatrixNorm::MaximumRowSumNormSparse(reinterpret_cast<MatrixSparse*>(xc->subtract(xp)));

            printf("k = %d -- delta_x = %.16f\n", k, delta_x);
//            printf("xc = \n%s\n", xc->toString().c_str());
//            printf("xp = \n%s\n", xp->toString().c_str());

            //copy xc in xp
            xp = reinterpret_cast<MatrixSparse*>(xc->clone());

            ++k;

        } while(delta_x >= epsilon && k <= 10000 && delta_x <= pow(10, 8));

        if (delta_x < epsilon) {
            return xc;
        }
        else {
            printf("The algorithm didn't converge on this matrix!\n");
            return nullptr;
        }
    }
    else {
        printf("We can't apply the method on a matrix with 0's on the diagonal!\n");
        return nullptr;
    }
}

bool MatrixSparse::isSymmetric() {
    return isSymmetric(this);
}

bool MatrixSparse::isSymmetric(Matrix *A) {

    //printf("A:\n%s\n", A->toString().c_str());

    for (int line = 0; line < A->getNoOfLines(); ++line) {
        for (auto p : ((MatrixSparse*) A)->getListElements(line)) {
            if (fabs(A->getElementAt(line, p.first) - A->getElementAt(p.first, line)) > epsilon) {
                printf("Different elements at: %d %d\n", line, p.first);
                return false;
            }
        }
    }

    return true;
}
