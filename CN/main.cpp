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
}

int main() {

	testNaiveMatrix();

	return 0;
}