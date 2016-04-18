#pragma once

# include "Matrix.h"

class MatrixNorm
{
public:
	/*
	The Maximum Column Sum Norm, or the maximum absolute column sum of x
	||x|| = max { |column 1 sum|, |column 2 sum|, ..., |column n sum| }

	@param x - the matrix for which we want to calculate the norm
	@return the maximum column sum norm
	*/
	static double MaximumColumnSumNorm(Matrix *x) {
		double norm = 0.0;

		for (int column = 0; column < x->getNoOfColumns(); ++column) {
			double columnSum = 0.0;

			for (int line = 0; line < x->getNoOfLines(); ++line) {
				columnSum += fabs(x->getElementAt(line, column));
			}

			if (columnSum > norm) {
				norm = columnSum;
			}
		}

		return norm;
	}

	/*
	The Maximum Row Sum Norm, or the maximum absolute row sum of x
	||x|| = max { |row 1 sum|, |row 2 sum|, ..., |row n sum| }

	@param x - the matrix for which we want to calculate the norm
	@return the maximum row sum norm
	*/
	static double MaximumRowSumNorm(Matrix *x) {
		double norm = 0.0;

		for (int line = 0; line < x->getNoOfLines(); ++line) {
			double rowSum = 0.0;

			for (int column = 0; column < x->getNoOfColumns(); ++column) {
				rowSum += fabs(x->getElementAt(line, column));
			}

			if (rowSum > norm) {
				norm = rowSum;
			}
		}

		return norm;
	}

	/*
	The Maximum Row Sum Norm, or the maximum absolute row sum of x
	||x|| = max { |row 1 sum|, |row 2 sum|, ..., |row n sum| }

	@param x - the matrix for which we want to calculate the norm
	@return the maximum row sum norm
	*/
	static double MaximumRowSumNormSparse(MatrixSparse *x) {
		double norm = 0.0;

		for (int line = 0; line < x->getNoOfLines(); ++line) {
			double rowSum = 0.0;

			for (auto element : x->getListElements(line)) {
				rowSum += fabs(element.second);
			}

			if (rowSum > norm) {
				norm = rowSum;
			}
		}

		return norm;
	}

	static double MaximumNorm(Matrix *x) {
		double norm = 0.0, sum;

		//check that we got a vector and not a matrix
		if (x == nullptr) {
			norm = -1;
		}
		else {
			//calculate the norm
			for (int i = 0; i < x->getNoOfLines(); ++i) {
				sum = 0.0;
				for (int j = 0; j < x->getNoOfColumns(); ++j) {
					sum += fabs(x->getElementAt(i, j));					
				}

				if (sum > norm)
					norm = sum;
			}
		}

		return norm;
	}
};

