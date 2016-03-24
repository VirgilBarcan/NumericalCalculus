#pragma once

# include "Matrix.h"

class VectorialNorm
{
public:
	/*
	The Euclidean Norm, or the distance from the origin to the point x in a n dimensional space
	||x|| = sqrt(x1^2 + x2^2 + ... + xn^2)

	@param x - the column vector for which we want to calculate the norm
	@return the Euclidean norm or -1 if x is not a vector
	*/
	static double EuclideanNorm(Matrix *x) {
		double norm = 0.0;

		//check that we got a vector and not a matrix
		if (x != nullptr && x->getNoOfColumns() == 1 && x->getNoOfLines() != 1) {
			norm = -1;
		}
		else {
			//calculate the norm
			for (int i = 0; i < x->getNoOfColumns(); ++i) {
				norm += pow(x->getElementAt(i, 0), 2);
			}

			norm = sqrt(norm);
		}

		return norm;
	}

	/*
	The Taxicab or Manhattan Norm, or the distance a taxi has to drive in a rectangular street grid to get from the origin to the point x
	||x|| = |x1| + |x2| + ... + |xn|
	This is also called the l1 distance

	@param x - the column vector for which we want to calculate the norm
	@return the Manhattan norm or -1 if x is not a vector
	*/
	static double ManhattanNorm(Matrix *x) {
		double norm = 0.0;

		//check that we got a vector and not a matrix
		if (x != nullptr && x->getNoOfColumns() != 1) {
			norm = -1;
		}
		else {
			//calculate the norm
			for (int i = 0; i < x->getNoOfColumns(); ++i) {
				norm += fabs(x->getElementAt(i, 0));
			}
		}

		return norm;
	}

	/*
	The Maximum Norm, or the maximum absolute value of all components of x
	||x|| = max { |x1|, |x2|, ..., |xn| }
	This is also called the infinity distance

	@param x - the column vector for which we want to calculate the norm
	@return the maximum norm or -1 if x is not a vector
	*/
	static double MaximumNorm(Matrix *x) {
		double norm = 0.0;

		//check that we got a vector and not a matrix
		if (x != nullptr && x->getNoOfColumns() != 1) {
			norm = -1;
		}
		else {
			//calculate the norm
			for (int i = 0; i < x->getNoOfColumns(); ++i) {
				if (fabs(x->getElementAt(0, i)) > norm)
					norm = fabs(x->getElementAt(i, 0));
			}
		}

		return norm;
	}

	/*
	The p-Norm generalizes the other norms: p=1 => taxicab norm, p=2 => Euclidean norm, p->infinity => maximum norm
	||x|| = (|x1|^p + |x2|^p + ... + |xn|^p)^(1/p)

	@param x - the column vector for which we want to calculate the norm
	@param p - the power
	@return the p norm or -1 if x is not a vector
	*/
	static double pNorm(Matrix *x, double p) {
		double norm = 0.0;

		//check that we got a vector and not a matrix
		if (x != nullptr && x->getNoOfColumns() != 1) {
			norm = -1;
		}
		else {
			//calculate the norm
			for (int i = 0; i < x->getNoOfColumns(); ++i) {
				norm += pow(fabs(x->getElementAt(i, 0)), p);
			}

			norm = pow(norm, 1 / p);
		}

		return norm;
	}
};

