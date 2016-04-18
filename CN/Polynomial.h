//
// Created by virgil on 18.04.2016.
//

# pragma once

# include <vector>

template <class T>
class Polynomial {
public:
    Polynomial();
    Polynomial(std::vector<T> coefficients);
    ~Polynomial();

    /**
     * Get the coefficients of the polynomial
     *
     * @return the coefficients of the polynomial
     */
    std::vector<T> getCoefficients();

    /**
     * Set the coefficients of the polynomial
     *
     * @param coefficients - the coefficients of the polynomial
     */
    void setCoefficients(std::vector<T> coefficients);

    /**
     * Get the degree of the polynomial
     *
     * @return the degree of the polynomial
     */
    unsigned long degree();

    /**
     * Compute the value of the polynomial at the given point
     *
     * @param x - the point where the polynomial will be evaluated
     * @return the value of the polynomial at the given point
     */
    T value(T x);

private:
    std::vector<T> coefficients;

    /**
     * Calculate the value of the polynomyal with the given coefficients at the given point using Horner's Method
     *
     * @param coefficients - the coefficients of the polynomyal to evaluate
     * @param x - the point where the polynomial will be evaluated
     * @return the value of the polynomial given by its coefficients at the wanted point
     */
    static T evaluate(std::vector<T> coefficients, double x);
};