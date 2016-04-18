//
// Created by virgil on 18.04.2016.
//

# include "Polynomial.h"

template <class T>
Polynomial<T>::Polynomial() {
    this->coefficients = {};
}

template <class T>
Polynomial<T>::Polynomial(std::vector<T> coefficients) {
    this->coefficients = coefficients;
}

template <class T>
Polynomial<T>::~Polynomial() {

}

template <class T>
std::vector<T> Polynomial<T>::getCoefficients() {
    return this->coefficients;
}

template <class T>
void Polynomial<T>::setCoefficients(std::vector<T> coefficients) {
    this->coefficients = coefficients;
}

template <class T>
unsigned long Polynomial<T>::degree() {
    return this->coefficients.size();
}

template <class T>
T Polynomial<T>::value(T x) {
    return evaluate(this->coefficients, x);
}

template <class T>
T Polynomial<T>::evaluate(std::vector<T> coefficients, double x) {
    //TODO
    return nullptr;
}