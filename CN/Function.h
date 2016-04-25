# pragma once

class Function
{
public:
	Function() {}
	~Function() {}

	/**
	This function returns the value of the function at the given point
	@param x - the point where the value of the function will be evaluated
	@return the value of the function evaluated in the given point
	*/
	virtual double evaluate(double x) = 0;

	/**
	This operator overloading helps us add two functions
	@param rhs - the right hand side part of plus
	@return a new function which is: this + rhs
	*/
	virtual Function* operator+(const Function& rhs) = 0;

	/**
	This operator overloading helps us multiply a function with a constant
	@param rhs - the right hand side part of product
	@return a new function which is: this * rhs
	*/
	virtual Function* operator*(const double& rhs) = 0;
};

