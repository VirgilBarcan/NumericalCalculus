#pragma once
#include "Function.h"
class SimpleFunction :
	public Function
{
public:
	SimpleFunction() {}
	~SimpleFunction() {}

	virtual double evaluate(double x) override = 0;

	virtual Function* operator+(const Function& rhs) override = 0;

	virtual Function* operator*(const double& rhs) override = 0;
};

