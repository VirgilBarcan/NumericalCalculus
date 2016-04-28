# pragma once
# include "SimpleFunction.h"
# include "CompositeFunction.h"
# include <cmath>

class ExponentialFunction :
	public SimpleFunction
{
public:
	ExponentialFunction() {
		this->base = 0;
		this->constant = 1;
	}

	ExponentialFunction(double base) {
		this->base = base;
		this->constant = 1;
	}

	ExponentialFunction(double constant, double base) {
		this->base = base;
		this->constant = constant;
	}

	~ExponentialFunction() {}

	void setBase(double base) {
		this->base = base;
	}

	void setConstant(double constant) {
		this->constant = constant;
	}

	double evaluate(double x) override {
		return constant * pow(base, x);
	}

	Function* operator+(const Function& rhs) {
		CompositeFunction *cf = new CompositeFunction();

		cf->addFunction(this);
		cf->addFunction(const_cast<Function *>(&rhs));

		return cf;
	}

	Function* operator*(const double& rhs) {
		CompositeFunction *cf = new CompositeFunction();

		ExponentialFunction *ef = new ExponentialFunction();
		ef->setBase(base); ef->setConstant(constant * rhs);

		cf->addFunction(ef);

		return cf;
	}

private:
	double base;
	double constant;
};