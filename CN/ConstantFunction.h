# pragma once
# include "SimpleFunction.h"
# include "CompositeFunction.h"

class ConstantFunction :
	public SimpleFunction
{
public:
	ConstantFunction() {
		this->constant = 0;
	}
	
	ConstantFunction(double value) {
		this->constant = value;
	}

	~ConstantFunction() {}

	void setValue(double value) {
		this->constant = value;
	}

	double evaluate(double x) override {
		return constant;
	}

	Function* operator+(const Function& rhs) {
		CompositeFunction *cf = new CompositeFunction();

		cf->addFunction(this);
		cf->addFunction(const_cast<Function *>(&rhs));

		return cf;
	}

	Function* operator*(const double& rhs) {
		CompositeFunction *cf = new CompositeFunction();

		cf->addFunction(new ConstantFunction(constant * rhs));

		return cf;
	}

private:
	double constant;
};