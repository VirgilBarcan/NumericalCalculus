# pragma once
# include "SimpleFunction.h"
# include <cmath>

class PowerFunction :
	public SimpleFunction
{
public:
	PowerFunction() {
		this->power = 0;
		this->constant = 1;
	}

	PowerFunction(double power) {
		this->power = power;
		this->constant = 1;
	}

	PowerFunction(double constant, double power) {
		this->power = power;
		this->constant = constant;
	}

	~PowerFunction() {}

	void setPower(double power) {
		this->power = power;
	}

	void setConstant(double constant) {
		this->constant = constant;
	}

	double evaluate(double x) override {
		return constant * pow(x, power);
	}

	Function* operator+(const Function& rhs) {
		CompositeFunction *cf = new CompositeFunction();

		cf->addFunction(this);
		cf->addFunction(const_cast<Function *>(&rhs));

		return cf;
	}

	Function* operator*(const double& rhs) {
		CompositeFunction *cf = new CompositeFunction();

		PowerFunction *pf = new PowerFunction();
		pf->setPower(power); pf->setConstant(constant * rhs);

		cf->addFunction(pf);

		return cf;
	}

private:
	double power;
	double constant;
};