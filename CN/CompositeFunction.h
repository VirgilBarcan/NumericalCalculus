# pragma once
# include "Function.h"
# include <vector>
# include <algorithm>

class CompositeFunction :
	public Function
{
public:
	CompositeFunction() {
		components = new std::vector<Function *>();
	}

	~CompositeFunction() {
		delete components;
	}

	void addFunction(Function *f) {
		components->push_back(f);
	}

	void removeFunction(Function *f) {
		components->erase(std::find(components->begin(), components->end(), f));
	}

	double evaluate(double x) override {
		double value = 0.0;
		
		for (auto it = components->begin(); it != components->end(); ++it) {
			value += (*it)->evaluate(x);
		}

		return value;
	}

	Function* operator+(const Function& rhs) {
		CompositeFunction *cf = new CompositeFunction();

		cf->addFunction(this);
		cf->addFunction(const_cast<Function *>(&rhs));

		return cf;
	}

	Function* operator*(const double& rhs) {
		CompositeFunction *cf = new CompositeFunction();

		for (auto it = components->begin(); it != components->end(); ++it) {
			cf->addFunction((*(*it)) * rhs);
		}

		return cf;
	}

private:
	std::vector<Function *> *components;
};

