#pragma once
#include "Node.h"

class Number2
{
	Node2* myNode;
public:
	//tape
	static std::vector<std::unique_ptr<Node2>> tape;


	inline void resetTape()
	{
		tape.clear();
	}

	inline Number2(double val) : myNode(new Leaf(val))
	{
		tape.push_back(std::unique_ptr<Node2>(myNode));
	}

	inline Number2(Node2* node) : myNode(node) {}

	Node2* getNode();
	void setVal(double val);
	double getVal();
	double& adjoint();
	void propogateAdjoints();

	



};



Number2 operator+(Number2 lhs, Number2 rhs);
Number2 operator-(Number2 lhs, Number2 rhs);
Number2 operator*(Number2 lhs, Number2 rhs);
Number2 operator/(Number2 lhs, Number2 rhs);


Number2 log(Number2 arg);
Number2 exp(Number2 arg);
Number2 sqrt(Number2 arg);
Number2 square(Number2 arg);
Number2 max(Number2 lhs, Number2 rhs);
Number2 CDF(Number2 arg);
