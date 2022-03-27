#pragma once
#include "Node.h"

class Number
{
	Node* myNode;
public:
	//tape
	static std::vector<std::unique_ptr<Node>> tape;


	inline void resetTape()
	{
		tape.clear();
	}

	inline Number(double val) : myNode(new Leaf(val))
	{
		tape.push_back(std::unique_ptr<Node>(myNode));
	}

	inline Number(Node* node) : myNode(node) {}

	Node* getNode();
	void setVal(double val);
	double getVal();
	double& adjoint();
	void propogateAdjoints();

	



};



Number operator+(Number lhs, Number rhs);
Number operator-(Number lhs, Number rhs);
Number operator*(Number lhs, Number rhs);
Number operator/(Number lhs, Number rhs);


Number log(Number arg);
Number exp(Number arg);
Number sqrt(Number arg);
Number square(Number arg);
Number max(Number lhs, Number rhs);
Number CDF(Number arg);
