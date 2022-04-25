#include "Number.h"


Node2*  Number2::  getNode()
{
	return myNode;
}


void Number2:: setVal(double val)
{
	dynamic_cast<Leaf*>(myNode)->setVal(val);
}



double Number2:: getVal()
{
	return dynamic_cast<Leaf*>(myNode)->getVal();
}



double& Number2:: adjoint()
{
	return myNode->adjoint();
}


void Number2:: propogateAdjoints()
{
	myNode->resetAdjoints();
	myNode->adjoint() = 1.0;

	auto it = tape.rbegin();
	while (it->get() != myNode)
	{
		++it;
	}

	while (it != tape.rend())
	{
		(*it)->propogateAdjoint();
		++it;
	}
}


//Operator overloading



Number2 operator+(Number2 lhs, Number2 rhs)
{
	Node2* n = new PlusNode2(lhs.getNode(), rhs.getNode());

	Number2::tape.push_back(std::unique_ptr<Node2>(n));

	return n;
}


Number2 operator-(Number2 lhs, Number2 rhs)
{
	Node2* n = new MinusNode2(lhs.getNode(), rhs.getNode());

	Number2::tape.push_back(std::unique_ptr<Node2>(n));

	return n;
}

Number2 operator*(Number2 lhs, Number2 rhs)
{

	Node2* n = new TimesNode2(lhs.getNode(), rhs.getNode());

	Number2::tape.push_back(std::unique_ptr<Node2>(n));

	return n;
}

Number2 log(Number2 arg)
{
	Node2* n = new LogNode2(arg.getNode());

	Number2::tape.push_back(std::unique_ptr<Node2>(n));

	return n;
}

Number2 exp(Number2 arg)
{
	Node2* n = new ExpNode2(arg.getNode());

	Number2::tape.push_back(std::unique_ptr<Node2>(n));

	return n;
}

Number2 sqrt(Number2 arg)
{
	Node2* n = new SquareRootNode2(arg.getNode());

	Number2::tape.push_back(std::unique_ptr<Node2>(n));

	return n;
}


Number2 square(Number2 arg)
{
	Node2* n = new SquareNode2(arg.getNode());

	Number2::tape.push_back(std::unique_ptr<Node2>(n));

	return n;
}


Number2 max(Number2 lhs, Number2 rhs)
{
	Node2* n = new MaxNode2(lhs.getNode(), rhs.getNode());

	Number2::tape.push_back(std::unique_ptr<Node2>(n));

	return n;
}


Number2 CDF(Number2 arg)
{
	Node2* n = new NormCDFNode2(arg.getNode());

	Number2::tape.push_back(std::unique_ptr<Node2>(n));

	return n;
}



Number2 operator/(Number2 lhs, Number2 rhs)
{
	Node2* n = new DivisionNode2(lhs.getNode(), rhs.getNode());

	Number2::tape.push_back(std::unique_ptr<Node2>(n));

	return n;
}


std::vector<std::unique_ptr<Node2>> Number2::tape;
