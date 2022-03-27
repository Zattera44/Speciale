#include "Number.h"


Node*  Number::  getNode()
{
	return myNode;
}


void Number:: setVal(double val)
{
	dynamic_cast<Leaf*>(myNode)->setVal(val);
}



double Number:: getVal()
{
	return dynamic_cast<Leaf*>(myNode)->getVal();
}



double& Number:: adjoint()
{
	return myNode->adjoint();
}


void Number:: propogateAdjoints()
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



Number operator+(Number lhs, Number rhs)
{
	Node* n = new PlusNode(lhs.getNode(), rhs.getNode());

	Number::tape.push_back(std::unique_ptr<Node>(n));

	return n;
}


Number operator-(Number lhs, Number rhs)
{
	Node* n = new MinusNode(lhs.getNode(), rhs.getNode());

	Number::tape.push_back(std::unique_ptr<Node>(n));

	return n;
}

Number operator*(Number lhs, Number rhs)
{

	Node* n = new TimesNode(lhs.getNode(), rhs.getNode());

	Number::tape.push_back(std::unique_ptr<Node>(n));

	return n;
}

Number log(Number arg)
{
	Node* n = new LogNode(arg.getNode());

	Number::tape.push_back(std::unique_ptr<Node>(n));

	return n;
}

Number exp(Number arg)
{
	Node* n = new ExpNode(arg.getNode());

	Number::tape.push_back(std::unique_ptr<Node>(n));

	return n;
}

Number sqrt(Number arg)
{
	Node* n = new SquareRootNode(arg.getNode());

	Number::tape.push_back(std::unique_ptr<Node>(n));

	return n;
}


Number square(Number arg)
{
	Node* n = new SquareNode(arg.getNode());

	Number::tape.push_back(std::unique_ptr<Node>(n));

	return n;
}


Number max(Number lhs, Number rhs)
{
	Node* n = new MaxNode(lhs.getNode(), rhs.getNode());

	Number::tape.push_back(std::unique_ptr<Node>(n));

	return n;
}


Number CDF(Number arg)
{
	Node* n = new NormCDFNode(arg.getNode());

	Number::tape.push_back(std::unique_ptr<Node>(n));

	return n;
}



Number operator/(Number lhs, Number rhs)
{
	Node* n = new DivisionNode(lhs.getNode(), rhs.getNode());

	Number::tape.push_back(std::unique_ptr<Node>(n));

	return n;
}


std::vector<std::unique_ptr<Node>> Number::tape;
