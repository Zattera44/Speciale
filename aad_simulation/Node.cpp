#include "Node.h"



void Node:: resetAdjoints()
{
	for (auto argument : myArguments)
	{
		argument->resetAdjoints();
	}
	myAdjoint = 0.0;
}



double Node:: result()
{
	return myResult;
}




//PlusNode

PlusNode:: PlusNode(Node* lhs, Node* rhs)
{
	myArguments.resize(2);
	myArguments[0] = lhs;
	myArguments[1] = rhs;

	myResult = lhs->result() + rhs->result();
}


void  PlusNode::  propogateAdjoint()
{
	myArguments[0]->adjoint() += myAdjoint;
	myArguments[1]->adjoint() += myAdjoint;
}



//MinusNode


MinusNode::MinusNode(Node* lhs, Node* rhs)
{
	myArguments.resize(2);
	myArguments[0] = lhs;
	myArguments[1] = rhs;

	myResult = lhs->result() - rhs->result();
}


void MinusNode:: propogateAdjoint()
{
	myArguments[0]->adjoint() += myAdjoint;
	myArguments[1]->adjoint() -= myAdjoint;
}


//TimesNode

TimesNode:: TimesNode(Node* lhs, Node* rhs)
{
	myArguments.resize(2);
	myArguments[0] = lhs;
	myArguments[1] = rhs;

	myResult = lhs->result() * rhs->result();
}


void TimesNode:: propogateAdjoint()
{
	myArguments[0]->adjoint() += myAdjoint * myArguments[1]->result();
	myArguments[1]->adjoint() += myAdjoint * myArguments[0]->result();
}

//Log Node

LogNode:: LogNode(Node* arg)
{
	myArguments.resize(1);
	myArguments[0] = arg;

	myResult = log(arg->result());
}


void LogNode:: propogateAdjoint()
{
	myArguments[0]->adjoint() += myAdjoint / myArguments[0]->result();
}

//Exp Node

ExpNode::ExpNode(Node* arg)
{
	myArguments.resize(1);
	myArguments[0] = arg;

	myResult = exp(arg->result());
}

void ExpNode::propogateAdjoint()
{
	myArguments[0]->adjoint() += myAdjoint * exp(myArguments[0]->result());
}

//SquareRoot Node

SquareRootNode:: SquareRootNode(Node* arg)
{
	myArguments.resize(1);
	myArguments[0] = arg;

	myResult = sqrt(arg->result());
}

void SquareRootNode:: propogateAdjoint()
{
	myArguments[0]->adjoint() += myAdjoint / (2 * sqrt(myArguments[0]->result()));
}

//Square Node

SquareNode:: SquareNode(Node* arg)
{
	myArguments.resize(1);
	myArguments[0] = arg;

	myResult = pow(arg->result(), 2);
}


void SquareNode:: propogateAdjoint()
{
	myArguments[0]->adjoint() += myAdjoint * 2 * myArguments[0]->result();
}


Leaf::Leaf(double val)
{
	myResult = val;
}


double Leaf:: getVal()
{
	return myValue;
}



void Leaf:: setVal(double val)
{
	myValue = val;
}


void Leaf:: evaluate()
{
	myResult = myValue;
}



//Max Node

MaxNode::MaxNode(Node* lhs, Node* rhs)
{
	myArguments.resize(2);
	myArguments[0] = lhs;
	myArguments[1] = rhs;

	myResult = (lhs->result() >= rhs->result()) ? lhs->result() : rhs->result();
}

void MaxNode::propogateAdjoint()
{
	myArguments[0]->adjoint() += myAdjoint * (myArguments[0]->result() >= myArguments[1]->result()) ? 1 : 0;
	myArguments[1]->adjoint() += myAdjoint *  (myArguments[0]->result() < myArguments[1]->result()) ? 1 : 0;
}



NormCDFNode :: NormCDFNode(Node* arg)
{
	myArguments.resize(1);
	myArguments[0] = arg;

	myResult = normCDF(arg->result());
}




void NormCDFNode:: propogateAdjoint()
{
	myArguments[0]->adjoint() += myAdjoint * normPDF(myArguments[0]->result());
}



DivisionNode :: DivisionNode(Node* lhs, Node* rhs)
{
	myArguments.resize(2);

	myArguments[0] = lhs;
	myArguments[1] = rhs;

	myResult = lhs->result() / rhs->result();
}



void DivisionNode::  propogateAdjoint()
{
	myArguments[0]->adjoint() += myAdjoint / myArguments[1]->result();
	myArguments[1]->adjoint() += -1 * myAdjoint * myArguments[0]->result() / pow(myArguments[1]->result(),2) ; 
}
