#include "Node.h"



void Node2:: resetAdjoints()
{
	for (auto argument : myArguments)
	{
		argument->resetAdjoints();
	}
	myAdjoint = 0.0;
}



double Node2:: result()
{
	return myResult;
}




//PlusNode2

PlusNode2 :: PlusNode2(Node2* lhs, Node2* rhs)
{
	myArguments.resize(2);
	myArguments[0] = lhs;
	myArguments[1] = rhs;

	myResult = lhs->result() + rhs->result();
}


void  PlusNode2::  propogateAdjoint()
{
	myArguments[0]->adjoint() += myAdjoint;
	myArguments[1]->adjoint() += myAdjoint;
}



//MinusNode2


MinusNode2::MinusNode2(Node2* lhs, Node2* rhs)
{
	myArguments.resize(2);
	myArguments[0] = lhs;
	myArguments[1] = rhs;

	myResult = lhs->result() - rhs->result();
}


void MinusNode2:: propogateAdjoint()
{
	myArguments[0]->adjoint() += myAdjoint;
	myArguments[1]->adjoint() -= myAdjoint;
}


//TimesNode2

TimesNode2:: TimesNode2(Node2* lhs, Node2* rhs)
{
	myArguments.resize(2);
	myArguments[0] = lhs;
	myArguments[1] = rhs;

	myResult = lhs->result() * rhs->result();
}


void TimesNode2:: propogateAdjoint()
{
	myArguments[0]->adjoint() += myAdjoint * myArguments[1]->result();
	myArguments[1]->adjoint() += myAdjoint * myArguments[0]->result();
}

//Log Node2

LogNode2:: LogNode2(Node2* arg)
{
	myArguments.resize(1);
	myArguments[0] = arg;

	myResult = log(arg->result());
}


void LogNode2:: propogateAdjoint()
{
	myArguments[0]->adjoint() += myAdjoint / myArguments[0]->result();
}

//Exp Node2

ExpNode2::ExpNode2(Node2* arg)
{
	myArguments.resize(1);
	myArguments[0] = arg;

	myResult = exp(arg->result());
}

void ExpNode2::propogateAdjoint()
{
	myArguments[0]->adjoint() += myAdjoint * exp(myArguments[0]->result());
}

//SquareRoot Node2

SquareRootNode2:: SquareRootNode2(Node2* arg)
{
	myArguments.resize(1);
	myArguments[0] = arg;

	myResult = sqrt(arg->result());
}

void SquareRootNode2:: propogateAdjoint()
{
	myArguments[0]->adjoint() += myAdjoint / (2 * sqrt(myArguments[0]->result()));
}

//Square Node2

SquareNode2:: SquareNode2(Node2* arg)
{
	myArguments.resize(1);
	myArguments[0] = arg;

	myResult = pow(arg->result(), 2);
}


void SquareNode2:: propogateAdjoint()
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



//Max Node2

MaxNode2::MaxNode2(Node2* lhs, Node2* rhs)
{
	myArguments.resize(2);
	myArguments[0] = lhs;
	myArguments[1] = rhs;

	myResult = (lhs->result() >= rhs->result()) ? lhs->result() : rhs->result();
}

void MaxNode2::propogateAdjoint()
{
	myArguments[0]->adjoint() += myAdjoint * (myArguments[0]->result() >= myArguments[1]->result()) ? 1 : 0;
	myArguments[1]->adjoint() += myAdjoint *  (myArguments[0]->result() < myArguments[1]->result()) ? 1 : 0;
}



NormCDFNode2 :: NormCDFNode2(Node2* arg)
{
	myArguments.resize(1);
	myArguments[0] = arg;

	myResult = normCDF(arg->result());
}




void NormCDFNode2:: propogateAdjoint()
{
	myArguments[0]->adjoint() += myAdjoint * normPDF(myArguments[0]->result());
}



DivisionNode2 :: DivisionNode2(Node2* lhs, Node2* rhs)
{
	myArguments.resize(2);

	myArguments[0] = lhs;
	myArguments[1] = rhs;

	myResult = lhs->result() / rhs->result();
}



void DivisionNode2::  propogateAdjoint()
{
	myArguments[0]->adjoint() += myAdjoint / myArguments[1]->result();
	myArguments[1]->adjoint() += -1 * myAdjoint * myArguments[0]->result() / pow(myArguments[1]->result(),2) ; 
}
