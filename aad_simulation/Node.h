#pragma once
#include <vector>
#include <memory>
#include "Utils.h"

class Node2
{
protected:
	std::vector<Node2*> myArguments;
	double myAdjoint = 0.0;
	double myResult;
public:
	 inline virtual ~Node2() {}
	 inline virtual void propogateAdjoint() = 0;
	 inline double& adjoint() {return myAdjoint;}
	 void resetAdjoints();
	 double result();
};



class PlusNode2 : public Node2
{

public:
	PlusNode2(Node2* lhs, Node2* rhs);
	void propogateAdjoint() override;
};


class MinusNode2 : public Node2
{
public:
	MinusNode2(Node2* lhs, Node2* rhs);
	void propogateAdjoint() override;
};


class TimesNode2 : public Node2
{
public:
	TimesNode2(Node2* lhs, Node2* rhs);
	void propogateAdjoint() override;
};

class LogNode2 : public Node2
{
public:
	LogNode2(Node2* arg);
	void propogateAdjoint() override;
};


class ExpNode2 : public Node2
{
public:
	ExpNode2(Node2* arg);
	void propogateAdjoint() override;
};


class SquareRootNode2 : public Node2
{
public:
	SquareRootNode2(Node2* arg);
	 void propogateAdjoint() override;
};


class SquareNode2 : public Node2
{
public:
	SquareNode2(Node2* arg);
	 void propogateAdjoint() override;
};



class MaxNode2 : public Node2
{
public:
	MaxNode2(Node2* lhs, Node2* rhs);
	void propogateAdjoint() override;
};


class NormCDFNode2 : public Node2
{
public:
	NormCDFNode2(Node2* arg);
	void propogateAdjoint() override;
};


class DivisionNode2 : public Node2
{
public:
	DivisionNode2(Node2* lhs, Node2* rhs);
	void propogateAdjoint() override;
};


class Leaf : public Node2
{
public:
	Leaf(double val);
	double getVal();
	void setVal(double val);
	void evaluate();
	inline void propogateAdjoint() override { }
private:
	double myValue;
};

