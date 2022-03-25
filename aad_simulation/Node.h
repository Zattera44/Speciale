#pragma once
#include <vector>
#include <memory>


class Node
{
protected:
	std::vector<Node*> myArguments;
	double myAdjoint = 0.0;
	double myResult;
public:
	 inline virtual ~Node() {}
	 inline virtual void propogateAdjoint() = 0;
	 inline double& adjoint() {return myAdjoint;}
	 void resetAdjoints();
	 double result();
};



class PlusNode : public Node
{

public:
	PlusNode(Node* lhs, Node* rhs);
	void propogateAdjoint() override;
};


class MinusNode : public Node
{
public:
	MinusNode(Node* lhs, Node* rhs);
	void propogateAdjoint() override;
};


class TimesNode : public Node
{
public:
	TimesNode(Node* lhs, Node* rhs);
	void propogateAdjoint() override;
};

class LogNode : public Node
{
public:
	LogNode(Node* arg);
	void propogateAdjoint() override;
};


class ExpNode : public Node
{
public:
	ExpNode(Node* arg);
	void propogateAdjoint() override;
};


class SquareRootNode : public Node
{
public:
	SquareRootNode(Node* arg);
	 void propogateAdjoint() override;
};


class SquareNode : public Node
{
public:
	SquareNode(Node* arg);
	 void propogateAdjoint() override;
};


class Leaf : public Node
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

