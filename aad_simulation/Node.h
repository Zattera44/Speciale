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
	virtual ~Node() {}
	virtual void propogateAdjoint() = 0;

	double& adjoint()
	{
		return myAdjoint;
	}


	void resetAdjoints()
	{
		for (auto argument : myArguments)
		{
			argument->resetAdjoints();
		}
		myAdjoint = 0.0;
	}


	double result()
	{
		return myResult;
	}

};





class PlusNode : public Node
{

public:
	PlusNode(Node* lhs, Node* rhs)
	{
		myArguments.resize(2);
		myArguments[0] = lhs;
		myArguments[1] = rhs;

		myResult = lhs->result() + rhs->result();
	}

	void propogateAdjoint() override
	{
		myArguments[0]->adjoint() += myAdjoint;
		myArguments[1]->adjoint() += myAdjoint;
	}


};

class TimesNode : public Node
{
public:

	TimesNode(Node* lhs, Node* rhs)
	{
		myArguments.resize(2);
		myArguments[0] = lhs;
		myArguments[1] = rhs;

		myResult = lhs->result() * rhs->result();
	}


	void propogateAdjoint() override
	{
		myArguments[0]->adjoint() += myAdjoint * myArguments[1]->result();
		myArguments[1]->adjoint() += myAdjoint * myArguments[0]->result();
	}
};

class LogNode : public Node
{
public:
	LogNode(Node* arg)
	{
		myArguments.resize(1);
		myArguments[0] = arg;

		myResult = log(arg->result());
	}

	void propogateAdjoint() override
	{
		myArguments[0]->adjoint() += myAdjoint / myArguments[0]->result();
	}

};


class Leaf : public Node
{
public:
	Leaf(double val) 
	{
		myResult = val;
	}

	double getVal()
	{
		return myValue;
	}

	void setVal(double val)
	{
		myValue = val;
	}


	void evaluate()
	{
		myResult = myValue;
	}


	void propogateAdjoint() override { }


private:
	double myValue;
};


class Number
{

	Node* myNode;

public:
	//tape
	static std::vector<std::unique_ptr<Node>> tape;

	Number(double val) : myNode(new Leaf(val))
	{
		tape.push_back(std::unique_ptr<Node>(myNode));
	}

	Number(Node* node) : myNode(node) {}

	Node* getNode()
	{
		return myNode;
	}

	void setVal(double val)
	{
		dynamic_cast<Leaf*>(myNode)->setVal(val);
	}
	
	double getVal()
	{
		return dynamic_cast<Leaf*>(myNode)->getVal();
	}


	double& adjoint()
	{
		return myNode->adjoint();
	}

	void propogateAdjoints()
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
};


std::vector<std::unique_ptr<Node>> Number::tape;


//Operator overriding

Number operator+(Number lhs, Number rhs)
{
	Node* n = new PlusNode(lhs.getNode(), rhs.getNode());

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

