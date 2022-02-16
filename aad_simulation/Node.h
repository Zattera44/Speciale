#pragma once
#include <vector>
#include <memory>

class Node
{
protected:
	std::vector<std::shared_ptr<Node>> myArguments;
	bool isVisited = false;
	unsigned myOrder = 0;
	double myAdjoint = 0.0;
	double myResult;
public:
	virtual ~Node() {}
	virtual void evaluate() = 0;
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


	void setOrder(unsigned order)
	{
		myOrder = order;
	}

	unsigned order()
	{
		return myOrder;
	}

	template <class V>
	void postorder(V visitFunc)
	{
		if (!isVisited)
		{
			for (auto argument : myArguments)
				argument->postorder(visitFunc);
			visitFunc(*this);
			isVisited = true;
	}
	}

	double result()
	{
		return myResult;
	}


	void resetProcessed()
	{
		for (auto argument : myArguments) argument->resetProcessed();
		isVisited = false;
	}

};





class PlusNode : public Node
{

public:
	PlusNode(std::shared_ptr<Node> lhs, std::shared_ptr<Node> rhs)
	{
		myArguments.resize(2);
		myArguments[0] = lhs;
		myArguments[1] = rhs;
	}


	void evaluate() override
	{
		myResult = myArguments[0]->result() + myArguments[1]->result();
	}


	void propogateAdjoint() override
	{
		std::cout << "Propagating node" << myOrder
			<< "adjoint = " << myAdjoint << std::endl;
		myArguments[0]->adjoint() += myAdjoint;
		myArguments[1]->adjoint() += myAdjoint;
	}


};

class TimesNode : public Node
{
public:

	TimesNode(std::shared_ptr<Node> lhs, std::shared_ptr<Node> rhs)
	{
		myArguments.resize(2);
		myArguments[0] = lhs;
		myArguments[1] = rhs;
	}


	void evaluate() override
	{
		myResult = myArguments[0]->result() * myArguments[1]->result();
	}

	void propogateAdjoint() override
	{
		std::cout << "Propagating node " << myOrder
			<< " adjoint = " << myAdjoint << std::endl;
		myArguments[0]->adjoint() += myAdjoint * myArguments[1]->result();
		myArguments[1]->adjoint() += myAdjoint * myArguments[0]->result();
	}

};

class LogNode : public Node
{
public:
	LogNode(std::shared_ptr<Node> arg)
	{
		myArguments.resize(1);
		myArguments[0] = arg;
	}

	void evaluate() override
	{
		myResult = log(myArguments[0]->result());
	}

	void propogateAdjoint() override
	{
		std::cout << "Propagating node " << myOrder
			<< " adjoint = " << myAdjoint << std::endl;
		myArguments[0]->adjoint() += myAdjoint / myArguments[0]->result();
	}

};


class Leaf : public Node
{
public:
	Leaf(double val) : myValue(val) {}

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


	void propogateAdjoint() override
	{
		std::cout << "Accumulating leaf" << myOrder
			<< " adjoint = " << myAdjoint << std::endl;
	}


private:
	double myValue;
};


class Number
{

	std::shared_ptr<Node> myNode;

public:

	Number(double val) : myNode(new Leaf(val)) {}
	Number(std::shared_ptr<Node> node) : myNode(node) {}

	std::shared_ptr<Node> getNode()
	{
		return myNode;
	}

	void setVal(double val)
	{
		std::dynamic_pointer_cast<Leaf>(myNode)->setVal(val);
	}
	
	double getVal()
	{
		return std::dynamic_pointer_cast<Leaf>(myNode)->getVal();
	}

	double evaluate()
	{
		myNode->resetProcessed();
		myNode->postorder([] (Node& n) { n.evaluate() ; } );
		return myNode->result();
	}


	void setOrder()
	{
		myNode->resetProcessed();
		unsigned order = 0;
		myNode->postorder(
			[&order](Node& n) mutable { n.setOrder(++order);
			 }
		);
	}

	void logResults()
	{
		myNode->resetProcessed();
		myNode->postorder([](Node& n) {
			std::cout << "Processed node "
				<< n.order() << " result = "
				<< n.result() << std::endl;
			; }

		);
	}

	double& adjoint()
	{
		return myNode->adjoint();
	}

	void propogateAdjoints()
	{
		myNode->resetAdjoints();
		myNode->adjoint() = 1.0;
	}



};

std::shared_ptr<Node> operator+(Number lhs, Number rhs)
{
	return std::shared_ptr<Node>(new PlusNode(lhs.getNode(), rhs.getNode()));
}


std::shared_ptr<Node> operator*(Number lhs, Number rhs)
{
	return std::shared_ptr<Node>(new TimesNode(lhs.getNode(), rhs.getNode()));
}

std::shared_ptr<Node> log(Number arg)
{
	return std::shared_ptr<Node>(new LogNode(arg.getNode()));
}

