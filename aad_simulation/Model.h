#pragma once

class  MarkovModel
{
public:
	virtual double increment(double spot, double dSpot) const = 0;
	virtual ~MarkovModel() {};
	virtual MarkovModel* clone() const = 0;
};


class BlackScholes : public MarkovModel
{
public:
	double increment(double spot, double dSpot);



private:
	double r;
	double sigma;
};