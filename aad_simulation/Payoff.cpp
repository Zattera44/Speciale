class PayOff
{
public:
	PayOff() {};
	virtual double operator()(double spot) const = 0;
	virtual ~PayOff() {};
private:
};


class EuropeanCall : public PayOff
{
public:
	EuropeanCall(double strike_, double TTM);
	virtual double operator()(double spot) const;
	virtual ~EuropeanCall() {};

private:
	double strike;
	double TTM;
};

class EuropeanPut : public PayOff
{
public:
	EuropeanPut(double strike_, double TTM);
	virtual double operator()(double spot) const;
	virtual ~EuropeanPut() {};

private:
	double strike;
	double TTM;
};