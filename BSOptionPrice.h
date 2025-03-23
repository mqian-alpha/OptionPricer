#include <cmath>

#define SMALL_SIGMA 1e-6
#define SMALL_T 1e-6

// Normal distribution function
double N(double x) {
    return 0.5 * erfc(-x / sqrt(2));
}

// Option contract class
class OptionContract {
public:
    OptionContract() {
        m_K = 0.0;
        m_T = 0.01;
        m_isCall = true;
    }
    OptionContract(double K, double T, bool isCall);
    double getStrike() const { return m_K; }
    double getMaturity() const { return m_T; }
    bool isCallOption() const { return m_isCall; }
private:
    double m_K; // strike price
    double m_T; // time to maturity
    bool m_isCall; // true for call option, false for put option
};

// Black-Scholes option pricing formula
class OptionSecurity {
    public:
        OptionSecurity(double S0, double sigma, double r, double q, const OptionContract& contract);
        double price() const;
        double price(double sigma) const;
        double delta() const;
        double gamma() const;
        double theta() const;
        double vega() const;
        double rho() const;
        double impliedVol(double targetPrice) const;
    private:
    // Private member variables for market data and option contract: 
    // m_S0: initial stock price
    // m_r: risk-free interest rate
    // m_sigma: volatility
    // m_q: dividend yield
    // m_contract: option contract
        double m_S0, m_r, m_sigma, m_q;
        OptionContract m_contract;
    };
    
    