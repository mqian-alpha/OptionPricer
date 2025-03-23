#include <cmath>
#include <iostream>
using namespace std;

#define SMALL_SIGMA 1e-6
#define SMALL_T 1e-6

// Normal distribution function
double N(double x) {
    return 0.5 * erfc(-x / sqrt(2));
}

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

OptionContract::OptionContract(double K, double T, bool isCall) {
    this->m_K = K;
    this->m_T = T > SMALL_T ? T : SMALL_T;
    this->m_isCall = isCall;
}


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
// Private member variables: 
// m_S0: initial stock price
// m_r: risk-free interest rate
// m_sigma: volatility
// m_q: dividend yield
// m_contract: option contract
    double m_S0, m_r, m_sigma, m_q;
    OptionContract m_contract;
};

OptionSecurity::OptionSecurity(double S0, double sigma, double r, double q, const OptionContract& contract) {
    this->m_S0 = S0;
    this->m_r = r;
    this->m_sigma = (sigma > SMALL_SIGMA) ? sigma : SMALL_SIGMA;
    this->m_q = q;
    this->m_contract = contract;
}   

// Implement the member functions
// price(): return the option price using the Black-Scholes formula
double OptionSecurity::price() const {
    double K = m_contract.getStrike();
    double T = m_contract.getMaturity();
    if(T <= SMALL_T) {
        return max(0.0, m_contract.isCallOption() ? m_S0 - K : K - m_S0);   
    }   
    bool isCall = m_contract.isCallOption();
    double d1 = (log(m_S0 / K) + (m_r - m_q + 0.5 * m_sigma * m_sigma) * T) / (m_sigma * sqrt(T));
    double d2 = d1 - m_sigma * sqrt(T);
    if (isCall) {
        return m_S0 * exp(-m_q * T) * N(d1) - K * exp(-m_r * T) * N(d2);
    } else {
        return K * exp(-m_r * T) * N(-d2) - m_S0 * exp(-m_q * T) * N(-d1);
    }
}

double OptionSecurity::price(double sigma) const {
    double K = m_contract.getStrike();
    double T = m_contract.getMaturity();
    bool isCall = m_contract.isCallOption();

    if(T <= SMALL_T) {
        return max(0.0, m_contract.isCallOption() ? m_S0 - K : K - m_S0);   
    }   

    double d1 = (log(m_S0 / K) + (m_r - m_q + 0.5 * sigma * sigma) * T) / (sigma * sqrt(T));
    double d2 = d1 - sigma * sqrt(T);
    if (isCall) {
        return m_S0 * exp(-m_q * T) * N(d1) - K * exp(-m_r * T) * N(d2);
    } else {
        return K * exp(-m_r * T) * N(-d2) - m_S0 * exp(-m_q * T) * N(-d1);
    }
}

// delta(): return the option delta defined as the derivative of the option price with respect to the stock price
double OptionSecurity::delta() const {
    double K = m_contract.getStrike();
    double T = m_contract.getMaturity();
    bool isCall = m_contract.isCallOption();
    if(T <= SMALL_T) {
        return isCall ? ((m_S0 > K) ? 1.0 : 0.0) : ((m_S0 < K) ? -1.0 : 0.0);
    }

    double d1 = (log(m_S0 / K) + (m_r - m_q + 0.5 * m_sigma * m_sigma) * T) / (m_sigma * sqrt(T));
    if (isCall) {
        return exp(-m_q * T) * N(d1);
    } else {
        return -exp(-m_q * T) * N(-d1);
    }
}

// gamma(): return the option gamma defined as the second derivative of the option price with respect to the stock price
double OptionSecurity::gamma() const {
    double K = m_contract.getStrike();
    double T = m_contract.getMaturity();
    if(T <= SMALL_T) {
        return 0.0;
    }

    double d1 = (log(m_S0 / K) + (m_r - m_q + 0.5 * m_sigma * m_sigma) * T) / (m_sigma * sqrt(T));
    return exp(-m_q * T) * N(d1) / (m_S0 * m_sigma * sqrt(T));
}

// theta(): return the option theta defined as the derivative of the option price with respect to the time to maturity
double OptionSecurity::theta() const {
    double K = m_contract.getStrike();
    double T = m_contract.getMaturity();
    bool isCall = m_contract.isCallOption();
    if(m_sigma <= SMALL_SIGMA) {
        return 0.0;
    }

    double d1 = (log(m_S0 / K) + (m_r - m_q + 0.5 * m_sigma * m_sigma) * T) / (m_sigma * sqrt(T));
    double d2 = d1 - m_sigma * sqrt(T);
    if (isCall) {
        return -m_S0 * exp(-m_q * T) * m_sigma * N(d1) / (2 * sqrt(T)) - m_r * K * exp(-m_r * T) * N(d2) + m_q * m_S0 * exp(-m_q * T) * N(d1);
    } else {
        return -m_S0 * exp(-m_q * T) * m_sigma * N(-d1) / (2 * sqrt(T)) + m_r * K * exp(-m_r * T) * N(-d2) - m_q * m_S0 * exp(-m_q * T) * N(-d1);
    }
}

// vega(): return the option vega defined as the derivative of the option price with respect to the volatility
double OptionSecurity::vega() const {
    double K = m_contract.getStrike();
    double T = m_contract.getMaturity();
    if(m_sigma <= SMALL_SIGMA) {
        return 0.0;
    }
    double d1 = (log(m_S0 / K) + (m_r - m_q + 0.5 * m_sigma * m_sigma) * T) / (m_sigma * sqrt(T));
    return m_S0 * exp(-m_q * T) * sqrt(T) * N(d1);
}

// rho(): return the option rho defined as the derivative of the option price with respect to the risk-free interest rate
double OptionSecurity::rho() const {
    double K = m_contract.getStrike();
    double T = m_contract.getMaturity();
    bool isCall = m_contract.isCallOption();
    double d2 = (log(m_S0 / K) + (m_r - m_q + 0.5 * m_sigma * m_sigma) * T) / (m_sigma * sqrt(T)) - m_sigma * sqrt(T);
    if (isCall) {
        return K * T * exp(-m_r * T) * N(d2);
    } else {
        return -K * T * exp(-m_r * T) * N(-d2);
    }
}

// impliedVol(): return the implied volatility given the target option price
double OptionSecurity::impliedVol(double targetPrice) const {
    double epsilon = 1e-6;
    double low = 0.0, high = 1.0;
    double mid = 0.5;
    while (high - low > epsilon) {
        double midPrice = price(mid);
        if (midPrice < targetPrice) {
            low = mid;
        } else {
            high = mid;
        }
        mid = 0.5 * (low + high);
    }
    return mid;
}

// Test the OptionPrice class
int main() {
    double S0 = 100.0, K = 100.0, r = 0.05, T = 1.0, sigma = 0.15, q = 0.0;
    bool isCall = true;
    OptionContract contract(K, T, isCall);
    OptionSecurity option(S0, sigma, r, q, contract);
    double price = option.price();
    double delta = option.delta();
    double gamma = option.gamma();
    double theta = option.theta();
    double vega = option.vega();
    double rho = option.rho();
    double impliedVol = option.impliedVol(price);
    cout << "Option price: " << price << endl
            << "Delta: " << delta << endl
            << "Gamma: " << gamma << endl
            << "Theta: " << theta << endl
            << "Vega: " << vega << endl
            << "Rho: " << rho << endl
            << "Implied volatility: " << impliedVol << endl;
    return 0;
}
