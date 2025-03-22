#include <cmath>
#include <iostream>
using namespace std;

// Normal distribution function
double N(double x) {
    return 0.5 * erfc(-x / sqrt(2));
}

// Black-Scholes option pricing formula
class OptionPrice {
public:
    OptionPrice(double S0, double K, double r, double T, double sigma, double q, bool isCall);
    double price() const;
    double delta() const;
    double gamma() const;
    double theta() const;
    double vega() const;
    double rho() const;
    double impliedVol(double targetPrice) const;
private:
// Private member variables: 
// S0: initial stock price
// K: strike price
// r: risk-free interest rate
// T: time to maturity
// sigma: volatility
// q: dividend yield
// isCall: true for call option, false for put option
    double S0, K, r, T, sigma, q;
    bool isCall;
};

OptionPrice::OptionPrice(double S0, double K, double r, double T, double sigma, double q, bool isCall) {
    this->S0 = S0;
    this->K = K;
    this->r = r;
    this->T = T;
    this->sigma = sigma;
    this->q = q;
    this->isCall = isCall;
}

// Implement the member functions
// price(): return the option price using the Black-Scholes formula
double OptionPrice::price() const {
    double d1 = (log(S0 / K) + (r - q + 0.5 * sigma * sigma) * T) / (sigma * sqrt(T));
    double d2 = d1 - sigma * sqrt(T);
    if (isCall) {
        return S0 * exp(-q * T) * N(d1) - K * exp(-r * T) * N(d2);
    } else {
        return K * exp(-r * T) * N(-d2) - S0 * exp(-q * T) * N(-d1);
    }
}

// delta(): return the option delta defined as the derivative of the option price with respect to the stock price
double OptionPrice::delta() const {
    double d1 = (log(S0 / K) + (r - q + 0.5 * sigma * sigma) * T) / (sigma * sqrt(T));
    if (isCall) {
        return exp(-q * T) * N(d1);
    } else {
        return -exp(-q * T) * N(-d1);
    }
}

// gamma(): return the option gamma defined as the second derivative of the option price with respect to the stock price
double OptionPrice::gamma() const {
    double d1 = (log(S0 / K) + (r - q + 0.5 * sigma * sigma) * T) / (sigma * sqrt(T));
    return exp(-q * T) * N(d1) / (S0 * sigma * sqrt(T));
}

// theta(): return the option theta defined as the derivative of the option price with respect to the time to maturity
double OptionPrice::theta() const {
    double d1 = (log(S0 / K) + (r - q + 0.5 * sigma * sigma) * T) / (sigma * sqrt(T));
    double d2 = d1 - sigma * sqrt(T);
    if (isCall) {
        return -S0 * exp(-q * T) * N(d1) * sigma / (2 * sqrt(T)) - r * K * exp(-r * T) * N(d2) + q * S0 * exp(-q * T) * N(d1);
    } else {
        return -S0 * exp(-q * T) * N(-d1) * sigma / (2 * sqrt(T)) + r * K * exp(-r * T) * N(-d2) - q * S0 * exp(-q * T) * N(-d1);
    }
}

// vega(): return the option vega defined as the derivative of the option price with respect to the volatility
double OptionPrice::vega() const {
    double d1 = (log(S0 / K) + (r - q + 0.5 * sigma * sigma) * T) / (sigma * sqrt(T));
    return S0 * exp(-q * T) * sqrt(T) * N(d1);
}

// rho(): return the option rho defined as the derivative of the option price with respect to the risk-free interest rate
double OptionPrice::rho() const {
    double d2 = (log(S0 / K) + (r - q + 0.5 * sigma * sigma) * T) / (sigma * sqrt(T)) - sigma * sqrt(T);
    if (isCall) {
        return K * T * exp(-r * T) * N(d2);
    } else {
        return -K * T * exp(-r * T) * N(-d2);
    }
}

// impliedVol(): return the implied volatility given the target option price
double OptionPrice::impliedVol(double targetPrice) const {
    double epsilon = 1e-6;
    double low = 0.0, high = 1.0;
    double mid = 0.5;
    while (high - low > epsilon) {
        double midPrice = OptionPrice(S0, K, r, T, mid, q, isCall).price();
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
    double S0 = 100.0, K = 100.0, r = 0.05, T = 1.0, sigma = 0.3, q = 0.0;
    bool isCall = true;
    OptionPrice option(S0, K, r, T, sigma, q, isCall);
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

