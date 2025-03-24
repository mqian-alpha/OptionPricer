#include "BSOptionPrice.h"

#include <iostream>

using namespace std;

// Test the OptionPrice class
int main() {
    // Create an option contract object
    double K = 100.0, T = 1.0;
    bool isCall = true;
    OptionContract contract(K, T, isCall);

    // Create an option security object with market data
    double S0 = 100.0, r = 0.05, sigma = 0.25, q = 0.0;
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
