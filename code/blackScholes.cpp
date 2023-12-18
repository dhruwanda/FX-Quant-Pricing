#include<iostream>
#include<cmath>

// Function to calculate the cumulative standard normal distribution
double normalCDF(double& x) {
   return ( 0.5 * erfc(-x * sqrt(0.5)) );
}

// Function to calculate the standard normal probability density function
double norm_pdf(const double& x){
    return ( (1/pow(2*M_PI, 0.5)) * exp(-0.5*x*x) );
}

// Function to calculate d1 in the Black-Scholes formula
double calc_d1(double& K, double& S, double& r, double& r_f, double& sigma, double& T) {
    return ( (log(S/K) + (r - r_f + sigma*sigma/2)*(T) )/(sigma * sqrt(T)) );
}

// Function to calculate d2 in the Black-Scholes formula
double calc_d2(double& d1, double& sigma, double& T) {
    return (d1 - sigma*sqrt(T));
}

// Function to calculate the Black-Scholes price of a European call option
double callPrice(double& S, double& K, double& r, double& r_f, double& sigma, double& T) {
    double d1 = calc_d1(K, S, r, r_f, sigma, T);
    double d2 = calc_d2(d1, sigma, T);
    return ( S * exp(-r_f*T) * normalCDF(d1) - K * exp(-r*T) * normalCDF(d2) );
}

// Function to calculate the Black-Scholes Delta of a European call option
// sensitivity of option price wrt FX/asset price
double calc_delta(double& S, double& K, double& r, double& r_f, double& sigma, double& T) {
    double d1 = calc_d1(K, S, r, r_f, sigma, T);
    return (exp(-r_f*T) * normalCDF(d1));
}

// Function to calculate the Black-Scholes Vega of a European call option
// sensitivity of optionPrice/portfolio wrt change in assumed volatility
double calc_vega(double& S, double& K, double& r, double& r_f, double& sigma, double& T) {
    double d1 = calc_d1(K, S, r, r_f, sigma, T);
    return (S*exp(-r_f * T) * norm_pdf(d1) * sqrt(T));
}


int main() {
    // Input parameters
    double sigma = 0.2;     // Volatility of the underlying Asset
    double r = 0.05;        // Risk-free rate (domestic)
    double r_f = 0.02;      // Foreign risk-free rate
    double K = 100.0;       // Strike price
    double S = 100;         // Current price (spot)
    double T = 1.0;         // 1 year
    // double S = 1.20;     // Current FX rate S_0 (e.g., USD/EUR)
    // double K = 1.30;     // Strike price
    // double r = 0.05;     // Risk-free rate (domestic)
    // double r_f = 0.02;   // Foreign risk-free rate
    // double sigma = 0.20; // Volatility of the underlying FX rate
    // double T = 1.0;      // Time to expiration in years

    double option_price = callPrice(S, K, r, r_f, sigma, T);
    double delta = calc_delta(S, K, r, r_f, sigma, T);
    double vega = calc_vega(S, K, r, r_f, sigma, T);
    
    std::cout << "The Black-Scholes price of the European call option is: " << option_price << std::endl;
    std::cout << "The Delta of the option is: " << delta << std::endl;
    std::cout << "The Vega of the option is: " << vega << std::endl;

    return 0;
}