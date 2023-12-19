#include<iostream>
#include<iomanip>
#include<cmath>
#include<vector>
#include<fstream>

// Finding the price of FX Derivative (call option) using Finite difference method

// Black Scholes parameters for FX Derivative
const double volatility = 0.2;
const double domestic_rate = 0.05;
const double foreign_rate = 0.02;
const double strike = 100.0;
const double maturity = 1.0; // 1 year

// Grid parameters
const int assetPriceSteps = 10;
const int timeSteps = 5;
const double maxAssetPrice = 200.0;
const double minAssetPrice = 0.0;
const double dS = (maxAssetPrice-minAssetPrice)/assetPriceSteps; //size of each step for S
const double dT = maturity/timeSteps; //size of each step for t

double round_up(double value, int decimal_places) {
    const double multiplier = std::pow(10.0, decimal_places);
    return std::ceil(value * multiplier) / multiplier;
}

// Finite difference method
std::vector<std::vector<double>> FiniteDiffMethod() {
    std::vector<double> assetPrices(assetPriceSteps+1);
    std::vector<double> optionPrice_old(assetPriceSteps+1);
    std::vector<double> optionPrice_new(assetPriceSteps+1);

    std::vector<double> timePoint(timeSteps+1);
    

    // set initial values from the terminal condition: at time T (tau=0) option value is the final payoff
    for(int i=0; i<=assetPriceSteps; i++) {
        assetPrices[i] = i*dS; // this is just a continuous incremented steps of assetprice (min to max)
        optionPrice_old[i] = std::max(0.0, assetPrices[i]*exp(-foreign_rate*maturity) - strike); // payoff for call option
    }
    
    std::vector<std::vector<double>> matrix (1, assetPrices); 
    // first row is asset prices

    matrix.push_back(optionPrice_old);
    // row for option prices at maturity
    
    // time stepping
    for(int j=1; j<=timeSteps; j++) {
        float tau = maturity - j*dT;  // time left to maturity
        std::cout << "tau: " << tau << std::endl;
        
        // setting the boundary conditions (dirichlet boundary condition)
        optionPrice_new[0] = 0; // option is worthless when asset price is 0
        optionPrice_new[assetPriceSteps] = std::max(0.0, 
                                                    maxAssetPrice* exp(-foreign_rate* (tau)) - 
                                                    strike* exp(-domestic_rate* (tau)) );
        
       for (int i=1; i<assetPriceSteps; i++) {
            // compute the finite differences in the PDE
            
            //approx derivatives with finite differences

            // delta: sensitivity of option price wrt asset price (first derivative of the option price)
            double delta = (optionPrice_old[i+1] -optionPrice_old[i-1]) /(2*dS);
            // gamma: sensitivity of delta wrt asset price (second derivative of option price)
            double gamma = (optionPrice_old[i+1] + optionPrice_old[i-1] - 2*optionPrice_old[i]) / (dS*dS);

            // apply the explicit finite difference approxiamtion to Black Scholes PDE
            optionPrice_new[i] = optionPrice_old[i] + 
                                2*dT *( (domestic_rate-foreign_rate)*assetPrices[i]*delta +
                                        (0.5*volatility*volatility*assetPrices[i]*assetPrices[i]*gamma) -
                                        domestic_rate*optionPrice_old[i] );
            
            // optionPrice_new[i] = round_up(optionPrice_new[i], 3);
        
        }
        // dT instead of (maturity-time)

        matrix.push_back(optionPrice_new);
        optionPrice_old = optionPrice_new;

    }
    return matrix;
}

int main() {

    std::vector<std::vector<double>> output_mat = FiniteDiffMethod();
    std::cout<< std::endl << "FINITE DIFF METHOD CALLED SUCCESSFULLY." <<std::endl;

    std::ofstream out("test_new.csv");
    for (auto& row : output_mat) {
        for (auto col : row)
            out << col <<',';
        out << '\n';
    }
    std::cout<< "FILE SAVED." << std::endl;
    return 0;
}

