#include<iostream>
#include<vector>
#include<random>
#include<cmath>
#include<fstream>

// initial derivative parameters 
const double S = 101.15;                        // stock price
const double K = 98.01;                         // strike price
const double vol = 0.0991;                      // volatility (%)
const double r = 0.01;                          // risk-free rate (%)
const int M = 1000;                             // number of simulations

const double market_value = 3.86;               // market price of option (TRUE VALUE)
const double T = (60.0)/365;                    // time in years
const int N = 10;                               // number of time steps
const double dt = T/N;                          // small change in time


std::random_device rd;
std::mt19937 gen(rd());

// values near the mean are the most likely
// standard deviation affects the dispersion of generated values from the mean
std::normal_distribution<> snv(0, 1);


// Function to generate a random path for the underlying asset price
std::vector<double> asset_price_path() {

    double St = S;
    std::vector<double> asset_prices(N+1, 0);
    asset_prices[0] = S;

    for (int j=1; j<=N; j++) {
        double nudt = (r - 0.5*vol*vol) *dt;        // constant drift
        double volsdt = vol * sqrt(dt);             // factor for random component
        double z = snv(gen);                        // SNV

        St = St * exp(nudt + volsdt*z);             // S at (t + delta_t)
        asset_prices[j] = St;
    }
    
    return asset_prices;
}

// Function to estimate the price of the FX derivative using Monte Carlo simulation
std::vector<std::vector<double>> MCpricing() {
    
    std::vector<std::vector<double>> prices_matrix;
    std::vector<double> prices_i(N+1);
    double ST = 0.0;
    double CT = 0.0;
    double sum_CT = 0.0;
    double sum_CT2 = 0.0;

    for (int i=0; i<M; i++) {
        prices_i = asset_price_path();

        ST = prices_i[N];                            // S at time of expiration T
        CT = std::max(0.0, ST - K);                  // payoff at expiration

        prices_i.push_back(0);
        prices_i.push_back(CT);

        sum_CT = sum_CT + CT;
        sum_CT2 = sum_CT2 + CT*CT;
        prices_matrix.push_back(prices_i);
    }
    
    // Compute Expectation and SE
    double mean_CT = sum_CT/M;
    double C0 = exp(-r*T)*mean_CT;
    
    double sigma = sqrt( (sum_CT2 - sum_CT*sum_CT/M)*exp(-2*r*T) / (M-1) );
    double SE = sigma/sqrt(M);
    
    std::cout<<"Call value is "<< C0 <<" with SE +/-"<< SE<< std::endl;
    return prices_matrix;
}


int main() {
    
    std::vector<std::vector<double>> option_price_mat = MCpricing();
    std::cout << "The prices of the FX option are estimated sucessfully." << std::endl;

    std::ofstream out("test_MC.csv");
    for (auto& row : option_price_mat) {
        for (auto col : row)
            out << col <<',';
        out << '\n';
    }
    std::cout<< "test_MC.csv FILE SAVED." << std::endl;


    return 0;
}