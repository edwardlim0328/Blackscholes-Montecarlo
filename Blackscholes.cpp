//Project Topic: Pricing European options using Black-Scholes formula and using Monte Carlo methods to price the greeks

#define _USE_MATH_DEFINES

#include <cmath>
#include <iostream>
#include <algorithm>

using namespace std;

class BlackScholes{
    public:
    double S; 
    double K;
    double Rate;
    double Vol;
    double divYield;
    double contractSize;
    double T;
    double x;
    int simulations;

    double D_1(){
      return (log(S/K) + ((Rate +(Vol*Vol/2))*(T)))/(Vol * sqrt(T));
      }
    double D_2(){
      return D_1() - (Vol * sqrt(T));
      }
    double norm_pdf(double x){
      return (1.0/(pow(2*M_PI,0.5)))*exp(-0.5*x*x);
      }
    double norm_cdf(double x){
      double k = 1.0/(1.0 + 0.2316419*x);
      double k_sum = k*(0.319381530 + k*(-0.356563782 + k*(1.781477937 + k*(-1.821255978 + 1.330274429*k))));
      if (x >= 0.0){
        return (1.0 - (1.0/(pow(2*M_PI,0.5)))*exp(-0.5*x*x) * k_sum);
        }
        else{
        return 1.0 - norm_cdf(-x);
      }
    }
    double CallValue(){
      return S*norm_cdf(D_1()) - (K*exp(-Rate*T)*norm_cdf(D_2()));
    }
    double PutValue(){
      return (exp(-Rate*T)*K*norm_cdf(-D_2())) - (S*norm_cdf(-D_1()));
    }
    double DeltaCall(){
      return norm_cdf(D_1());
    }
    double GammaCall(){
      return norm_cdf(D_1())/(S*Vol*sqrt(T));
    }
    double VegaCall(){
      return S*norm_pdf(D_1())*sqrt(T);
   }
    double ThetaCall(){
      return -(S*norm_pdf(D_1())*Vol)/(2*sqrt(T)) - Rate*K*exp(-Rate*T)*norm_cdf(D_2());
    }
    double RhoCall(){
      return K*T*exp(-Rate*T)*norm_cdf(D_2());
    }
    double DeltaPut(){
      return norm_cdf(D_1()) - 1;
    }
    double GammaPut(){
      return GammaCall();
    }
    double VegaPut(){
      return VegaCall();
    }
    double ThetaPut(){
      return -(S*norm_pdf(D_1())*Vol)/(2*sqrt(T)) + Rate*K*exp(-Rate*T)*norm_cdf(-D_2());
    }
    double RhoPut(){
      return -T*K*exp(-Rate*T)*norm_cdf(-D_2());
    }    
    double gaussian_box_muller() {
      double x = 0.0;
      double y = 0.0;
      double euclid_sq = 0.0;
        do {
            x = 2.0 * rand() / static_cast<double>(RAND_MAX)-1;
            y = 2.0 * rand() / static_cast<double>(RAND_MAX)-1;
            euclid_sq = x*x + y*y;
    } while(euclid_sq >= 1.0);
         return x*sqrt(-2*log(euclid_sq)/euclid_sq);
    }
};

class Montecarlo: public BlackScholes{
    public:
    double monteCarloCallPrice() {
    double S_adjust = S*exp(T*(Rate-0.5*Vol*Vol));
    double S_cur = 0.0;
    double payoff_sum = 0.0;
    for (int i=0; i<simulations; i++) {
        double gauss_bm = gaussian_box_muller();
        S_cur = S_adjust * exp(sqrt(Vol*Vol*T)*gauss_bm);
        payoff_sum += std::max(S_cur - K, 0.0);
  }
    return (payoff_sum / static_cast<double>(simulations)) * exp(-Rate*T);
}
    double monteCarloPutPrice() {
    double S_adjust = S *exp(T*(Rate-0.5*Vol*Vol));
    double S_cur = 0.0;
    double payoff_sum = 0.0;

  for (int i=0; i<simulations; i++) {
    double gauss_bm = gaussian_box_muller();
    S_cur = S_adjust * exp(sqrt(Vol*Vol*T)*gauss_bm);
    payoff_sum += std::max(K - S_cur, 0.0);
  }

  return (payoff_sum / static_cast<double>(simulations)) * exp(-Rate*T);
}
    
};

int main() {
  BlackScholes option1;
  Montecarlo option2;
  
  option1.S = 100;
  option1.K = 100;
  option1.Rate = 0.05;
  option1.Vol = 0.2;
  option1.T = 1;
    
  
  option2.S = 100;
  option2.K = 100;
  option2.Rate = 0.05;
  option2.Vol = 0.2;
  option2.T = 1;
  
  double K;
  double S;
  double T;
  double Rate;
  double Vol;
  double divYield;
  double contractSize;
  double simulations;
  

  cout << "What is the underlying stock price? " << option1.S << endl;
  cout << "What is the strike price? " << option1.K << endl;
  cout << "What is the volatility of the option? " << option1.Vol << endl;
  cout << "What is the risk free interest rate? " << option1.Rate << endl;
  cout << "What is the dividend yield? ";
  cin >> option1.divYield;
  cout << "How long until maturity? ";
  cin >> option1.T;
  cout << "How many contracts do you have? ";
  cin >> option1.contractSize;
  cout << "How many simulations would you like to run? ";
  cin >> simulations;
  cout << '\n';
  cout<<"Price of the Call option is: "<< option1.CallValue() <<endl;
  cout<<"Delta of the call is: "<<option1.DeltaCall()<<endl;
  cout<<"Gamma of the call is: "<<option1.GammaCall()<<endl;
  cout<<"Vega of the call is: "<<option1.VegaCall()<<endl;
  cout<<"Theta of the call is: "<<option1.ThetaCall()<<endl;
  cout<<"Rho of the call is: "<<option1.RhoCall()<< endl; 
  cout << '\n';
  cout<<"Price of the Put option is: "<< option1.PutValue() <<endl;
  cout<<"Delta of the Put is: "<<option1.DeltaPut()<<endl;
  cout<<"Gamma of the Put is: "<<option1.GammaPut()<<endl;
  cout<<"Vega of the Put is: "<<option1.VegaPut()<<endl;
  cout<<"Theta of the Put is: "<<option1.ThetaPut()<<endl;
  cout<<"Rho of the Put is: "<<option1.RhoPut()<<endl;  
  cout << '\n';
  cout<<"Monte Carlo Call Price: "<<option2.monteCarloCallPrice()<<endl;
  cout<<"Monte Carlo Put Price: "<<option2.monteCarloPutPrice()<<endl;

  return 0;
}
