#include <Rcpp.h>
using namespace Rcpp;


// [[Rcpp::export]]
NumericVector reulermultinom(double size, NumericVector rate, double dt) {
  int ncol = rate.size();
  NumericVector trans(ncol); // transition events
  double p = sum(rate); //total event rate
  double tmpp = p;
  double tmpsize = R::rbinom(size, (1 - exp(-tmpp*dt))); // total number of events
  for (int k = 0; k < (ncol-1); k++) {
    double tr = R::rbinom(tmpsize, rate(k)/p);
    trans(k) = tr;
    tmpsize -= tr;
    tmpp -= rate(k);
  }
  trans(ncol-1) = tmpsize;
  return(trans);
}

// [[Rcpp::export]]
List seir_2grp_tauleap(List params) {
  // chained operations are tricky in cpp
  // pull out list w/in list into its own object
  List init = params["init"];
  // use Rcpp as() function to "cast" R vector to cpp scalar
  int nsteps = as<int>(params["nsteps"]);
  // ints (IntegerVector), since rpois returns double,
  NumericVector SS1( nsteps, init["S1"] );
  NumericVector EE1_1(nsteps, init["E1_1"] );
  NumericVector EE1_2(nsteps, init["E1_2"] );
  NumericVector II1(nsteps, init["I1"] );
  NumericVector AA1(nsteps, init["A1"] );//asymptomatic
  NumericVector XX1(nsteps, init["X1"] );//isolated
  NumericVector RR1(nsteps, init["R1"] );
  NumericVector CI1(nsteps, init["CI1"] ); // cumulative
  
  NumericVector SS2(nsteps, init["S2"] );
  NumericVector EE2_1(nsteps, init["E2_1"] );
  NumericVector EE2_2(nsteps, init["E2_2"] );
  NumericVector II2(nsteps, init["I2"] );
  NumericVector AA2(nsteps, init["A2"] );
  NumericVector XX2(nsteps, init["X2"] );
  NumericVector RR2(nsteps, init["R2"] );
  NumericVector CI2(nsteps, init["CI2"] ); // cumulative
  // fill time w/zeros
  NumericVector time(nsteps);
  
  // pull out params for easy reading 
  NumericVector cm = params["cm"];//contact matrix
  double delta = params["delta"]; // 1 / incubation period
  double gamma = params["gamma"]; // 1 / recovery period
  double epsilon = params["epsilon"]; // 1 / latent period
  double rate_E2_I = 1 / (1/delta - 1/epsilon); // 1 / rate from pre_symptomatic infectious period to infectious period
  
  double tau = params["tau"]; // tim step size
  double fa = params["frac_asymp"]; // fraction asymptomatic
  // double rrxa = params["rel_rate_isol_asymp"]; //relative isolation rate for asymptomatic patients
  double time_intervention_start = params["time_intervention_start"];
  double time_intervention_stop = params["time_intervention_stop"];
  double dur_intervention = time_intervention_stop - time_intervention_start;
  double max_rate_x = params[ "max_rate_isol"]; // maximum isolation rate
  double final_R01 = params["final_R01"]; 
  double final_R02 = params["final_R02"];
  double R01 = params["R01"];
  double R02 = params["R02"];
  
  // Rprintf("the value of beta[1] : %f \n", beta1);
  // Rprintf("the value of beta[2] : %f \n", beta2);
  
  // Calculate the number of events for each step, update state vectors
  for (int istep = 0; istep < (nsteps-1); istep++) {
    // pull out this step's scalars for easier reading
    // and to avoid compiler headaches
    double iS1 = SS1[istep]; //  rpois returns double
    double iS2 = SS2[istep];
    double iE1_1 = EE1_1[istep];
    double iE1_2 = EE1_2[istep];
    double iE2_1 = EE2_1[istep];
    double iE2_2 = EE2_2[istep];
    double iI1 = II1[istep];
    double iI2 = II2[istep];
    double iA1 = AA1[istep];
    double iA2 = AA2[istep];
    double iX1 = XX1[istep];
    double iX2 = XX2[istep];
    double iR1 = RR1[istep];
    double iR2 = RR2[istep];
    double iCI1 = CI1[istep];
    double iCI2 = CI2[istep];
    
    /////////////////////////
    // State Equations
    /////////////////////////
    
    // R::rpois always returns a single value
    // to return multiple (e.g. Integer/NumericVector, 
    // use Rcpp::rpois(int ndraw, param) and friends
    // Prevent negative states
    double rate_x = params["rate_isol"]; //initialize
    //relative isolation rate for asymptomatic patients
    double rrxa = params["rel_rate_isol_asymp"]; 
    double beta1 = R01 / (1/rate_E2_I + fa/(gamma + rate_x*rrxa) + (1-fa)/(gamma+rate_x));
    double beta2 = R02 / (1/rate_E2_I + fa/(gamma + rate_x*rrxa) + (1-fa)/(gamma+rate_x));
    double final_beta1 = final_R01 / (1/rate_E2_I + fa/(gamma+rate_x*rrxa) + (1-fa)/(gamma+rate_x));
    double final_beta2 = final_R02 / (1/rate_E2_I + fa/(gamma+rate_x*rrxa) + (1-fa)/(gamma+rate_x));
  
    double time_diff = (istep+1)*tau - time_intervention_start;
    
    if (0 < time_diff && time_diff <= dur_intervention) {
      rate_x = rate_x + time_diff*(max_rate_x - rate_x)/dur_intervention;
      beta1 = beta1 - time_diff*(beta1 - final_beta1)/dur_intervention;
      beta2 = beta2 - time_diff*(beta2 - final_beta2)/dur_intervention;
      // rrxa = time_diff/dur_intervention;
      rrxa = rrxa + time_diff *(1 - rrxa)/dur_intervention;
    }
    
    else if (dur_intervention < time_diff) {
      beta1 = final_beta1;
      beta2 = final_beta2;
      rrxa = 1;
      rate_x = max_rate_x;
    }
    
    if ((istep+1)*tau < (time_intervention_start-3)) {// case identification occurs from Feb 18
      rate_x = 0;
    };
    
    double I1 = iE1_2 + iI1 + iA1; // asymptomatic individuals have the same transmissibility
    double I2 = iE2_2 + iI2 + iA2;
    double N1 = iS1 + iE1_1 + iE1_2 + iI1 + iA1 + iR1 + iX1;
    double N2 = iS2 + iE2_1 + iE2_2 + iI2 + iA2 + iR2 + iX2;
    double patch1 = beta1 * (I1*cm[0] + I2*cm[2]) / (N1*cm[0] + N2*cm[2]);
    double patch2 = beta2 * (I1*cm[1] + I2*cm[3]) / (N1*cm[1] + N2*cm[3]);
    double foi1 = (cm[0] * patch1 + cm[1] * patch2);
    double foi2 = (cm[2] * patch1 + cm[3] * patch2);
    
    // Rprintf("the value of new infect : %f, %f \n", new_inf1, new_inf2);
    // Rprintf("the value of contact: %f, %f,  %f, %f \n", cm[0], cm[1],cm[2], cm[3] );
     
    double newcase1 = R::rbinom(iS1, (1 - exp(-foi1 * tau)));
    double newcase2 = R::rbinom(iS2, (1 - exp(-foi2 * tau)));
    double E1_1toE1_2 = R::rbinom(iE1_1, (1 - exp(- epsilon * tau))); //
    NumericVector E1_2_rates = NumericVector::create(rate_E2_I*(1-fa), rate_E2_I*fa);
    NumericVector fromE1_2 = reulermultinom(iE1_2, E1_2_rates, tau);
    double E1_2toI1 = fromE1_2(0); //
    double E1_2toA1 = fromE1_2(1); //
    
    NumericVector I1_rates = NumericVector::create(gamma, rate_x);
    NumericVector fromI1 = reulermultinom(iI1, I1_rates, tau);
    double I1toR1 = fromI1(0); //
    double I1toX1 = fromI1(1); //
    
    NumericVector A1_rates = NumericVector::create(gamma, rate_x * rrxa);
    NumericVector fromA1 = reulermultinom(iA1, I1_rates, tau);
    double A1toR1 = fromA1(0); //
    double A1toX1 = fromA1(1); //

    double E2_1toE2_2 = R::rbinom(iE2_1, (1 - exp(- epsilon * tau))); //
    NumericVector fromE2_2 = reulermultinom(iE2_2, E1_2_rates, tau);

    double E2_2toI2 = fromE2_2(0); //
    double E2_2toA2 = fromE2_2(1); //
    
    NumericVector fromI2 = reulermultinom(iI2, I1_rates, tau);
    double I2toR2 = fromI2(0); //
    double I2toX2 = fromI2(1); //
    NumericVector fromA2 = reulermultinom(iA2, A1_rates, tau);
    double A2toR2 = fromA2(0); //
    double A2toX2 = fromA2(1); //
    
    // Calculate the change in each state variable
    double dS1 = - newcase1;
    double dE1_1 = + newcase1 - E1_1toE1_2;
    double dE1_2 = E1_1toE1_2 - E1_2toI1 - E1_2toA1;
    double dI1 = E1_2toI1 - I1toR1 - I1toX1;
    double dA1 = E1_2toA1 - A1toR1 - A1toX1;
    double dR1 = I1toR1 + A1toR1;
    double dX1 = I1toX1 + A1toX1;
    
    double dS2 = - newcase2;
    double dE2_1 = newcase2 - E2_1toE2_2;
    double dE2_2 = E2_1toE2_2 - E2_2toI2 - E2_2toA2;
    double dI2 = E2_2toI2 - I2toR2 - I2toX2;
    double dA2 = E2_2toA2 - A2toR2 - A2toX2;
    double dR2 = I2toR2 + A2toR2;
    double dX2 = I2toX2 + A2toX2;
    
    // Update next timestep
    // high-risk population
    SS1[istep+1] = iS1 + dS1;
    EE1_1[istep+1] = iE1_1 + dE1_1;
    EE1_2[istep+1] = iE1_2 + dE1_2;
    II1[istep+1] = iI1 + dI1;
    AA1[istep+1] = iA1 + dA1;
    RR1[istep+1] = iR1 + dR1;
    XX1[istep+1] = iX1 + dX1;
    
    SS2[istep+1] = iS2 + dS2;
    EE2_1[istep+1] = iE2_1 + dE2_1;
    EE2_2[istep+1] = iE2_2 + dE2_2;
    II2[istep+1] = iI2 + dI2;
    AA2[istep+1] = iA2 + dA2;
    RR2[istep+1] = iR2 + dR2;
    XX2[istep+1] = iX2 + dX2;
    
    // cumulative incidence
    CI1[istep+1] = iCI1 + newcase1;
    CI2[istep+1] = iCI2 + newcase2;
    // time in fractional years (ie units parameters are given in)
    time[istep+1] = (istep+1)*tau;
  }
  // Return results as data.frame
  DataFrame sim = DataFrame::create(
    Named("time") = time,
    Named("S1") = SS1,
    Named("E1") = EE1_1 + EE1_2,
    Named("I1") = II1,
    Named("A1") = AA1,
    Named("R1") = RR1,
    Named("X1") = XX1,
    Named("CI1") = CI1,
    
    Named("S2") = SS2,
    Named("E2") = EE2_1 + EE2_2,
    Named("I2") = II2,
    Named("A2") = AA2,
    Named("R2") = RR2,
    Named("X2") = XX2,
    Named("CI2") = CI2
  );
  
  return sim;
};

