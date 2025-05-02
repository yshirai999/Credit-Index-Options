# A Levy driven OU process to evaluate credit index options

# Summary 

- The goal is to build a model to price options on credit indexes

- The intensity process is assumed common for each name in the index, and its dynamics is given by an OU process driven by a double gamma process (i.e. a gamma time changed gamma process)

- The threshold that determines default for each underlying name is an exponential random variable, with the collection of all such variate (one for each name) being independent 

- The dynamics of the risk free rate (here taken to be the Treasury yield) is given by an OU process driven by a gamma process

- For detail description, see the accompanying paper at: https://arxiv.org/abs/2301.05332
  
## Data folder

- Data on credit index options for this project was provided by Morgan Stanley

- A small sample of this dataset and its structure is provided in the Data folder

- Data on Treasuries was downloaded from https://Treasury.gov

## Double Gamma Intensity

Below is a summary of the main MATLAB files and scripts in this folder

- 

## Results

Below is a summary of the main figures in the plots folder

- CalibrationResult.pdf: shows the fitting of this model to options data as of January 2 2020
  - Calibration was performed to OTM options with ATM strike spread at 44.5 bps using the IG33 series index
  - A receiver option is OTM if the spot spread is below the option's strike spread, and the opposite holds for payer options
  - True prices are represented by stars, whereas model prices correspond to circles
  - Different parameters were fitted for each maturity considered (6 maturities in total)
  - The model fits data better for shorter maturities and it seems generally unable to price both payer and receiver options

- Kurt1, Skew1 and Vol1L: depict time series of kurtosis, skewness and volatility of the spread as implied by the model and as extracted from options data
  - For options data extraction, the Carr-Madan formula was used
  - These moments are computed under the risk neutral probability that uses as numeraire the option annuity
  - The option annuity is a a measure of the number of underlying issuers for which a credit event has not occurred yet
  - It can also be seen as an hypothetical bond issued at par by the average firm in the index pool
