# A double gamma driven OU process to evaluate credit index options

## Summary

- The goal is to build a model to price options on credit indexes

- The intensity process is assumed common for each name in the index, and its dynamics is given by an OU process driven by a double gamma process (i.e. a gamma time changed gamma process)

- The threshold that determines default for each underlying name is an exponential random variable, with the collection of all such variate (one for each name) being independent

- The dynamics of the risk free rate (here taken to be the Treasury yield) is given by an OU process driven by a gamma process

- For detail description, see the accompanying paper at: [https://arxiv.org/abs/2301.05332](https://arxiv.org/abs/2301.05332)

## About the Author

**Yoshihiro Shirai** is a Pearson Fellow at the University of Washington with expertise in applied mathematics, quantitative modeling, and macro strategy research.

- 🔗 [LinkedIn](https://www.linkedin.com/in/yoshihiro-shirai/)
- 📖 [Google Scholar](https://scholar.google.com/citations?user=...)
- 💻 [GitHub](https://github.com/yshirai999)
- 🌐 [Personal Website](https://www.yoshihiroshirai.com)

  
## Data folder

- Data on credit index options for this project was provided by Morgan Stanley

- A small sample of this dataset and its structure is provided in the Data folder

- Data on Treasuries was downloaded from [Treasury.gov](https://Treasury.gov)

## Double Gamma Intensity

Below is a summary of the main MATLAB files and scripts in this folder

- Calibration folder:
  - CalibrationRecPay1_2020IG33.m and CalibrationRecPay.m: perfom the following two tasks for, respectively, the shortest maturity traded each day between 1/2/2020 and 6/30/2020 and for each of the maturities traded on January 2 2020
    - Calibration to option market data, achieved as follows:
      - The function CDXOMCCalibrationRecPay.m, takes as input option market data (as 'Data') and initial parameters and outputs the MSE of model prices to market prices
        - Model options prices are computed via MC simulation of risk free rate and default intensity processes at the option maturity (or the start of the underlying CDX contract)
        - The price of the contract at option maturity is computed analytically by the function CDXPIDE based on the formula at the bottom of page 18 of the accompanying paper at: [https://arxiv.org/abs/2301.05332](https://arxiv.org/abs/2301.05332)
          - This function incluldes the front end protection, which is computed by the EFEP.m function
      - The current value of the intensity process is not directly observable, but can be retrieved by equating to zero the value of the forward contract at the strike for which put-call parity hold
        - This is done in line 65 of the CDXOMCCalibrationRecPay.m
      - The calibration uses MATLAB's fminsearch function, which implements Nelder Mead algorithm (lines 192 and 225 respectively in the two scripts)
  - At each day/maturities considered, the two scripts also compute model and market implied spread statistics under the annuity measure, i.e. the martingale measure corresponding to taking the annuity as numeraire
    - In both cases, this is achieved using formulas (6.2), (6.3) and (6.4) in the accompanying paper
      - These formulas require one to compute the index annuity and the forward credit spread, i.e. the expectation of the spread at option maturity under the annuity measure
        - The relationship between these two is given by the relation (6.5) in the accompanying paper
        - The annuity is a a measure of the number of underlying issuers for which a credit event has not occurred yet
        - It can also be seen as an hypothetical bond issued at par by the average firm in the index pool
        - Therefore, the annuity is in typically estimated by the fixed income team (see e.g. Table 5.1 of inputs provided by JPM for <https://onlinelibrary.wiley.com/doi/epdf/10.1111/j.1467-9965.2010.00444.x>)
        - In our case this annuity was not provided together with the data, so we approximated the forward spread with the spot spread, and then computed the annuity using the above mentioned formula (6.5)
      - For the model implied statistics, one can assess the goodness of our approximation since we have a model for the annuity measure and so (:
        - The annuity can be computed, which is done in the Annuity.m function
        - The spot spread can also be computed according to formula (3.3) in the accompanying paper, and this calculations is performed in the CDXspread.m function
        - This check is performed in lines 374-398 of the CalibrationRecPay.m script

  - Density folder:
    - MultiIntGammaFFT.m is the most important script in this folder
    - It computes the density of the vector $(r,\lambda)$ using Fourier inversion and the formula (4.12) in the accompanying paper for the characteristic exponent of the vector $(r,\lambda)$ and its integrated process
    - Similarly, univariate densities for $r$ and $\lambda$ are computed in the scripts IntGammaFFT_lambda.m and IntGammaFFT_r.m

  - FDM folder:
    - CDXO.m implements the finite difference method illustraded in the paper's appendix to solve the PIDE (4.16)
    - CDXOMC.m computes the PIDE's solution using Montecarlo
    - The error between the pricing methods is computed in the script CDXO_Test_Strike.m using the parameters calibrated to market prices on 1/2/2020 for the shortest maturity
    - The finite difference metho is also tested for the price of the forward start CDXO against the analytical formula implemented in the CDXPIDE.m function

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
