````

----------------------------------------------------------------
MODEL: treat ~ ldh
-----------------------------------------------------------------
                 Generalized Linear Model Regression Results
==============================================================================
Dep. Variable:                  treat   No. Observations:                   76
Model:                            GLM   Df Residuals:                       74
Model Family:                Binomial   Df Model:                            1
Link Function:                  logit   Scale:                          1.0000
Method:                          IRLS   Log-Likelihood:                -44.527
Date:                Mon, 06 Aug 2018   Deviance:                       89.054
Time:                        15:12:37   Pearson chi2:                     76.0
No. Iterations:                     4   Covariance Type:             nonrobust
==============================================================================
                 coef    std err          z      P>|z|      [0.025      0.975]
------------------------------------------------------------------------------
Intercept      1.2528      0.359      3.494      0.000       0.550       1.956
ldh           -1.7123      0.514     -3.329      0.001      -2.720      -0.704
==============================================================================

----------------------------------------------------------------
MODEL: treat ~ 1
-----------------------------------------------------------------
                 Generalized Linear Model Regression Results
==============================================================================
Dep. Variable:                  treat   No. Observations:                   76
Model:                            GLM   Df Residuals:                       75
Model Family:                Binomial   Df Model:                            0
Link Function:                  logit   Scale:                          1.0000
Method:                          IRLS   Log-Likelihood:                -50.527
Date:                Mon, 06 Aug 2018   Deviance:                       101.05
Time:                        15:12:37   Pearson chi2:                     76.0
No. Iterations:                     4   Covariance Type:             nonrobust
==============================================================================
                 coef    std err          z      P>|z|      [0.025      0.975]
------------------------------------------------------------------------------
Intercept      0.4829      0.236      2.045      0.041       0.020       0.946
==============================================================================
76.0
C:\Users\zivic\AppData\Local\Programs\Python\Python36\lib\site-packages\lifelines\utils\__init__.py:250: RuntimeWarning: It looks like your weights are not integers, possibly prospenity scores then?
It's important to know that the naive variance estimates of the coefficients are biased. Instead use Monte Carlo to
estimate the variances. See paper "Variance estimation when using inverse probability of treatment weighting (IPTW) with survival analysis"
or "Adjusted Kaplan-Meier estimator and log-rank test with inverse probability of treatment weighting for survival data."

  """, RuntimeWarning)
````