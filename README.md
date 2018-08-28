# zEpid Replications
This repository contains replications of various publications using the *zEpid* package. When I find
a publication (with publicly available data) that uses a method implemented in *zEpid*, I will 
attempt to replicate the corresponding analysis.

The purpose of the replication process is; (1) verify functions within *zEpid*, (2) demonstrate 
the utility of *zEpid*, and (3) go through use-cases of *zEpid* as a way to user-test the library

Current studies replicated are;
1) Keil AP, Edwards JK, Richardson DR, Naimi AI, Cole SR. The parametric G-formula for time-to-event data: towards 
intuition with a worked example. Epidemiology (Cambridge, Mass). 2014;25(6):889-897. doi:10.1097/EDE.0000000000000160.
2) Cole SR, Hernan MA. Adjusted survival curves with inverse probability weights. Comput Methods Programs Biomed. 
2004;75(1):45-9.

If you have any recommendations on studies you would like to see replicated with *zEpid*, you 
can request either through GitHub or gmail (zepidpy)

#### Version Requirements
The following code will only run with *zEpid* version 0.3.0+. In this version, both data sets 
are available to load directly

## G-computation algorithm
The following are replications of publications that use the g-formula

### Keil et al
Using a dataset of bone marrow transplant recepients (n = 137), the marginal effect of preventing graph-versus-host 
disease compared to the natural course. SAS code for the implementation is available at 
https://github.com/alexpkeil1/Gformula-tutorial

zEpid-replications/Keil_2014 contains: 

	gformula_keil.ipynb 	(code and results)


## IPW
The following are replications of publications that use IPW (any type)

### Cole & Hernan
Cole and Hernan demonstrate adjusting survival curves via inverse probability of treatment weights (IPTW) using 
Ewing's sarcoma dataset (n = 76). 

zEpid-replications/Cole_2014 contains: 

	ipw_cole.py 	    (code and results)



