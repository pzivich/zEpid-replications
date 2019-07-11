DATA q1; *Data input;
	INPUT treat $ outcome $ count;
	DATALINES;
P F 25
P U 25
T F 25
T U 25
;
RUN;

DATA inc; *Data input;
	INPUT treat $ outcome $ count t;
	DATALINES;
P F 5 100
P U 5 100
T F 5 100
T U 5 100
;
RUN;
DATA inc;
	SET inc;
	logt = LOG(t);
RUN;

PROC GENMOD DATA=inc DESC;
	CLASS treat;
	MODEL outcome = treat / OFFSET=logt LINK=log DIST=poisson;
	WEIGHT count;
RUN;



DATA caries2;
   INPUT treatment failure months;
   nmonths=log(months);
   DATALINES;
0 14 400
1 6  100
;
RUN;
PROC GENMOD DATA=caries2;
	*CLASS treatment(REF='con');
	MODEL failure = treatment / DIST=poisson LINK=log OFFSET=nmonths TYPE3;
	ESTIMATE 'I(0)' intercept 1 treatment 0 / EXP;
	ESTIMATE 'IR' treatment 1 / EXP;
    store out=insmodel; 
RUN;

proc plm source=insmodel;
    score data=caries2 out=PredRates pred stderr lclm uclm / nooffset ilink;
run;
   
proc print data=PredRates noobs label;
run;

proc nlmixed data=caries2;
    lambda = exp(b0 + b1*(treatment=1) + logt);
    model failure ~ poisson(lambda);
    estimate "Rate Difference" exp(b0+b1)-exp(b0) df=1e8;
run;

DATA input;
	SET input;
	logt = LOG(t);
RUN;
PROC GENMOD DATA=input;
	MODEL dead = art / LINK=log DIST=poisson OFFSET=logt;
	ESTIMATE 'IR' art 1 / EXP;
RUN;

PROC MEANS DATA=input N SUM;
	VAR dead t;
	CLASS art;
RUN;

DATA caries2;
   INPUT treatment failure months;
   nmonths=log(months);
   DATALINES;
0 77 23236.53
1 10 4077.67
;
RUN;
PROC GENMOD DATA=caries2;
	*CLASS treatment(REF='con');
	MODEL failure = treatment / DIST=poisson LINK=identity OFFSET=nmonths TYPE3;
	ESTIMATE 'IR' treatment 1 / EXP;
    store out=insmodel; 
RUN;

proc plm source=insmodel;
    score data=caries2 out=PredRates pred stderr lclm uclm / nooffset ilink;
run;
   
proc print data=PredRates noobs label;
run;

proc nlmixed data=caries2;
    lambda = exp(b0 + b1*(treatment=1) + logt);
    model failure ~ poisson(lambda);
    estimate "Rate Difference" exp(b0+b1)-exp(b0) df=1e8;
run;



DATA q1; *Data input;
	INPUT treat $ outcome $ count;
	DATALINES;
P F 40
P U 10
T F 15
T U 35
;
RUN;
PROC FREQ DATA=q1;
	TABLE treat*outcome / RISKDIFF;
	WEIGHT count;
RUN;


/* Testing IPW */
data ipw2;
	set ipw;
	where dead ^= .;
run;
proc logistic data=ipw2 descending; 
	model art = male age0 age_rs1 age_rs2 cd40 cd4_rs1 cd4_rs2 dvl0; 
	output out=denom p=d; 
run;
proc logistic data=ipw2 descending; 
	model art=; 
	output out=num p=n; 
run;
data iptw; 
	merge ipw2 denom num; 
	by id; 
	if art=1 then do;
		uw=1/d; *unstabilized weight;
		sw=n/d;
		swe=1; *unstabilized weight for exposed participants to estimate effect of exposure in the exposed;
		swu=((1-d)/d) /** ((1-n)/n)*/;*unstabilized weight for unexposed participants to estimate the effect of exposure in the unexposed;
	end;
	else if art=0 then do;
		uw=1/(1-d);
		sw=(1-n)/(1-d);
		swe=(d/(1-d)) /** (n/(1-n))*/; *unstabilized weight for exposed participants to estimate effect of exposure in the exposed;
		swu=1;  *unstabilized weight for unexposed participants to estimate effect of exposure in the unexposed;
	end;
run;

PROC MEANS DATA=iptw SUM N;
	VAR uw sw swe swu;
RUN;
PROC GENMOD DATA=iptw DESCENDING;
	CLASS id;
	MODEL dead = art / LINK=identity DIST=binomial;
	WEIGHT uw;
	REPEATED SUBJECT=id / TYPE=ind;
	ESTIMATE 'art' art 1;
RUN;

PROC GENMOD DATA=iptw DESCENDING;
	CLASS id;
	MODEL dead = art / LINK=identity DIST=binomial;
	WEIGHT sw;
	REPEATED SUBJECT=id / TYPE=ind;
	ESTIMATE 'art' art 1;
RUN;

PROC GENMOD DATA=iptw DESCENDING;
	CLASS id;
	MODEL dead = art / LINK=identity DIST=binomial;
	WEIGHT swe;
	REPEATED SUBJECT=id / TYPE=ind;
	ESTIMATE 'art' art 1;
RUN;

PROC GENMOD DATA=iptw DESCENDING;
	CLASS id;
	MODEL dead = art / LINK=identity DIST=binomial;
	WEIGHT swu;
	REPEATED SUBJECT=id / TYPE=ind;
	ESTIMATE 'art' art 1;
RUN;

*IPCW;
DATA long;
	SET long;
	IF enter=59 THEN drop=0;
RUN;

proc logistic data=long;
	model drop= enter enter_q enter_c;
	output out=n p=n; 
	run;
proc logistic data=long;
	model drop= enter enter_q enter_c male age0 age0_q age0_c dvl0 cd40 cd40_q cd40_c dvl cd4 cd4_q cd4_c;
	output out=d p=d; 
run;
data longw; 
	merge long n d; 
	by id enter;
	retain num den;
	if first.id then do; 
		num=1; 
		den=1; 
	end;
	num=num*n;
	den=den*d;
	w=num/den;
run;
proc means data=longw; 
	var w; 
run;


PROC FREQ DATA=long;
	TABLE dead*drop enter enter*drop;
RUN;


PROC FREQ DATA=mdat;
	TABLE dead;
RUN;


DATA mdat;
	SET mdat;
	IF dead = . THEN miss = 1;
		ELSE miss = 0;
RUN;

proc logistic data=mdat;
	model miss = ;
	output out=n p=n; 
	run;
proc logistic data=mdat;
	model miss = male age0 dvl0 cd40;
	output out=d p=d; 
run;
data mweight; 
	merge mdat n d; 
	sw=n/d;
	w = 1/d;
run;
proc means data=mweight; 
	var sw w; 
run;
