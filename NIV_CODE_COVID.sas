/*************************************PART1*************************************/
/* Load the dataset */
proc import datafile="/home/DATA.csv" 
    out=covid_data 
    dbms=csv 
    replace;
    getnames=yes;
run;

/* Print the dataset */
proc print data=covid_data(obs=10);
run;

TITLE 'Descriptive statistics  ';
/* Descriptive statistics */
proc means data=covid_data n mean std min p25 median p75 max;
    var weight age duration;
run;

/* Frequency table for binary variables */
proc freq data=covid_data;
    tables covidsever diabetes sex;
run;
/*************************************PART1- A*************************************/
/* Separate the data into males and females */
proc sql;
    create table males as 
    select weight 
    from covid_data 
    where sex = 0;
    
    create table females as 
    select weight 
    from covid_data 
    where sex = 1;
quit;
TITLE 'RESAMPLING METHOD ';
/* Randomization test for independent samples */
/*	 H_0 : The mean weight of females is the same as the mean weight of males.
	 H_1 : The mean weight of females is different from the mean weight of males.*/
	
PROC IML;
/* Load the data */
use males;
read all var {weight} into M;
use females;
read all var {weight} into F;

/* Calculate the observed difference in means */
obsdiff = mean(F) - mean(M);
print "Observed Difference", obsdiff;

/* Combine the data */
alldata = F // M;                     
N1 = nrow(F);  /* Number of females */
N2 = nrow(M);  /* Number of males */
N = N1 + N2;   /* Total number of observations */
NRepl = 9999;  /* Number of permutations */
nulldist = j(NRepl, 1);  /* Allocate vector to hold permuted differences */

/* Perform permutation resampling */
do k = 1 to NRepl;
   x = sample(alldata, N, "WOR");  /* Permuting the data without replacement */
   nulldist[k] = mean(x[1:N1]) - mean(x[(N2+1):N]);  /* Mean difference */
end;

/* Generate a histogram of the null distribution */
title "Histogram of Null Distribution";
refline = "refline " + char(obsdiff) + " / axis=x lineattrs=(color=red);";
call Histogram(nulldist) other=refline;

/* Calculate the p-value */
pval = (1 + sum(abs(nulldist) >= abs(obsdiff))) / (NRepl+1);
print pval;
quit;

/*The p-value associated with the observed difference is approximately 0.5851
A high p-value indicates that there is no significant evidence to reject the null hypothesis.*/

/* Bootstrap to estimate the difference in mean weights */
proc iml;
/* Load the data */
use males;
read all var {weight} into M;
use females;
read all var {weight} into F;

/* Calculate the observed difference in means */
obsdiff = mean(F) - mean(M);
print "Observed Difference", obsdiff;

/* Set the parameters for bootstrap resampling */
N1 = nrow(F);
N2 = nrow(M);
NRepl = 10000;  /* Number of bootstrap samples */
call randseed(12345);  /* Set random seed */

/* Allocate vector to hold bootstrap mean differences */
nulldist = j(NRepl, 1);

/* Perform bootstrap resampling */
do k = 1 to NRepl;
    x1 = sample(F, N1);  /* Bootstrap sample for females */
    x2 = sample(M, N2);  /* Bootstrap sample for males */
    nulldist[k] = mean(x1[1:N1]) - mean(x2[1:N2]);  /* Mean difference */
end;

/* Calculate the 2.5th and 97.5th percentiles for the confidence interval */
p = {0.025, 0.975};
call qntl(q, nulldist, p);  /* Compute quantiles */
print q;

/* Generate a histogram of the bootstrap distribution */
title "Bootstrap Distribution of Mean Weight Difference (Females - Males)";
refline = "refline " + char(q[1:2]) + " / axis=x lineattrs=(color=blue);";
call Histogram(nulldist) other=refline;

quit;

/*The observed difference in means (Females - Males) is approximately -0.461 kg. 
This indicates that, on average, females weigh slightly less than males in this dataset.
The p-value is 0.5851, so we fail to reject  H_0 . 
We do not have enough statistical evidence to assert that the mean weight is different between females and males. 
i.e., when  H_0  is true it would be unlikely to observe a difference of -0.461 kg 
and we would expect that most of the mean weight differences fluctuate between -2.092624 and 1.1684276 kg approximately. 
The 95% confidence interval for the mean difference, given by the 2.5th and 97.5th percentiles of the bootstrap distribution, is approximately (-2.092624, 1.1684276). 
Since this interval includes 0, it supports the conclusion that there is no significant difference in mean weight between females and males.*/


TITLE 'Bayesian analysis ';
/* Bayesian analysis using PROC MCMC */
ods graphics on;
proc mcmc data=covid_data nmc=20000 thin=10 nbi=5000 seed=1 monitor=(b0 b1 sigma2) 
    outpost=weight_mcmc;
    /*Parameters: (b0, b1) Regression coefficients (sigma2) variance.*/
    parms b0 b1 sigma2;
    /*Priors: Diffuse normal priors for b0 and b1*/
    prior b0 ~ normal(0, var=1000000);
    prior b1 ~ normal(0, var=1000000);
    /*inverse gamma prior for sigma2*/
    prior sigma2 ~ igamma(shape=0.01, scale=0.01);
    /*likelihood*/
    mu = b0 + b1*sex;
    model weight ~ normal(mu, var=sigma2);
run;

/* Generate predictive distributions */
data topred;
    input sex;
    datalines;
    1
    0
;
run;

data weight_mcmc_pred;
    set weight_mcmc;
    if _N_ = 1 then do;
        array sex_vals[2] _temporary_ (1 0);
    end;
    array weight_pred[2];
    do i = 1 to 2;
        weight_pred[i] = b0 + b1*sex_vals[i];
    end;
    output;
    keep weight_pred1 weight_pred2 b0 b1 sigma2;
run;

/* 95% Credible intervals for predictions */
proc iml;
    varNames = {"weight_pred1" "weight_pred2"};
    use weight_mcmc_pred;
    read all var varNames into X;
    close weight_mcmc_pred;

    Prob = {2.5, 97.5} / 100;   /* prob in (0,1) */
    call qntl(CIs, X, Prob);
    print Prob CIs[c=varNames];
run;

/* Plot the predictive distribution */
proc sgplot data=weight_mcmc_pred;
    histogram weight_pred1 / transparency=0.5 fillattrs=(color=blue) legendlabel="Females";
    histogram weight_pred2 / transparency=0.5 fillattrs=(color=red) legendlabel="Males";
    title "Predictive Distribution of Weights for Females and Males";
    xaxis label="Predicted Weight";
    yaxis label="Frequency";
run;
/*The overlapping credible intervals indicate that 
 there is no significant difference in the mean weights between females and males.
 The means of these distributions suggest that 
 the average weight for males is slightly higher than that for females, 
 but the credible intervals overlap considerably, 
 indicating that the difference is not statistically significant.*/


/*************************************PART1 B*************************************/
TITLE 'Linear Regression Analysis';
/* Fit ordinary least squares (OLS) regression model */
proc reg data=research_task;
    model duration = weight age;
    title "ordinary least squares (OLS) regression model";
    
run;
TITLE 'Bayesian analysis ';
proc mcmc data=research_task outpost=posterior nmc=20000 thin=10 seed=42;
    parms beta0 0 beta1 0 beta2 0 sigma2 1;
    prior beta0 beta1 beta2 ~ normal(0, var=1000);
    prior sigma2 ~ igamma(0.01, scale=0.01);

    model duration ~ normal(beta0 + beta1*weight + beta2*age, var=sigma2);
    title "Bayesian MCMC for Predicting Duration";
run;


/*main analysis*/
/* Logistic Model */
 PROC LOGISTIC DATA=stand_covid descending ;
 CLASS diabetes sex / param=ref ref=first ;
MODEL covidsever = weight duration age sex diabetes age|diabetes / outroc=rocdata;
ROC; ROCCONTRAST;
OUTPUT OUT=predicted P=pred; /* Save the predicted probabilities */
RUN;

DATA CUTOFF;
 SET ROCDATA;
_SPECIF_ = (1 - _1MSPEC_);
 LOGIT=LOG(_PROB_/(1-_PROB_));
 CUT_POINT=(LOGIT + 1.914)/14.744;
 J = _SENSIT_ + _SPECIF_ - 1;
 D = SQRT ((1-_SENSIT_)**2 + (1-_SPECIF_)**2);
 DIFF= ABS (_SENSIT_ - _SPECIF_);
RUN;

/* Max value of J is 0.9151794 */
PROC MEANS data=cutoff MAX;
VAR J;
RUN;

/* Confusion matrix */
DATA predicted;
    SET predicted;
    predicted_class = (pred >= 0.20357); /* You can adjust the threshold as needed */
RUN;

/* Generating Confusion Matrix */
PROC FREQ DATA=predicted;
    TABLES covidsever*predicted_class / NOCOL NOROW NOPERCENT CHISQ;
RUN;

/*GAMs Model for Continuous Predictors*/
PROC GAM DATA = covid plots = components(clm commonaxes);
MODEL covidsever (EVENT = '1') = spline(weight, df = 3) spline(duration, df = 3) spline(age, df = 3)   / dist = binary;
RUN;
 
/* Load the dataset */
proc import datafile="/home/DATA2.csv" 
    out=gp 
    dbms=csv 
    replace;
    getnames=yes;
run; 
/* Poisson Model */
PROC GENMOD DATA=gp PLOTS=all;
CLASS gps;
MODEL severeCases = gps / dist = poisson link = log;
RUN;

/* Negative Binomial Model */
PROC GENMOD DATA=gp PLOTS=all;
CLASS gps;
MODEL severeCases = gps / dist = negbin link = log;
RUN;

/* Quasi-Poisson Model */
PROC GENMOD DATA=gp PLOTS=all;
CLASS gps;
MODEL severeCases = gps / dist = poisson link = log;
SCALE=PEARSON;
RUN;