{smcl}
{* *! version 1.0.0  22nov2023}{...}
{hi:help malllows_m}{...}
{right:{browse "https://github.com/MatthiasHofstede/mallows_m"}}
{hline}

{title:Title}
{pstd}{hi:mallows_m} {hline 2} Mallows M estimator


{title:Syntax}

{p 8 16 2}
{opt mallows_m} {depvar} [{indepvars}] {ifin} [{it:{help mallows_m##weight:weight}}]
   [{cmd:,} {it:options}]

{synoptset 20 tabbed}{...}
{synopthdr}
{synoptline}
{syntab:Model}
{synopt :{opt hu:ber}} Huber psi function; this is default{p_end}
{synopt :{opt tu:key}} Tukey psi function{p_end}
{synopt :{opt ls}} initial least squares fit; this is default {p_end}
{synopt :{opt lts}} initial unweighted least-trimmed squares fit{p_end}
{synopt :{opt case}} weights are case weights; this is default{p_end}
{synopt :{opt inv_var}} weights are inverses of variances{p_end}
{synopt :{opt mad}} scale estimation using re-scaled MAD of the residuals; this is default{p_end}
{synopt :{opt prop:osal2}} scale estimation using Hubers proposal2{p_end}
{synopt :{opt k2(#)}} tuning constant for Huber proposal2 scale estimation; default is
    {cmd:k2(1.345)}{p_end}

{syntab:Optimisation}
{synopt :{opt iter:ate(#)}} the limit on the number of IRLS iterations; default is
    {cmd:set maxiter} {p_end}
{synopt :{opt tol:erance(#)}} tolerance for the Mallows-M estimator; default is
    {cmd:acc(1e-10)} {p_end}
{synopt :{opt relax}}do not return error if convergence is not reached
    {p_end}
{synoptline}


{marker weight}{...}
{p 4 6 2}
Only {cmd:iweight}s are allowed; see {help weight}. More information on how weights are treated can be found in the options {opt inv_var} and {opt case}


{marker description}{...}
{title:Description}

{pstd}
{cmd:mallows_m} fits a linear model by robust regression using a M estimator of Mallows type.


{title:Dependencies}

{pstd}
    {cmd:mallows_m} requires {cmd:robreg}; see {net "describe robreg, from(http://fmwww.bc.edu/repec/bocode/r)":ssc describe robreg}.


{marker options}{...}
{title:Options}

{dlgtab:Model}

{phang}
	{opt huber} specifies that Huber's loss function should be used as a weighting function. This is the default weighting function. Huber's loss function is defined as H(u) = min(1, 1.345 / abs(u)). 

{phang}
	{opt tukey} specifies that the Tukeys Bisweight function should be used as a weighting function. The definition of this function is T(u) = (1−min(1,∣u/4.685∣)^2)^2.

{phang}
	{opt ls} finds an initial starting point using a least squares fit. This is the default initialisation method. 

{phang}
	{opt lts} finds an initial starting point using an unweighted least-trimmed squares fit with 200 samples.

{phang}
	{opt inv_var} indicates that the weights are treated as inverses of the variances. So, a weight of two means this error is half as variable. 
	
{phang}
	{opt case} indicates that the weights are treated as case weights. So, a weight of 2 means there are two of this instance. 

{phang}
	{opt mad} specifies that the scale is estimated using re-scaled Median Absolute Deviation of the residuals. This is the default scale estimation method. 

{phang}
	{opt proposal2} indicates that the scale is estimated using Huber's proposal 2. 

{phang}
	{opt k2(#)} sets the tuning constant used in Huber's proposal 2 scale estimation. The default is {opt k2(1.345)}.


{dlgtab:optimize}

{phang}
	{opt tolerance(#)} sets the tolerance for the Mallows-M estimator. The default is {opt tolerance}(1e-10).

{phang}
	{opt iterate(#)} sets the maximum number of iterations for the Mallows-M estimator.  If convergence is not reached within {opt iterate()} iterations, the algorithm stops and
        returns an error. The default is set by {opt set maxiter}.

{phang}
    {opt relax} causes the algorithm to return the current results
    instead of returning error if convergence is not reached within
    {cmd:iterate()} iterations.


{title:Example}
{pstd}Setup{p_end}
{phang2}{cmd:. sysuse auto}{p_end}

{pstd}get Mallows M estimator{p_end}
{phang2}{cmd:. mallows_m mpg weight foreign}{p_end}

{pstd}try different initialisations{p_end}
{phang2}{cmd:. mallows_m mpg weight foreign, lts}{p_end}
{phang2}{cmd:. mallows_m mpg weight foreign, lts proposal2}{p_end}
{phang2}{cmd:. mallows_m mpg weight foreign, proposal2 tukey}{p_end}

{pstd}change the tuning constant{p_end}
{phang2}{cmd:. mallows_m mpg weight foreign, proposal2 tukey k2(1.5)}{p_end}



{marker results}{...}
{title:Stored results}

{synoptset 20 tabbed}{...}
{p2col 7 20 24 2: Scalars}{p_end}
{synopt:{cmd:e(N)}}number of observations{p_end}
{synopt:{cmd:e(numbIt)}}number of iterations used by the algorithm{p_end}
{synopt:{cmd:e(scale)}}the robust scale estimate used{p_end}
{synopt:{cmd:e(tolerance)}}tolerance for Mallows-M estimation{p_end}
{synopt:{cmd:e(iterate)}}maximum number of iterations of Mallows-M estimation{p_end}

{synoptset 20 tabbed}{...}
{p2col 7 20 24 2: Macros}{p_end}
{synopt:{cmd:e(cmdline)}}command as typed{p_end}
{synopt:{cmd:e(cmd)}}{cmd:mallows_m}{p_end}
{synopt:{cmd:e(depvar)}}name of dependent variable{p_end}
{synopt:{cmd:e(exog)}}name of exogenous variables{p_end}
{synopt:{cmd:e(weight)}}name of weight variable{p_end}
{synopt:{cmd:e(scale_est)}}scale estimation method{p_end}
{synopt:{cmd:e(wt_method)}}the weighting method{p_end}
{synopt:{cmd:e(init)}}method used to initialise the algorithm{p_end}
{synopt:{cmd:e(psi)}}name of psi weighting function{p_end}

{synoptset 20 tabbed}{...}
{p2col 7 20 24 2: Matrices}{p_end}
{synopt:{cmd:e(fitted)}}the fitted values of the estimation{p_end}
{synopt:{cmd:e(residuals)}}the residuals of the estimation{p_end}
{synopt:{cmd:e(coef)}}coefficients of the Mallows M estimator{p_end}

{title:Authors}

{pstd}
    Matthias Hofstede (Erasmus University Rotterdam),

{pstd}
    Support: hofstede@ese.eur.nl

{pstd}
    Thanks for citing this software as follows:

{pmore}
    Hofstede M. (2024). mallows_m: Stata module providing Mallows M estimation. Available from
    {browse "https://github.com/MatthiasHofstede/mallows_m"}.


