version 17     			// version #.# fixes the version of Stata

program define mallows_m, eclass


	syntax varlist [if] [in] [iw/] [,  ///
					HUber /// use huber psi function; this is default		
					TUkey /// use tukey psi function
					mad /// use MAD scale est; this is default
					PROPosal2 /// use huber proposal2 scale estimationn
					inv_var /// weights are inverse of the variances, so a weight of two means this error is half as variable
					case  /// weights are case weights (weight of 2 means there are two of these)
					ls /// init with least squares; this is default
					lts /// init with an unweighted least-trimmed squares fit with 200 samples
					k2(real 1.345) /// tuning constant used for Huber proposal 2 scale estimation.
					ITERate(int `c(maxiter)') /// the limit on the number of IWLS iterations.
					TOLerance(real 1e-10) /// the accuracy for the stopping criterion.
					relax /// use tukey psi function
					]
	tempname x y weights xx yy w wt res resid coef n1 theta gamma scale result minimum done weighted_resid som som2 convi weight2 x1 y1 it wt_scalar
	local number: word count `varlist'  
	
	preserve 
	marksample touse 
	qui keep if `touse'
	
	
	//means no weight specified, so wt is initialized as a scalar
	if "`weight'" == "" {
		scalar `wt_scalar' = 0		
	}
	else {	
		scalar `wt_scalar' = 1	 
		mkmat `exp', matrix(`weights')
	}
	
	local stop = `number'-1

	local depvar1: word 1 of `varlist'    
	forval i = 1/`stop' {
	   local val = `i'+ 1
	   local x`i' : word `val' of `varlist' 	
	}
	local var_list1""
	forval i = 1/`stop' {
		local var_list1 "`var_list1' `x`i''"
	}
	mkmat `var_list1', matrix(`x') 	
	mkmat `depvar1', matrix(`y')

	mat `xx' =  `x'
	mat `yy' = `y'
	local dim = rowsof(`y')
	mat `w' = J(`dim',1,1)

	
	if ("`lts'" != "" & "`ls'" != "") {
		di as error `"you can only specify one initialization method"'
		exit 198		
	}	
	if ("`huber'" != "" & "`tukey'" != "") {
		di as error `"you can only specify one psi function"'
		exit 198		
	}
	if ("`mad'" != "" & "`proposal2'" != "") {
		di as error `"you can only specify one scale estimation method"'
		exit 198		
	}
	if ("`inv_var'" != "" & "`case'" != "") {
		di as error `"you can only specify one weight method"'
		exit 198		
	}

	
	//set defaults
	if "`lts'" != "" {
			local init "lts"
	} 
	else {
		local init "ls"
	}
	
	if  "`case'" != ""{
			local wt_method "case"
	} 
	else {
		local wt_method "inv_var"
	}
	
	if "`proposal2'" != "" {
			local scale_est "proposal2"
	} 
	else {
		local scale_est "mad"
	}

	
	if `wt_scalar' != 0 {
		if `"`wt_method'"' == "inv_var" || `"`wt_method'"' == "" { 
			matrix `x1' = .
			matrix `y1' = .
			mata: operation(st_matrix("`weights'"), st_matrix("`y'"), st_matrix("`x'"), "`y1'", "`x1'")
			matrix `x' = `x1'
			matrix `y' = `y1'
			scalar `wt_scalar' = 0
		}
		else if `"`wt_method'"' == "case" {
			
			mat `wt' = `weights'
			mat `w' = `weights'	
		}
	}
	svmat `w', name(`weight2')
	svmat `x', name(`x')
	svmat `y', name(`y')
	

	if `"`init'"' == "ls" {
		qui regress `y'1 `x'* [aweight = `weight2'1], nocons
		predict `res', residuals
		mkmat `res', matrix(`resid')
		
	}
	else if `"`init'"' == "lts" {
		
		qui robreg lqs `y'1 `x'*, nocons nsamp(200)
		predict `res', residuals
		mkmat `res', matrix(`resid')	
	}
	
	if `wt_scalar' == 0 {
		scalar `n1' = `dim' - colsof(`x')
		
	}
	else {
		tempname sum
		scalar `sum' = 0
		forvalues i = 1/`dim' {
			scalar `sum' = `sum' + `wt'[`i',1]	
		}
		scalar `n1' = `sum' - colsof(`x')
	}

	scalar `done' = 0
	scalar `theta' = 2 * normal(`k2') - 1
    scalar `gamma' =  `theta' + `k2'^2 * (1 - `theta') - 2 * `k2' *  normalden(`k2')

	if `wt_scalar' == 0 {
		tempname abs median
		gen `abs' = abs(`res')
		qui sum `abs', d
		scalar `median' = r(p50)
		scalar `scale' = `median' *  1.4826 
		
	}
	else {
		
		matrix `result' = .
		mata: wmad(st_matrix("`resid'"), st_matrix("`wt'"), "`result'")
		scalar `scale' = `result'[1,1]	
	}

	
	local it = 0
	forval i = 1/`iterate'  {

		local it = `it' + 1
		
		tempname old_res
		mat `old_res' = `resid'
		
		if `"`scale_est'"' == "mad"  {
			
			if `wt_scalar' == 0 {
				
				drop `abs' 
				tempname abs 
				gen `abs' = abs(`res')
				qui sum `abs', d
				scalar `median' = r(p50)
				scalar `scale' = `median' / 0.6745
			}
			else {
				
				matrix `result' = .
				mata: wmad(st_matrix("`resid'"), st_matrix("`wt'"), "`result'")
				scalar `scale' = `result'[1,1]	
			}
			
		}
		else {
			
			if `wt_scalar' == 0 {
				
				scalar `minimum' = 0
				forval i = 1/`dim'  {
					scalar `minimum' = `minimum' + min((`k2' * `scale')^2, `resid'[`i',1]^2)
				}
				scalar `scale' = sqrt(`minimum'/(`n1' * `gamma'))
				
			}
			else {
			
				scalar `minimum' = 0
				forval i = 1/`dim'  {
					scalar `minimum' = `minimum' + min((`k2' * `scale')^2, `resid'[`i',1]^2) * `wt'[`i',1]
				}
				
				scalar `scale' = sqrt(`minimum'/(`n1' * `gamma'))
			}
	
		}
				
		if (`scale' == 0) {
			scalar `done' = 1
			continue, break
		}		
		 
		 mat `weighted_resid' = `resid'/`scale'
		 if "`tukey'" != "" {
		 	mata: psi_tukey(st_matrix("`weighted_resid'"), "`w'")
			
		 }
		 else {	
			mata: psi_huber(st_matrix("`weighted_resid'"), "`w'")	
		 } 
		
		 
		 
		 if `wt_scalar' != 0 {
		 	
			forvalues i = 1/`dim' {
					matrix `w'[`i',1] = `weights'[`i',1] * `w'[`i',1]
			}
			
		 }
		 
		drop `res' 
		tempname res w_new
		svmat `w', name(`w_new')
		qui regress `y'1 `x'* [aweight = `w_new'1], nocons
		predict `res', residuals
		mkmat `res', matrix(`resid')
		mat `coef' = e(b)'
		
		
		mat `convi' = .
		mata: irls_delta(st_matrix("`old_res'"), st_matrix("`resid'"), "`convi'")
		scalar `convi' = `convi'[1,1]
		
		if (`convi' <= `tolerance') {
			scalar `done' = 1
			continue, break
		}		
		
		drop `w_new'1 
		
	}
			
	if (!`done') {
		if "`relax'" != "" {
		}
		else{
			di as error `"Mallows M failed to converge in `maxit' iterations"'
			eret clears
			exit 430
		}
	}
	
	tempname old_res final_res fitted
	mat `fitted' = `xx' * `coef'
	mat `final_res' = `yy' -  `fitted'
	mat `coef' = e(b)
	ereturn clear
	
	ereturn scalar N = _N
	ereturn scalar numbIt = `it'
	ereturn scalar scale = `scale'
	ereturn scalar tolerane = `tolerance'
	ereturn scalar iterate = `iterate'

	ereturn local depvar `"`depvar1'"'
	ereturn local exog `"`var_list1'"'
	ereturn local weight `"`exp'"'
	ereturn local scale_est `scale_est'
	ereturn local wt_method `wt_method'
	ereturn local init `init'
	if "`tukey'" != "" {
		ereturn local psi `tukey'
	}
	else {
		ereturn local psi `huber'
	}
	ereturn local cmdline `"mallows_m `0'"'
    ereturn local cmd "mallows_m"


	
	mat colnames `fitted' = "fitted "
	ereturn matrix fitted = `fitted'
	mat colnames `final_res' = residuals
	ereturn matrix residuals = `final_res'
	mat colnames `coef' =`var_list1'
	ereturn matrix coef = `coef'
	
	di "Coefficients:"
	matlist e(coef)
	
	di "Scale estimate: " `scale'
	
	
	restore
end
	
		
mata:

void operation(vector weights, vector y, matrix x, string scalar yresult, string scalar xresult) 
{
	 facc = sqrt(weights)
     y = y :* facc
	 x = x :* facc
	
	st_matrix(xresult,x)
	st_matrix(yresult,y)
	
	
}

void psi_huber(vector u, string scalar w)
{
	k = 1.345
	n = rows(u)
	step = J(n,1,1)

	 for (j = 1; j <= n; j++) {
	 	
		element = k/abs(u[j,1])
		
		if (element < 1) {
			
			step[j,1] = element
		}
	 }
	 
	
	 st_matrix(w, step)

}

void psi_tukey(vector u,  string scalar w) 
{
	
	c = 4.685
	n = rows(u)
	min = J(n,1,1)
	
	for (j = 1; j <= n; j++) {
		
		
		element = abs(u[j,1] / c)

		if (element < 1) {
			
			min[j,1] = element^2
		}
		
	 }
	 step = (J(n,1,1)- min):^2
	
	 st_matrix(w, step)
	
}

void irls_delta(vector old, vector neww, string scalar result)
{
 
	numerator = sum((old - neww):^2)
	som = sum(old:^2)
	
	if (som > 1e-20) {
		denum = som
	} else {
		denum = 1e-20
		
	}
	
	
	
	outcome = sqrt(numerator/denum)
	
	st_matrix(result, outcome)
   
}
	
void wmad(vector x, vector w, string scalar result)
{
	x = abs(x)
	o = order(x, 1)
	x = x[o, .]
	w = w[o,.]
	p = runningsum(w)/sum(w)
	
	n = 0
	for (i = 1; i <= rows(p); i++) {
		
		if (p[i] < 0.5) {
			n++
		}
	}
	
	if (p[n+1] > 0.5) {
		outcome = x[n + 1]/0.6745
	} else {
		outcome = (x[n + 1] + x[n + 2])/(2*0.6745)
		
	}
	st_matrix(result, outcome)
}

end

