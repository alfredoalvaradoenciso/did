/*
	This simulated example illustrates how to estimate causal effects with event studies using a range of methods
	and plot the coefficients & confidence intervals using the event_plot command.
	
	Date: 28/05/2021
	Author: Kirill Borusyak (UCL), k.borusyak@ucl.ac.uk
	
	You'll need the following commands:
		- did_imputation (Borusyak et al. 2021): currently available at https://github.com/borusyak/did_imputation
		- did_multiplegt (de Chaisemartin and D'Haultfoeuille 2020): available on SSC
		- eventstudyinteract (Sun and Abraham 2020): currently available at https://github.com/lsun20/EventStudyInteract (SSC?)
		- csdid (Callaway and Sant'Anna 2020): currently available at https://friosavila.github.io/playingwithstata/main_csdid.html

*/

/*
	This code written by Kirill Borusyak has been adapted by David Burgherr
	(LSE, d.m.burgherr@lse.ac.uk) for the Metrics Reading Group meeting
	organized by Chloe East on 08/07/2021.
	
	Thanks to Kirill for kindly allowing us to post the adapted version of
	the code.
*/

/*
// Install commands

	*Event plot command (provided by Borusyak, Jaravel & Spiess, 2021)
	ssc install event_plot, replace

	*Borusyak, Jaravel & Spiess (2021)
	ssc install did_imputation, replace

	*Callaway & Sant'Anna (2020) (v1.5 written by Fernando Rios-Avila @friosavila)
	// ado files can be downloaded from:
	// https://friosavila.github.io/playingwithstata/main_csdid.html

	*de Chaisemartin & D'Haultfoeuille (2020)
	ssc install did_multiplegt_old, replace

	*Sun & Abraham (2020)
	github install lsun20/eventstudyinteract, replace
*/

// Generate a complete panel of 300 units observed in 15 periods

	clear all
	timer clear
	set seed 10
	global T = 15
	global I = 300

	set obs `=$I*$T'
	gen i = int((_n-1)/$T )+1 		// unit id
	gen t = mod((_n-1),$T )+1		// calendar period
	tsset i t

	
// Randomly generate treatment rollout periods uniformly across Ei=10..16
// (note that periods t>=16 would not be useful since all units are treated by then)

	*Period when unit is first treated
	gen Ei = ceil(runiform()*7)+$T -6 if t==1
	bys i (t): replace Ei = Ei[1]
	
	*Relative time, i.e. number of periods since treated (could be missing if never treated)
	gen K = t-Ei
	
	*Treatment indicator
	gen D = K>=0 & Ei!=.

	
// Generate the outcome with parallel trends and heterogeneous treatment effects

	*Heterogeneous treatment effects (in this case vary over calendar periods)
	gen tau = cond(D==1, (t-12.5), 0)
	
	*Error term
	gen eps = rnormal()
	
	*Outcome (FEs play no role since all methods control for them)
	gen Y = i + 3*t + tau*D + eps


// did_imputation of Borusyak et al. (2021)

	*Estimation
	did_imputation Y i t Ei, allhorizons pretrends(5)
	// Y:	outcome variable
	// i:	unit id variable
	// t:	time period variable
	// Ei:	variable for unit-specific treatment date (never-treated: Ei == missing)

	// allhorizons: include all non-negative horizons available
	// pretrends(): number of pre-treatment coefficients to be estimated
	// standard errors are clustered at unit level by default

	*Plotting
	event_plot, default_look graph_opt(xtitle("Periods since the event") ytitle("Average causal effect") ///
		title("Borusyak et al. (2021) imputation estimator") xlabel(-5(1)5) name(BJS))

	*Storing estimates for later
	estimates store bjs
	
	
// did_multiplegt of de Chaisemartin and D'Haultfoeuille (2020)

	*Estimation
	did_multiplegt_old Y i t D, robust_dynamic dynamic(5) placebo(5) longdiff_placebo breps(100) cluster(i)
	// Y:	outcome variable
	// i:	unit id variable
	// t:	time period variable
	// D:	treatment variable
	
	// robust_dynamic: uses estimator from dCdH (2021) on DID with intertemporal effects
	// dynamic(): number of dynamic treatment effects to be estimated (can only be used with robust_dynamic)
	// placebo(): number of placebo estimates to be estimated
	// longdiff_placebo: estimates placebo effects using long differences (comparable to dynamic TE estimates)
	// breps(): number of bootstrap iterations for computation of standard errors
	// cluster(i): computes standard errors using block bootstrap at level specified
	
	// Note that, according to the help file, by default "placebos are
	// first-difference estimators, while dynamic effects are long-difference
	// estimators, so they are not really comparable." Thus, we should not plot
	// and compare them in the same graph, if the "longdiff_placebo" option is
	// not specified.

	*Plotting
	event_plot e(estimates)#e(variances), default_look graph_opt(xtitle("Periods since the event") ///
		ytitle("Average causal effect") title("de Chaisemartin and D'Haultfoeuille (2020)") xlabel(-5(1)5) ///
		name(dCdH)) stub_lag(Effect_#) stub_lead(Placebo_#) together
	// bmat#vmat: name of point estimate matrix and name of variance-covariance matrix
	// stub_lag((prefix#postfix): name of lag coefficients in estimation output
	// stub_lead(prefix#postfix): name of lead coefficients in estimation output
	// graph_opt(): twoway options for graph overall
	// together: show leads and lags as one line

	*Storing estimates for later
	matrix dcdh_b = e(estimates)
	matrix dcdh_v = e(variances)


// csdid of Callaway and Sant'Anna (2020) (v1.5 written by Fernando Rios-Avila @friosavila)

	*Preparation
	gen gvar = cond(Ei>15, 0, Ei) // group variable as required for the csdid command

	*Estimation
	csdid Y, ivar(i) time(t) gvar(gvar) agg(event)
	// Y: 		outcome variable
	// ivar():	unit id variable
	// time():	time period variable
	// gvar():	variable for unit-specific treatment date (never treated: gvar == 0)
	// (defines "group" in CS jargon)
	
	// agg(): aggregation to use
	// wboot: Wild Bootstrap standard errors (default is asymptotic normal)
	// cluster(): should in principle be possible but throws an error (for me at least)
	// by default uses never treated units as control group (could specify "notyet")
	
	// Note that this command is work in progress. As such, it may subject to
	// ongoing changes. For example, the wboot option currently seems to
	// throw an error if specified. Further, the confidence intervals are not
	// yet correct (not uniform as in CS).
	
	// Also, note that Nick Huntington-Klein provides a Stata package that
	// acts as a wrapper for the "did" package by CS in R (via rcall, i.e. 
	// need to have R installed). It is available on his github page:
	// https://github.com/NickCH-K/did

	*Plotting
	event_plot e(b)#e(V), default_look graph_opt(xtitle("Periods since the event") ///
		ytitle("Average causal effect") xlabel(-14(1)5) title("Callaway and Sant'Anna (2020)") name(CS)) ///
		stub_lag(T+#) stub_lead(T-#) together

	*Storing estimates for later
	matrix cs_b = e(b)
	matrix cs_v = e(V)


// eventstudyinteract of Sun and Abraham (2020)

	*Preparation
	sum Ei
	gen lastcohort = Ei==r(max) // dummy for the latest- or never-treated cohort
	forvalues l = 0/5 {
		gen L`l'event = K==`l'
	}
	forvalues l = 1/14 {
		gen F`l'event = K==-`l'
	}
	drop F1event // normalize K=-1 (and also K=-15) to zero

	*Estimation
	eventstudyinteract Y L*event F*event, vce(cluster i) absorb(i t) cohort(Ei) control_cohort(lastcohort)
	// Y: outcome variable
	// L*event: lags to include
	// F*event: leads to include
	// vce(): options for variance-covariance matrix (cluster SE)
	// absorb(): absorb unit and time fixed effects
	// cohort(): variable for unit-specific treatment date (never-treated: Ei == missing)
	// control_cohort(): indicator variable for control cohort (either latest-treated or never-treated units)

	*Plotting
	event_plot e(b_iw)#e(V_iw), default_look graph_opt(xtitle("Periods since the event") ///
		ytitle("Average causal effect") xlabel(-14(1)5) title("Sun and Abraham (2020)") name(SA)) ///
		stub_lag(L#event) stub_lead(F#event) together

	*Storing estimates for later
	matrix sa_b = e(b_iw)
	matrix sa_v = e(V_iw)


// TWFE OLS estimation

	*Estimation
	reghdfe Y F*event L*event, absorb(i t) vce(cluster i)

	*Plotting
	event_plot, default_look stub_lag(L#event) stub_lead(F#event) together ///
		graph_opt(xtitle("Days since the event") ytitle("OLS coefficients") xlabel(-14(1)5) ///
		title("OLS") name(OLS))

	*Saving estimates for later
	estimates store ols
	
	
// Construct vector of true average treatment effects by number of periods since treatment

	matrix btrue = J(1,6,.)
	matrix colnames btrue = tau0 tau1 tau2 tau3 tau4 tau5
	qui forvalues h = 0/5 {
		sum tau if K==`h'
		matrix btrue[1,`h'+1]=r(mean)
	}


// Combine all plots using the stored estimates

event_plot btrue# bjs dcdh_b#dcdh_v cs_b#cs_v sa_b#sa_v ols, ///
	stub_lag(tau# tau# Effect_# T+# L#event L#event) stub_lead(pre# pre# Placebo_# T-# F#event F#event) ///
	plottype(scatter) ciplottype(rcap) ///
	together perturb(-0.325(0.13)0.325) trimlead(5) noautolegend ///
	graph_opt(title("Event study estimators in a simulated panel (300 units, 15 periods)", size(medlarge)) ///
		xtitle("Periods since the event") ytitle("Average causal effect") xlabel(-5(1)5) ylabel(0(1)3) ///
		legend(order(1 "True value" 2 "Borusyak et al." 4 "de Chaisemartin-D'Haultfoeuille" ///
				6 "Callaway-Sant'Anna" 8 "Sun-Abraham" 10 "OLS") rows(3) region(style(none))) ///
	/// the following lines replace default_look with something more elaborate
		xline(-0.5, lcolor(gs8) lpattern(dash)) yline(0, lcolor(gs8)) graphregion(color(white)) bgcolor(white) ylabel(, angle(horizontal)) ///
	) ///
	lag_opt1(msymbol(+) color(cranberry)) lag_ci_opt1(color(cranberry)) ///
	lag_opt2(msymbol(O) color(cranberry)) lag_ci_opt2(color(cranberry)) ///
	lag_opt3(msymbol(Dh) color(navy)) lag_ci_opt3(color(navy)) ///
	lag_opt4(msymbol(Th) color(forest_green)) lag_ci_opt4(color(forest_green)) ///
	lag_opt5(msymbol(Sh) color(dkorange)) lag_ci_opt5(color(dkorange)) ///
	lag_opt6(msymbol(Oh) color(purple)) lag_ci_opt6(color(purple)) 
graph export "five_estimators_example.png", replace
