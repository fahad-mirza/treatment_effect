* Install packages needed for this visualization
ssc install palettes, replace
ssc install colrspace, replace
ssc install schemepack, replace


* Generating Sample Data
* Requires a baseline variable, and a treatment variable (one group)
sysuse auto, clear

generate treat = 0 if _n <= _N/5
replace treat = 1 if missing(treat)


tab treat, gen(treatment)
replace treatment1 = . if treatment1 != 1

* Generating sample baseline price
generate price1 = price - 3000, after(price)


********************************************************************************
* Note:
* No need to edit from this point onwards unless you want to modify the look
* of the visual and possibly change the way data is displayed or fix errors
********************************************************************************

* Regression with first treatment arm
quietly regress price treatment2 trunk, robust
matrix A1 = r(table)'

	* Control mean figures
	matrix define C1a = J(1, 9, .)
	ci means price if treatment2 == 0
	matrix C1a[1,1] = r(mean)
	matrix C1a[1,5] = r(lb)
	matrix C1a[1,6] = r(ub)
	matrix colnames C1a = b e t pvalue ll ul df crit eform
	matrix rownames C1a = control

	matrix define C1b = J(1, 9, .)
	ci means price1 if treatment2 == 0
	matrix C1b[1,1] = r(mean)
	matrix C1b[1,5] = r(lb)
	matrix C1b[1,6] = r(ub)
	matrix colnames C1b = b e t pvalue ll ul df crit eform
	matrix rownames C1b = controlbase

	
* Append Regression results with control mean	
matrix G1 = A1\C1a\C1b

* Generating column vector of 1 to denote treatment arm 1
local rows1 : rowsof G1
matrix T1 = J(`rows1', 1, 1)
matrix colnames T1 = treatarm

* Combining columnwise the regression results and treatment arm
matrix combined = G1,T1

matrix list combined

********************************************************************************
* Saving matrix as a dataset

* get original row names of matrix (and row count)
local rownames : rowfullnames combined
local c : word count `rownames'

* get original column names of matrix and substitute out _cons
local names : colfullnames combined
local newnames : subinstr local names "_cons" "cons", word

* rename columns of matrix
matrix colnames combined = `newnames'

* convert to dataset
clear
svmat combined, names(col)

* add matrix row names to dataset
gen rownames = ""
forvalues i = 1/`c' {
    replace rownames = "`:word `i' of `rownames''" in `i'
}
order rownames

********************************************************************************

* Plotting
generate treatment = treatarm

replace treatment = 1 if treatment == 1 & rownames == "controlbase"
replace treatment = 5 if treatment == 1 & rownames == "control"

replace treatment = 5 if treatment == 1 & rownames == "treatment2"

* Generate full coefficient
sort treatarm, stable

by treatarm: generate srno = _n
reshape wide rownames b se t pvalue ll ul df crit eform treatment, i(srno) j(treatarm)

generate control1 = b1 if rownames1 == "control", after(b1)
ereplace control1 = max(control1)


replace b1 = b1 + control1 if rownames1 == "treatment2"
replace ll1 = ll1 + control1 if rownames1 == "treatment2"
replace ul1 = ul1 + control1 if rownames1 == "treatment2"


* Intercept point calculation
* Since our intervention point is X=3 on the horizontal axis, we need Y value
* We will use Y = MX + C
* First we need slope: M = (y2 - y1)/(x2 - x1)

local X = 3

local bl : display `=_N'
local ctrl : display (`=_N' - 1)

local M = (b1[`ctrl'] - b1[`bl']) / (treatment1[`ctrl'] - treatment1[`bl'])
local C = b1[`ctrl'] - (`M' * treatment1[`ctrl'])
local Y = (`M' * `X') + `C'

local newobs : display `=_N' + 1
set obs `newobs'

replace treatment1 = `X' if _n == _N
replace b1 = `Y' if _n == _N
replace rownames1 = "intercept" if _n == _N


* Plot

if (b1[`ctrl'] > b1[`bl']) & (b1[`ctrl'] > b1[1]) {
	local min : display b1[`bl'] - abs(b1[`bl']*0.25)
	local max : display b1[`ctrl'] + abs(b1[`ctrl']*0.25)
	local textpos: display `min' + abs(b1[`bl']*(1/8))
}
else if (b1[`ctrl'] > b1[`bl']) & (b1[`ctrl'] < b1[1]) {
	local min : display b1[`bl'] - abs(b1[`bl']*0.25)
	local max : display b1[`ctrl'] + abs(b1[`ctrl']*0.25)
	local textpos: display `min' + abs(b1[`bl']*(1/8))	
}
else if (b1[`ctrl'] < b1[`bl']) & (b1[`ctrl'] > b1[1]) {
	local min : display b1[1] - abs(b1[1]*0.25)
	local max : display b1[`bl'] + abs(b1[`bl']*0.25)
	local textpos: display `min' + abs(b1[1]*(1/8))		
}
else {
	local min : display b1[`ctrl'] - abs(b1[`ctrl']*0.25)
	local max : display b1[`bl'] + abs(b1[`bl']*0.25)
	local textpos: display `min' + abs(b1[`ctrl']*(1/8))

}


local trt : display b1[1]
local ctr : display b1[`ctrl']

local teff : display round(abs(`ctr' - `trt'), .01)

if `trt' > `ctr' {
local val : display `trt' - abs(`teff'/2)
}
else {
local val : display `trt' + abs(`teff'/2)	
}

* Text Box values
local scale = `max'-`min'
local topb : display `textpos' + abs(`scale'*0.02)
local botb : display `textpos' - abs(`scale'*0.02)

*		(scatteri `ctr' 5.05 `ctr' 5.1, recast(line) lcolor(gs4))
*		(scatteri `trt' 5.05 `trt' 5.1, recast(line) lcolor(gs4))

#delimit ;
twoway 	(connected b1 treatment1 if inlist(rownames1, "treatment2", "intercept"), ms(i))
		(connected b1 treatment1 if inlist(rownames1, "control", "controlbase"), ms(i))
		(pcarrowi `ctr' 5.1 `trt' 5.1 , lcolor(gs4) mcolor(gs4))
		(scatteri `min' 1 `max' 1, recast(line) lcolor(dkgreen%40))
		(scatteri `min' 3 `max' 3, recast(line) lcolor(gs3%50) lpattern(dash))
		(scatteri `min' 5 `max' 5, recast(line) lcolor(red%40))
		(scatteri `val' 5.35 "`teff'", ms(i) mlabpos(0))
		(scatteri `botb' 0.7 `botb' 1.3 `topb' 1.3 `topb' 0.7, recast(area) color(white) lwidth(0.2) lcolor(dkgreen%40))
		(scatteri `textpos' 1 "Baseline", ms(i) mlabpos(0) mlabsize(2.25))
		(scatteri `botb' 2.7 `botb' 3.3 `topb' 3.3 `topb' 2.7, recast(area) color(white) lwidth(0.2) lcolor(gs3%50))
		(scatteri `textpos' 3 "Intervention", ms(i) mlabpos(0) mlabsize(2.25))
		(scatteri `botb' 4.7 `botb' 5.3 `topb' 5.3 `topb' 4.7, recast(area) color(white) lwidth(0.2) lcolor(red%40))
		(scatteri `textpos' 5 "Endline", ms(i) mlabpos(0) mlabsize(2.25))
		,
		title("{bf} Treatment Effect Visualization",  pos(11) margin(b+3 t=-1) size(*.6)) 
		subtitle("Single group automation of graph using sysuse auto",  pos(11) margin(b+6 t=-3 l=1.75) size(*.6)) 
		xlabel(none)
		xtitle("")
		ylabel(none)
		ytitle("")
		yscale(range(`min' `max') off)
		legend(order(2 1) label(1 "Treatment") label(2 "Control") size(2.25) pos(4))
		scheme(white_tableau)
		
		;
#delimit cr
