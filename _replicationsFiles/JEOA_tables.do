*******
** This do-file builds all tables of 
**			ECONOMIC EFFECTS OF DEMOGRAPHIC DIVIDEND IN BRAZILIAN REGIONS		
*******

** Main folder
eststo clear

use "JEOA_Data.dta", replace

gen reg=substr(string(mesocode), 1, 1)
destring reg, replace

label var y "$\ln y_{it}$"
label var r_lyw "$\ln \tilde{y}_{it-1}$"
label var r_gk "$\Delta \ln k_{it}$"
label var r_gl "$\Delta \ln L_{it}$" 
label var r_gn "$\Delta \ln N_{it}$"
label var r_ge "$\Delta s_{it}$"

label var r_ly "$\ln y_{it-1}$"
label var r_lw "$\ln wa_{it-1}$"
label var r_lp "$\ln p_{it-1}$"
label var r_le "$ s_{it-1}$"

label var r_nr "$ nr_{it}$"
label var r_nm "$ nm_{it}$"
label var r_br "$ cbr_{it}$"
label var r_dr "$ cdr_{it}$"

gen dy70 = 1 if year == 1970
gen dy80 = 1 if year == 1980
gen dy91 = 1 if year == 1991
gen dy00 = 1 if year == 2000
recode dy70 dy80 dy91 dy00 (. = 0)

bysort microcode (year): gen Ddy80 = dy80 - dy80[_n-1]
bysort microcode (year): gen Ddy91 = dy91 - dy91[_n-1]
bysort microcode (year): gen Ddy00 = dy00 - dy00[_n-1]

egen decade = group(year) 
egen ccode = group(mesocode)
drop if reg==1 | reg == 5

xtset microcode decade

xtreg y r_ly r_lp r_lw r_gk r_gl r_gn dy*, fe cluster(ccode)

*** TESTS

local level (_b[r_lw] = -_b[r_ly]+1) (_b[r_lp] = -_b[r_ly]+1)
local growth1 (_b[r_gl] = 1) (_b[r_gn] = -1)
local growth2 (_b[r_gl] = 1) (_b[r_nr] = -1) (_b[r_nm] = -1)
local growth3 (_b[r_gl] = 1) (_b[r_br] = -1) (_b[r_dr] = 1) (_b[r_nm] = -1)

qui forvalues i = 1/3 {

	global dep y
	global pre r_ly r_lp r_lw
	global preeduc r_le
	global endeduc r_ge

	if `i' == 1{
		*** CLS (2014) 
		global end r_gk r_gl r_gn 
		global endctrl r_gk r_gl r_gn
		*global endctrl r_gk r_gl r_br r_dr r_nm
	}
	if `i' == 2{
		*** CLS (2014) + NRI and NM 
		global end r_gk r_gl r_nr r_nm 
		global endctrl r_gk r_gl r_nr r_nm
		*global endctrl r_gk r_gl r_br r_dr r_nm
	}
	
	if `i' == 3{
		*** CLS (2014) + CBR, CDR and NM 
		global end r_gk r_gl r_br r_dr r_nm 
		global endctrl r_gk r_gl r_br r_dr r_nm
		
	}
	
	* BODY
	
	xtabond2 $dep $pre $end dy*, gmm($pre) ///
	gmm($endctrl, lag(2 .)) iv(dy*, eq(level)) robust small nodiffsargan cluster(ccode) two 
	eststo reg`i'_1
	qui estadd local col "No"
	estadd local mode "Sys."
	estadd local pred "Pred."
	estadd scalar mr = e(N)/3
	estadd scalar inst = e(j)
	estadd scalar hjp = e(hansenp)
	estadd scalar sp = e(sarganp)
	
	qui xtabond2 $dep $pre $end $endeduc dy*, gmm($pre) ///
	gmm($endctrl $endeduc, lag(2 .)) iv(dy*, eq(level)) robust small nodiffsargan cluster(ccode) two 
	eststo reg`i'_2
	estadd local col "No"
	estadd local mode "Sys."
	estadd local pred "Pred."
	estadd scalar mr = e(N)/3
	estadd scalar inst = e(j)
	estadd scalar hjp = e(hansenp)
	estadd scalar sp = e(sarganp)
	
	qui xtabond2 $dep $pre $end $endeduc $preeduc dy*, gmm($pre $preeduc) ///
	gmm($endctrl $endeduc, lag(2 .)) iv(dy*, eq(level)) robust small nodiffsargan cluster(ccode) two 
	eststo reg`i'_3
	estadd local col "No"
	estadd local mode "Sys."
	estadd local pred "Pred."
	estadd scalar mr = e(N)/3
	estadd scalar inst = e(j)
	estadd scalar hjp = e(hansenp)
	estadd scalar sp = e(sarganp)
	
	qui xtreg $dep $pre $end dy*, fe cluster(ccode)
	eststo reg`i'_4
	estadd scalar mr = e(N)/3
	estadd local mode "WG"
	
	qui xtreg $dep $pre $end $endeduc dy*, fe cluster(ccode)
	eststo reg`i'_5
	estadd scalar mr = e(N)/3
	estadd local mode "WG"
	
	qui xtreg $dep $pre $end $endeduc $preeduc dy*, fe cluster(ccode)
	eststo reg`i'_6
	estadd scalar mr = e(N)/3
	estadd local mode "WG"


	*** APPENDIX

	* INSTRUMENTS

	forvalues d =2/4{
		forvalues lag = 0/`=`d'-1'{
			foreach v in $pre $preeduc {
				sort microcode decade
				qui gen pZ`d'L`lag'L`v' = L`lag'.`v' if decade == `d'
				qui gen pZ`d'L`lag'D`v' = D.L`lag'.`v' if decade == `d'
			}
		}
		forvalues lag = 1/`=`d'-1'{
			foreach v in $endctrl $endeduc {
				sort microcode decade
				qui gen eZ`d'L`lag'L`v' = L`lag'.`v' if decade == `d'
				qui gen eZ`d'L`lag'D`v' = D.L`lag'.`v' if decade == `d'
			}
		}
	}
	qui recode eZ* pZ* (. = 0)

	* COLLAPSED INSTRUMENTS

	forvalues lag = 0/3{
		foreach v in $pre $preeduc  {
			qui gen pCL`lag'L`v'=L`lag'.`v'
			qui gen pCL`lag'D`v'=D.L`lag'.`v'
		}
	}
	forvalues lag = 1/3{
		foreach v in $endctrl $endeduc {
			qui gen eCL`lag'L`v'=L`lag'.`v'
			qui gen eCL`lag'D`v'=D.L`lag'.`v'
		}
	}
	qui recode eC* pC* (. = 0)

	global Elevel eZ*L2L* eZ*L3L* 
	global Plevel pZ*L1L* pZ*L2L* pZ*L3L*
	global Ediff  eZ*L1D* eZ*L2D* eZ*L3D* 
	global Pdiff  pZ*L0D* pZ*L1D* pZ*L2D* pZ*L3D*

	global EPlevel pZ*L2L* pZ*L3L*
	global EPdiff  pZ*L1D* pZ*L2D* pZ*L3D*

	global cElevel eCL2L* eCL3L*
	global cPlevel pCL1L* pCL2L* pCL3L*
	global cEdiff  eCL1D* eCL2D* eCL3D*
	global cPdiff  pCL0D* pCL1D* pCL2D* pCL3D*

	qui reg $dep $pre $preeduc $end $endeduc dy*, vce(cluster ccode)
	eststo rob_reg`i'_1
	estadd scalar mr = e(N)/3
	estadd local mode "OLS"

	***FE
	qui xtreg $dep $pre $end $endeduc $preeduc dy*, fe cluster(ccode)
	eststo rob_reg`i'_2
	estadd scalar mr = e(N)/3
	estadd local mode "WG"

	*Diff.

	qui ivreg2 D.($dep) (D.($pre $end $endeduc $preeduc) =  $Elevel $Plevel) Ddy*, cluster(ccode) noconstant robust
	scalar kplm0 = e(idp)
	scalar kpf0 = e(widstat)
	qui xtabond2 $dep $pre $end $endeduc $preeduc dy*, gmm($pre $preeduc) ///
	gmm($endctrl $endeduc, lag(2 .)) iv(dy*, eq(level)) robust small nodiffsargan cluster(ccode) two nolevel
	eststo rob_reg`i'_3
	estadd local col "No"
	estadd local mode "Diff."
	estadd local pred "Pred."
	estadd scalar mr = e(N)/2
	estadd scalar inst = e(j)
	estadd scalar hjp = e(hansenp)
	estadd scalar sp = e(sarganp)
	estadd scalar kplm = kplm0
	estadd scalar kpf = kpf0

	qui ivreg2 D.($dep) (D.($pre $end $endeduc $preeduc) = $cElevel $cPlevel) Ddy*, cluster(ccode) noconstant robust
	scalar kplm0 = e(idp)
	scalar kpf0 = e(widstat)
	qui xtabond2 $dep $pre $end $endeduc $preeduc dy*, gmm($pre $preeduc, collapse) ///
	gmm($endctrl $endeduc, collapse lag(2 .)) iv(dy*, eq(level)) robust small nodiffsargan cluster(ccode) two nolevel
	eststo rob_reg`i'_4
	estadd local col "Yes"
	estadd local mode "Diff."
	estadd local pred "Pred."
	estadd scalar mr = e(N)/2
	estadd scalar inst = e(j)
	estadd scalar hjp = e(hansenp)
	estadd scalar sp = e(sarganp)
	estadd scalar kplm = kplm0
	estadd scalar kpf = kpf0

	*** Sys.

	qui ivreg2 $dep ($pre $end $endeduc $preeduc = $Elevel $Plevel $Ediff $Pdiff) dy*, robust cluster(ccode)
	scalar kplm0 = e(idp)
	scalar kpf0 = e(widstat)
	qui xtabond2 $dep $pre $end $endeduc $preeduc dy*, gmm($pre $preeduc) ///
	gmm($endctrl $endeduc, lag(2 .)) iv(dy*, eq(level)) robust small nodiffsargan cluster(ccode) two 
	eststo rob_reg`i'_5
	estadd local col "No"
	estadd local mode "Sys."
	estadd local pred "Pred."
	estadd scalar mr = e(N)/3
	estadd scalar inst = e(j)
	estadd scalar hjp = e(hansenp)
	estadd scalar sp = e(sarganp)
	estadd scalar kplm = kplm0
	estadd scalar kpf = kpf0

	qui ivreg2 $dep ($pre $end $endeduc $preeduc = $cElevel $cPlevel $cEdiff $cPdiff) dy*, cluster(ccode) robust
	scalar kplm0 = e(idp)
	scalar kpf0 = e(widstat)
	qui xtabond2 $dep $pre $end $endeduc $preeduc dy*, gmm($pre $preeduc, collapse) ///
	gmm($endctrl $endeduc, collapse lag(2 .)) iv(dy*, eq(level)) iv(dy*, eq(level)) robust nodiffsargan small cluster(ccode) two 
	eststo rob_reg`i'_6
	estadd local col "Yes"
	estadd local mode "Sys."
	estadd local pred "Pred."
	estadd scalar mr = e(N)/3
	estadd scalar inst = e(j)
	estadd scalar hjp = e(hansenp)
	estadd scalar sp = e(sarganp)
	estadd scalar kplm = kplm0
	estadd scalar kpf = kpf0

	qui ivreg2 $dep ($pre $end $endeduc $preeduc = $Elevel $EPlevel $Ediff $EPdiff) dy*, cluster(ccode) robust
	scalar kplm0 = e(idp)
	scalar kpf0 = e(widstat)
	qui xtabond2 $dep $pre $end $endeduc $preeduc dy*, gmm($pre $preeduc, lag(2 .)) ///
	gmm($endctrl $endeduc, lag(2 .)) iv(dy*, eq(level)) robust small nodiffsargan cluster(ccode) two
	eststo rob_reg`i'_7
	estadd local col "No"
	estadd local mode "Sys."
	estadd local pred "End."
	estadd scalar mr = e(N)/3
	estadd scalar inst = e(j)
	estadd scalar hjp = e(hansenp)
	estadd scalar sp = e(sarganp)
	estadd scalar kplm = kplm0
	estadd scalar kpf = kpf0
	
	qui ivreg2 $dep ($pre $end $endeduc $preeduc = $Elevel $cPlevel $Ediff $cPdiff) dy*, cluster(ccode) robust
	scalar kplm0 = e(idp)
	scalar kpf0 = e(widstat)
	qui xtabond2 $dep $pre $end $endeduc $preeduc dy*, gmm($pre $preeduc, lag(1 1)) ///
	gmm($endctrl $endeduc, lag(2 .)) iv(dy*, eq(level)) robust small nodiffsargan cluster(ccode) two
	eststo rob_reg`i'_8
	estadd local col "No"
	estadd local mode "Sys."
	estadd local pred "Lag 1"
	estadd scalar mr = e(N)/3
	estadd scalar inst = e(j)
	estadd scalar hjp = e(hansenp)
	estadd scalar sp = e(sarganp)
	estadd scalar kplm = kplm0
	estadd scalar kpf = kpf0
	
	forvalues j=1/6{
		est restore reg`i'_`j'
		qui test `growth`i''
		estadd scalar test1 = r(p)
		qui test `level'
		estadd scalar test2 = r(p)
		qui test `growth`i'' `level'
		estadd scalar test3 = r(p)
	}

	forvalues j=1/8{
		est restore rob_reg`i'_`j'
		qui test `growth`i''
		estadd scalar test1 = r(p)
		qui test `level'
		estadd scalar test2 = r(p)
		qui test `growth`i'' `level'
		estadd scalar test3 = r(p)
	}
	

	drop eZ* pZ* eC* pC*
}



* TABLE 1
esttab reg1*, b(2) se(2) replace label nodep nomtitle mgroups("Sys-GMM" "Within Groups", pattern(1 0 0 1 0 0) prefix(\multicolumn{@span}{c}{) suffix(}) span erepeat(\cmidrule(lr){@span})) drop(r_ly r_gk dy* _cons*) compress substitute(% \%) ///
	scalars("mr Micro-Regions" ///
			"hjp Hansen-J (p-value)" ///
			"test1 Growth (p-value)" ///
			"test2 Level (p-value)" ///
			"test3 All (p-value)") /// 
	sfmt(%9.0g %8.2f %8.2f %8.2f %8.2f) nonotes nogaps

* TABLE 2
esttab reg2*, b(2) se(2) replace label nodep nomtitle mgroups("Sys-GMM" "Within Groups", pattern(1 0 0 1 0 0) prefix(\multicolumn{@span}{c}{) suffix(}) span erepeat(\cmidrule(lr){@span})) order(r_lp r_lw r_gl r_nm r_nr r_ge r_le) drop(r_ly r_gk dy* _cons*) compress substitute(% \%) ///
	scalars("mr Micro-Regions" ///
			"hjp Hansen-J (p-value)" ///
			"test1 Growth (p-value)" ///
			"test2 Level (p-value)" ///
			"test3 All (p-value)") /// 
	sfmt(%9.0g %8.2f %8.2f %8.2f %8.2f) nonotes nogaps	

* TABLE 3
esttab reg3*, b(2) se(2) replace label nodep nomtitle mgroups("Sys-GMM" "Within Groups", pattern(1 0 0 1 0 0) prefix(\multicolumn{@span}{c}{) suffix(}) span erepeat(\cmidrule(lr){@span})) order(r_lp r_lw r_gl r_nm r_br r_dr r_ge r_le) drop(r_ly r_gk dy* _cons*) compress substitute(% \%) ///
	scalars("mr Micro-Regions" ///
			"hjp Hansen-J (p-value)" ///
			"test1 Growth (p-value)" ///
			"test2 Level (p-value)" ///
			"test3 All (p-value)") /// 
	sfmt(%9.0g %8.2f %8.2f %8.2f %8.2f) nonotes nogaps		

*** APPENDIX

* TABLE 4
esttab rob_reg1*, b(2) se(2) replace label nodep nomtitle drop(dy* _cons*) compress substitute(% \%) ///
	scalars("mr Microregions" ///
			"test1 Growth (p-value)" ///
			"test2 Level (p-value)" ///
			"test3 All (p-value)" ///	
			"mode Method" ///
			"pred Lagged Var." ///			
			"inst Instruments" ///
			"col Collapsed" ///
			"hjp Hansen-J p-value" ///
			"kplm Kleibergen-Paap LM" ///
			"kpf Kleibergen-Paap rk F") ///
	sfmt(%9.0g %8.2f %8.2f %8.2f a1 %9.0g a1 a1 %8.2f %8.2f) nonotes nogaps

* TABLE 5
esttab rob_reg2*, b(2) se(2) replace label nodep nomtitle drop(dy* _cons*) compress substitute(% \%) ///
	scalars("mr Microregions" ///
			"test1 Growth (p-value)" ///
			"test2 Level (p-value)" ///
			"test3 All (p-value)" ///	
			"mode Method" ///
			"pred Lagged Var." ///			
			"inst Instruments" ///
			"col Collapsed" ///
			"hjp Hansen-J p-value" ///
			"kplm Kleibergen-Paap LM" ///
			"kpf Kleibergen-Paap rk F") ///
	sfmt(%9.0g %8.2f %8.2f %8.2f a1 %9.0g a1 a1 %8.2f %8.2f) nonotes nogaps
 
* TABLE 6
esttab rob_reg3*, b(2) se(2) replace label nodep nomtitle drop(dy* _cons*) compress substitute(% \%) ///
	scalars("mr Microregions" ///
			"test1 Growth (p-value)" ///
			"test2 Level (p-value)" ///
			"test3 All (p-value)" ///	
			"mode Method" ///
			"pred Lagged Var." ///			
			"inst Instruments" ///
			"col Collapsed" ///
			"hjp Hansen-J p-value" ///
			"kplm Kleibergen-Paap LM" ///
			"kpf Kleibergen-Paap rk F") ///
	sfmt(%9.0g %8.2f %8.2f %8.2f a1 %9.0g a1 a1 %8.2f %8.2f) nonotes nogaps

