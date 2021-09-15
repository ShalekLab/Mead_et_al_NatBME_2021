R QUALITY scripts

AS OF: 9/29/17
Pipeline.R: Calls functions and loads experimental data to execute quality control tests
	Experiment_Details.xlsx must have:
		- Assays: assays performed ("CTG", "LYZ", "CTO")
		- Cont: names of controls in the compound list ("POS", "POSL")
		- C.Std:  number of standard deviations to define outliers for the control population (3)
		- Dates: dates the experiments were performed (MMDD)
		- Plates: plates labels to be listed in data frame 
      			must be "mouse.plate#"" e.g. c("27.P1","27.P2","28.P1","28.P2")
		- File: name of excel file containing raw data to import
		- Comp: headers of the compound list excel sheet (varied)
		- H.Quad: headers of the quadrant stamp excel sheet (varied)			
		- H.POS: headers of the positive control quadrant stamp excel sheet (varied)
		- Headers: headers of the raw data excel sheet for assay results ("Mouse_plate", consistent with raw data naming format)
		- QCtest: List of quality control tests (CTO, CTG, LYZ.NS.NC, LYZ.S.NS)
		- QCval: B score limit for corresponding QCtest (1 -> fairly moderate, set higher for stronger effect size) 

	Quadrant labeling requires all ==> Dose Cell.Type S/NS-Replicate Number 


Get Experiment Details
"Assays"  	- vector of assays performed ("CTG", "LYZ", "CTO")
"Cont"    	- vector of names of controls in the compound list ("POS", "POSL")
"dates"   	- vector of dates the experiments were performed (MMDD)
                                                             
Import.R: Reads Quadrants, Compound lists, raw data arrays and saves as RDS

Convert.R: Converts arrays to tables, adds LC assay and saves as RDS with columns [row, col, Quad, Plate, Treatment, Index, Assays]
	- Saves B scored table
	- Saves NC subtracted table
	- Saves Raw data table


Quality.R: 
	- Removes Control outliers based on C.Std number of standard deviations
	- Finds LOESS result on Raw data
	- Filters tables for plates that fail the Quality control test for the NC subtracted data
		- CTG compares NC AND POS
		- CTO compares POSL and POS
		- LYZ compares NC and NonStim, NonStim and Stim
	- Saves RDS and CSV of filtered tables
	- Plots all plates in one PNG file



arrangement of the CTG.arr.rds: dim1&2: rows & columns of raw plates, dim3: named plates

arrangement of the LYZ.arr.rds: dim1&2: rows & columns of raw plates, dim3: named plates

arrangement of the CTO.arr.rds: dim1&2: rows & columns of raw plates, dim3: named plates

arrangement of the Comp_layout.rds: dim1&2: rows & columns of compound plate, dim3: named plates

arrangement of the Fulltbl.rds: data table of all expts B-score, NC-subtracted, values (all NONE and NC have been removed)

Revisions:

Quality.R and Quality_fn.R: 
Version 2.5 	- Finds Z score of S and NS population separately
Version 2.6  	- Does not remove outliers from Treatment population
Version 3 	  - Removed quality control for CTG and CTO
		- Created new script Optional_Quality_fn.R with old Ctrl_test() and B_Plots() code
Version 4     - Incorporated cell type for quadrants
