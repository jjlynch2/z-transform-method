# SE Compare
A C++ Tool for Pair-wise Skeletal Element Osteometric Sorting

# Synopsis
This is a c++ tool for pair-wise skeletal element osteometric sorting. This tool implements several methods for pair-wise skeletal element osteometric sorting. Methods implemented include: 

* T-test approach (Byrd and Adams)
* Mean t-test approach (Lynch)
* Absolute value t-test approach (Lynch) 
* Unweighted Z-transform method (Warnke-Sommer)
* Z-transform method weighted with the inverse standard error (Warnke-Sommer)
* Z-transform method weighted with the standardized effect size (Warnke-Sommer).

This tool generates p-values for skeletal element pairs using a reference set. Leave-one-out testing can be applied to the reference sets to obtain p-values for known skeletal element pairs for validation purposes. 

# Installation
## Installation from Source
### Requirements

This tool requires r packages for some of its statistical functions. 
* Rcpp
* RInside

To install these packages using R, use the following commands in R: 
```
> install.packages("Rcpp")
> install.packages("RInside")
```
Once the packages are installed, `cd` into the main directory of the SE Compare tool and run `make` to build the executable.

The resulting binary **SE_Compare** will be located in the **build** directory.

## Binary
Precompiled binary versions of SE Compare are available for Linux and Mac OS. Requires R. (See above).  

```
SE_Compare_x86_64-pc-linux-gnu
```
```
SE_Compare_Mac_OS
```

# Input File Format
The SE Compare tool takes tab delimeted files as input. These files must contain the following columns in order:
`Skeletal Element ID`    `Side`    `Element Type`    `Measure 1`    `Measure 2`    `...`    `Measure N`     

The `Skeletal Element ID` column can be alphanumeric values. The `Side` column denotes whether the skeletal element is **Right** or **Left**. 

The `Element Type` column describes the skeletal element type; ie. Radius, Fibula, etc.  

The `Measurement Columns: 1-N` contain numeric skeletal element measurement values. 

# Example Usage

Perform Leave-One-Out Crossvalidation on reference population set using unweighted Z method

`./SE_Compare --rTrain fileR --lTrain fileL --uweightedZ TRUE` 

Use t-test method to obtain p-values for test data using a reference population set

`./SE_Compare --rTrain fileR --lTrain fileL --rClass file2R --lClass file2L --tTest TRUE`

# Sample Data

Sample data for leave-one-out validation and assigning p-values to test data sets are located in the Sample_Data directory.

The data are a subset of the DPAA CIL Osteometric Reference Data Set used for the Z-transform method by Sommer et al 2018.  Use of the data should reference both the DPAA CIL Osteometric Reference Set and the associated Sommer paper.

