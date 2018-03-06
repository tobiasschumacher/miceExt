# miceExt: Extension Package to mice

This package extends and builds on the mice package by adding a functionality to perform 
multivariate predictive mean matching on imputed data as well as new functionalities to
perform predictive mean matching on factor variables.



## Installation

The `miceExt` package can be installed from CRAN as follows:  
	`install.packages("miceExt")`

	
	
## Overview

Overall, miceExt provides three funtions, namely

* `mice.post.matching()`,
* `mice.binarize()`,
* `mice.factorize()`,

out of which the first function post-processes results of the mice()-algorithm by performing 
multivariate predictive mean matching on a user-defined set of column tuples, and results in 
imputations that are always equal to already-observed values, which annihilates the chance of 
getting unrealistic output values. 
The latter two functions provide a new option to impute categorical data by even extending the
functionality of `mice.post.matching()`. The function `mice.binarize()` transforms categorical
attributes of a given data frame into a binary dummy representation, which results in an
exclusively numerical data set that mice can handle well. Inconsistencies within the imputed
dummy columns can then be handled by `mice.post.matching()`, and `mice.factorize()` finally 
serves the purpose of retransforming the imputed binary data into the corresponding original
categories, resulting in a proper imputation of the given categorical data.



## Examples

### 1  Post-processing of imputated data by multivariate PMM

In this example, we work on a modification of the  `mammalsleep` data set from mice, `mammal_data`,
which is included in the miceExt-package and which has identical missing data patterns on the column
tuples (`ps`,`sws`) and (`mls`,`gt`). We want to post-process the imputations gained from after running
`mice()` on this data by performing multivariate PMM on these tuples. This procedure works in two simple
steps:

1. Run `mice()` on data set `mammal_data` and obtain a mids object to post-process:  
	 `mids_mammal <- mice(mammal_data)`


2. Run `mice.post.matching()`. As column argument `blocks` has not been specified, it will
   automatically detect the column tuples with identical missing data patterns and then 
   impute on these:  
	 `post_mammal <- mice.post.matching(mids_mammal)`

Now we can look into the resulting imputations via `post_mammal$midsobj$imp` or analyze the results via 
the `with()` function.


### 2   Imputation of categorical data

In this example, we want to impute the categorical columns `gen` and `phb` in the data set `boys` that is
included in the mice-package with the functionalities of the package. This works in three main steps:

1. Binarize the factor columns in boys the we want to impute on. By default, `mice.binarize()` will automatically identify all factor columns with missing values and binarize them.  
	`boys_bin <- mice.binarize(boys)`

2. Run `mice()` on binarized data and post-process the result with `mice.post.matching()`, as it is very likely that `mice()` imputed multiple ones among one set of dummy variables:  
	 ```
	 # run mice, note that we need to grab boys_bin$data and also use boys_bin$pred_matrix as predictor matrix for mice() 
	 # to obtain cleaner models
	 mids_boys <- mice(boys_bin$data, predictorMatrix = boys_bin$pred_matrix)  
	 
	 # post_process with mice.post.matching, use output weights from mice.binarize() to avoid imbalanced imputations
	 post_boys <- mice.post.matching(mids_boys, weights = boys_bin$weights)
	 ```

3. Retransform the resulting imputations back into categorical format:  
	 `res_boys <- mice.factorize(post_boys$midsobj, boys_bin$par_list)`

Also in this case, we can analyze the resulting imputed dataset via the `with()` function. If, e.g., we want to take
a closer look at the distribution of the values of `gen`, we can use:  
	`with(res_boys, table(gen))`