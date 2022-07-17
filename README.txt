README
-This is an initial commit of the code used to generate the results in the paper “Novel methods for estimating the instantaneous and overall COVID-19 case fatality risk among care home residents in England” by Overton, C.E. et al. 
-The code wrapper is very specific to the structure of the data used in this project.
-The model functions however are easily generalisable to other settings.
-I intend to make the model functions more generalisable when I have time, but the key parts are easy to extract from the code.

The key model functions are:
-estimate_CFR_backwards: Estimates the backward CFR for specific level - daily, weekly, overall.
-estimate_CFR_forwards: Estimates the forward CFR for specific level - daily, weekly, overall.
-estimate_delays_variable: Estimates the delay distribution for each month from line-list data.


