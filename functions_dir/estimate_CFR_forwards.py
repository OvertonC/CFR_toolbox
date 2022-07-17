def estimate_CFR_forwards(death_type,new_df2,dateList,parameters,CFR_type,variable_option):
    #### Define binomial based likelihood functions
    import numpy as np
    import pandas as pd
    from scipy import integrate
    from scipy import stats as sp
    import matplotlib.pyplot as plt
    import math
    from scipy.optimize import minimize
    from scipy.stats import weibull_min
    from scipy.special import gamma
    #### load required packages end
    from tqdm.notebook import tqdm
    from datetime import datetime, date
    from datetime import timedelta
    def loglikelihood(x):
        # calculates the inverse of the loglikelihood function for observing data I
        # when the distribution is truncated at time T, given shape a and scale b 
        event1=round(I[0]); # first events
        event2=round(I[1]); # second events
        p = x
        L = sp.binom.pmf(event1,event2,p);

        logL=((np.log(L))**(-1)); # calculate inverse loglikelihood
        return(logL)
    def function_upper(x):
        a = loglikelihood(x)**(-1) - loglikelihood(output_data.x)**(-1)
        b = abs(3.84/2 + a)
        return(b + np.sign(output_data.x - x))

    def function_lower(x):
        a = loglikelihood(x)**(-1) 
        c = loglikelihood(output_data.x)**(-1)
        b = abs(3.84/2 + a - c)
        return(b + np.sign(x - output_data.x))

    def function_75(x):
        a = loglikelihood(x)**(-1) - loglikelihood(output_data.x)**(-1)
        b = abs(0.455/2 + a)
        return(b + np.sign(output_data.x - x))

    def function_25(x):
        a = loglikelihood(x)**(-1) 
        c = loglikelihood(output_data.x)**(-1)
        b = abs(0.455/2 + a - c)
        return(b + np.sign(x - output_data.x))


    #### Estimate the case fatality ratio
    if death_type != "60":
        death_cut = 28
    else:
        death_cut = 28


     #### Calculate daily CFR

     ## add naive method - 13 days
    cases = new_df2["Cases"][:]
    deaths = new_df2["Deaths"][:]
    C = np.zeros(len(dateList))
    g = np.zeros(len(dateList))
    output = np.zeros(len(dateList))
    CFR_daily = np.zeros(len(dateList))
    output_lower = np.zeros(len(dateList))
    output_upper = np.zeros(len(dateList))
    output_25 = np.zeros(len(dateList))
    output_75 = np.zeros(len(dateList))
    numerator = np.zeros(len(dateList))
    denominator = np.zeros(len(dateList))
    adj = np.zeros(len(dateList))
    
    if death_type == "Confirmed":
        a = parameters[0,0]
        b = parameters[0,1] + 1.5
    else:
        a = parameters[0,0]
        b = parameters[0,1]
    month_index = ["2020-01","2020-02","2020-03","2020-04","2020-05","2020-06",
                   "2020-07","2020-08","2020-09","2020-10","2020-11","2020-12",
                   "2021-01","2021-02","2021-03","2021-04","2021-05","2021-06",
                   "2021-07","2021-08","2021-09","2021-10","2021-11","2021-12",
                  "2022-01"]
    index_overall = len(month_index)
    if CFR_type == "daily": 
######### Forwards
        adj = np.zeros(len(dateList)+death_cut)
        adj2 = np.zeros(len(dateList)+death_cut)
        adjust = np.zeros(len(dateList)+death_cut)
        for i in tqdm(range(0,len(dateList)-death_cut)):
            adjust = np.zeros(len(dateList)+death_cut)            
            for j in range(i,i+death_cut):
                index_cur = pd.to_datetime(dateList[j]).to_period('M').strftime("%Y-%m")
                index_cur = month_index.index(index_cur)
                if variable_option == 0 :
                    index_cur = index_overall
                if np.isnan(parameters[index_cur,0]):
                    index_cur = "2021-01"
                    index_cur = month_index.index(index_cur)
                    index_cur = index_overall
                if death_type == "Confirmed":
                    a = parameters[index_cur,0]
                    b = parameters[index_cur,1] 
                else:
                    a = parameters[index_cur,0]
                    b = parameters[index_cur,1]
                C = np.zeros(len(dateList))
                g = np.zeros(len(dateList))
                for k in range(max(0,j-death_cut),j):
                    C[k] = cases[k]
                    g[k] = sp.gamma.cdf(j-k,a,scale = b) - sp.gamma.cdf(j-k-1,a,scale = b)
                    if (j-k) > death_cut:
                        g[k] = 0
                adj[j] = np.dot(C,g)
                adj2[j] = deaths[j]*(sp.gamma.cdf(j-i+1,a,scale = b) - sp.gamma.cdf(j-i,a,scale = b))
                adjust[j] = adj2[j]/adj[j]
            CFR_daily_est = np.sum(adjust)
            CFR_daily[i] = CFR_daily_est
            I = np.array([cases[i]*CFR_daily_est,cases[i]])
 #             if 0 in I:
 #                 output[i] = float('nan')
 #                 output_lower[i] = float('nan')
 #                 output_upper[i] = float('nan')
 #             else:

            x0 = round(I[0])/round(I[1])
            output_data = minimize(loglikelihood, x0,method='nelder-mead',
                               options={'xatol': 1e-8, 'disp': False}); # minimize inverse loglikelihood
            output[i] = output_data.x
            x0L = 0.9*output_data.x
            x0R = 1.1*output_data.x
            output_upper_data = minimize(function_upper, x0R,method='nelder-mead',
                                       options={'xatol': 1e-8, 'disp': False}); # minimize inverse loglikelihood
            output_upper[i] = output_upper_data.x
            output_lower_data = minimize(function_lower, x0L,method='nelder-mead',
                                       options={'xatol': 1e-8, 'disp': False}); # minimize inverse loglikelihood
            output_lower[i] = output_lower_data.x
            numerator[i] = I[0]
            denominator[i] = I[1]
        dateList = range(0,len(dateList))

 #######################
       # CFR_average = CFR_daily[:].rolling(window=14).mean()
         
 #########
    elif CFR_type == "weekly":
        adj = np.zeros(len(dateList)+death_cut)
        adj2 = np.zeros(len(dateList)+death_cut)
        adjust = np.zeros(len(dateList)+death_cut)
        num_weeks = np.floor((len(dateList)-death_cut)/7).astype(int)
        output = np.zeros(num_weeks)
        output_lower = np.zeros(num_weeks)
        output_upper = np.zeros(num_weeks)
        output_25 = np.zeros(num_weeks)
        output_75 = np.zeros(num_weeks)
        numerator = np.zeros(num_weeks)
        denominator = np.zeros(num_weeks)
        a = parameters[index_overall,0]
        b = parameters[index_overall,1]
        for i in tqdm(range(0,num_weeks)):
            CFR_daily_est_temp = np.zeros(7)
            denom_temp = np.zeros(7)
            for tau in range(0,7):
                adjust = np.zeros(len(dateList)+death_cut)   
                start_point = 7*(i) + tau
                for j in range(start_point,start_point+death_cut):
                    index_cur = pd.to_datetime(dateList[j]).to_period('M').strftime("%Y-%m")
                    index_cur = month_index.index(index_cur)
                    if variable_option == 0 :
                        index_cur = index_overall
                    if np.isnan(parameters[index_cur,0]):
                        index_cur = "2021-01"
                        index_cur = month_index.index(index_cur)
                        index_cur = index_overall
                    if death_type == "Confirmed":
                        a = parameters[index_cur,0]
                        b = parameters[index_cur,1] 
                    else:
                        a = parameters[index_cur,0]
                        b = parameters[index_cur,1]
                    ### summand calculation
                    C = np.zeros(len(dateList))
                    g = np.zeros(len(dateList)) 
                    for k in range(max(0,j-death_cut),j):
                        C[k] = cases[k]
                        g[k] = sp.gamma.cdf(j-k,a,scale = b) - sp.gamma.cdf(j-k-1,a,scale = b)
                        if (j-k) > death_cut:
                            g[k] = 0
                    adj[j] = np.dot(C,g)
                    adj2[j] = cases[start_point]*deaths[j]*(sp.gamma.cdf(j-start_point+1,a,scale = b) - sp.gamma.cdf(j-start_point,a,scale = b))
                    adjust[j] = adj2[j]/adj[j]
                CFR_daily_est_temp[tau] = np.sum(adjust)
                denom_temp[tau] = cases[start_point]
            CFR_daily_est = np.sum(CFR_daily_est_temp)
            denom = np.sum(denom_temp)
            CFR_daily[i] = CFR_daily_est/denom
            I = np.array([cases[i]*CFR_daily_est/denom,cases[i]])  
            I = np.array([CFR_daily_est,denom])
            # I = np.array([CFR_daily_est/denom,cases[i]])    
            x0 = round(I[0])/round(I[1])
            output_data = minimize(loglikelihood, x0,method='nelder-mead',
                              options={'xatol': 1e-8, 'disp': False}); # minimize inverse loglikelihood
            output[i] = output_data.x
            x0L = 0.9*output_data.x
            x0R = 1.1*output_data.x
            output_upper_data = minimize(function_upper, x0R,method='nelder-mead',
                                      options={'xatol': 1e-8, 'disp': False}); # minimize inverse loglikelihood
            output_upper[i] = output_upper_data.x
            output_lower_data = minimize(function_lower, x0L,method='nelder-mead',
                                      options={'xatol': 1e-8, 'disp': False}); # minimize inverse loglikelihood
            output_lower[i] = output_lower_data.x

            output_upper_data = minimize(function_75, x0R,method='nelder-mead',
                                      options={'xatol': 1e-8, 'disp': False}); # minimize inverse loglikelihood
            output_75[i] = output_upper_data.x
            output_lower_data = minimize(function_25, x0L,method='nelder-mead',
                                      options={'xatol': 1e-8, 'disp': False}); # minimize inverse loglikelihood
            output_25[i] = output_lower_data.x
            numerator[i] = I[0]
            denominator[i] = I[1]

        #CFR_daily = deaths/adj
        dateList = range(0,num_weeks)
      #  CFR_average = CFR_daily[:].rolling(window=14).mean()
        ###########
        output_array = np.array((dateList,output_lower,output_25,output,output_75,output_upper,numerator,denominator))
        output_array = pd.DataFrame(output_array).T
 ###########
    elif  CFR_type == "overall":
       adj = np.zeros(len(dateList)+death_cut)
       adj2 = np.zeros(len(dateList)+death_cut)
       adjust = np.zeros(len(dateList)+death_cut)
       a = parameters[index_overall,0]
       b = parameters[index_overall,1]
       #for i in tqdm(range(0,len(dateList)-death_cut)):
       for i in tqdm(range(0,len(dateList) - death_cut)):
           adjust = np.zeros(len(dateList)+death_cut)  
           num_2_temp = np.zeros(len(dateList))
           for y in range(0,i+death_cut):
               C = np.zeros(len(dateList))
               g = np.zeros(len(dateList))
               num_1_temp = np.zeros(len(dateList))
               for z in range(0,i):
                   num_1_temp[z] = deaths[y]*(sp.gamma.cdf(y-z+1,a,scale = b))*cases[z]
               denom_1_temp = np.zeros(len(dateList))
               for x in range(0,y):
                   denom_1_temp[x] = cases[x]*(sp.gamma.cdf(y-x+1,a,scale = b))
               num_1 = np.nansum(num_1_temp)
               denom_1 = np.nansum(denom_1_temp)
               num_2_temp[y] = num_1/denom_1
           num_2 = np.nansum(num_2_temp)
           denom_2_temp = np.zeros(len(dateList))
           for w in range(0,i):
               denom_2_temp[w] = cases[w]
           denom_2 = np.nansum(denom_2_temp)
           CFR_daily_est = num_2/denom_2
           CFR_daily[i] = CFR_daily_est
           cases_adj = np.sum(cases[0:i+1])
           I = np.array([denom_2*CFR_daily_est,denom_2])
#             if 0 in I:
#                 output[i] = float('nan')
#                 output_lower[i] = float('nan')
#                 output_upper[i] = float('nan')
#             else:
           x0 = round(I[0])/round(I[1])
           output_data = minimize(loglikelihood, x0,method='nelder-mead',
                              options={'xatol': 1e-8, 'disp': False}); # minimize inverse loglikelihood
           output[i] = output_data.x
           x0L = 0.9*output_data.x
           x0R = 1.1*output_data.x
           output_upper_data = minimize(function_upper, x0R,method='nelder-mead',
                                      options={'xatol': 1e-8, 'disp': False}); # minimize inverse loglikelihood
           output_upper[i] = output_upper_data.x
           output_lower_data = minimize(function_lower, x0L,method='nelder-mead',
                                      options={'xatol': 1e-8, 'disp': False}); # minimize inverse loglikelihood
           output_lower[i] = output_lower_data.x
           numerator[i] = I[0]
           denominator[i] = I[1]
           
       dateList = range(0,len(dateList))        ###########        
    output_array = np.array((dateList,output_lower,output_25,output,output_75,output_upper,numerator,denominator))
    output_array = pd.DataFrame(output_array).T
    #array_test.to_csv("array_test.csv")

    return(CFR_daily,output_array)

