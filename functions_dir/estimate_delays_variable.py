def estimate_delays_variable(df_delay,death_type):
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
    df_delay = df_delay[(df_delay["death_indicator"] == 1)]
    df_delay = df_delay.drop_duplicates()

    nboot = 0
    if death_type == "28":
        death_cut = 28
    else:
        death_cut = 60    

    ##### define functions start 
    def loglikelihood(I,a,b,T):
        # calculates the inverse of the loglikelihood function for observing data I
        # when the distribution is truncated at time T, given shape a and scale b 
        event1=I[:,0]; # first events
        event2=I[:,1]; # second events
        L = sp.gamma.pdf(event2-event1,b,0,a)/sp.gamma.cdf(T-event1,b,0,scale = a);
        #L = sp.weibull_min.pdf(event2-event1,b,0,a)/sp.weibull_min.cdf(T-event1,b,0,scale = a);
        logL=(sum(np.log(L))**(-1)); # calculate inverse loglikelihood
        return(logL);

    def function(x):
        # turning the loglikelihood into a function of the underlying parameters
        # a and b 
        a=x[0];
        b=x[1];
        return(loglikelihood(I, a, b, T))
    #### define functions end


    #### input parameters start
    #T = 43977; # truncation date, currently based on excel date conversion
    x0=[1,1]; # initial value for minimisation
    #### input parameters end
    index = 0
    I = df_delay[["specimen_date","date_of_death"]].copy()
    I = I[["specimen_date","date_of_death"]]
    I['specimen_date'] = (I["specimen_date"] - pd.to_datetime('2020-01-01')).dt.days
    I['date_of_death'] = (I["date_of_death"] - pd.to_datetime('2020-01-01')).dt.days
    I["diff"] = I["date_of_death"] - I["specimen_date"]
    I = I[(I["diff"] >= 0) & (I["diff"] <= death_cut)]
    I.loc[I["diff"] ==0 ,"specimen_date"] = I.loc[I["diff"] == 0,"date_of_death"] - 0.5 
    I["diff"] = I["date_of_death"] - I["specimen_date"]
    I = np.array(I)

    T = I.max() + 0.5
    month_index = ["2020-01","2020-02","2020-03","2020-04","2020-05","2020-06",
                   "2020-07","2020-08","2020-09","2020-10","2020-11","2020-12",
                   "2021-01","2021-02","2021-03","2021-04","2021-05","2021-06",
                   "2021-07","2021-08","2021-09","2021-10","2021-11","2021-12",
                  "2022-01"]
    
    a_point = np.zeros(len(month_index)+1)
    b_point = np.zeros(len(month_index)+1)     
    parameters = np.zeros((len(month_index)+1,3))
    
    for i in range(0,len(month_index)):

    #### load data files start
        I = df_delay[["specimen_date","date_of_death"]].copy()
        I = I[I["specimen_date"].dt.to_period('M') == month_index[i]]

        # I = I[I["specimen_date"] > start_date_est]
        # I = I[I["specimen_date"] < end_date_est]
    
    #     I = I[I["admission_week"] == j]
        I = I[["specimen_date","date_of_death"]]
        I['specimen_date'] = (I["specimen_date"] - pd.to_datetime('2020-01-01')).dt.days
        I['date_of_death'] = (I["date_of_death"] - pd.to_datetime('2020-01-01')).dt.days
        I["diff"] = I["date_of_death"] - I["specimen_date"]
        I = I[(I["diff"] >= 0) & (I["diff"] <= death_cut)]
        I.loc[I["diff"] ==0 ,"specimen_date"] = I.loc[I["diff"] == 0,"date_of_death"] - 0.5 
        I["diff"] = I["date_of_death"] - I["specimen_date"]
        I = np.array(I)

        if len(I) < 30:
            a_point[index] = float('nan')
            b_point[index] = float('nan')
            for i in range(0,nboot):
                a_boot[index,i] = float('nan')
                b_boot[index,i] = float('nan') 
            parameters[index,0] = float("nan")
            parameters[index,1] = float("nan")
            parameters[index,2] = len(I)
            index = index + 1
        else:
            #### load data files end

    #         T = 404.5
            #### find MLE start
            output = minimize(function, x0,method='nelder-mead',
                              options={'xatol': 1e-8, 'disp': False}); # minimize inverse loglikelihood
            N = len(I)
            a = output.x[0] # MLE scale parameter
            b = output.x[1] # MLE shape parameter
            a_point[index] = a
            b_point[index] = b
            mean=a*gamma(1+1/b); # MLE mean
            variance=a**2*(gamma(1+2/b)-(gamma(1+1/b))**2); # MLE variance
            SD=math.sqrt(variance); # MLE standard deviation
    
            N = len(I)
            a_boot = np.zeros(nboot)
            b_boot = np.zeros(nboot)
    
            for i in range(0,nboot):
                vect = np.array((sp.gamma.rvs(b, scale = a, size = N)))
                col_1 = np.zeros((N,1))
                col_2 = col_1.T + vect
                I = np.concatenate((col_1,vect[:,None]),axis = 1)
                output = minimize(function, x0,method='nelder-mead'); # minimize inverse loglikelihood
                a_boot[index,i] = output.x[0] # MLE scale parameter
                b_boot[index,i] = output.x[1] # MLE shape parameter   
    
            #### find MLE end
            parameters[index,0] = b
            parameters[index,1] = a
            parameters[index,2] = N
            index = index + 1
    I = df_delay[["specimen_date","date_of_death"]].copy()

    # I = I[I["specimen_date"] > start_date_est]
    # I = I[I["specimen_date"] < end_date_est]

#     I = I[I["admission_week"] == j]
    I = I[["specimen_date","date_of_death"]]
    I['specimen_date'] = (I["specimen_date"] - pd.to_datetime('2020-01-01')).dt.days
    I['date_of_death'] = (I["date_of_death"] - pd.to_datetime('2020-01-01')).dt.days
    I["diff"] = I["date_of_death"] - I["specimen_date"]
    I = I[(I["diff"] >= 0) & (I["diff"] <= death_cut)]
    I.loc[I["diff"] ==0 ,"specimen_date"] = I.loc[I["diff"] == 0,"date_of_death"] - 0.5 
    I["diff"] = I["date_of_death"] - I["specimen_date"]
    I = np.array(I)

    if len(I) < 30:
        a_point[index] = float('nan')
        b_point[index] = float('nan')
        for i in range(0,nboot):
            a_boot[index,i] = float('nan')
            b_boot[index,i] = float('nan') 
        parameters[index,0] = float("nan")
        parameters[index,1] = float("nan")
        parameters[index,2] = len(I)
        index = index + 1
    else:
        #### load data files end

#         T = 404.5
        #### find MLE start
        output = minimize(function, x0,method='nelder-mead',
                          options={'xatol': 1e-8, 'disp': False}); # minimize inverse loglikelihood
        N = len(I)
        a = output.x[0] # MLE scale parameter
        b = output.x[1] # MLE shape parameter
        a_point[index] = a
        b_point[index] = b
        mean=a*gamma(1+1/b); # MLE mean
        variance=a**2*(gamma(1+2/b)-(gamma(1+1/b))**2); # MLE variance
        SD=math.sqrt(variance); # MLE standard deviation

        N = len(I)
        a_boot = np.zeros(nboot)
        b_boot = np.zeros(nboot)

        for i in range(0,nboot):
            vect = np.array((sp.gamma.rvs(b, scale = a, size = N)))
            col_1 = np.zeros((N,1))
            col_2 = col_1.T + vect
            I = np.concatenate((col_1,vect[:,None]),axis = 1)
            output = minimize(function, x0,method='nelder-mead'); # minimize inverse loglikelihood
            a_boot[index,i] = output.x[0] # MLE scale parameter
            b_boot[index,i] = output.x[1] # MLE shape parameter   

        #### find MLE end
        parameters[index,0] = b
        parameters[index,1] = a
        parameters[index,2] = N
    return(parameters)