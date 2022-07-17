def generate_single_df(death_type,dateList,tests_df,deaths_df,deaths_df_2,scenario_name_short,region):
    #### Generate continuous time series to join data frames 
    #date_start = date_end - timedelta(days = 77)
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
    # numdays = (date_end-date_start).days
    # dateList = []
    # # for x in range (0, numdays):
    # #     dateList.append(date_start + timedelta(days = x))

    if region == 1:
        region = "England"
    
    if death_type != "Confirmed":
        df1 = deaths_df
        df2 = tests_df
        df2["Date"] =  pd.to_datetime(df2["Date"])
        df1["Date"] =  pd.to_datetime(df1["Date"],format='%d/%m/%Y')

        dates_df = pd.DataFrame({"Date": dateList})

        new_df = pd.merge(dates_df,  
                             df2,  
                             on ='Date',  
                             how ='left') 
        new_df2 = pd.merge(new_df,  
                             df1,  
                             on ='Date',  
                             how ='left') 
        new_df2 = new_df2.fillna(0)
        new_df2.to_csv('output_files/region_data/Data'+scenario_name_short+region+'.csv',index = False)
    else:
        df1 = deaths_df_2
        df2 = tests_df
        df2["Date"] =  pd.to_datetime(df2["Date"])
        df1["Date"] =  pd.to_datetime(df1["Date"],format='%d/%m/%Y')

        dates_df = pd.DataFrame({"Date": dateList})
        new_df = pd.merge(dates_df,  
                             df2,  
                             on ='Date',  
                             how ='left') 
        new_df2 = pd.merge(new_df,  
                             df1,  
                             on ='Date',  
                             how ='left') 
        new_df2 = new_df2.fillna(0)
        new_df2.to_csv('output_files/region_data/Data'+scenario_name_short+region+'.csv',index = False)
    return(new_df2)
