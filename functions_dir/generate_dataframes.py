def generate_dataframes(death_type,care_type,df_PHE,df_CQC,region,age_type):
    ### Generate deaths data frame from CQC notification time series
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
    #### Format PHE data time series
    print(region)
    df = df_PHE[df_PHE['PHEC_name'] == region]
    
    # df = df_PHE[df_PHE["PHEC_name"].isin(np.array(["Bolton", "Bury", "Oldham", "Rochdale", "Stockport",
    #                                              "Tameside", "Trafford", "Wigan", "Manchester","Salford"]))]
    df["date_of_death"] =  pd.to_datetime(df["date_of_death"])
    df["specimen_date"] =  pd.to_datetime(df["specimen_date"])#,format='%d/%m/%Y')
    ## look at removing dom care....

    df["death_indicator"] = (df['date_of_death'] - df['specimen_date']).dt.days
    
    if death_type == "28":
        df.loc[df["death_indicator"]<=28,"death_indicator"] = 1
        df.loc[df["death_indicator"]>28,"death_indicator"] = 0
    elif death_type == "60":
        df.loc[df["death_indicator"]<=60,"death_indicator"] = 1
        df.loc[df["death_indicator"]>60,"death_indicator"] = 0
    else:
        print("What type of death?")
        df.loc[df["death_indicator"]<=28,"death_indicator"] = 1
        df.loc[df["death_indicator"]>28,"death_indicator"] = 0
    if age_type == 1:
        age_band = range(65,1000)
    elif age_type == 2:
        age_band = range(65,75)
    elif age_type == 3:
        age_band = range(75,85)
    elif age_type == 4:
        age_band = range(85,1000)
    if care_type == "All":
        df = df[(df['age'] >= age_band[0]) & (df['age'] < age_band[-1]) & ((df['primaryservicetype_x'] == 'Care home service with nursing')
              | (df['primaryservicetype_x'] == 'Care home service without nursing'))]
        #df = df[(df['age'] >= 65)]
    elif care_type == "With":
        df = df[(df['age'] >= 65) & (df['primaryservicetype_x'] == 'Care home service with nursing')]
    elif care_type == "Without":
        df = df[(df['age'] >= 65) & (df['primaryservicetype_x'] == 'Care home service without nursing')]
    else:
        print("What type of care?")
        raise KeyboardInterrupt

    if (care_type != "All"):
        df = df[["FINALID", "date_of_death","primaryservicetype_x","specimen_date","PHEC_name","UTLA_name","death_indicator"]]
    else:
        df = df[["FINALID", "date_of_death","specimen_date","PHEC_name","UTLA_name","death_indicator"]]
        
    df = df.drop_duplicates()
    #df_processed = df

    df['Unnamed: 0'] = 1
    # cur_df = df.groupby(['date_of_death','PHEC_name'])['death_indicator'] ### PHEC region
    cur_df = df.groupby(['date_of_death'])['death_indicator'] ### PHEC region
    deaths_df = cur_df.aggregate(np.sum).reset_index()
    # deaths_df.columns = ["Date","Location","Deaths"]
    deaths_df.columns = ["Date","Deaths"]



    df['Unnamed: 0'] = 1
    # cur_df = df.groupby(['specimen_date','PHEC_name'])['Unnamed: 0'] ### PHEC region
    cur_df = df.groupby(['specimen_date'])['Unnamed: 0'] ### PHEC region
    tests_df = cur_df.aggregate(np.sum).reset_index()
    # tests_df.columns = ["Date","Location","Cases"]
    tests_df.columns = ["Date","Cases"]

    df_delay = df
    
    #### Format CQC deaths time series
    if region == "Yorkshire and Humber":
        region = "Yorkshire and The Humber"
    df = df_CQC[df_CQC['region'] == region]

    if care_type == "All":
        df = df[(df['coviddeathtype'] == 'Confirmed') & (df['primaryservicetype'] != 'Domiciliary care service')]
    elif care_type == "With":
        df = df[(df['coviddeathtype'] == 'Confirmed') & (df['primaryservicetype'] == 'Care home service with nursing')]
    elif care_type == "Without":
        df = df[(df['coviddeathtype'] == 'Confirmed') & (df['primaryservicetype'] == 'Care home service without nursing')]
    else:
        print("What type of care?")
        raise KeyboardInterrupt
        

    df['coviddeathtype'] = 1
    cur_df = df.groupby(['raiseddate','region'])['deaths'] ### use laname or region
    out_df = cur_df.aggregate(np.sum).reset_index()
    out_df.columns = ["Date","Location","Deaths"]
    deaths_df_2 = out_df 
    
    return(df_delay,tests_df,deaths_df,deaths_df_2)