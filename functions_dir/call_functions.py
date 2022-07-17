#### Call functions - no input required

def call_functions(backwards,saved,aggregated,PHE_filename,CQC_filename,care_type_filename,death_type,care_type,age_type,date_start,scenario_name_short,CFR_type,variable_option):
    # Process/load data
    import numpy as np
    import pandas as pd
    from scipy import integrate
    from scipy import stats as sp
    import matplotlib.pyplot as plt
    import math
    from scipy.optimize import minimize
    from scipy.stats import weibull_min
    from scipy.special import gamma
    from tqdm.notebook import tqdm
    from datetime import datetime, date
    from datetime import timedelta
    #### load required packages end
    
    #### Import required functions
    from functions_dir.process_data import process_data
    from functions_dir.generate_dataframes import generate_dataframes
    from functions_dir.generate_single_df import generate_single_df
    from functions_dir.estimate_delays import estimate_delays
    from functions_dir.estimate_delays_variable import estimate_delays_variable
    
    from functions_dir.estimate_CFR_forwards import estimate_CFR_forwards
    from functions_dir.estimate_CFR_backwards import estimate_CFR_backwards
    from functions_dir.initialise_variables import initialise_variables
    
    df_PHE,df_CQC,dateList = process_data(saved,aggregated,PHE_filename,CQC_filename,care_type_filename,date_start)

    ### Generate unique location ids for looping over regions of interest.
    location_ids = df_PHE['PHEC_name'].unique()
    location_ids = location_ids[~pd.isnull(location_ids)]
    forwards = 1 - backwards
    ### Loop model over regions of interest
    for i in tqdm(range(0,len(location_ids))):
        # Generate deaths and tests data frames
        df_delay,tests_df,deaths_df,deaths_df_2 = generate_dataframes(death_type,care_type,df_PHE,df_CQC,location_ids[i],age_type)

        # Merge to a single data frame
        new_df2 = generate_single_df(death_type,dateList,tests_df,deaths_df,deaths_df_2,scenario_name_short,location_ids[i])

        # Calculate test to death delay distributions
        parameters = estimate_delays_variable(df_delay,death_type)

        # Estimate CFR, using point or binomial methods - will display warning messages, can be ignored
        if backwards:
            CFR_daily,output_array = estimate_CFR_backwards(death_type,new_df2,dateList,parameters,CFR_type,variable_option)
        elif forwards:
            CFR_daily,output_array = estimate_CFR_forwards(death_type,new_df2,dateList,parameters,CFR_type,variable_option)
        else:
             print("Which approach")

        # Save output
        if aggregated == 0:
            output_array.to_csv('output_files/EN_data/output_'+scenario_name_short+'.csv',index = False)
        elif aggregated == 1:
            output_array.to_csv('output_files/region_data/output_'+scenario_name_short+location_ids[i]+'_Regions.csv',index = False)
        elif aggregated == 2:
            output_array.to_csv('output_files/la_data/output_'+scenario_name_short+location_ids[i]+'_Regions.csv',index = False)  
        else:
            print("What level of aggregation?")
    #pd.DataFrame(parameters).to_csv("carehome_death_delays"+str(age_type)+".csv")
    
    return(output_array,parameters)


