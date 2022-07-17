def process_data(saved,aggregated,PHE_filename,CQC_filename,care_type_filename,date_start):
    ## Assigns type of care to each entry
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
    
    #print("loaded")
    if saved == 0:
        print("processing")
        df = pd.read_csv(PHE_filename)  
        df2 = pd.read_csv(care_type_filename)
        df2 = df2[["CQC Location (for office use only","Postcode","Service types"]]
        #df = df["coviddeathtype" == "Confirmed"]
        #df1 = df1[df1["result"]=="positive"]
        df_caretype = df2.drop_duplicates()
        df_caretype.columns = ['locationid','Postcode','primaryservicetype']

        for i in tqdm(range(0,len(df_caretype))):
            temp_df = df_caretype.iloc[i]
            if ("nursing" in temp_df['primaryservicetype']) | ("Nursing" in temp_df['primaryservicetype']):
                df_caretype.iloc[i]['primaryservicetype'] = "Care home service with nursing"
            elif ("residential" in temp_df['primaryservicetype']) | ("Residential" in temp_df['primaryservicetype']):
                df_caretype.iloc[i]['primaryservicetype'] = "Care home service without nursing"
            else:
                df_caretype.iloc[i]['primaryservicetype'] = "Other"


        #df_caretype1 = df_caretype[['locationid','primaryservicetype']]
        df_caretype1 = df_caretype[['Postcode','primaryservicetype']]

        # df_match = df1[['all_location_ids','patient_postcode']].copy()
        # df_match.columns = ['locationid','Postcode']

        df_match = df.copy()
        #df_match_final = df_match.rename(columns = {'all_location_ids':'locationid','postcode':'Postcode'}).copy()
        df_match_final = df_match.rename(columns = {'postcode':'Postcode'}).copy()
        #df_match.rename(columns = {'patient_postcode':'Postcode'})
        #df_match_final = df_match.copy()
        
        # vlookup_df = pd.merge(df_match_final,  
        #                      df_caretype1,  
        #                      on =['locationid'],  
        #                      how ='left')
        
        vlookup_df = pd.merge(df_match_final,  
                             df_caretype1,  
                             on =['Postcode'],  
                             how ='left')

        vlookup_df = pd.merge(vlookup_df,  
                             df_caretype,  
                             on =['Postcode'],  
                             how ='left')
        #vlookup_df
        output_df = vlookup_df.copy()

        for i in tqdm(range(0,len(output_df))):
            temp_df = vlookup_df.iloc[i]
            if pd.isna(temp_df["primaryservicetype_x"]):
                output_df.loc[i,'primaryservicetype_x'] = vlookup_df.loc[i,"primaryservicetype_y"]

        ########################################################################################################    

        df_PHE = output_df.copy()   
        print("saving")
        df_PHE.to_csv("processed data location")

    else:
        print("loading PHE data, please wait...")
        df_PHE = pd.read_csv('"processed data location")
    print("loading CQC data, please wait...")    
    df_CQC = pd.read_csv(CQC_filename)
    date_begin = datetime.fromisoformat('1900-01-01') ### Start date for analysis
    adjust = np.array([timedelta(days=i) for i in df_CQC["raiseddate"]])
    df_CQC["raiseddate"] = date_begin + adjust
    if aggregated == 0:
        df_PHE["PHEC_name"] = 1
        df_CQC["region"] = 1
    elif aggregated == 1:
        df_PHE["PHEC_name"] = df_PHE["PHEC_name"]
        df_CQC["region"] = df_CQC["region"]
    elif aggregated == 2:
        df_PHE["PHEC_name"] = df_PHE["UTLA_name"]
        df_CQC["region"] = df_CQC["laname"]
        
        #### Generate time series for dates
    date_end = pd.to_datetime(df_PHE["specimen_date"].max())
    numdays = (date_end-date_start).days
    dateList = []
    for x in range (0, numdays):
        dateList.append(date_start + timedelta(days = x))
    return(df_PHE,df_CQC,dateList)