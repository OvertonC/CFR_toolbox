#### Initialise variables and scenarios - no input required
def initialise_variables(backwards,death_type,care_type,CFR_type,age_type,variable_option,date_start):
    #### Assign method variable
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
    forwards = 1 - backwards

    if variable_option:
        variable_type = "variable_lag"
    else:
        variable_type = "constant_lag"

    if backwards:
        method = "backwards"
    elif forwards:
        method = "forwards"
    else:
        print("What method?")


    #### Define scenario names for output files
    if care_type  == "All":
        scenario_name = ": over 65 - " + death_type + " day deaths - All" + method + CFR_type + str(age_type) + str(variable_type)
        scenario_name_short = "over_65_" + death_type + "day_all" + method + CFR_type + str(age_type) + str(variable_type)
        scenario_color = 'c'
    elif care_type =="With":
        scenario_name = ": over 65 - " + death_type + " day deaths - With nursing" + method + CFR_type + str(age_type) + str(variable_type)
        scenario_name_short = "over_65_" + death_type + "day_with" + method + CFR_type + str(age_type) + str(variable_type)
        scenario_color = 'c'
    elif care_type == "Without":
        scenario_name = ": over 65 - " + death_type + " day deaths - Without nursing" + method + CFR_type + str(age_type) + str(variable_type)
        scenario_name_short = "over_65_" + death_type + "day_without" + method + CFR_type + str(age_type) + str(variable_type)
        scenario_color = 'c'
    else:
        print("What type of care?")
        raise KeyboardInterrupt


    return list((scenario_name,scenario_name_short))