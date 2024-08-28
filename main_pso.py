##################
"""
Python Modelule for river bed characteristic calibration for SFR2 package using PSO.

Author: Ranveer Kumar
IIT (BHU)

"""

# from data_management_tools import *
# Import modules
import sys
import numpy as np
import pyswarms as ps
import pickle
import subprocess 
# from wells_clustering import *  #
import pandas as pd
import os
import tables


import warnings
warnings.filterwarnings('ignore')  # Corrected to suppress warnings

###################################
### Global variables
filename = 'GM51c_SFR2'
model  ='MODFLOW-NWT'
environment_file_path = os.getcwd()+"/"+filename+"_MODFLOW/"+filename+".h5"


####################################
def go_to_the_root():
    # This function return to the parent location of the program in the folder architecture
    os.chdir(os.path.realpath(os.path.dirname(__file__)))
    # os.chdir("../")
    cwd = os.getcwd()
    return cwd

def mf_simulator2():
    """This function run the Modflow model with the correct project."""
    parent_path = go_to_the_root()
    
    modflow_nwt_path = f"{parent_path}/{filename}_MODFLOW/"+model+"_h5_dbl.exe"
    mfn_file_name = filename+".mfn"  

    # Change working directory to the folder containing MODFLOW NWT files
    os.chdir(f"{parent_path}/{filename}_MODFLOW")

    # Call MODFLOW NWT using subprocess with the MFN file name as input
    with subprocess.Popen([modflow_nwt_path, mfn_file_name], stdout=subprocess.PIPE, universal_newlines=True) as process:
        
        for line in process.stdout:
            print(line, end='')
            
            
    # print(line)
    if 'successfully' in str(line):
        sucess= 1
    else:
        sucess = 0
    go_to_the_root()

    return sucess

def modify_SFR2(parameter, new_values):
    """
    parametrs: 
    
    'ROUGHCH', 'ROUGHBK','HCOND1', 'THICKM1', 'HCOND2', 'THICKM2'

    new_values: 
    
    1d array with number of element equal to number of segments

    """
    pars = ['ROUGHCH', 'ROUGHBK','HCOND1', 'THICKM1', 'HCOND2', 'THICKM2']
    pars_index = [8, 9, 14, 15, 19, 20 ]

    go_to_the_root()

    with tables.open_file(environment_file_path, "r+") as f:
        sfr2_prop = f.root["Stream (SFR2)"]["14. Segment Property"][:,:,:]
        stress_p = sfr2_prop.shape[2]
        # check the data
        if len(new_values)!= sfr2_prop.shape[1]:
            print(f'length of value array ({len(new_values)}) is not equalt to number of segments ({sfr2_prop.shape[1]})')
            return None
        else:
            pass
        
        try:
            pidx = [pars_index[i] for i in range(len(pars)) if pars[i]==parameter][0]
        except:
            print(f' invalid parameter name : {parameter} ::: chosse from: [ROUGHCH, ROUGHBK,HCOND1, THICKM1, HCOND2, THICKM2]')

        f.root["Stream (SFR2)"]["14. Segment Property"][pidx,:,:] = np.tile(new_values, (stress_p, 1)).T

        print(f'parameter: {parameter} modified in .H5 file.')

    return None

def SFR2_info():
    go_to_the_root()

    with tables.open_file(environment_file_path, "r+") as f:
        n_segment = f.root["Stream (SFR2)"]["13. Number of Segments"][0]
        n_reach = f.root["Stream (SFR2)"]["00. Number of BCs"][0]
    return n_segment, n_reach


def read_gage(gage):
    """
    gage: gage file names

        example: 'stream_gage_1'
    
    """
    go_to_the_root()

    time_series = pd.read_csv(os.path.join(os.getcwd(), filename+'_MODFLOW', 'model_ts.csv'))
    # print(time_series)
    file_path = os.path.join(os.getcwd(), filename+'_MODFLOW', filename+f'_{gage}.gag')

    with open(file_path, 'r+') as f:
        lines = f.readlines()

    column_names = lines[1].strip().split()[1:]
    column_names  =[cn.split('"')[0] for cn in column_names]
    
    df = pd.read_csv(file_path, sep=r'\s+', skiprows=2, header=None, names=column_names)
    df['Time'] = time_series['Time']

    return df



def objective(df_obs, df_sim, obs_col):
    """Match observations to simulated data by time and calculate RMSE."""

    # print(df_obs.columns)
    # print('*'*20)
    # print(df_sim.columns)
    # print('*'*20)

    time_col_obs='Time'
    time_col_sim='Time'
    # obs_col=obs_column
    sim_col='Flow'
    
    # Check if time columns are numeric
    if pd.api.types.is_numeric_dtype(df_obs[time_col_obs]) and pd.api.types.is_numeric_dtype(df_sim[time_col_sim]):
        # Merge directly on time column for numeric time
        df_merged = pd.merge(df_obs, df_sim, left_on=time_col_obs, right_on=time_col_sim, how='inner')
    else:
        # If time columns are not numeric, convert to datetime and then use nearest match
        df_obs[time_col_obs] = pd.to_datetime(df_obs[time_col_obs])
        df_sim[time_col_sim] = pd.to_datetime(df_sim[time_col_sim])
        
        # Use nearest match for datetime data
        df_merged = pd.merge_asof(df_obs.sort_values(time_col_obs),
                                  df_sim.sort_values(time_col_sim),
                                  left_on=time_col_obs,
                                  right_on=time_col_sim,
                                  direction='nearest')
    
    # Calculate RMSE
    rmse = np.sqrt(np.nanmean((df_merged[obs_col] - df_merged[sim_col]) ** 2))
    
    return rmse  #, df_merged



def cost_function(X):
    """
    X = 2d array of paramters by pso, each element contains parameter values [A1, A2, A3....An, B1, B2, B3, ... Bn]
    
    
    """
    ###
    # pars = ['HCOND1', 'HCOND2', 'THICKM1', 'THICKM2']

    ###

    # np.save('slave_1_obj.npy', X)
    rmse = []
    for x in X:
        # get number of segments
        seg_n, _ = SFR2_info()
        
        # number of parameters
        # n_pars = int(len(x)/(seg_n))

        ## extract par values
        hcond1 =  x[:seg_n]
        hcond2 =  x[:seg_n]
        thickm1 =  x[seg_n:]
        thickm2 =  x[seg_n:]

        # write parameters value to SFR2 package
        modify_SFR2('HCOND1', hcond1)
        modify_SFR2('HCOND2', hcond2)
        modify_SFR2('THICKM1', thickm1)
        modify_SFR2('THICKM2', thickm2)

        # Simulate the model
        sucess_code = mf_simulator2()

        # read gage_file
        obs_data = pd.read_csv(os.path.join(os.getcwd(), filename+'_MODFLOW', 'observed_gage.csv'))
        # print(obs_data)
        # print('#'*20)
        gage_names  =obs_data.columns[1:]
        # print(gage_names)
        # print('#'*20)
        with open('model_conv_logs.txt', 'a') as f:
            f.write(f'pars: {x}  :: sucess: {sucess_code} \n')
    
        
        if sucess_code==1: # number of stress period
            rmse_arr = []
            for gage in gage_names:
                gage_data = read_gage(gage) ## should have a column name "Time" with time series or datetime
                obs  =obs_data[['Time', gage]] ## should have a column name "Time" with time series or datetime
                rmse_arr.append(objective(obs, gage_data, gage))
        
            rmse.append(sum(rmse_arr))
            print('Iteration RMSE: ', sum(rmse_arr))
        else:
            rmse.append(1.0e8) ## penalty for not convervesion of gw model
            print('penalty RMSE: ', 1.0e8)
        # print('Cost: ', sum(rmse_arr))

    return np.array(rmse)



# Set-up hyperparameters
options = {'c1': 0.5, 'c2': 0.5, 'w': 0.3}

# Parameter bounds
seg_n, _ = SFR2_info()  # Ensure this function is defined elsewhere
n_slave= 15


hk_lower_bound = np.ones(seg_n) * 0.01
hk_upper_bound = np.ones(seg_n) * 0.3

thickness_lower_bound = np.ones(seg_n) * 0.1
thickness_upper_bound = np.ones(seg_n) * 0.7
##################################

par_lower_bnd = np.append(hk_lower_bound, thickness_lower_bound)  # hk then thickness as per cost function
par_upper_bnd = np.append(hk_upper_bound, thickness_upper_bound)

# Call instance of PSO with convergence criteria
optimizer = ps.single.GlobalBestPSO(
    n_particles=2,
    dimensions=seg_n * 2,
    options=options,
    bounds=(par_lower_bnd, par_upper_bnd),
    #ftol = 100, # tolerance for convergence
)

# Define a stopping criterion based on the tolerance for change in the cost function value
# ftol = 100  # Stop if the change in the cost function is less than this value

# Perform optimization with the convergence criterion
cost, pos = optimizer.optimize(
    cost_function,  # Ensure 'cost_function' is defined
    iters=1 # Maximum number of iterations
           # number of parallel process
      )

# Extract the relevant data from the optimizer
optimizer_data = {
    'cost_history': optimizer.cost_history,
    'position': pos,  # Final position of the particles
    'cost': cost  # Final cost value
}

# Save the data to disk
with open('optimizer_data.pkl', 'wb') as file:
    pickle.dump(optimizer_data, file)
