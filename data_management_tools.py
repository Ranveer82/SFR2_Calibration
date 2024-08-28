"""
Author Ranveer Kumar 
Latest update: 12/01/2023

Functions usefull to read information out of a GMS project and compute fitness functions for optimisation

commented
"""
import global_var
from modules import * 
import subprocess 
# from wells_clustering import *  #
import pandas as pd
import os
import flopy


# Used functions#######################
def go_to_the_root():
    # This function return to the parent location of the program in the folder architecture
    os.chdir(os.path.realpath(os.path.dirname(__file__)))
    # os.chdir("../")
    cwd = os.getcwd()
    return cwd


##### new functions  ######    

def mf_simulator():
    """This function run the Modflow model with the correct project."""
    parent_path = go_to_the_root()
    os.chdir(f"{parent_path}/{global_var.filename}_MODFLOW")
    filename_ = global_var.filename+".mfn"
    modflow_model = global_var.model+"_h5.exe"
    os.system(f"{modflow_model} {filename_}")
    go_to_the_root()

def mf_simulator2():
    """This function run the Modflow model with the correct project."""
    parent_path = go_to_the_root()
    
    modflow_nwt_path = f"{parent_path}/{global_var.filename}_MODFLOW/"+global_var.model+"_h5_dbl.exe"
    mfn_file_name = global_var.filename+".mfn"  

    # Change working directory to the folder containing MODFLOW NWT files
    os.chdir(f"{parent_path}/{global_var.filename}_MODFLOW")

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


def get_boundary_cells():
    """
    This function calulates the cell id of the boundary of the model by using High pass filter(3*3)
    dependency: from scipy import ndimage

    Argument:
    --------
    ibound = 2D array of GMS ibound dataset

    row = total number of rows in the model

    col = total number of columns in the model

    Returns:
    -------
    a list of cell ids of all boundary cells in the model

    """
    dis = get_discretization()
    row = dis['nrow']
    col = dis['ncol']
    layer = dis['nlay']

    go_to_the_root()
    boundary_cells_id = []

    with tables.open_file(global_var.environment_file_path, "r+")as f:

        for l in range(layer):
            
            ibound = f.root['Arrays/ibound'+str(l+1)]
            

            ibound2d = np.reshape(ibound,(row,col))
            ibound_index = np.arange(1,row*col+1,1).reshape(row,col)

            kernel = np.array([[-1,-1,-1],[-1,8,-1],[-1,-1,-1]])
            highpass_3 = ndimage.convolve(ibound2d, kernel)

            ## filtering the boundary cells
            for i in range(row):
                for j in range(col):
                    if highpass_3[i,j]>0:
                        boundary_cells_id.append(ibound_index[i,j]+len(ibound)*l)
                    else:
                        pass


    # ibound2d = np.reshape(ibound,(row,col))
    # ibound_index = np.arange(1,row*col+1,1).reshape(row,col)

    # kernel = np.array([[-1,-1,-1],[-1,8,-1],[-1,-1,-1]])
    # highpass_3 = ndimage.convolve(ibound2d, kernel)
   
        # print(boundary_cells_id)
    return boundary_cells_id

def get_boundary_wells():
    """
    This function compares the boundary cells with the well ids to find the list of well ids which are on the boundary 
    as SFB (specified flow boundaries).
    """
    go_to_the_root()
    b_cells = get_boundary_cells()
    with tables.open_file(global_var.environment_file_path, "r+") as f:

        well_IDs = f.root["Well"]["02. Cell IDs"]
        common_cells = np.intersect1d(b_cells,well_IDs)

    return common_cells

def get_pwell_ids():
    """
    this function will give the list of the pumping wells id in the model
    """
    go_to_the_root()
    with tables.open_file(global_var.environment_file_path, "r+") as f:

        well_IDs = f.root["Well"]["02. Cell IDs"]

        Total_BC = f.root["Well/00. Number of BCs"]
        BC_ids = well_IDs[0:Total_BC[0]]
        b_wells = get_boundary_wells()
        pumping_well_ids = well_IDs[len(b_wells):Total_BC[0]]
    
    return pumping_well_ids #, BC_ids

def get_MNW2_ids():
    """
    this function will give the list of the pumping wells id in the model
    """
    go_to_the_root()
    h5_data = tables.open_file(global_var.environment_file_path, "r+")
    well_IDs = h5_data.root["MNW2"]["02. Cell IDs"]

    Total_BC = h5_data.root["MNW2/00. Number of BCs"]
    BC_ids = well_IDs[0:Total_BC[0]]
    b_wells = get_boundary_wells()
    pumping_well_ids = well_IDs[len(b_wells):Total_BC[0]]
    h5_data.close()

    return pumping_well_ids #, BC_ids


def get_active_cell_id(layer=1):
    """
    argument: layer number

    output: the index of active cells in that layer
    """
    go_to_the_root()
    ib = "ibound"+str(layer)
    with tables.open_file(global_var.environment_file_path, "r+") as f:

        ibound = f.root["Arrays"][ib]
        active = [i+1 for i in range(len(ibound)) if ibound[i]==1]
    
    return active

def get_stream_cell_id():
    """
    this function will give the list of the stream cells id in the model
    """
    go_to_the_root()

    with tables.open_file(global_var.environment_file_path, "r+") as f:
    
        Stream_IDs = f.root["Stream"]["02. Cell IDs"]
        Total_BC = f.root["Stream/00. Number of BCs"]
        BC_ids = Stream_IDs[0:Total_BC[0]]
    
    return BC_ids

def get_drain_cell_id():
    """
    this function will give the list of the drain cells id in the model
    """
    go_to_the_root()
    h5_data = tables.open_file(global_var.environment_file_path, "r+")
    Drain_IDs = h5_data.root["Drain"]["02. Cell IDs"]

    Total_BC = h5_data.root["Drain/00. Number of BCs"]
    BC_ids = Drain_IDs[0:Total_BC[0]]
    h5_data.close()

    return BC_ids

def get_river_cell_id():
    """
    this function will give the list of the river cells id in the model
    """
    go_to_the_root()
    h5_data = tables.open_file(global_var.environment_file_path, "r+")
    Stream_IDs = h5_data.root["River"]["02. Cell IDs"]

    Total_BC = h5_data.root["River/00. Number of BCs"]
    BC_ids = Stream_IDs[0:Total_BC[0]]
    h5_data.close()
    return BC_ids


def read_arrays(array_type, layer=1):
    """ 
    This function will read the given array data from the model properties 

    Arguments: 
            array_type: 
            -----------
            # horizontal anisotropy:- HANI
            # Hydraulic conductivity:- HK
            # Specific storage: SS
            # Specific Yield: SY
            # Starting head: StartHead
            # Vertical anisotropy: VANI
            # Bottom elevation: bot
            # top elevation: top
            -------------
            row :- total number of rows in the model
            -------------
            col :- total number of columns in the model
            ------------
            layer : number of layer
    Output:
        data1d:    1D array of data (for manipulating)
        data2d:    2D array data of the array_type.(for representation)
    """
    go_to_the_root()
    dis = get_discretization()
    hdf5_file = tables.open_file(global_var.environment_file_path, "r+")
    
    row = dis['nrow']
    col = dis['ncol']
    data1d = np.array(hdf5_file.root['Arrays'][array_type+str(layer)])
    data2d = np.array(data1d).reshape(row,col)
    hdf5_file.close()
    return data1d, data2d

def write_arrays(array_w, array_type):
    """
    this function will write arrays values in the active cells only.

    argument:

    array_w = 1d array data in currect sequence to be written to the active cells. 
              (the data will be for all the cells, but only be written to the active cells)
    
    array_type = the model propertie to be written. (string)

    layer =  layer number to be updated 

    """
    
    go_to_the_root()
    arr_path = 'Arrays/'+array_type
    hdf5_file = tables.open_file(global_var.environment_file_path, "r+")

    hdf5_file.root[arr_path][:] = array_w

    # array_initial_values = hdf5_file.root[arr_path]
    # print(array_initial_values)
    # active_cells = get_active_cell_id(layer)
    # print(active_cells)
    # array_updated_values = np.array(array_initial_values)
    # for i in active_cells:
    #     array_updated_values[i] = array_w[i]

    # array_initial_values[...] = array_updated_values
    # # hdf5_file.root['Arrays'][arr_name] = array_updated_values
    print(f"Array data for {array_type} in has been updated.")
    hdf5_file.close()

    return None 

def head_DESC_index(filecontent):
    """
    This function find the index of the starting of bytes b'           HEAD' (16 char or bytes), 
    the chain after which the file is begining.

    the data structure is as follows in order:
                KSTP: the time step number, an integer, 4 bytes.
                KPER: the stress period number, an integer, 4 bytes.
                PERTIM: the time in the current stress period, a real number, either 4 or 8 bytes.
                TOTIM: the total elapsed time, a real number, either 4 or 8 bytes.
                DESC: a description of the array, 16 ANSI characters, 16 bytes. ### our identifier ##
                NCOL: the number of columns in the array, an integer, 4 bytes.
                NROW: the number of rows in the array, an integer, 4 bytes.
                ILAY: the layer number, an integer, 4 bytes.

                Next is the list of data array of size NROW*NCOL of real numbers

            Input: filecontent---> the bytes data as read by the binary file as bytes.
            output: list of indexes where we encounter DESC (the array type)
    """
    list_index = []
    for i in range(len(filecontent)):
        if filecontent[i:i+16] == b'            HEAD':
            list_index.append(i)
    print('DESC Index list created...')
    return list_index

def read_mf_barray(grid_bytes):
    """
    This function read an array produce by MODFLOW in bytes

    input: grid_bytes: the data array bytes for one stress period and for each layer of length col*row*4 (for 4 bytes data)

    returns: 1d data array 
    """
    dis = get_discretization()
    nrow = dis['nrow']
    ncol = dis['ncol']
    nlay = dis['nlay']

    if len(grid_bytes)/4 != nrow*ncol:
        print(len(grid_bytes)/4)
        print(global_var.filename)
        raise Exception(f"{global_var.filename}The number of cells in the grid does not match with the number bytes detected.")
    grid = []
    i = 0
    while i < len(grid_bytes):
        grid.append(struct.unpack('f', grid_bytes[i:i+4])[0]) # unpacking each pack of bytes and adding it to the array called grid 
        i += 4 # It is encoding on 4 bytes
    # print('*')
    return grid

def read_head():
    """
    This function will read the head out of a .hed file.
    MODFLOW 2000
    MODFLOW 2005
    MODFLOW NWT

    Output: 
    
    2D array [[head layer 1 step 1]
                [head layer 2 step 1]
                [head layer 1 step 2]
                ....................]
    
    
    3D array  [[[head, head.....ncol]
                        [..................]     ----> layer 1 time step 1
                        [......        nrow]]

                        [[head, head.....ncol]
                        [..................]     ----> layer 2 time step 1
                        [......        nrow]].... stp*ilay]

    """
    go_to_the_root()
    dis = get_discretization()
    row = dis['nrow']
    col = dis['ncol']
    layer = dis['nlay']

    file = open(global_var.heads_file_path, 'rb')
    filecontent = file.read()
    file.close()
    head_kstp_ilay_3d = []
    head_kstp_ilay_2d = []  # Head at each time step and wach layer 
    indexes = head_DESC_index(filecontent)
    for i in indexes:
        array_1d = read_mf_barray(filecontent[i+24: i+24+4*col*row])
        array_2d = np.reshape(array_1d, (row, col))
        head_kstp_ilay_3d.append(array_2d)
        head_kstp_ilay_2d.append(array_1d)
    print('Head data has been read...!!')
    return head_kstp_ilay_2d,head_kstp_ilay_3d

def drw_DESC_index(filecontent):
    """
    This function find the index of the starting of bytes b'       DRAWDOWN' (16 char or bytes), 
    the chain after which the file is begining.

    the data structure is as follows in order:
                KSTP: the time step number, an integer, 4 bytes.
                KPER: the stress period number, an integer, 4 bytes.
                PERTIM: the time in the current stress period, a real number, either 4 or 8 bytes.
                TOTIM: the total elapsed time, a real number, either 4 or 8 bytes.
                DESC: a description of the array, 16 ANSI characters, 16 bytes. ### our identifier ##
                NCOL: the number of columns in the array, an integer, 4 bytes.
                NROW: the number of rows in the array, an integer, 4 bytes.
                ILAY: the layer number, an integer, 4 bytes.

                Next is the list of data array of size NROW*NCOL of real numbers

            Input: filecontent---> the bytes data as read by the binary file as bytes.
            output: list of indexes where we encounter DESC (the array type) '        DRAWDOWN'
    """
    list_index = []
    for i in range(len(filecontent)):
        if filecontent[i:i+16] == b'        DRAWDOWN':
            list_index.append(i)
    return list_index

def read_drawdown():
    """
    This function will read the drawdown out of a .drw file.
    MODFLOW 2000
    MODFLOW 2005
    MODFLOW NWT
    Output: 3D array  [[[drawdown, drawdown.....ncol]
                        [..................]     ----> layer 1 time step 1
                        [......        nrow]]

                        [[drawdown, drawdown.....ncol]
                        [..................]     ----> layer 2 time step 1
                        [......        nrow]].... stp*ilay]

    """
    go_to_the_root()
    dis = get_discretization()
    row = dis['nrow']
    col = dis['ncol']
    layer = dis['nlay']

    file = open(global_var.drawdown_file_path, 'rb')
    filecontent = file.read()
    file.close()
    drw_kstp_ilay2d = []
    drw_kstp_ilay3d = []

    indexes = drw_DESC_index(filecontent)
    for i in indexes:
        array1d = read_mf_barray(filecontent[i+24: i+24+4*col*row])
        array2d = np.reshape(array1d, (row, col))
        drw_kstp_ilay3d.append(array2d)
        drw_kstp_ilay2d.append(array1d)


    return drw_kstp_ilay2d,drw_kstp_ilay3d


def h5_backup():
    """
    This function will create a backup of the original h5 file.

    """
    go_to_the_root()
    h5_data =  tables.open_file(global_var.environment_file_path, "r")
    backup_cell_id = np.array(h5_data.root["Well"]["02. Cell IDs"])
    backup_well_property = np.array(h5_data.root["Well/07. Property"])
    backup_no_BC = int(h5_data.root["Well/00. Number of BCs"][0])
    h5_data.close()
    return backup_cell_id,backup_well_property,backup_no_BC



def pure_h5(well_id, well_prty, no_ofBC):
    """
    This function will refesh the model to initial stage.

    """
    # well_id, well_prty, no_ofBC = h5_backup()

    ##.wel file
    fn = global_var.filename+".wel"
    file_path = os.getcwd()+"/"+global_var.filename+"_MODFLOW/"+fn
    with open(file_path, 'r') as file:
        data = file.read()
        # c_BC_no = int(data[13:17])
    data = data.replace(data[13:17],str(no_ofBC))
    with open(file_path,'w') as file:
        file.write(data)
    print('.wel file refreshed')

    ## h5 file
    h5_data =  tables.open_file(global_var.environment_file_path, "r+")
    h5_data.root["Well"]["02. Cell IDs"][...] = well_id
    h5_data.root["Well/07. Property"][...] = well_prty
    h5_data.root["Well/00. Number of BCs"][...] = no_ofBC
    h5_data.close()
    print('H5 file refreshed')

    return None



def write_well(data_array, cell_id):
    """
    This function will write the wells in the H5 input file. 

    Arguments:

    data_array---> 1d array of well discharges for  each time step (m cols).

    cell_id----> cell ids to be updated in the model as wells.

    """
    go_to_the_root()
    h5_data = tables.open_file(global_var.environment_file_path, "r+")
    well_IDs = h5_data.root["Well"]["02. Cell IDs"]
    well_property = h5_data.root["Well/07. Property"]

    time_step = len(well_property[0][0])

    if len(data_array)!=time_step:
        raise Exception("data of well discharge does match with the time steps")
    # elif len(data_array[0])!=time_step:
    #     raise Exception("data of well discharge does match with the time steps")
    else:
        pass

    # Total_BC = int(h5_data.root["Well/00. Number of BCs"][0])
    # no_of_BC =Total_BC+1

    # ## updating no of BC in .wel file
    # fn = global_var.filename+".wel"
    # file_path = os.getcwd()+"/"+global_var.filename+"_MODFLOW/"+fn
    # with open(file_path, 'r') as file:
    #     data = file.read()
    # data = data.replace(str(Total_BC),str(no_of_BC))
    # with open(file_path,'w') as file:
    #     file.write(data)
    # print('.wel file updated')

    # ## updated data in no of BC in h5
    # h5_data.root["Well/00. Number of BCs"][...] = no_of_BC

  
    # ## Updating the cell ids
    # print(f'##### total bc {Total_BC}')
    # well_IDs[Total_BC] = cell_id
    
    ## updating the well properties
    # for wells in range(len(data_array)):
    h5_data.root["Well"]["07. Property"][0, 0, :] =np.array(data_array)

    h5_data.close()
    print('Well properties updated')
    return None

def write_MNW2(data_array, cell_id, well_properties, skin_properties):
    """
    This function will write the wells in the H5 input file. 

    Arguments:

    data_array---> 2D array of well discharges for n wells (n rows). the each row will have data for one well for each time step (m cols).

    cell_id----> n number of cell ids to be updated in the model as wells.
    well_properties ;---- list of well properties as
                    [Ztop,Zbotm,ROW,COL]
    Skin_properties:  list of well and skin properties
                    [Rw,Rskin,Kskin]

    """
    go_to_the_root()
    h5_data = tables.open_file(global_var.environment_file_path, "r+")
    well_IDs = h5_data.root["MNW2"]["02. Cell IDs"]
    well_property = h5_data.root["MNW2/07. Property"]

    time_step = len(well_property[0][0])

    if len(data_array)!=len(cell_id):
        raise Exception("No of wells mismatch in data and cell id")
    elif len(data_array[0])!=time_step:
        raise Exception("data of well discharge does match with the time steps")
    else:
        pass

    Total_BC = int(h5_data.root["MNW2/00. Number of BCs"][0])
    no_of_BC =Total_BC #+len(cell_id) # uncheck this len(cell_id) to update or add new wells, 
                                        #for updating an existing well this should be off

    ## updating no of BC in .MNW2 file #####################################
    fn = global_var.filename+".mnw2"
    file_path = os.getcwd()+"/"+global_var.filename+"_MODFLOW/"+fn

    with open(file_path, 'r') as file:
        lines = file.readlines()

    # for i, line in enumerate(lines):
    #     if '2c.' in line:
    #         # Replace Rw, Rskin, Kskin values
    #         lines[i] = f"{' '.join(map(str, skin_properties))}\n"
    #     elif '2d2.' in line:
    #         # Replace Ztop, Zbotm, ROW, COL values
    #         lines[i] = f"{' '.join(map(str, well_properties))}\n"
    
    for i, line in enumerate(lines):
        if i==5:
            # Replace Rw, Rskin, Kskin values
            lines[i] = f"{' '.join(map(str, skin_properties))}\n"
        elif i==6:
            # Replace Ztop, Zbotm, ROW, COL values
            lines[i] = f"{' '.join(map(str, well_properties))}\n"

    with open(file_path, 'w') as file:
        file.writelines(lines)
    print('.MNW2 file updated')

    ########################################################################

    ## updated data in no of BC in h5
    h5_data.root["Well/00. Number of BCs"][...] = no_of_BC

  
    ## Updating the cell ids

    # well_IDs[Total_BC:no_of_BC] = cell_id   # change the index to Total_BC:no_of_BC fro adding other wells
    
    ## updating the well properties
    for wells in range(len(data_array)):
        h5_data.root["MNW2/07. Property"][1,wells,:] = data_array[wells] # replasce the index to Total_BC+wells for adding other wells

    h5_data.close()
    print('MNW2 properties updated')
    return None


def write_well_array(data_array, cell_id):
    """
    This function will write the wells in the H5 input file. 

    Arguments:

    data_array---> 2D array of well discharges for n wells (n rows). the each row will have data for one well for each time step (m cols).

    cell_id----> n number of cell ids to be updated in the model as wells.

    """
    go_to_the_root()
    h5_data = tables.open_file(global_var.environment_file_path, "r+")
    well_IDs = h5_data.root["Well"]["02. Cell IDs"]
    well_property = h5_data.root["Well/07. Property"]

    time_step = len(well_property[0][0])

    if len(data_array)!=len(cell_id):
        raise Exception("No of wells mismatch in data and cell id")
    elif len(data_array[0])!=time_step:
        raise Exception("data of well discharge does match with the time steps")
    else:
        pass

    Total_BC = int(h5_data.root["Well/00. Number of BCs"][0])
    no_of_BC =Total_BC+len(cell_id)

    ## updating no of BC in .wel file
    fn = global_var.filename+".wel"
    file_path = os.getcwd()+"/"+global_var.filename+"_MODFLOW/"+fn
    with open(file_path, 'r') as file:
        data = file.read()
    data = data.replace(str(Total_BC),str(no_of_BC))
    with open(file_path,'w') as file:
        file.write(data)
    print('.wel file updated')

    ## updated data in no of BC in h5
    h5_data.root["Well/00. Number of BCs"][...] = no_of_BC

  
    ## Updating the cell ids

    well_IDs[Total_BC:no_of_BC] = cell_id
    
    ## updating the well properties
    for wells in range(len(data_array)):
        h5_data.root["Well"]["07. Property"][0,(Total_BC+wells),:] =data_array[wells]

    h5_data.close()
    print('Well properties updated')
    return None

def MAR_output(Initial_head, well_id, threshold_change = 0.05, skip_timestep=9):
    """
    This function will calculate the average influence area after a MAR, averaged over the stress periods.
    the area will be calulated as--->  Number of cells exceding x meter change in head* the area of each cell

    Initial_head:  2D array of initial heads with each array as head data for each layer in each stp

    well_id: well id in any layer (cell_id of the well)(index)

    skip_timestep =  Time step to be skipped for average change calculations (if warmup period is used)

    Possible Outputs:
            1.  1d array of change in water table over the stress period from the initial water table (wt_change)
            2.  Average area of influence due to MAR (Inf_area_avg)
            3.  Final head -initial head for current simulation  (final_head_change)
            4.  1d array of influence area variation  (Inf_area_change)
            5.  1d array of differenece in head observed with MAR structure (head_change)
    """
    discretization = get_discretization()
    row, col = discretization['nrow'], discretization['ncol']

    layer_index = int((well_id-1)/(row*col))
    index_in_array = (well_id-1)-(layer_index*row*col)
    # print('##################')
    # print(index_in_array)
    head_array, _ = read_head()
    diff_head = np.array(head_array)-np.array(head_array[layer_index])
    # final_head_change = head_array[-1][well_id]-head_array[0][well_id]
    # wt_change = []
    final_head = []
    head_change = []
    Inf_area_change = []
    t=layer_index  # layer number to get data from
    while t < len(diff_head):
        cell_count = [c for c in diff_head[t] if c>=threshold_change]
        Inf_area_change.append(len(cell_count))
        head_change.append(diff_head[t][index_in_array])
        final_head.append(head_array[t][index_in_array])
        # wt_change.append(head_array[t]-head_array[0])

        t+=discretization['nlay']
    
    # Inf_area_avg = np.mean(np.array(Inf_area_change[skip_timestep:]))
    
    # average_wt_change = np.mean(np.array(head_change[skip_timestep:]))

    return   Inf_area_change, head_change, final_head #Inf_area_avg,  final_head_change, average_wt_change

def get_ijk(wellid):
    dis = get_discretization()
    nrow =dis['nrow']
    ncol = dis['ncol']
    lay = (wellid//(ncol*nrow))+1
    ## cell ids
    row_index = -1
    col_index = -1
    total_cells =np.arange(1, nrow*ncol+1).reshape(nrow, ncol)
    for i in range(len(total_cells)):
        if wellid in total_cells[i]:
            row_index = i+1
            col_index = list(total_cells[i]).index(wellid)
            break

    if row_index == -1 and col_index == -1:
        print("wellid not available")

    return row_index, col_index+1, lay

def head_ts(well_id):
    """


    """
    go_to_the_root()
    row, col, lay = get_ijk(well_id)
    # head_file = global_var.heads_file_path
    head_obj = flopy.utils.HeadFile(global_var.heads_file_path)
    head_ts = head_obj.get_ts((lay-1, row-1, col-1))

    return  head_ts

def MAR_storage(Initial_head, well_id):
    """
    This function will calculate the average influence area after a MAR, averaged over the stress periods.
    the area will be calulated as--->  Number of cells exceding x meter change in head* the area of each cell

    Initial_head:  2D array of initial heads with each array as head data for each layer in each stp

    well_id: well id in layer 1 (index)

    skip_timestep =  Time step to be skipped for average change calculations (if warmup period is used)

    Possible Outputs:
            1.  1d array of change in water table over the stress period from the initial water table (wt_change)
            2.  Average area of influence due to MAR (Inf_area_avg)
            3.  Final head -initial head for current simulation  (final_head_change)
            4.  1d array of influence area variation  (Inf_area_change)
            5.  1d array of differenece in head observed with MAR structure (head_change)
    """
    head_array, _ = read_head()
    diff_head = np.array(head_array)-np.array(Initial_head)
    # final_head_change = head_array[-1][well_id]-head_array[0][well_id]
    # wt_change = []
    total_stored_volume = np.sum(diff_head[-2])*(global_var.cell_size**2)
    head_change = []
    Inf_area = []
    t=0
    while t < len(diff_head):
        # cell_count = [c for c in diff_head[t] if c>=threshold_change]
        # Inf_area_change.append(len(cell_count))
        head_change.append(diff_head[t][well_id])
        # final_head.append(head_array[t][well_id])
        # wt_change.append(head_array[t]-head_array[0])

        t+=global_var.glo_n_layers
    
    # Inf_area_avg = np.mean(np.array(Inf_area_change[skip_timestep:]))
    
    # average_wt_change = np.mean(np.array(head_change[skip_timestep:]))

    return   total_stored_volume, head_change[-1] #, final_head #Inf_area_avg,  final_head_change, average_wt_change

def Initiate_model():
    """
    This function will run the model for the first time without changes and generate save the required results
    """

    go_to_the_root()
    mf_simulator()

    ### the initial values can be activated or deacativated based on use case###

    head2d, head3d = read_head()
    drawdown2d,drawdown3d = read_drawdown()

    return head2d,drawdown2d


def get_discretization():
    """
    this function will read the .dis file to get the descritzation of the model
    """
    go_to_the_root()

    fn = global_var.filename+".dis"
    file_path = os.getcwd()+"/"+global_var.filename+"_MODFLOW/"+fn
    with open(file_path, 'r') as file:
        data = file.readlines()
        # c_BC_no = int(data[13:17])
    temp = data[4].split(' ')
   
    dis_data = {'nlay':int(temp[0]), 
                'nrow':int(temp[1]),
                'ncol':int(temp[2]),
                'nper':int(temp[3]),
                't_unit':int(temp[4]),
                'l_unit':int(temp[5]),
                'cell_x_size':float(data[7].split(' ')[0]),
                'cell_y_size':float(data[9].split(' ')[0])}
    return dis_data
##############################################
############### SFR2 Calibration #############


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

    with tables.open_file(global_var.environment_file_path, "r+") as f:
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

    with tables.open_file(global_var.environment_file_path, "r+") as f:
        n_segment = f.root["Stream (SFR2)"]["13. Number of Segments"][0]
        n_reach = f.root["Stream (SFR2)"]["00. Number of BCs"][0]
    return n_segment, n_reach


def read_gage(gage):
    """
    gage: gage file names

        example: 'stream_gage_1'
    
    """
    go_to_the_root()

    time_series = pd.read_csv(os.path.join(os.getcwd(), global_var.filename+'_MODFLOW', 'model_ts.csv'))
    # print(time_series)
    file_path = os.path.join(os.getcwd(), global_var.filename+'_MODFLOW', global_var.filename+f'_{gage}.gag')

    with open(file_path, 'r+') as f:
        lines = f.readlines()

    column_names = lines[1].strip().split()[1:]
    column_names  =[cn.split('"')[0] for cn in column_names]
    
    df = pd.read_csv(file_path, sep=r'\s+', skiprows=2, header=None, names=column_names)
    df['Time'] = time_series['Time']

    return df



