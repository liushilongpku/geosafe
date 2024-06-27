"""
*******************************************************************
*******************************************************************
*******************************************************************
FileName:   GeoSafe.py
Author:     Enyi Yu
Email:      18801228090@163.com
Date:       January 20, 2024
Affiliation:Peking University
Description:
This script does determing whether the rock mass will undergo shear
and tensile failure.

Copyright (C) 2024 Enyi Yu
*******************************************************************
*******************************************************************
*******************************************************************
"""

import math
import re
import numpy as np
import vtk
import os
import glob
import numbers
from vtk.util import numpy_support
from vtk.util.numpy_support import numpy_to_vtk

# 定义一个函数，获取指定文件夹内所有的路径和文件名，并分别形成两个列表
def get_folder_contents(folder_path, paths=[], filenames=[], last_folder_name=[]):
    contents = os.listdir(folder_path)
    for content in contents:
        full_path = os.path.join(folder_path, content)
        if os.path.isdir(full_path):
            get_folder_contents(full_path, paths, filenames)
        elif full_path.endswith("vtu"):
            paths.append(full_path)
            filenames.append(content)
            last_backslash_index = full_path.rfind('\\')
            second_last_backslash_index = full_path.rfind('\\', 0, last_backslash_index)
            last_folder_name.append(full_path[second_last_backslash_index + 1: last_backslash_index])
    return paths, filenames, last_folder_name

# 假设函数用于提取时刻信息
def extract_time_from_path(path):
    # 定义一个正则表达式来匹配路径中的6位数字时间信息
    time_pattern = re.compile(r'\\(\d{6})\\')
    match = time_pattern.search(path)
    if match:
        # 如果找到匹配项，则返回匹配的时间字符串
        return match.group(1)
    else:
        # 如果没有找到匹配项，返回None
        return None
    
# 定义一个自定义排序函数
def custom_sort(item):
    if '000001' in item:
        return 0  # 如果路径中包含 '000001'，则给予最高优先级
    else:
        return 1  # 否则，优先级较低
    
def checknameindict(filename, dicname):
    if filename in dicname:
        value = dicname[filename]
        # 检查值是否为数字（整数或浮点数）
        if isinstance(value, numbers.Number):
            return value
        # 如果值是可以索引的类型（如列表或元组），返回前两个元素
        elif isinstance(value, (list, tuple)) and len(value) >= 2:
            return value[0], value[1]
        else:
            # 如果值既不是数字也不是可索引的类型，或者是可索引的类型但长度小于2
            # 这里可以根据需要返回一个特定的值或抛出异常
            return None  # 或根据需要进行其他处理
    else:
        print(filename + ' does not exist in the dictionary.')

# 计算安全系数kafang for non-fault
def caclu_safenumber_kafang(sigma_1, sigma_3, c, phi):
    kafang = 1-(sigma_1-sigma_3)/(2*c* math.cos(phi)+(sigma_1+sigma_3)*math.sin(phi))
    return kafang

# 计算破坏状态
def does_fault_been_broken(sigma_nn, tauu, c, phi):
    signalofbroken = 1 - (abs(tauu) / (abs(sigma_nn)*math.tan(phi)+c))
    return signalofbroken

# 计算斜面上的应力状态
def cauclatestressonfaces(tensor):
    n_i_times_n_j = np.outer(n, n)
    sigma_n = np.sum(tensor * n_i_times_n_j)
    T_N = np.sum(tensor * n, axis=1)
    tau_n = math.sqrt(np.sum(T_N**2) - sigma_n**2)
    return [sigma_n, tau_n]

# 计算破坏状态
def cacul_broken_out_fault(tensor, c, phi):
    sigma_oct = np.trace(tensor)/3
    tau_oct = np.sqrt((tensor[0][0]-tensor[1][1])**2+(tensor[1][1]-tensor[2][2])**2+(tensor[2][2]-tensor[0][0])**2+6*(tensor[0][1]**2+tensor[0][2]**2+tensor[1][2]**2))/3
    signalofincapandbedbroken  = does_fault_been_broken(sigma_oct, tau_oct, c, phi)
    return signalofincapandbedbroken

#计算断层之外的应力
def cacul_eight_out_fault(tensor):
    sigma_oct = np.trace(tensor)/3
    tau_oct = np.sqrt((tensor[0][0]-tensor[1][1])**2+(tensor[1][1]-tensor[2][2])**2+(tensor[2][2]-tensor[0][0])**2+6*(tensor[0][1]**2+tensor[0][2]**2+tensor[1][2]**2))/3
    return [sigma_oct, tau_oct]

def search_stress_names(vtk_file_path):
    cell_data_names = []
    # 获取并遍历单元格数据，搜索包含"stress"的数组
    cell_data = data.GetCellData()
    for i in range(cell_data.GetNumberOfArrays()):
        array_name = cell_data.GetArrayName(i)
        if "stress" in array_name.lower():  # 不区分大小写
            cell_data_names.append(array_name)  
    return cell_data_names

# 计算应力状态
def calculate_stress(tensor):
    # 计算给定张量的特征值
    eigenvalues, _ = np.linalg.eig(tensor)
    # 对特征值进行排序
    sorted_indices = np.argsort(eigenvalues)[::-1]
    eigenvalues = eigenvalues[sorted_indices]
    # 提取前三个主应力
    principal_stresses = eigenvalues[:3]
    # 计算平均应力
    mean_stress = np.mean(principal_stresses)
    # 计算最大剪切应力
    max_shear_stress = np.sqrt(((principal_stresses[0] - principal_stresses[1])**2 + 
                                (principal_stresses[1] - principal_stresses[2])**2 + 
                                (principal_stresses[2] - principal_stresses[0])**2) / 3)
    # 最大主应力、中间主应力和最小主应力
    max_principal_stress = principal_stresses[0]
    mid_principal_stress = principal_stresses[1]
    min_principal_stress = principal_stresses[2]

    return mean_stress, max_shear_stress, max_principal_stress, mid_principal_stress, min_principal_stress

def remove_paths_with_wellRegion(list_of_paths):
    # 使用列表推导式来过滤掉包含特定字符串的路径
    pathsfilter = [path for path in list_of_paths if "wellRegion" not in path]
    return pathsfilter
"""
*******************************************************************
数据输入段
*******************************************************************
"""
# input and start！！！
folder_path = r"F:\multiwellblackoil\测试用例\测试数据\vtkOutput"
paths, filenames ,last_folder_name= get_folder_contents(folder_path)

#注入速率，kg/s
q = 23.15
#储层厚度， m
h = 200 
#断层平面的法向向量 60°
#n = [0.866,0.0,0.5]
n = [0.49497996, -0.68802381, -0.530677]
#粘聚力与内摩擦角字典
cAndPhi={"segment":(5000000,math.pi/6),"fault":(500000,math.pi*10/180),"reservoir":(500000,math.pi*15/180),"domain":(5000000,math.pi*2/9),
         "cap":(1000000,math.pi*20/180),"base":(1000000,math.pi*20/180),"wellRegion1":(1000000,math.pi*20/180),"wellRegion2":(1000000,math.pi*20/180),"wellRegion3":(1000000,math.pi*20/180),"wellRegion4":(1000000,math.pi*20/180),"wellRegion22":(1000000,math.pi*20/180)}
#cAndPhi={"segment":(5000000,math.pi/6),"fault":(0,math.pi*1/18),"reservoir":(3000000,math.pi*2/9),"domain":(5000000,math.pi*2/9),
#         "cap":(5000000,math.pi/6),"base":(5000000,math.pi/6)}
#最大拉伸强度字典
X_T_d = {"segment":10000000,"fault":500000,"reservoir":1000000,"domain":1000000,"cap":1000000,"base":1000000,"wellRegion1":1000000,"wellRegion2":1000000,"wellRegion3":1000000,"wellRegion4":1000000,"wellRegion22":1000000}

# allowable maximum tensile strain字典: amts
amts_d = {"segment":0.003,"fault":0.003,"reservoir":0.003,"domain":0.003,"cap":0.003,"base":0.003,"wellRegion1":0.003,"wellRegion2":0.003,"wellRegion3":0.003,"wellRegion4":0.003,"wellRegion22":0.003}

# E, miu
EandMiu={"segment":(10000000,0.2),"fault":(10000000000,0.15),"reservoir":(10000000000,0.15),"domain":(10000000,0.2),"cap":(15000000000,0.2),"base":(15000000000,0.2),"wellRegion1":(15000000000,0.2),"wellRegion2":(15000000000,0.2),"wellRegion3":(15000000000,0.2),"wellRegion4":(15000000000,0.2),"wellRegion22":(15000000000,0.2)}

# K, fai
KandFai={"fault":(1e-15,0.2),"reservoir":(50e-15,0.25),"cap":(0.1e-15,0.1),"base":(0.1e-15,0.1),"wellRegion1":(0.1e-15,0.1),"wellRegion2":(0.1e-15,0.1),"wellRegion3":(0.1e-15,0.1),"wellRegion4":(0.1e-15,0.1),"wellRegion22":(0.1e-15,0.1)}


"""
*******************************************************************
主程序段
*******************************************************************
"""
# 在主程序段之前，对 paths, filenames, last_folder_name 进行排序
# 确保包含 '000001' 的路径排在最前面



# 使用自定义排序函数对 paths 进行排序
sorted_indices = sorted(range(len(paths)), key=lambda i: custom_sort(paths[i]))
# 根据排序后的索引重新排列 paths, filenames, last_folder_name
paths = [paths[i] for i in sorted_indices]
filenames = [filenames[i] for i in sorted_indices]
last_folder_name = [last_folder_name[i] for i in sorted_indices]

# 现在 paths 列表中包含 '000001' 的路径将排在最前面，可以按照新的顺序进行处理

paths = remove_paths_with_wellRegion(paths)

for i in range(len(paths)):
    foldername = last_folder_name[i]
    c , phi = checknameindict(foldername,cAndPhi)
    X_T = checknameindict(foldername,X_T_d)
    amts = checknameindict(foldername,amts_d)
    E , miu = checknameindict(foldername,EandMiu)
    K , fai = checknameindict(foldername,KandFai)

    time_info = extract_time_from_path(paths[i])

    reader = vtk.vtkXMLUnstructuredGridReader()
    reader.SetFileName(paths[i])
    reader.Update()
    data = reader.GetOutput()
    cell_data = data.GetCellData()

    # 创建储层数据的数组
    mean_stress_list = vtk.vtkDoubleArray()
    mean_stress_list.SetName("MeanStress_MPa")
    max_shear_stress_list = vtk.vtkDoubleArray()
    max_shear_stress_list.SetName("ShearStress_MPa")
    max_principal_stress_list = vtk.vtkDoubleArray()
    max_principal_stress_list.SetName("MaxPrincipalStress_MPa")
    mid_principal_stress_list = vtk.vtkDoubleArray()
    mid_principal_stress_list.SetName("MidPrincipalStress_MPa")
    min_principal_stress_list = vtk.vtkDoubleArray()
    min_principal_stress_list.SetName("MinPrincipalStress_MPa")
    slip_list = vtk.vtkDoubleArray()
    slip_list.SetName("sliptendency")
    eipsilon_list = vtk.vtkDoubleArray()
    eipsilon_list.SetName("MaxTensileStressCriterion")
    thelta_list = vtk.vtkDoubleArray()
    thelta_list.SetName("MaxTensileStrainCriterion")
    destructivestate_list = vtk.vtkDoubleArray()
    destructivestate_list.SetName("DestructiveState_if_not_fault")
    destructivestate_if_fault_list = vtk.vtkDoubleArray()
    destructivestate_if_fault_list.SetName("DestructiveState_if_is_fault")
    sigma_n_list = vtk.vtkDoubleArray()
    sigma_n_list.SetName("NormalStress_if_not_fault")
    sigma_n_if_fault_list = vtk.vtkDoubleArray()
    sigma_n_if_fault_list.SetName("NormalStress_if_is_fault")
    tau_n_list = vtk.vtkDoubleArray()
    tau_n_list.SetName("ShearStress_if_not_fault")
    tau_n_if_fault_list = vtk.vtkDoubleArray()
    tau_n_if_fault_list.SetName("ShearStress_if_is_fault")


    # 将pressure单位转化为MPa
    array_pressure = numpy_support.vtk_to_numpy(data.GetCellData().GetArray("pressure"))
    array_pressure = array_pressure / 1000000.0
    new_array_pressure = numpy_support.numpy_to_vtk(array_pressure, deep=True, array_type=vtk.VTK_DOUBLE)
    new_array_pressure.SetName("pressure_MPa")
    cell_data.AddArray(new_array_pressure)

    # 将deltaPressure单位转化为MPa
    array_deltaPressure = numpy_support.vtk_to_numpy(data.GetCellData().GetArray("deltaPressure"))
    array_deltaPressure = array_deltaPressure / 1000000.0
    new_array_deltaPressure = numpy_support.numpy_to_vtk(array_deltaPressure, deep=True, array_type=vtk.VTK_DOUBLE)
    new_array_deltaPressure.SetName("deltaPressure_MPa")
    cell_data.AddArray(new_array_deltaPressure)

    # 读取stress数据并计算应力
    cell_data_names = search_stress_names(paths[i])
    array_stress = numpy_support.vtk_to_numpy(data.GetCellData().GetArray(cell_data_names[0]))

    for stress_values in array_stress:    
        row1 = [stress_values[0], stress_values[5], stress_values[4]]
        row2 = [stress_values[5], stress_values[1], stress_values[3]]
        row3 = [stress_values[4], stress_values[3], stress_values[2]]
        stress_tensor = [row1, row2, row3]
        mean_stress, max_shear_stress, max_principal_stress, mid_principal_stress, min_principal_stress = calculate_stress(stress_tensor)
        # 平均有效应力、最大剪切应力、最大主应力、中间主应力和最小主应力
        mean_stress_list.InsertNextValue(mean_stress / 1000000)
        max_shear_stress_list.InsertNextValue(max_shear_stress / 1000000)
        max_principal_stress_list.InsertNextValue(max_principal_stress / 1000000)
        mid_principal_stress_list.InsertNextValue(mid_principal_stress / 1000000)
        min_principal_stress_list.InsertNextValue(min_principal_stress / 1000000)
        # 判断最大拉应力风险
        eipsilon = 1 - (max_principal_stress/X_T)
        eipsilon_list.InsertNextValue(eipsilon)
        # 判断最大拉应变风险
        #thelta = 1 - (-max_principal_stress / E + miu * (min_principal_stress + mid_principal_stress) / E ) / amts
        thelta = (-max_principal_stress / E + miu * (min_principal_stress + mid_principal_stress) / E )
        thelta_list.InsertNextValue(thelta)
        # 斜面上的破坏风险，针对断层，但实际所有单元都计算了
        stressonfaces = cauclatestressonfaces(np.array(stress_tensor))
        #计算断层之外的应力
        stressoutfault = cacul_eight_out_fault(np.array(stress_tensor))
        # 计算安全系数kafang
        slip = abs(stressonfaces[1])/abs(stressonfaces[0])
        slip_list.InsertNextValue(slip)
        """
        if foldername == "fault":
            #计算破坏状态
            destructivestate_list.InsertNextValue(does_fault_been_broken(stressonfaces[0],stressonfaces[1],c,phi))
            sigma_n_list.InsertNextValue(stressonfaces[0] / 1000000)
            tau_n_list.InsertNextValue(stressonfaces[1] / 1000000)   
        else:
            destructivestate_list.InsertNextValue(cacul_broken_out_fault(stress_tensor,c,phi))
            sigma_n_list.InsertNextValue(stressoutfault[0] / 1000000)
            tau_n_list.InsertNextValue(stressoutfault[1] / 1000000) 
        """
        #断层判断
        destructivestate_if_fault_list.InsertNextValue(does_fault_been_broken(stressonfaces[0], stressonfaces[1], c, phi))
        sigma_n_if_fault_list.InsertNextValue(stressonfaces[0] / 1000000)
        tau_n_if_fault_list.InsertNextValue(stressonfaces[1] / 1000000)
        #非断层判断
        destructivestate_list.InsertNextValue(cacul_broken_out_fault(stress_tensor, c, phi))
        sigma_n_list.InsertNextValue(stressoutfault[0] / 1000000)
        tau_n_list.InsertNextValue(stressoutfault[1] / 1000000)


        # 将计算结果添加到数据中
    cell_data.AddArray(mean_stress_list)
    cell_data.AddArray(max_shear_stress_list)
    cell_data.AddArray(max_principal_stress_list)
    cell_data.AddArray(mid_principal_stress_list)
    cell_data.AddArray(min_principal_stress_list)
    cell_data.AddArray(slip_list)
    cell_data.AddArray(eipsilon_list)
    cell_data.AddArray(thelta_list)
    cell_data.AddArray(sigma_n_list)
    cell_data.AddArray(tau_n_list)
    cell_data.AddArray(destructivestate_list)
    cell_data.AddArray(sigma_n_if_fault_list)
    cell_data.AddArray(tau_n_if_fault_list)
    cell_data.AddArray(destructivestate_if_fault_list)


    array_deltaPressure = numpy_support.vtk_to_numpy(data.GetCellData().GetArray("deltaPressure"))
    array_density = numpy_support.vtk_to_numpy(data.GetCellData().GetArray("fluid_phaseDensity"))
    array_viscosity = numpy_support.vtk_to_numpy(data.GetCellData().GetArray("fluid_phaseViscosity"))
    array_fraction = numpy_support.vtk_to_numpy(data.GetCellData().GetArray("fluid_phaseFraction"))
    array_pressure = numpy_support.vtk_to_numpy(data.GetCellData().GetArray("pressure"))
    q_v = q / array_density[:,0]
    μ = array_viscosity[:,0]
    B = 203110 * array_pressure ** (-1.03176)

    if time_info == '000001':
        # 对于0000001时刻的文件，计算并保存const
        sigma = numpy_support.vtk_to_numpy(data.GetCellData().GetArray("NormalStress_if_not_fault"))
        tau = numpy_support.vtk_to_numpy(data.GetCellData().GetArray("ShearStress_if_not_fault"))
        #计算破坏状态
        p_d = abs(sigma) - (abs(tau) / math.sin(phi) - c / (math.tan(phi)))
        #p_d = K * h * p_d / (q_v * B * μ)
        new_array_p_d = numpy_support.numpy_to_vtk(p_d, deep=True, array_type=vtk.VTK_DOUBLE)
        cell_data.AddArray(new_array_p_d)    
        
    else: 
        # 计算无量纲化的delta P   
        #p_d = K * h * array_deltaPressure / (q_v * B * μ)
        p_d = array_deltaPressure 
        new_array_p_d = numpy_support.numpy_to_vtk(p_d, deep=True, array_type=vtk.VTK_DOUBLE)
        cell_data.AddArray(new_array_p_d)

    new_array_p_d.SetName("deltaPressure_Dimensionless")

    # 保存数据到新文件路径
    writer = vtk.vtkXMLUnstructuredGridWriter()
    writer.SetFileName(paths[i])
    writer.SetInputData(data)
    writer.Write()
    print(f"文件 {paths[i]} 的计算已完成。")
