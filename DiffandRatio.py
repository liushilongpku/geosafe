import vtk
from vtk.util import numpy_support
import numpy as np
import glob
import os



folder_path = r"F:\multiwellblackoil\测试用例\测试数据\vtkOutput"
# 注意可能需要更改11行和19行的路径
# 获取所有的 vtu 文件，并按照文件名排序
vtu_files_1 = sorted(glob.glob(os.path.join(folder_path, r'000001\*\Level0\*\*.vtu'), recursive=True))

# 获取所有的同级文件夹
folder_paths = sorted([os.path.join(folder_path, folder) for folder in os.listdir(folder_path) if os.path.isdir(os.path.join(folder_path, folder))])


for folder in folder_paths:
    vtu_files_2 = sorted(glob.glob(os.path.join(folder, r'*\Level0\*\*.vtu'), recursive=True))
    for vtu_file_1, vtu_file_2 in zip(vtu_files_1, vtu_files_2):
        print(vtu_file_1, vtu_file_2)
        # 读取时刻1的文件
        reader = vtk.vtkXMLUnstructuredGridReader()
        reader.SetFileName(vtu_file_1)
        reader.Update()
        data_1 = reader.GetOutput()

        # 读取时刻2的文件
        reader = vtk.vtkXMLUnstructuredGridReader()
        reader.SetFileName(vtu_file_2)
        reader.Update()
        data_2 = reader.GetOutput()

        # 对每个需要处理的数组进行处理
        for array_name in ['MeanStress_MPa', 'ShearStress_MPa', 'pressure', "deltaPressure_Dimensionless"]:
            # 如果数组不存在，跳过这个数组
            if data_1.GetCellData().GetArray(array_name) is None or data_2.GetCellData().GetArray(array_name) is None:
                continue

            array_1 = numpy_support.vtk_to_numpy(data_1.GetCellData().GetArray(array_name))
            array_2 = numpy_support.vtk_to_numpy(data_2.GetCellData().GetArray(array_name))

            # 对于pressure数组，使用不同的计算方法
            if array_name == 'pressure':
                new_array = (array_2 - array_1)/1000000
                array_ratio = (array_2 - array_1) / array_1
            elif array_name == 'deltaPressure_Dimensionless':
                new_array = 1- (array_2 / array_1)
                array_ratio = array_1
            else:
                new_array = -(array_2 - array_1)
                array_ratio = -(array_2 - array_1) / array_1

            # 创建一个新的 VTK 数组来保存新的数组值和数组比率
            new_array_vtk = numpy_support.numpy_to_vtk(new_array, deep=True, array_type=vtk.VTK_DOUBLE)
            new_array_vtk.SetName(array_name + '_diff')
            array_ratio_vtk = numpy_support.numpy_to_vtk(array_ratio, deep=True, array_type=vtk.VTK_DOUBLE)
            array_ratio_vtk.SetName(array_name + '_ratio')

            # 将新的数组值和数组比率添加到数据中
            data_2.GetCellData().AddArray(new_array_vtk)
            data_2.GetCellData().AddArray(array_ratio_vtk)


        # 写入新的 vtu 文件
        writer = vtk.vtkXMLUnstructuredGridWriter()
        writer.SetFileName(vtu_file_2) # 这里直接使用原文件名，将会覆盖原文件
        writer.SetInputData(data_2)
        writer.Write()

        print(f'文件 {vtu_file_2} 的计算和保存已完成。')



"""
import vtk
from vtk.util import numpy_support
import numpy as np
import glob
import os

# 获取所有的 vtu 文件，并按照文件名排序
vtu_files_1 = sorted(glob.glob(r'000000\mesh\Level0\*\*.vtu', recursive=True))

# 获取所有的同级文件夹
folders = sorted([folder for folder in os.listdir('.') if os.path.isdir(folder)])

# 读取文件
reader = vtk.vtkXMLUnstructuredGridReader()
for folder in folders:
    vtu_files_2 = sorted(glob.glob(os.path.join(folder, r'mesh\Level0\*\*.vtu'), recursive=True))
    for vtu_file_1, vtu_file_2 in zip(vtu_files_1, vtu_files_2):
        # 读取时刻1的文件
        reader.SetFileName(vtu_file_1)
        reader.Update()
        data_1 = reader.GetOutput()

        # 如果 'mean_stress' 数组不存在，跳过这个文件
        if data_1.GetCellData().GetArray('mean_stress') is None:
            continue

        stress_1 = numpy_support.vtk_to_numpy(data_1.GetCellData().GetArray('mean_stress'))

        # 读取时刻2的文件
        reader.SetFileName(vtu_file_2)
        reader.Update()
        data_2 = reader.GetOutput()

        # 如果 'mean_stress' 数组不存在，跳过这个文件
        if data_2.GetCellData().GetArray('mean_stress') is None:
            continue

        stress_2 = numpy_support.vtk_to_numpy(data_2.GetCellData().GetArray('mean_stress'))

        # 计算新的应力值
        new_stress = -(stress_2 - stress_1)/1000000
        stress_ratio = -(stress_2 - stress_1) / stress_1

        # 创建一个新的 VTK 数组来保存新的应力值和应力比率
        new_stress_vtk = numpy_support.numpy_to_vtk(new_stress, deep=True, array_type=vtk.VTK_FLOAT)
        new_stress_vtk.SetName('mean_stress_diff')
        stress_ratio_vtk = numpy_support.numpy_to_vtk(stress_ratio, deep=True, array_type=vtk.VTK_FLOAT)
        stress_ratio_vtk.SetName('mean_stress_ratio')

        # 将新的应力值和应力比率添加到数据中
        data_2.GetCellData().AddArray(new_stress_vtk)
        data_2.GetCellData().AddArray(stress_ratio_vtk)

        # 写入新的 vtu 文件
        writer = vtk.vtkXMLUnstructuredGridWriter()
        writer.SetFileName(vtu_file_2) # 这里直接使用原文件名，将会覆盖原文件
        writer.SetInputData(data_2)
        writer.Write()
"""