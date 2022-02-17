import os
import re
      
def renaming(target):
    return target.replace(',', ' ,').replace('VTK', 'COBRA').replace('Vtk', 'Cobra').replace('.regression', '_1.regression')
    
def replace_inside_file(file_path):
    with open(file_path, 'r', errors = 'ignore')  as file :
      filedata = file.read()
    filedata = renaming(filedata)
    with open(file_path, 'w') as file:
      file.write(filedata)
      
def rename_files(target_dir, option):
    for path, dirs,files in os.walk(target_dir):
        for dir in dirs:
            old_dir_path = os.path.join(path,dir)
            if '.git' in old_dir_path:
                continue
            new_dir_name = renaming(dir)
            new_dir_path = os.path.join(path,new_dir_name)
            if new_dir_name!=dir:
                os.rename(old_dir_path, new_dir_path)
            rename_files(new_dir_path,option)
            
        if option:
            continue
        
        for file in files:
            new_file = renaming(file)
            old_file_path = os.path.join(path, file)
            if 'rename.py' not in old_file_path:
                replace_inside_file(old_file_path)
                
            new_file_path = os.path.join(path, new_file)
            if old_file_path!=new_file_path:
                os.rename(old_file_path, new_file_path)
           
rename_files(".",option=True)
rename_files(".",option=False)