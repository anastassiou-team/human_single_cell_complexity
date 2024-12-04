import sys
json_fil_name = sys.argv[1]


import glob, os
filelist = glob.glob('./Required_Files/json/*_stage_1_seed_1_hof_0_*.json')
filelist = filelist + (glob.glob('./Required_Files/json/*_stage_1_seed_2_hof_0_*.json'))


print(json_fil_name + ' / ', str(len(filelist)))

json_fil_name = filelist[int(json_fil_name)]
print(json_fil_name)

json_fil_name =[ os.path.basename(json_fil_name)]




hof_number = ['1250', '1562', '1953', '2441', '3052', '3815', '4768', '5960', 
            '7451', '9312', '11640', '14550', '18188', '22735', '28419', '35523',
            '44304', '55400', '69250', '86562', '108203', '135254', '168442', '200000']


# hof_number = ['25000','50000', '100000','200000']
core_usage = 24


commands = []
x = 0
from subprocess import Popen
for each_cell in json_fil_name:
    for each_hof in hof_number:
        commands.append("python syn_testing.py '" +each_cell + "' " + each_hof)
        x = x+1
        if x % core_usage == 0:
            procs = [ Popen(i ,cwd = "./" ,shell=True) for i in commands ]
            for p in procs:
                p.wait()
            x = 0
            commands = []

procs = [ Popen(i ,cwd = "./" ,shell=True) for i in commands ]
for p in procs:
   p.wait()