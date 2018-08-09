import os
import argparse
import subprocess
import time

def CheckRunning(shlist):
    storage="/storage/plzen1/home/holcovam/ScanTools/"
    for e in range(0, len(shlist)):
        os.chmod(shlist[e], 0o755)
        proc = ('qsub '+ storage + shlist[e])
        p1 = subprocess.Popen(proc, shell=True)
        sts1 = os.waitpid(p1.pid, 1)[1]

def FSC2(input_dir, num_reps=50, min_sims=100000, max_ecm=20, calc_CI=False, numcores=1, scratch_mb='200', time_scratch="01:50:00", mem="200", print1=False, overwrite="None", fsc2_path="/storage/plzen1/home/holcovam/programs/fsc26_linux64/fsc26"):
    """This method parallelises job submission of fastsimcoal2, but requires a very specific set up of input files.  The output of '.generateFSC2input' should be a folder that contains the multi-dimensional SFS.  Place this folder in a new folder that will be the FSC2_Data_Parent_Directory.  This directory should also contain one or more template (.tpl) and estimates (.est) files whose format can be found in the fastsimcoal2 documentation.  For each sub-directory containing input data, this method will re-format and rename the .tpl and .est files to reflect the necessary information in the sub-directory multi-dimensional SFS and then submit these jobs to the cluster.  I've tried to make the code as general as possible, but this is one method that will likely require the user to read and understand the code in order to get things working well for them.  Also, a major potential source of errors is in the correct formatting of the .tpl and .est files, so it is worthwhile to ensure that these are correct (by running FSC2 on a subset of your sub-directories) before launching full-scale"""
    Data_Files = []
    tpl_files = []
    est_files = []
    CI_Data_Files = []
    shlist = []

    if input_dir.endswith("/") is False:
        input_dir += "/"

    for path in os.listdir(input_dir):
        if os.path.isdir(input_dir + path) and path.startswith("FSC2input"):
            samp_name = path.split("_")[1]
            #folder_name = samp_name
            if samp_name + "_DSFS.obs" in os.listdir(input_dir + path):
                for i in range(0, num_reps):
                    new_file = open(input_dir + path + "/" + samp_name + str(i) + "_DSFS.obs", 'w')
                    with open(input_dir + path + "/" + samp_name + "_DSFS.obs") as data_file:
                        for line in data_file:
                            new_file.write(line)
                        new_file.close()
                    Data_Files.append(input_dir + path + "/" + samp_name + str(i) + "_DSFS.obs")
            else:
                print("Did not find input data file for: ", samp_name)
            if calc_CI == "True":
                num_files = 0
                for file in os.listdir(input_dir + path):
                    if file.endswith("_DSFS.obs") and file.split("_")[-2].split(".")[-1][0:3] == "rep" and file != samp_name + "_DSFS.obs":
                        for i in range(0, num_reps):
                            new_file = open(input_dir + path + "/" + samp_name + file.split("_")[-2].split(".")[-1].split("_")[0]+ "_" + str(i) + "_DSFS.obs", 'w')
                            with open(input_dir + path + "/" + file) as data_file:
                                for line in data_file:
                                    new_file.write(line)
                                new_file.close()
                            CI_Data_Files.append(input_dir + path + "/" + samp_name + file.split("_")[-2].split(".")[-1].split("_")[0]+ "_" + str(i) + "_DSFS.obs")
                            num_files += 1
                if len(CI_Data_Files) < 1:
                    print("Did not find bootstrap replicates for: ", samp_name)
                else:
                    print("Found ", num_files, " replicate dsfs files for CI calculation for ", samp_name)
        if path.endswith(".tpl"):
            tpl_files.append(path)
            est_files.append(path.split(".")[0])
    if len(tpl_files) == 0:
        print("Did not find any tpl files!! Aborting!!")
    else:
        if calc_CI == "True":
            Data_Files = CI_Data_Files
        for file in Data_Files:
            name = file.split("_DSFS")[0]
            samp_name = name.split("/")[-1]
            folder_name = samp_name [0:11]
            for tpl in tpl_files:
                tpl_name = tpl.split(".tpl")[0]
                if os.path.isdir(name + "_" + tpl_name) is False or overwrite == "hard":
                    new_tpl = open(name + "_" + tpl_name + ".tpl", 'w')
                    new_data = open(name + "_" + tpl_name + "_DSFS.obs", 'w')

                    with open(file, 'r') as data:
                        for i, line in enumerate(data):
                            if i == 1:
                                pop_info = line.strip("\n").strip("\t").split("\t")
                                pop_num = int(pop_info[0])
                                samp_nums = pop_info[-pop_num:]
                            new_data.write(line)
                    with open(input_dir + tpl, 'r') as template:
                        samp_num_lines = pop_num + 4
                        for i, line in enumerate(template):
                            if i < samp_num_lines:
                                new_tpl.write(line)
                            elif i == samp_num_lines:
                                for num in samp_nums:
                                    new_tpl.write(num + "\n")
                            elif i >= samp_num_lines + len(samp_nums):
                                new_tpl.write(line)
                    new_est = open(name + "_" + tpl_name + ".est", 'w')
                    try:
                        with open(input_dir + tpl_name + ".est") as est:
                            for line in est:
                                new_est.write(line)
                    except FileNotFoundError:
                        print("Did not find est file for: ", tpl)
                    #folder_name = samp_name ''.join(i for i in s if not i.isdigit())
                    shname = name + "_" + tpl_name + ".sh"
                    shfile5 = open(shname, 'w')
                    shfile5.write('#!/bin/bash -e\n' +
                                  '#PBS -N '+samp_name+'\n' +
                                  '#PBS -l walltime='+str(time_scratch)+'\n' +
                                  '#PBS -l select=1:ncpus='+str(numcores)+':mem='+str(mem)+'mb:scratch_local='+str(scratch_mb)+'mb\n' +
                                  '#PBS -m abe\n' +
                                  '#PBS -j oe\n\n' +
                                  'module add python-3.4.1-gcc\n'+
                                  'module add python34-modules-gcc\n'+
                                  'trap \'clean_scratch\' TERM EXIT\n'+
                                  'if [ ! -d "$SCRATCHDIR" ] ; then echo "Scratch not created!" 1>&2; exit 1; fi \n' +
                                  'DATADIR="/storage/plzen1/home/holcovam/ScanTools"\n' +
                                  'cp $DATADIR/'+ input_dir + "FSC2input_" + folder_name+ "/" + samp_name + "_" + tpl_name + '* $SCRATCHDIR || exit 1\n'+
                                  'cp '+fsc2_path+' $SCRATCHDIR || exit 1\n'+
                                  'cd $SCRATCHDIR || exit 2\n' +
                                  'echo data loaded at `date`\n\n' +
                                  'chmod +x fsc26 \n' +
                                  #'ls -l \n' +
                                  './fsc26 -t ' + samp_name + "_" + tpl_name + '.tpl -e ' + samp_name + "_" + tpl_name + '.est -n ' + str(min_sims) + ' -u -d -q -L ' + str(max_ecm) + ' -M \n' +                                     
                                  'rm seed.txt \n'+
                                  'rm fsc26\n'+
                                  'rm *DSFS.obs\n'+
                                  'rm *.sh\n'+
                                  'rm *.tpl \n'+
                                  'rm *.est \n'+
                                  #'ls -l \n' +
                                  'cp $SCRATCHDIR/*.par $DATADIR/'+ input_dir + "FSC2input_" + folder_name+' || exit 1\n'+
                                  'rm *.par \n'+
                                  'cp -r $SCRATCHDIR/* $DATADIR/'+input_dir+' || export CLEAN_SCRATCH=false\n'+
                                  'printf "\\nFinished\\n\\n"\n')
                    shfile5.close()
                    shlist.append(shname)

############IF PROBLEM WITH EXCESS OF NONCONVERGED CHAINS, COPY /home/majda/alpine/fastsimcoal2/afterWPSG/scripts/notConverged.py here ###################

                else:
                    print("Output for " + samp_name + "_" + tpl_name + " already exists.  Use hard_overwrite = True to overwrite.")
    return shlist

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('-i', type=str, metavar='input_directory', required=True, help='input file created with recode012.py')
    parser.add_argument('-reps', type=int, metavar='Number_Replicates', required=True, help='Number of replicates for each scenario')
    parser.add_argument('-minsims', type=int, metavar='minimum_simulations_for_esfs', required=True, help='')
    parser.add_argument('-max_ecm', type=int, metavar='number_populations', required=True, help='')
    parser.add_argument('-mb', type=int, metavar='scratch_mb', required=True, help='memory on scratch')
    parser.add_argument('-mem', type=int, metavar='amount of memory', required=True, help='integer designating amount of memory to use for each job in megabites')
    parser.add_argument('-t', type=str, metavar='time_scratch', required=True, help='')
    parser.add_argument('-nc', type=int, metavar='number_of_cores', required=True, help='number of cores to request for each job')
    parser.add_argument('-ci', type=str, metavar='Calculate_Confidence_intervals', required=True, help='Must have requested bootstrap replicate data sets from generateFSC2input in order to calculate confidence intervals.')
    parser.add_argument('-Ov', type=str, metavar='Overwrite_type', required=True, help='set to hard if you want to overwrite .bestlhoods files')
    parser.add_argument('-fsc2path', type=str, metavar='absolute_path_to_fsc2_executable', required=True, help='')
    parser.add_argument('-print1', type=str, metavar='print_shell_scripts?', required=True, help='boolean designating whether shell scripts should be submitted or printed')

    args = parser.parse_args()
    shell_script_list = FSC2(input_dir=args.i, num_reps=args.reps, min_sims=args.minsims, max_ecm=args.max_ecm, calc_CI=args.ci, numcores=args.nc, mem=args.mem, time_scratch=args.t,scratch_mb=args.mb, print1=args.print1, overwrite=args.Ov, fsc2_path=args.fsc2path)
    print(shell_script_list)
    CheckRunning(shell_script_list)


