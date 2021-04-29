#!/usr/bin/env python35
import os
import subprocess
import pandas
import math
import datetime
import time
import glob


class scantools:

    def __init__(self, WorkingDir, encoding="-9"):
        if WorkingDir.endswith("/") is False:
            WorkingDir += "/"
        if os.path.exists(WorkingDir) is False:
            os.mkdir(WorkingDir)
            os.mkdir(WorkingDir + "OandE/")
        if os.path.exists(WorkingDir + "OandE/") is False:
            os.mkdir(WorkingDir + "OandE/")
        if encoding != "-9":
            POP_file = pandas.read_csv(WorkingDir + "PopKey.csv", header=0, encoding=encoding)
        else:
            POP_file = pandas.read_csv(WorkingDir + "PopKey.csv", header=0)
        POP_names = list(POP_file.columns.values)[1:]
        sample_names = list(POP_file['Samples'])
        samps = {}
        samp_nums = {}
        for pop in POP_names:
            pop_list = []
            include_index = list(POP_file[pop])
            for i, sample in enumerate(sample_names):
                if include_index[i] == 1:
                    pop_list.append(sample)
            samps[pop] = pop_list
            samp_nums[pop] = len(pop_list)


        # Determine number of individuals to downsample all populations to
        min_ind = min([sum(list(POP_file[pop])) for pop in POP_names])
        self.pops = POP_names
        self.samps = samps
        self.samp_nums = samp_nums
        self.min_ind = min_ind
        self.dir = WorkingDir
        self.oande = WorkingDir + "OandE/"
        self.code_dir = os.getcwd()
        self.log_file = open(WorkingDir + "LogFile.txt", 'a+')
        time = 'New Instance at: {:%Y-%m-%d %H:%M:%S}'.format(datetime.datetime.now()) + "\n"
        self.log_file.write(time)
        self.split_dirs = []
        for path in os.listdir(self.dir):
            if path.split("_")[0] == "VCF":
                self.split_dirs.append(self.dir + path)


    def removePop(self, pops_to_be_removed):
        '''Purpose: remove population from all object and recalculate min_ind'''
        if isinstance(pops_to_be_removed, list) is False:
            pops = [pops_to_be_removed]
        for popname in pops:
            popname = str(popname)
            if popname in self.pops:
                self.pops.remove(popname)
                self.samps.pop(popname, None)
                self.samp_nums.pop(popname, None)
                self.log_file.write("Removed Pop: " + popname + "\n")
            else:
                print("Population does not exist")
        min_ind = min([self.samp_nums[pop] for pop in self.pops])
        self.min_ind = min_ind

    def removeInds(self, ind_list):
        if isinstance(ind_list, list) is False:
            ind_list = [ind_list]
        for pop in self.samps:
            for indname in ind_list:
                indname = str(indname)
                if indname in self.samps[pop]:
                    self.samps[pop].remove(indname)
                    self.samp_nums[pop] -= 1
                    self.log_file.write("Removed Ind: " + indname + "\n")
        min_ind = min([self.samp_nums[pop] for pop in self.pops])
        self.min_ind = min_ind



    def combinePops(self, pops, popname):
        new_samps = []
        for pop in pops:
            for samp in self.samps[pop]:
                new_samps.append(samp)
        self.pops.append(popname)
        self.samps[popname] = new_samps
        self.samp_nums[popname] = len(new_samps)
        self.log_file.write("Combined Pops: " + str(pops) + " as " + popname + "\n")

    def splitVCFs(self, vcf_dir, min_dp, mffg, ref_path="/storage/pruhonice1-ibot/home/holcovam/references/lyrataV2/", ref_name="alygenomes.fasta", gatk_path="$GATK/GenomeAnalysisTK.jar", repolarization_key="repolarized.lookupKey.perSpeciesThreshold.txt", pops='all', mem=16, time_scratch='4:00:00', ncpu=4, scratch_path="$SCRATCHDIR",print1=True, overwrite=False, scratch_gb="10", keep_intermediates=False, use_scratch=True):
        '''Purpose:  Find all vcfs in vcf_dir and split them by population according to samples associated with said population.
                    Then, take only biallelic snps and convert vcf to table containing scaff, pos, ac, an, dp, and genotype fields.
                    Finally, concatenate all per-scaffold tables to one giant table. Resulting files will be put into ~/Working_Dir/VCFs/
            Notes: mffg is maximum fraction of filtered genotypes.  Number of actual genotypes allowed will be rounded up.
                    If you want to print the batch scripts, you must set print1 and overwrite to True'''

        if vcf_dir.endswith("/") is False:
            vcf_dir += "/"
        vcf_dir_name = vcf_dir.split("/")[-2]
        if use_scratch is True:
            outdir = "VCF_" + str(vcf_dir_name) + "_DP" + str(min_dp) + ".M" + str(mffg) + "/"
        else:
            outdir = self.dir + "VCF_" + str(vcf_dir_name) + "_DP" + str(min_dp) + ".M" + str(mffg) + "/"
        self.vcf_dir = vcf_dir
        if outdir not in self.split_dirs:
            self.split_dirs.append(outdir)

        ref_spec = ref_name.split(".")[0]

        if os.path.exists(outdir) is False:
            os.mkdir(outdir)
        elif overwrite is True:
            print("Overwriting files in existing VCF directory")

        if pops == 'all':
            pops = self.pops

        for pop in pops:

            if (os.path.exists(outdir + pop + '.table.repol.txt') is True or os.path.exists(outdir + pop + '.table.recode.txt') is True) and overwrite is not True:
                print("Found file for pop = " + pop + '.  Set overwrite = True to overwrite files.')
            else:
                # Add samples to list for each population according to PF file
                sample_string1 = ""
                for samp in self.samps[pop]:
                    sample_string1 += " -sn " + samp 
                joblist = []

                mfg = int(math.ceil(float(self.samp_nums[pop]) * float(mffg)))

                vcf_list = []
                vcf_basenames = []
                for file in os.listdir(vcf_dir):
                    if file[-6:] == 'vcf.gz':
                        vcf_list.append(file)
                        vcf_basenames.append(file[:-7])
                    elif file[-3:] == 'vcf':
                        vcf_list.append(file)
                        vcf_basenames.append(file[:-4])
                print(vcf_list)
                    # Select single population and biallelic SNPs for each scaffold and convert to variants table
                    
                shfile1 = open(pop + vcf_dir_name + '.sh', 'w')
                shfile1.write('#!/bin/bash\n' +
                             '#PBS -N '+pop + vcf_dir_name +'\n' +
                             '#PBS -l walltime='+time_scratch+'\n' +
                             '#PBS -l select=1:ncpus='+ncpu+':mem='+mem+'gb:scratch_local='+scratch_gb+'gb\n' +
                             '#PBS -j oe\n\n' +
                             'module add gatk-3.7 \n'+
                             'module add python-3.4.1-gcc \n'+
                             'module add parallel-20160622 \n'+
                             'trap \'clean_scratch\' TERM EXIT\n'+
                             'if [ ! -d "$SCRATCHDIR" ] ; then echo "Scratch not created!" 1>&2; exit 1; fi \n' +
                             'DATADIR="/storage/pruhonice1-ibot/home/holcovam/ScanTools"\n' +
                             'cp '+ref_path+ref_spec+'* $SCRATCHDIR || exit 1\n' +
                             'cp $DATADIR/'+ vcf_dir +'*vcf.gz* $SCRATCHDIR || exit 1\n'+
                             'cd $SCRATCHDIR || exit 2\n' +
                             'echo data all scaffolds present in the vcf_dir loaded at `date`\n' +
                             'ls *vcf.gz | parallel -j '+str(int(ncpu)-2)+' "java -Xmx' + str(int(int(mem)/(int(ncpu)-2))) + 'g -jar ' + gatk_path + ' -T SelectVariants -R ' + ref_name + ' -V {} '  + sample_string1 + ' -o {.}.' + pop + '.pop.vcf"\n' +
                             'ls *pop.vcf | parallel -j '+str(int(ncpu)-2)+' "java -Xmx' + str(int(int(mem)/(int(ncpu)-2))) + 'g -jar ' + gatk_path + ' -T VariantFiltration -R ' + ref_name + ' -V {} --genotypeFilterExpression \\"DP < ' + str(min_dp) + '\\" --genotypeFilterName \\"DP\\" -o {.}.dp1.vcf"\n') # some bug in this part, well yea " were missing?
                if keep_intermediates is False: shfile1.write('ls *pop.vcf| parallel -j '+str(int(ncpu)-2)+' "rm {} {}.idx"\n')
                shfile1.write('ls *dp1.vcf| parallel -j '+str(int(ncpu)-2)+' "java -Xmx' + str(int(int(mem)/(int(ncpu)-2))) + 'g -jar ' + gatk_path + ' -T VariantFiltration -R ' + ref_name + ' -V {} --setFilteredGtToNocall -o {.}.nc.vcf"\n')
                if keep_intermediates is False: shfile1.write('ls *dp1.vcf| parallel -j '+str(int(ncpu)-2)+' "rm {} {}.idx"\n')
                shfile1.write('ls *nc.vcf |  parallel -j '+str(int(ncpu)-2)+' "java -Xmx' + str(int(int(mem)/(int(ncpu)-2))) + 'g -jar ' + gatk_path + ' -T SelectVariants -R ' + ref_name + ' -V {} --maxNOCALLnumber ' + str(mfg) + ' -o {.}.bi.vcf"\n')
                if keep_intermediates is False: shfile1.write('ls *nc.vcf| parallel -j '+str(int(ncpu)-2)+' "rm {} {}.idx"\n')
                shfile1.write('ls *bi.vcf| parallel -j '+str(int(ncpu)-2)+' "java -Xmx' + str(int(int(mem)/(int(ncpu)-2))) + 'g -jar ' + gatk_path + ' -T VariantsToTable -R ' + ref_name + ' -V {} -F CHROM -F POS -F REF -F AN -F DP -GF GT -o {.}_raw.table"\n') # it will contain the long string have to rename it
                if keep_intermediates is False: shfile1.write('ls *bi.vcf| parallel -j '+str(int(ncpu)-2)+' "rm {} {}.idx"\n' +
                        'echo filtering done at `date` now continue with combining the scaffolds\n\n\n')
                
                shfile1.write('cp $DATADIR/recode012.py $SCRATCHDIR || exit 1\n'+
                        'cp $DATADIR/repol.py $SCRATCHDIR || exit 1\n'+
                        'cp $DATADIR/'+repolarization_key+' $SCRATCHDIR || exit 1\n'+
                        'cat *_raw.table | tail -n+2 > '+ pop +'.table\n'+
                        'python3 recode012.py -i ' + pop + '.table -pop ' + pop + ' -o $SCRATCHDIR/\n'+
                        'python3 repol.py -i ' + pop + '.table.recode.txt -o ' + pop + ' -r ' + repolarization_key + '\n')

                shfile1.write('rm '+ref_spec+'*\n'+
                            'rm *vcf.gz* \n'+
                            'rm *_raw.table\n'+
                            'rm recode012.py\n'+
                            'rm repol.py\n'+
                            'rm '+repolarization_key+'\n'+
                            'rm '+ pop + '.table\n'+
                            'cp $SCRATCHDIR/* $DATADIR/'+outdir+' || export CLEAN_SCRATCH=false\n'+
                            'printf "\\nFinished\\n\\n"\n')
                shfile1.close()

                if print1 is False:  # send job to the MetaCentrum
                   cmd1 = ('qsub ' + pop + vcf_dir_name + '.sh')
                   p1 = subprocess.Popen(cmd1, shell=True)
                   sts1 = os.waitpid(p1.pid, 0)[1]
                   joblist.append(p1.pid)

                else:
                   file1 = open(pop + vcf_dir_name + '.sh', 'r')
                   data1 = file1.read()
                   print(data1)
                   
            if print1 is False:
                self.log_file.write("###  Split VCFs  ###\n" +
                                    "VCF Directory: " + vcf_dir + "\n" +
                                    "Reference Path: " + ref_path + "\n" +
                                    "Repolarization Key: " + repolarization_key + "\n" +
                                    "Output Directory: " + outdir + "\n" +
                                    "Min Depth Per Individual: " + str(min_dp) + "\n" +
                                    "Max Fraction of Filtered Genotypes: " + str(mffg) + "\n" +
                                    "Populations: " + str(pops) + "\n")


    def splitVCFsNorepol(self, vcf_dir, min_dp, mffg, ref_path="/storage/pruhonice1-ibot/home/holcovam/references/lyrataV2/", ref_name="alygenomes.fasta", gatk_path="$GATK/GenomeAnalysisTK.jar", pops='all', mem=16, time_scratch='4:00:00', ncpu=4, scratch_path="$SCRATCHDIR",print1=True, overwrite=False, scratch_gb="10", keep_intermediates=False, use_scratch=True):
        '''Purpose:  Find all vcfs in vcf_dir and split them by population according to samples associated with said population.
                    Then, take only biallelic snps and convert vcf to table containing scaff, pos, ac, an, dp, and genotype fields.
                    Finally, concatenate all per-scaffold tables to one giant table. Resulting files will be put into ~/Working_Dir/VCFs/
            Notes: mffg is maximum fraction of filtered genotypes.  Number of actual genotypes allowed will be rounded up.
                    If you want to print the batch scripts, you must set print1 and overwrite to True'''

        if vcf_dir.endswith("/") is False:
            vcf_dir += "/"
        vcf_dir_name = vcf_dir.split("/")[-2]
        if use_scratch is True:
            outdir = "VCF_" + str(vcf_dir_name) + "_DP" + str(min_dp) + ".M" + str(mffg) + "/"
        else:
            outdir = self.dir + "VCF_" + str(vcf_dir_name) + "_DP" + str(min_dp) + ".M" + str(mffg) + "/"
        self.vcf_dir = vcf_dir
        if outdir not in self.split_dirs:
            self.split_dirs.append(outdir)

        ref_spec = ref_name.split(".")[0]

        if os.path.exists(outdir) is False:
            os.mkdir(outdir)
        elif overwrite is True:
            print("Overwriting files in existing VCF directory")

        if pops == 'all':
            pops = self.pops

        for pop in pops:

            if (os.path.exists(outdir + pop + '.table.recode.txt') is True) and overwrite is not True:
                print("Found file for pop = " + pop + '.  Set overwrite = True to overwrite files.')
            else:
                # Add samples to list for each population according to PF file
                sample_string1 = ""
                for samp in self.samps[pop]:
                    sample_string1 += " -sn " + samp 
                joblist = []

                mfg = int(math.ceil(float(self.samp_nums[pop]) * float(mffg)))

                vcf_list = []
                vcf_basenames = []
                for file in os.listdir(vcf_dir):
                    if file[-6:] == 'vcf.gz':
                        vcf_list.append(file)
                        vcf_basenames.append(file[:-7])
                    elif file[-3:] == 'vcf':
                        vcf_list.append(file)
                        vcf_basenames.append(file[:-4])
                print(vcf_list)
                    # Select single population and biallelic SNPs for each scaffold and convert to variants table
                    
                shfile1 = open(pop + vcf_dir_name + '.sh', 'w')
                shfile1.write('#!/bin/bash\n' +
                             '#PBS -N '+pop + vcf_dir_name +'\n' +
                             '#PBS -l walltime='+time_scratch+'\n' +
                             '#PBS -l select=1:ncpus='+ncpu+':mem='+mem+'gb:scratch_local='+scratch_gb+'gb\n' +
                             '#PBS -j oe\n\n' +
                             'module add gatk-3.7 \n'+
                             'module add python-3.4.1-gcc \n'+
                             'module add parallel-20160622 \n'+
                             'trap \'clean_scratch\' TERM EXIT\n'+
                             'if [ ! -d "$SCRATCHDIR" ] ; then echo "Scratch not created!" 1>&2; exit 1; fi \n' +
                             'DATADIR="/storage/pruhonice1-ibot/home/holcovam/ScanTools"\n' +
                             'cp '+ref_path+ref_spec+'* $SCRATCHDIR || exit 1\n' +
                             'cp $DATADIR/'+ vcf_dir +'*vcf.gz* $SCRATCHDIR || exit 1\n'+
                             'cd $SCRATCHDIR || exit 2\n' +
                             'echo data all scaffolds present in the vcf_dir loaded at `date`\n' +
                             'ls *vcf.gz | parallel -j '+str(int(ncpu)-2)+' "java -Xmx' + str(int(int(mem)/(int(ncpu)-2))) + 'g -jar ' + gatk_path + ' -T SelectVariants -R ' + ref_name + ' -V {} '  + sample_string1 + ' -o {.}.' + pop + '.pop.vcf"\n' +
                             'ls *pop.vcf | parallel -j '+str(int(ncpu)-2)+' "java -Xmx' + str(int(int(mem)/(int(ncpu)-2))) + 'g -jar ' + gatk_path + ' -T VariantFiltration -R ' + ref_name + ' -V {} --genotypeFilterExpression \\"DP < ' + str(min_dp) + '\\" --genotypeFilterName \\"DP\\" -o {.}.dp1.vcf"\n') # some bug in this part, well yea " were missing?
                if keep_intermediates is False: shfile1.write('ls *pop.vcf| parallel -j '+str(int(ncpu)-2)+' "rm {} {}.idx"\n')
                shfile1.write('ls *dp1.vcf| parallel -j '+str(int(ncpu)-2)+' "java -Xmx' + str(int(int(mem)/(int(ncpu)-2))) + 'g -jar ' + gatk_path + ' -T VariantFiltration -R ' + ref_name + ' -V {} --setFilteredGtToNocall -o {.}.nc.vcf"\n')
                if keep_intermediates is False: shfile1.write('ls *dp1.vcf| parallel -j '+str(int(ncpu)-2)+' "rm {} {}.idx"\n')
                shfile1.write('ls *nc.vcf |  parallel -j '+str(int(ncpu)-2)+' "java -Xmx' + str(int(int(mem)/(int(ncpu)-2))) + 'g -jar ' + gatk_path + ' -T SelectVariants -R ' + ref_name + ' -V {} --maxNOCALLnumber ' + str(mfg) + ' -o {.}.bi.vcf"\n')
                if keep_intermediates is False: shfile1.write('ls *nc.vcf| parallel -j '+str(int(ncpu)-2)+' "rm {} {}.idx"\n')
                shfile1.write('ls *bi.vcf| parallel -j '+str(int(ncpu)-2)+' "java -Xmx' + str(int(int(mem)/(int(ncpu)-2))) + 'g -jar ' + gatk_path + ' -T VariantsToTable -R ' + ref_name + ' -V {} -F CHROM -F POS -F REF -F AN -F DP -GF GT -o {.}_raw.table"\n') # it will contain the long string have to rename it
                if keep_intermediates is False: shfile1.write('ls *bi.vcf| parallel -j '+str(int(ncpu)-2)+' "rm {} {}.idx"\n' +
                        'echo filtering done at `date` now continue with combining the scaffolds\n\n\n')
                
                shfile1.write('cp $DATADIR/recode012.py $SCRATCHDIR || exit 1\n'+
                        'cat *_raw.table | tail -n+2 > '+ pop +'.table\n'+
                        'python3 recode012.py -i ' + pop + '.table -pop ' + pop + ' -o $SCRATCHDIR/\n')

                shfile1.write('rm '+ref_spec+'*\n'+
                            'rm *vcf.gz* \n'+
                            'rm *_raw.table\n'+
                            'rm recode012.py\n'+
                            'rm '+ pop + '.table\n'+
                            'cp $SCRATCHDIR/* $DATADIR/'+outdir+' || export CLEAN_SCRATCH=false\n'+
                            'printf "\\nFinished\\n\\n"\n')
                shfile1.close()

                if print1 is False:  # send job to the MetaCentrum
                   cmd1 = ('qsub ' + pop + vcf_dir_name + '.sh')
                   p1 = subprocess.Popen(cmd1, shell=True)
                   sts1 = os.waitpid(p1.pid, 0)[1]
                   joblist.append(p1.pid)

                else:
                   file1 = open(pop + vcf_dir_name + '.sh', 'r')
                   data1 = file1.read()
                   print(data1)
                   
            if print1 is False:
                self.log_file.write("###  Split VCFs  ###\n" +
                                    "VCF Directory: " + vcf_dir + "\n" +
                                    "Reference Path: " + ref_path + "\n" +
                                    "Output Directory: " + outdir + "\n" +
                                    "Min Depth Per Individual: " + str(min_dp) + "\n" +
                                    "Max Fraction of Filtered Genotypes: " + str(mffg) + "\n" +
                                    "Populations: " + str(pops) + "\n")



    def splitVCFsTreeMix(self, vcf_dir, min_dp, mffg, ref_path="/storage/pruhonice1-ibot/home/holcovam/references/lyrataV2/", ref_name="alygenomes.fasta", gatk_path="$GATK/GenomeAnalysisTK.jar", pops='all', mem=16, time_scratch='4:00:00', ncpu=4, scratch_path="$SCRATCHDIR",print1=True, overwrite=False, scratch_gb="10", keep_intermediates=False, use_scratch=True):
        '''Purpose:  Find all vcfs in vcf_dir and split them by population according to samples associated with said population.
                    Then, take only biallelic snps and convert vcf to table containing scaff, pos, ac, an, dp, and genotype fields.
                    Finally, concatenate all per-scaffold tables to one giant table. Resulting files will be put into ~/Working_Dir/VCFs/
            Notes: mffg is maximum fraction of filtered genotypes.  Number of actual genotypes allowed will be rounded up.
                    If you want to print the batch scripts, you must set print1 and overwrite to True'''

        if vcf_dir.endswith("/") is False:
            vcf_dir += "/"
        vcf_dir_name = vcf_dir.split("/")[-2]
        if use_scratch is True:
            outdir = "VCF_" + str(vcf_dir_name) + "_DP" + str(min_dp) + ".M" + str(mffg) + "/"
        else:
            outdir = self.dir + "VCF_" + str(vcf_dir_name) + "_DP" + str(min_dp) + ".M" + str(mffg) + "/"
        self.vcf_dir = vcf_dir
        if outdir not in self.split_dirs:
            self.split_dirs.append(outdir)

        ref_spec = ref_name.split(".")[0]

        if os.path.exists(outdir) is False:
            os.mkdir(outdir)
        elif overwrite is True:
            print("Overwriting files in existing VCF directory")

        if pops == 'all':
            pops = self.pops

        for pop in pops:

            if (os.path.exists(outdir + pop + '.table.repol.txt') is True or os.path.exists(outdir + pop + '.table.recode.txt') is True) and overwrite is not True:
                print("Found file for pop = " + pop + '.  Set overwrite = True to overwrite files.')
            else:
                # Add samples to list for each population according to PF file
                sample_string1 = ""
                for samp in self.samps[pop]:
                    sample_string1 += " -sn " + samp 
                joblist = []

                mfg = int(math.ceil(float(self.samp_nums[pop]) * float(mffg)))

                vcf_list = []
                vcf_basenames = []
                for file in os.listdir(vcf_dir):
                    if file[-6:] == 'vcf.gz':
                        vcf_list.append(file)
                        vcf_basenames.append(file[:-7])
                    elif file[-3:] == 'vcf':
                        vcf_list.append(file)
                        vcf_basenames.append(file[:-4])
                print(vcf_list)
                    # Select single population and biallelic SNPs for each scaffold and convert to variants table
                    
                shfile1 = open(pop + vcf_dir_name + '.sh', 'w')
                shfile1.write('#!/bin/bash\n' +
                             '#PBS -N '+pop + vcf_dir_name +'\n' +
                             '#PBS -l walltime='+time_scratch+'\n' +
                             '#PBS -l select=1:ncpus='+ncpu+':mem='+mem+'gb:scratch_local='+scratch_gb+'gb\n' +
                             '#PBS -j oe\n\n' +
                             'module add gatk-3.7 \n'+
                             'module add python-3.4.1-gcc \n'+
                             'module add parallel-20160622 \n'+
                             'trap \'clean_scratch\' TERM EXIT\n'+
                             'if [ ! -d "$SCRATCHDIR" ] ; then echo "Scratch not created!" 1>&2; exit 1; fi \n' +
                             'DATADIR="/storage/pruhonice1-ibot/home/holcovam/ScanTools"\n' +
                             'cp '+ref_path+ref_spec+'* $SCRATCHDIR || exit 1\n' +
                             'cp $DATADIR/'+ vcf_dir +'*vcf.gz* $SCRATCHDIR || exit 1\n'+
                             'cd $SCRATCHDIR || exit 2\n' +
                             'echo data all scaffolds present in the vcf_dir loaded at `date`\n' +
                             'ls *vcf.gz | parallel -j '+str(int(ncpu)-2)+' "java -Xmx' + str(int(int(mem)/(int(ncpu)-2))) + 'g -jar ' + gatk_path + ' -T SelectVariants -R ' + ref_name + ' -V {} '  + sample_string1 + ' -o {.}.' + pop + '.pop.vcf"\n'+
                             'ls *pop.vcf | parallel -j '+str(int(ncpu)-2)+' "java -Xmx' + str(int(int(mem)/(int(ncpu)-2))) + 'g -jar ' + gatk_path + ' -T VariantFiltration -R ' + ref_name + ' -V {} --genotypeFilterExpression \\"DP < ' + str(min_dp) + '\\" --genotypeFilterName \\"DP\\" -o {.}.dp1.vcf"\n') # some bug in this part, well yea " were missing?
                if keep_intermediates is False: shfile1.write('ls *pop.vcf| parallel -j '+str(int(ncpu)-2)+' "rm {} {}.idx"\n')
                shfile1.write('ls *dp1.vcf| parallel -j '+str(int(ncpu)-2)+' "java -Xmx' + str(int(int(mem)/(int(ncpu)-2))) + 'g -jar ' + gatk_path + ' -T VariantFiltration -R ' + ref_name + ' -V {} --setFilteredGtToNocall -o {.}.nc.vcf"\n')
                if keep_intermediates is False: shfile1.write('ls *dp1.vcf| parallel -j '+str(int(ncpu)-2)+' "rm {} {}.idx"\n')
                shfile1.write('ls *nc.vcf |  parallel -j '+str(int(ncpu)-2)+' "java -Xmx' + str(int(int(mem)/(int(ncpu)-2))) + 'g -jar ' + gatk_path + ' -T SelectVariants -R ' + ref_name + ' -V {} --maxNOCALLnumber ' + str(mfg) + ' -o {.}.bi.vcf"\n')
                if keep_intermediates is False: shfile1.write('ls *nc.vcf| parallel -j '+str(int(ncpu)-2)+' "rm {} {}.idx"\n')
                shfile1.write('ls *bi.vcf| parallel -j '+str(int(ncpu)-2)+' "java -Xmx' + str(int(int(mem)/(int(ncpu)-2))) + 'g -jar ' + gatk_path + ' -T VariantsToTable -R ' + ref_name + ' -V {} -F CHROM -F POS -F AC -F AN -o {.}_raw.table"\n') # it will contain the long string have to rename it
                if keep_intermediates is False: shfile1.write('ls *bi.vcf| parallel -j '+str(int(ncpu)-2)+' "rm {} {}.idx"\n' +
                        'echo filtering done at `date` now continue with combining the scaffolds\n\n\n')
                
                shfile1.write('cat *_raw.table | grep -v "CHROM" > '+ pop +'_tm.table\n')
                #        'cp $DATADIR/repol.py $SCRATCHDIR || exit 1\n'+
                #        'cp $DATADIR/'+repolarization_key+' $SCRATCHDIR || exit 1\n'+
                 #       'cat *_raw.table | tail -n+2 > '+ pop +'.table\n'+
                #        'python3 recode012.py -i ' + pop + '.table -pop ' + pop + ' -o $SCRATCHDIR/\n'+
                #        'python3 repol.py -i ' + pop + '.table.recode.txt -o ' + pop + ' -r ' + repolarization_key + '\n')

                shfile1.write('rm '+ref_spec+'*\n'+
                          #  'rm *vcf.gz* \n'+
                            'rm *_raw.table\n'+
                          #  'rm recode012.py\n'+
                        #    'rm repol.py\n'+
                      #      'rm '+repolarization_key+'\n'+
                        #    'rm '+ pop + '.table\n'+
                            'cp $SCRATCHDIR/* $DATADIR/'+outdir+' || export CLEAN_SCRATCH=false\n'+
                            'printf "\\nFinished\\n\\n"\n')
                shfile1.close()

                if print1 is False:  # send job to the MetaCentrum
                   cmd1 = ('qsub ' + pop + vcf_dir_name + '.sh')
                   p1 = subprocess.Popen(cmd1, shell=True)
                   sts1 = os.waitpid(p1.pid, 0)[1]
                   joblist.append(p1.pid)

                else:
                   file1 = open(pop + vcf_dir_name + '.sh', 'r')
                   data1 = file1.read()
                   print(data1)
                   
            if print1 is False:
                self.log_file.write("###  Split VCFs  ###\n" +
                                    "VCF Directory: " + vcf_dir + "\n" +
                                    "Reference Path: " + ref_path + "\n" +
                                #    "Repolarization Key: " + repolarization_key + "\n" +
                                    "Output Directory: " + outdir + "\n" +
                                    "Min Depth Per Individual: " + str(min_dp) + "\n" +
                                    "Max Fraction of Filtered Genotypes: " + str(mffg) + "\n" +
                                    "Populations: " + str(pops) + "\n")
                                    
#MAJDA: protein evolution
    def splitVCFsAnn(self, vcf_dir, min_dp, mffg, ref_path="/storage/pruhonice1-ibot/home/holcovam/references/lyrataV2/", ref_name="alygenomes.fasta", gatk_path="$GATK/GenomeAnalysisTK.jar", repolarization_key="repolarized.lookupKey.perSpeciesThreshold.txt", pops='all', mem=16, time_scratch='4:00:00', ncpu=4, scratch_path="$SCRATCHDIR",print1=True, overwrite=False, scratch_gb="10", keep_intermediates=False, use_scratch=True):
        '''Purpose:  Find all vcfs in vcf_dir and split them by population according to samples associated with said population.
                    Then, take only biallelic snps and convert vcf to table containing scaff, pos, ac, an, dp, and genotype fields.
                    Finally, concatenate all per-scaffold tables to one giant table. Resulting files will be put into ~/Working_Dir/VCFs/
            Notes: mffg is maximum fraction of filtered genotypes.  Number of actual genotypes allowed will be rounded up.
                    If you want to print the batch scripts, you must set print1 and overwrite to True'''

        if vcf_dir.endswith("/") is False:
            vcf_dir += "/"
        vcf_dir_name = vcf_dir.split("/")[-2]
        if use_scratch is True:
            outdir = "VCF_" + str(vcf_dir_name) + "_DP" + str(min_dp) + ".M" + str(mffg) + "/"
        else:
            outdir = self.dir + "VCF_" + str(vcf_dir_name) + "_DP" + str(min_dp) + ".M" + str(mffg) + "/"
        self.vcf_dir = vcf_dir
        if outdir not in self.split_dirs:
            self.split_dirs.append(outdir)

        ref_spec = ref_name.split(".")[0]

        if os.path.exists(outdir) is False:
            os.mkdir(outdir)
        elif overwrite is True:
            print("Overwriting files in existing VCF directory")

        if pops == 'all':
            pops = self.pops
            print(pops)

        for pop in pops:

            if (os.path.exists(outdir + pop + '.table.repol.txt') is True or os.path.exists(outdir + pop + '.table.recode.txt') is True) and overwrite is not True:
                print("Found file for pop = " + pop + '.  Set overwrite = True to overwrite files.')
            else:
                # Add samples to list for each population according to PF file
                sample_string1 = ""
                for samp in self.samps[pop]:
                    sample_string1 += " -sn " + samp 
                joblist = []

                mfg = int(math.ceil(float(self.samp_nums[pop]) * float(mffg)))

                vcf_list = []
                vcf_basenames = []
                for file in os.listdir(vcf_dir):
                    if file[-6:] == 'vcf.gz':
                        vcf_list.append(file)
                        vcf_basenames.append(file[:-7])
                    elif file[-3:] == 'vcf':
                        vcf_list.append(file)
                        vcf_basenames.append(file[:-4])
                print(vcf_list)         
                # Select single population and biallelic SNPs for each scaffold and convert to variants table
                
                shfile1 = open(pop + vcf_dir_name + '.sh', 'w')
                shfile1.write('#!/bin/bash\n' +
                             '#PBS -N '+pop + vcf_dir_name +'\n' +
                             '#PBS -l walltime='+time_scratch+'\n' +
                             '#PBS -l select=1:ncpus='+ncpu+':mem='+mem+'gb:scratch_local='+scratch_gb+'gb\n' +
                             '#PBS -j oe\n\n' +
                             'module add gatk-3.7 \n'+
                             'module add python-3.4.1-gcc \n'+
                             'module add parallel-20160622 \n'+
                             'trap \'clean_scratch\' TERM EXIT\n'+
                             'if [ ! -d "$SCRATCHDIR" ] ; then echo "Scratch not created!" 1>&2; exit 1; fi \n' +
                             'DATADIR="/storage/pruhonice1-ibot/home/holcovam/ScanTools"\n' +
                             'cp '+ref_path+ref_spec+'* $SCRATCHDIR || exit 1\n' +
                             'cp $DATADIR/'+ vcf_dir +'*vcf* $SCRATCHDIR || exit 1\n'+###!!!!!!!!!!   add .gz
                             'cd $SCRATCHDIR || exit 2\n' +
                             'echo data all scaffolds present in the vcf_dir loaded at `date`\n' + ##ADD *VCF.GZ!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                             'ls *.vcf* | parallel -j '+str(int(ncpu)-2)+' "java -Xmx' + str(int(int(mem)/(int(ncpu)-2))) + 'g -jar ' + gatk_path + ' -T SelectVariants -R ' + ref_name + ' -V {} '  + sample_string1 + ' -o {.}.' + pop + '.pop.vcf"\n' +
                             'ls *pop.vcf | parallel -j '+str(int(ncpu)-2)+' "java -Xmx' + str(int(int(mem)/(int(ncpu)-2))) + 'g -jar ' + gatk_path + ' -T VariantFiltration -R ' + ref_name + ' -V {} --genotypeFilterExpression \\"DP < ' + str(min_dp) + '\\" --genotypeFilterName \\"DP\\" -o {.}.dp1.vcf"\n') # some bug in this part, well yea " were missing?
                if keep_intermediates is False: shfile1.write('ls *pop.vcf| parallel -j '+str(int(ncpu)-2)+' "rm {} {}.idx"\n')
                shfile1.write('ls *dp1.vcf| parallel -j '+str(int(ncpu)-2)+' "java -Xmx' + str(int(int(mem)/(int(ncpu)-2))) + 'g -jar ' + gatk_path + ' -T VariantFiltration -R ' + ref_name + ' -V {} --setFilteredGtToNocall -o {.}.nc.vcf"\n')
                if keep_intermediates is False: shfile1.write('ls *dp1.vcf| parallel -j '+str(int(ncpu)-2)+' "rm {} {}.idx"\n')
                shfile1.write('ls *nc.vcf |  parallel -j '+str(int(ncpu)-2)+' "java -Xmx' + str(int(int(mem)/(int(ncpu)-2))) + 'g -jar ' + gatk_path + ' -T SelectVariants -R ' + ref_name + ' -V {} --maxNOCALLnumber ' + str(mfg) + ' -o {.}.bi.vcf"\n')
                if keep_intermediates is False: shfile1.write('ls *nc.vcf| parallel -j '+str(int(ncpu)-2)+' "rm {} {}.idx"\n')
                shfile1.write('ls *bi.vcf| parallel -j '+str(int(ncpu)-2)+' "java -Xmx' + str(int(int(mem)/(int(ncpu)-2))) + 'g -jar ' + gatk_path + ' -T VariantsToTable -R ' + ref_name + ' -V {} -F CHROM -F POS -F REF -F AN -F DP -GF GT -F ANN -o {.}_raw.table"\n') # it will contain the long string have to rename it
                if keep_intermediates is False: shfile1.write('ls *bi.vcf| parallel -j '+str(int(ncpu)-2)+' "rm {} {}.idx"\n' +
                        'echo filtering done at `date` now continue with combining the scaffolds\n\n\n')

                shfile1.write('cp $DATADIR/recode012ann.py $SCRATCHDIR || exit 1\n'+
                        'cp $DATADIR/repolann.py $SCRATCHDIR || exit 1\n'+
                        'cp $DATADIR/'+repolarization_key+' $SCRATCHDIR || exit 1\n'+
                        'cat *_raw.table | tail -n+2 > '+ pop +'.table\n'+
                      #  'python3 recode012ann.py -i ' + pop + '.table -pop ' + pop + ' -o $SCRATCHDIR/\n'+ ###You might want to uncomment this!!!
                        'python3 recode012ann.py -i ' + pop + '.table -pop ' + pop + ' -o $SCRATCHDIR/\n')

                    #    'python3 repolann.py -i ' + pop + '.table.recode.txt -o ' + pop + ' -r ' + repolarization_key + '\n') ###You might want to uncomment this!!!


                shfile1.write('rm '+ref_spec+'*\n'+
                            'rm *.vcf* \n'+ ####ADD .gz
                            'rm *_raw.table\n'+
                            'rm recode012ann.py\n'+
                            'rm repolann.py\n'+
                            'rm '+repolarization_key+'\n'+
                            'rm '+ pop + '.table\n'+
                            'cp $SCRATCHDIR/* $DATADIR/'+outdir+' || export CLEAN_SCRATCH=false\n'+
                            'printf "\\nFinished\\n\\n"\n')
                shfile1.close()

                if print1 is False:  # send job to the MetaCentrum
                   cmd1 = ('qsub ' + pop + vcf_dir_name + '.sh')
                   p1 = subprocess.Popen(cmd1, shell=True)
                   sts1 = os.waitpid(p1.pid, 0)[1]
                   joblist.append(p1.pid)

                else:
                   file1 = open(pop + vcf_dir_name + '.sh', 'r')
                   data1 = file1.read()
                   print(data1)
                   
            if print1 is False:
                self.log_file.write("###  Split VCFs  ###\n" +
                                    "VCF Directory: " + vcf_dir + "\n" +
                                    "Reference Path: " + ref_path + "\n" +
                                    "Repolarization Key: " + repolarization_key + "\n" +
                                    "Output Directory: " + outdir + "\n" +
                                    "Min Depth Per Individual: " + str(min_dp) + "\n" +
                                    "Max Fraction of Filtered Genotypes: " + str(mffg) + "\n" +
                                    "Populations: " + str(pops) + "\n")
                   
    def getPloidies(self, recode_dir, use_repol=True):
        '''Purpose: Create new methods of scantools object containing ploidy of each population (.ploidies) as well as a list of dips (.dips) and tetraploid populations (.tets)
           Notes: Can only be executed after recode has been executed on vcfs'''

        print("Be sure that 'recode' scripts have all finished")
        if recode_dir.endswith("/") is False:
            recode_dir += "/"
        ploidies = {}
        dips = []
        tets = []
        if os.path.exists(recode_dir) is True:
            for pop in self.pops:
                try:
                    if use_repol is True:
                        tmp = open(recode_dir + pop + '.table.repol.txt', 'r')
                    else:    
                        tmp = open(recode_dir + pop + '.table.recode.txt', 'r')
                    line = tmp.readline()
                    ploidy = line.split("\t")[1]
                    ploidies[pop] = ploidy
                    if ploidy == "4.0":
                        tets.append(pop)
                    elif ploidy == "2.0":
                        dips.append(pop)
                    else:
                        print("Ploidy level not recognized")
                except (FileNotFoundError, IndexError):
                    print("Error determining ploidy for population: ", pop)
            self.ploidies = ploidies
            self.dips = dips
            self.tets = tets
        else:
            print("recode_dir does not exist")


    def calcwpm(self, recode_dir, window_size, min_snps, pops="all", print1=False, mem=16, ncpu=1, sampind="-99", scratch_gb=2, use_repol=True, time_scratch="4:00:00", overwrite=False):
        '''Purpose: Calculate within population metrics including: allele frequency, expected heterozygosity, Wattersons theta, Pi, ThetaH, ThetaL and neutrality tests: D, normalized H, E
           Notes:  Currently, all populations are downsampled to same number of individuals.  By default, this minimum individuals across populations minus 1 to allow for some missing data
                    It is worth considering whether downsampling should be based on number of individuals or number of alleles.
                    Results are held ~/Working_Dir/Recoded/ in series of files ending in _WPM.txt.  These can be concatenated using concatWPM'''

        if use_repol is True:
            suffix = '.table.repol.txt'
        else:
            suffix = '.table.recode.txt'

        if sampind == "-99":
            sind = self.min_ind - 1
        else:
            sind = sampind

        if pops == "all":
            pops = self.pops
        if recode_dir.endswith("/") is False:
            recode_dir += "/"

        dir_name = recode_dir.split("/")[-2]

        if os.path.exists(recode_dir) is True and sind > 3: 

            for pop in pops:

                prefix = pop + ".WS" + str(window_size / 1000) + "k_MS" + str(min_snps) + "_" + str(sind) + "ind"
                if os.path.exists(recode_dir + prefix + '_WPM.txt') is True and overwrite is not True:
                    print("Output file for pop " + pop + ' already exists.  Set overwrite = True to overwrite.  Aborting.')
                else:
                    if os.path.exists(recode_dir + pop + suffix) is True:
                        shfile3 = open(dir_name + "." + pop + '.sh', 'w')

                        shfile3.write('#!/bin/bash -e\n' +
                                      '#PBS -N ' +pop+ '\n' +
                                      '#PBS -l walltime='+str(time_scratch)+'\n' +
                                      '#PBS -l select=1:ncpus='+str(ncpu)+':mem='+str(mem)+'gb:scratch_local='+str(scratch_gb)+'gb\n' +
                                      '#PBS -m abe\n' +
                                      '#PBS -j oe\n\n' +
                                      'module add python-3.4.1-gcc\n'+
                                      'module add python34-modules-gcc\n'+
                                      'trap \'clean_scratch\' TERM EXIT\n'+
                                      'if [ ! -d "$SCRATCHDIR" ] ; then echo "Scratch not created!" 1>&2; exit 1; fi \n' +
                                      'DATADIR="/storage/pruhonice1-ibot/home/holcovam/ScanTools"\n' +
                                      'cp $DATADIR/wpm.py $SCRATCHDIR || exit 1\n'+
                                      'cp $DATADIR/'+ recode_dir + pop + suffix +' $SCRATCHDIR || exit 1\n'+
                                      'cd $SCRATCHDIR || exit 2\n' +
                                      'echo data loaded at `date`\n' +
                                      'python3 wpm.py -i ' + pop + suffix + ' -o $SCRATCHDIR/ -p ' + prefix + ' -sampind ' + str(sind) + ' -ws ' + str(window_size) + ' -ms ' + str(min_snps) + '\n'+
                                      'rm wpm.py\n'+
                                      'rm ' + pop + suffix + '\n'+
                                      'cp $SCRATCHDIR/* $DATADIR/'+recode_dir+' || export CLEAN_SCRATCH=false\n'+
                                      'printf "\\nFinished\\n\\n"\n')
                                      

                        shfile3.close()

                        if print1 is False:
                            cmd3 = ('qsub ' + dir_name + "." + pop + '.sh')
                            p3 = subprocess.Popen(cmd3, shell=True)
                            sts3 = os.waitpid(p3.pid, 0)[1]

                        else:
                            file3 = open(dir_name + "." + pop + '.sh', 'r')
                            data3 = file3.read()
                            print(data3)

                        os.remove(dir_name + "." + pop + '.sh')
                    else:
                        print("Did not find input files for: ", pop)

            if print1 is False:
                self.log_file.write("###  Calculate Within-Population-Metrics  ###\n" +
                                    "Input Directory: " + recode_dir + "\n" +
                                    "Window Size: " + str(window_size) + "\n" +
                                    "Minimum SNPs in a window: " + str(min_snps) + "\n" +
                                    "Number of individuals to be downsampled to: " + str(sampind) + "\n" +
                                    "Use repolarized data: " + str(use_repol) + "\n" +
                                    "Populations: " + str(pops) + "\n")

        elif sind <= 3:
            print("Number of individuals to be used/downsampled to is <= 3.  Unable to calculate within-population-metrics on so few individuals.")
        else:
            print("Did not find recode_dir.  Must run splitVCFs followed by recode before able to calculate within population metrics")

    def calcbpm(self, recode_dir, pops, output_name, window_size, min_snps, mem=16, ncpu=1, use_repol=True, keep_intermediates=False, time_scratch="1:00:00",scratch_gb=2, print1=False):
        '''Purpose:  Calculate between population metrics including: Dxy, Fst (using Weir and Cockerham 1984), and Rho (Ronfort et al. 1998)
           Notes: User provides a list of populations to be included.  For pairwise estimates, simply provide two populations
                    Calculations are done for windows of a given bp size.  User also must specify the minimum number of snps in a window
                    for calculations to be made'''

        if use_repol is True:
            suffix = '.table.repol.txt'
        else:
            suffix = '.table.recode.txt'

        if recode_dir.endswith("/") is False:
            recode_dir += "/"

        output_name += "_WS" + str(window_size) + "_MS" + str(min_snps)
        if os.path.exists(recode_dir) is True and len(pops) > 1:
            pop_num = 0
            file_string = ""
            file_string_noRecode = ""
            for pop in pops:
                try:
                    a = open(recode_dir + pop + suffix, 'r')
                    a.close()
                    file_string += "$DATADIR/"+ recode_dir + pop + suffix + " "
                    file_string_noRecode += pop + suffix + " "
                    pop_num += 1
                except IOError:
                    print("Did not find input file for pop ", pop)
            if len(pops) != pop_num:
                print("Did not find all input files!!  Aborting.")
                os.remove(recode_dir + output_name + '.concat.txt')
            else:
                shfile3 = open(recode_dir + output_name + '.bpm.sh', 'w')

                shfile3.write('#!/bin/bash -e\n' +
                              '#PBS -N '+output_name +'\n' +
                              '#PBS -l walltime='+str(time_scratch)+'\n' +
                              '#PBS -l select=1:ncpus='+str(ncpu)+':mem='+str(mem)+'gb:scratch_local='+str(scratch_gb)+'gb\n' +
                              '#PBS -m abe\n' +
                              '#PBS -j oe\n\n' +
                              'module add python-3.4.1-gcc\n'+
                              'module add python34-modules-gcc\n'+
                              'trap \'clean_scratch\' TERM EXIT\n'+
                              'if [ ! -d "$SCRATCHDIR" ] ; then echo "Scratch not created!" 1>&2; exit 1; fi \n' +
                              'DATADIR="/storage/pruhonice1-ibot/home/holcovam/ScanTools"\n' +
                              'cp $DATADIR/bpm.py $SCRATCHDIR || exit 1\n'+
                              'cp '+ file_string.split(" ")[0] +' $SCRATCHDIR || exit 1\n'+
                              'cp '+ file_string.split(" ")[1] +' $SCRATCHDIR || exit 1\n'+
                              'cd $SCRATCHDIR || exit 2\n' +
                              'echo data loaded at `date`\n' +
                              'sort -k3,3 -k4,4n -m ' + file_string_noRecode + '> ' + output_name + '.concat.txt\n' +
                              'python3 bpm.py -i ' + output_name + '.concat.txt -o $SCRATCHDIR/ -prefix ' + output_name + ' -ws ' + str(window_size) + ' -ms ' + str(min_snps) + ' -np ' + str(pop_num) + '\n'+
                              'rm bpm.py\n'+
                              'rm ' + file_string_noRecode.split(" ")[0] + '\n'+
                              'rm ' + file_string_noRecode.split(" ")[1] + '\n')

                if keep_intermediates is False:
                    shfile3.write('rm ' + output_name + '.concat.txt\n')
                    
                shfile3.write('cp $SCRATCHDIR/* $DATADIR/'+recode_dir+' || export CLEAN_SCRATCH=false\n'+
                              'printf "\\nFinished\\n\\n"\n')
                shfile3.close()

                if print1 is False:
                    cmd3 = ('qsub ' + recode_dir + output_name + '.bpm.sh')
                    p3 = subprocess.Popen(cmd3, shell=True)
                    sts3 = os.waitpid(p3.pid, 0)[1]

                    self.log_file.write("###  Calculate Between-Population-Metrics  ###\n" +
                                        "Input Directory: " + recode_dir + "\n" +
                                        "Window Size: " + str(window_size) + "\n" +
                                        "Minimum SNPs in a window: " + str(min_snps) + "\n" +
                                        "Use repolarized data: " + str(use_repol) + "\n" +
                                        "Populations: " + str(pops) + "\n")

                else:
                    file3 = open(recode_dir + output_name + '.bpm.sh', 'r')
                    data3 = file3.read()
                    print(data3)

                os.remove(recode_dir + output_name + '.bpm.sh')
        elif len(pops) < 2:
            print("'pops' argument must be a list of strings specifiying two or more population names as they appear in input file prefixes.  len(pops) was < 2")
        else:
            print("Did not find recode_dir.  Must run splitVCFs followed by recode before able to calculate between population metrics")

#MAJDA
#calculate Fst,... output it with alcode, ann, aas
#maybe add 4d info? - i.e. like repol.py?? 
    def calcbpmAnn(self, recode_dir, pops, output_name, window_size, min_snps, mem=16, ncpu=1, use_repol=True, keep_intermediates=False, time_scratch="1:00:00",scratch_gb=2, print1=False):
        '''Purpose:  Calculate between population metrics including: Dxy, Fst (using Weir and Cockerham 1984), and Rho (Ronfort et al. 1998)
           Notes: User provides a list of populations to be included.  For pairwise estimates, simply provide two populations
                    Calculations are done for windows of a given bp size.  User also must specify the minimum number of snps in a window
                    for calculations to be made'''

        if use_repol is True:
            suffix = '.table.repol.txt'
        else:
            suffix = '.table.recode.txt'

        if recode_dir.endswith("/") is False:
            recode_dir += "/"

        output_name += "_WS" + str(window_size) + "_MS" + str(min_snps)
        if os.path.exists(recode_dir) is True and len(pops) > 1:
            pop_num = 0
            file_string = ""
            file_string_noRecode = ""
            for pop in pops:
                try:
                    a = open(recode_dir + pop + suffix, 'r')
                    a.close()
                    file_string += "$DATADIR/"+ recode_dir + pop + suffix + " "
                    file_string_noRecode += pop + suffix + " "
                    pop_num += 1
                except IOError:
                    print("Did not find input file for pop ", pop)
            if len(pops) != pop_num:
                print("Did not find all input files!!  Aborting.")
                os.remove(recode_dir + output_name + '.concat.txt')
            else:
                shfile3 = open(recode_dir + output_name + '.bpm.sh', 'w')

                shfile3.write('#!/bin/bash -e\n' +
                              '#PBS -N '+output_name +'\n' +
                              '#PBS -l walltime='+str(time_scratch)+'\n' +
                              '#PBS -l select=1:ncpus='+str(ncpu)+':mem='+str(mem)+'gb:scratch_local='+str(scratch_gb)+'gb\n' +
                              '#PBS -m abe\n' +
                              '#PBS -j oe\n\n' +
                              'module add python-3.4.1-gcc\n'+
                              'module add python34-modules-gcc\n'+
                              'trap \'clean_scratch\' TERM EXIT\n'+
                              'if [ ! -d "$SCRATCHDIR" ] ; then echo "Scratch not created!" 1>&2; exit 1; fi \n' +
                              'DATADIR="/storage/pruhonice1-ibot/home/holcovam/ScanTools"\n' +
                              'cp $DATADIR/bpmann.py $SCRATCHDIR || exit 1\n'+
                              'cp '+ file_string.split(" ")[0] +' $SCRATCHDIR || exit 1\n'+
                              'cp '+ file_string.split(" ")[1] +' $SCRATCHDIR || exit 1\n'+
                              'cd $SCRATCHDIR || exit 2\n' +
                              'echo data loaded at `date`\n' +
                              'sort -k3,3 -k4,4n -m ' + file_string_noRecode + '> ' + output_name + '.concat.txt\n' +
                              'python3 bpmann.py -i ' + output_name + '.concat.txt -o $SCRATCHDIR/ -prefix ' + output_name + ' -ws ' + str(window_size) + ' -ms ' + str(min_snps) + ' -np ' + str(pop_num) + '\n'+
                              'rm bpmann.py\n'+
                              'rm ' + file_string_noRecode.split(" ")[0] + '\n'+
                              'rm ' + file_string_noRecode.split(" ")[1] + '\n')

                if keep_intermediates is False:
                    shfile3.write('rm ' + output_name + '.concat.txt\n')
                    
                shfile3.write('cp $SCRATCHDIR/* $DATADIR/'+recode_dir+' || export CLEAN_SCRATCH=false\n'+
                              'printf "\\nFinished\\n\\n"\n')
                shfile3.close()

                if print1 is False:
                    cmd3 = ('qsub ' + recode_dir + output_name + '.bpm.sh')
                    p3 = subprocess.Popen(cmd3, shell=True)
                    sts3 = os.waitpid(p3.pid, 0)[1]

                    self.log_file.write("###  Calculate Between-Population-Metrics  ###\n" +
                                        "Input Directory: " + recode_dir + "\n" +
                                        "Window Size: " + str(window_size) + "\n" +
                                        "Minimum SNPs in a window: " + str(min_snps) + "\n" +
                                        "Use repolarized data: " + str(use_repol) + "\n" +
                                        "Populations: " + str(pops) + "\n")

                else:
                    file3 = open(recode_dir + output_name + '.bpm.sh', 'r')
                    data3 = file3.read()
                    print(data3)

                os.remove(recode_dir + output_name + '.bpm.sh')
        elif len(pops) < 2:
            print("'pops' argument must be a list of strings specifiying two or more population names as they appear in input file prefixes.  len(pops) was < 2")
        else:
            print("Did not find recode_dir.  Must run splitVCFs followed by recode before able to calculate between population metrics")


    def calcPairwisebpm(self, recode_dir, pops, window_size, min_snps, print1=False, mem=16, ncpu=1, use_repol=True, keep_intermediates=False, time_scratch="1:00:00", overwrite=False, scratch_gb=1):
        '''Purpose:  Calculate between population metrics including: Dxy, Fst (using Weir and Cockerham 1984), and Rho (Ronfort et al. 1998)
           Notes: User provides a list of populations to be included.  For pairwise estimates, simply provide two populations
                    Calculations are done for windows of a given bp size.  User also must specify the minimum number of snps (min_snps) in a window
                    for calculations to be made'''

        if use_repol is True:
            suffix = '.table.repol.txt'
        else:
            suffix = '.table.recode.txt'

        if recode_dir.endswith("/") is False:
            recode_dir += "/"

        if os.path.exists(recode_dir) is True and len(pops) > 1:
            # Concatenate input files and sort them
            for i, pop1 in enumerate(pops):  # Add data from all populations to single, huge list
                for pop2 in pops[i + 1:]:
                    output_name = pop1 + pop2 + "_WS" + str(window_size) + "_MS" + str(min_snps)
                    skip = False
                    try:
                        a = open(recode_dir + pop1 + suffix, 'r')
                        a.close()
                    except IOError:
                        print("Did not find input file for pop ", pop1)
                        skip = True
                    try:
                        a = open(recode_dir + pop2 + suffix, 'r')
                        a.close()
                    except IOError:
                        print("Did not find input file for pop ", pop2)
                        skip = True
                    if skip is True:
                        print("Did not find all input files!!  Aborting pairwise bpm for contrast: ", output_name)
                    elif os.path.exists(recode_dir + output_name + '_BPM.txt') and overwrite is False:
                        print(recode_dir + output_name + '_BPM.txt already exists.  Set "overwrite" to True if you want to overwrite.')
                    else:
                        shfile3 = open(recode_dir + output_name + '.bpm.sh', 'w')

                        shfile3.write('#!/bin/bash -e\n' +
                                      '#PBS -N '+output_name +'\n' +
                                      '#PBS -l walltime='+str(time_scratch)+'\n' +
                                      '#PBS -l select=1:ncpus='+str(ncpu)+':mem='+str(mem)+'gb:scratch_local='+str(scratch_gb)+'gb\n' +
                                      '#PBS -m abe\n' +
                                      '#PBS -j oe\n\n' +
                                      'module add python-3.4.1-gcc\n'+
                                      'module add python34-modules-gcc\n'+
                                      'trap \'clean_scratch\' TERM EXIT\n'+
                                      'if [ ! -d "$SCRATCHDIR" ] ; then echo "Scratch not created!" 1>&2; exit 1; fi \n' +
                                      'DATADIR="/storage/pruhonice1-ibot/home/holcovam/ScanTools"\n' +
                                      'cp $DATADIR/bpm.py $SCRATCHDIR || exit 1\n'+
                                      'cp $DATADIR/'+ recode_dir + pop1 + suffix +' $SCRATCHDIR || exit 1\n'+
                                      'cp $DATADIR/'+ recode_dir + pop2 + suffix +' $SCRATCHDIR || exit 1\n'+
                                      'cd $SCRATCHDIR || exit 2\n' +
                                      'echo data loaded at `date`\n' +
                                      'sort -k3,3 -k4,4n -m '+ pop1 + suffix + " "+ pop2 + suffix + " > " + output_name + '.concat.txt\n' +
                                      'python3 bpm.py -i ' + output_name + '.concat.txt -o $SCRATCHDIR/ -prefix ' + output_name + ' -ws ' + str(window_size) + ' -ms ' + str(min_snps) + ' -np 2\n'+
                                      'rm bpm.py\n'+
                                      'rm '+ pop1 + suffix +'\n'+
                                      'rm '+ pop2 + suffix +'\n')
                        if keep_intermediates is False:
                            shfile3.write('rm ' + output_name + '.concat.txt\n')
                        shfile3.write('cp $SCRATCHDIR/* $DATADIR/'+recode_dir+' || export CLEAN_SCRATCH=false\n'+
                                      'printf "\\nFinished\\n\\n"\n')
                        shfile3.close()

                        if print1 is False:
                            cmd3 = ('qsub ' + recode_dir + output_name + '.bpm.sh')
                            p3 = subprocess.Popen(cmd3, shell=True)
                            sts3 = os.waitpid(p3.pid, 0)[1]

                        else:
                            file3 = open(recode_dir + output_name + '.bpm.sh', 'r')
                            data3 = file3.read()
                            print(data3)

                        os.remove(recode_dir + output_name + '.bpm.sh')

            if print1 is False:
                self.log_file.write("###  Calculate PAIRWISE Between-Population-Metrics  ###\n" +
                                    "Input Directory: " + recode_dir + "\n" +
                                    "Window Size: " + str(window_size) + "\n" +
                                    "Minimum SNPs in a window: " + str(min_snps) + "\n" +
                                    "Use repolarized data: " + str(use_repol) + "\n" +
                                    "Populations: " + str(pops) + "\n")
        elif len(pops) < 2:
            print("'pops' argument must be a list of strings specifiying two or more population names as they appear in input file prefixes.  len(pops) was < 2")
        else:
            print("Did not find recode_dir.  Must run splitVCFs followed by recode before able to calculate between population metrics")

###Majda
    def calcPairwisebpmAnn(self, recode_dir, pops, window_size, min_snps, print1=False, mem=16, ncpu=1, use_repol=True, keep_intermediates=False, time_scratch="1:00:00", overwrite=False, scratch_gb=1):
        '''Purpose:  Calculate between population metrics including: Dxy, Fst (using Weir and Cockerham 1984), and Rho (Ronfort et al. 1998)
           Notes: User provides a list of populations to be included.  For pairwise estimates, simply provide two populations
                    Calculations are done for windows of a given bp size.  User also must specify the minimum number of snps (min_snps) in a window
                    for calculations to be made'''

        if use_repol is True:
            suffix = '.table.repol.txt'
        else:
            suffix = '.table.recode.txt'

        if recode_dir.endswith("/") is False:
            recode_dir += "/"

        if os.path.exists(recode_dir) is True and len(pops) > 1:
            # Concatenate input files and sort them
            for i, pop1 in enumerate(pops):  # Add data from all populations to single, huge list
                for pop2 in pops[i + 1:]:
                    output_name = pop1 + pop2 + "_WS" + str(window_size) + "_MS" + str(min_snps)
                    skip = False
                    try:
                        a = open(recode_dir + pop1 + suffix, 'r')
                        a.close()
                    except IOError:
                        print("Did not find input file for pop ", pop1)
                        skip = True
                    try:
                        a = open(recode_dir + pop2 + suffix, 'r')
                        a.close()
                    except IOError:
                        print("Did not find input file for pop ", pop2)
                        skip = True
                    if skip is True:
                        print("Did not find all input files!!  Aborting pairwise bpm for contrast: ", output_name)
                    elif os.path.exists(recode_dir + output_name + '_BPM.txt') and overwrite is False:
                        print(recode_dir + output_name + '_BPM.txt already exists.  Set "overwrite" to True if you want to overwrite.')
                    else:
                        shfile3 = open(recode_dir + output_name + '.bpm.sh', 'w')

                        shfile3.write('#!/bin/bash -e\n' +
                                      '#PBS -N '+output_name +'\n' +
                                      '#PBS -l walltime='+str(time_scratch)+'\n' +
                                      '#PBS -l select=1:ncpus='+str(ncpu)+':mem='+str(mem)+'gb:scratch_local='+str(scratch_gb)+'gb\n' +
                                      '#PBS -m abe\n' +
                                      '#PBS -j oe\n\n' +
                                      'module add python-3.4.1-gcc\n'+
                                      'module add python34-modules-gcc\n'+
                                      'trap \'clean_scratch\' TERM EXIT\n'+
                                      'if [ ! -d "$SCRATCHDIR" ] ; then echo "Scratch not created!" 1>&2; exit 1; fi \n' +
                                      'DATADIR="/storage/pruhonice1-ibot/home/holcovam/ScanTools"\n' +
                                      'cp $DATADIR/bpmann.py $SCRATCHDIR || exit 1\n'+
                                      'cp $DATADIR/'+ recode_dir + pop1 + suffix +' $SCRATCHDIR || exit 1\n'+
                                      'cp $DATADIR/'+ recode_dir + pop2 + suffix +' $SCRATCHDIR || exit 1\n'+
                                      'cd $SCRATCHDIR || exit 2\n' +
                                      'echo data loaded at `date`\n' +
                                      'sort -k3,3 -k4,4n -m '+ pop1 + suffix + " "+ pop2 + suffix + " > " + output_name + '.concat.txt\n' +
                                      'python3 bpmann.py -i ' + output_name + '.concat.txt -o $SCRATCHDIR/ -prefix ' + output_name + ' -ws ' + str(window_size) + ' -ms ' + str(min_snps) + ' -np 2\n'+
                                      'rm bpmann.py\n'+
                                      'rm '+ pop1 + suffix +'\n'+
                                      'rm '+ pop2 + suffix +'\n')
                        if keep_intermediates is False:
                            shfile3.write('rm ' + output_name + '.concat.txt\n')
                        shfile3.write('cp $SCRATCHDIR/* $DATADIR/'+recode_dir+' || export CLEAN_SCRATCH=false\n'+
                                      'printf "\\nFinished\\n\\n"\n')
                        shfile3.close()

                        if print1 is False:
                            cmd3 = ('qsub ' + recode_dir + output_name + '.bpm.sh')
                            p3 = subprocess.Popen(cmd3, shell=True)
                            sts3 = os.waitpid(p3.pid, 0)[1]

                        else:
                            file3 = open(recode_dir + output_name + '.bpm.sh', 'r')
                            data3 = file3.read()
                            print(data3)

                        os.remove(recode_dir + output_name + '.bpm.sh')

            if print1 is False:
                self.log_file.write("###  Calculate PAIRWISE Between-Population-Metrics  ###\n" +
                                    "Input Directory: " + recode_dir + "\n" +
                                    "Window Size: " + str(window_size) + "\n" +
                                    "Minimum SNPs in a window: " + str(min_snps) + "\n" +
                                    "Use repolarized data: " + str(use_repol) + "\n" +
                                    "Populations: " + str(pops) + "\n")
        elif len(pops) < 2:
            print("'pops' argument must be a list of strings specifiying two or more population names as they appear in input file prefixes.  len(pops) was < 2")
        else:
            print("Did not find recode_dir.  Must run splitVCFs followed by recode before able to calculate between population metrics")

    #MAJDA
    def outlierAASs(self, pops, recode_dir, outlier_quantile, print1, time_scratch, ncpu, mem, scratch_gb, alcode_path, min_afd):

#contrasts = str("WCATET", "SECWCA", "SECTET", "PANWCA", "PANTET", "PANSEC", "PANDIN", "PANBAL", "DINWCA", "DINTET", "DINSEC", "DINBAL", "CROWCA", "CROTET", "CROSEC", "CROPAN", "CRODIN", "CROBAL", "BALWCA", "BALTET", "BALSEC"), recode_dir = "VCF_synNon_DP8.M0.5", outlier_quantile = .99, print1 = True, alcode_path = "../references/lyrataV2/ALcodesAll.txt"
#fish://holcovam@nympha.metacentrum.cz/auto/pruhonice1-ibot/home/holcovam/ScanTools/VCF_all300_DP8.M0.5/TETDIP_WS1_MS1_BPM.txt
        if recode_dir.endswith("/") is False:
            recode_dir += "/"
        
        if os.path.exists(recode_dir) is True and len(pops) > 0:
            for pop in tuple(pops):
                   
                rfile = open(recode_dir + 'outlierAASs'+str(outlier_quantile)+pop+'.r', 'w')
                rfile.write('con<-"'+pop+'"\n'+
                            'msum<-matrix(nrow = 1, ncol = 15,dimnames =list(c(con),c("Smin","S1Q","Smedian","Smean","S3Q","Smax","Ssd", "Nmin","N1Q","Nmedian","Nmean","N3Q","Nmax","Nsd","Squantile_'+str(outlier_quantile)+'")))\n'+
                            'al <- read.table("'+alcode_path.split("/")[-1]+'",header = F)\n'+
                            'mat <- matrix(nrow = nrow(al), ncol = 1,dimnames = list(as.character((al[,1])),c(con)))\n'+
                            'inp<-head(read.table(paste(con,"synNon_WS1_MS1_BPM.txt",sep=""),header=T,fill = T),-1)\n'+
                            's<-subset(inp,inp$ann %in% "synonymous_variant" & inp$AFD >= '+str(min_afd)+')\n'+
                            'n<-subset(inp,inp$ann %in% "missense_variant")\n'+
                            'sq <- quantile(s$FstH,probs=c('+str(outlier_quantile)+'))\n'+
                            'msum[1,1:6]<-summary (s$FstH)\n'+
                            'msum[1,7]<-sd(s$FstH)\n'+
                            'msum[1,8:13]<-summary (n$FstH)\n'+
                            'msum[1,14]<-sd(n$FstH)\n'+
                            'msum[1,15]<-sq\n'+
                            'for (g in readLines("'+alcode_path.split("/")[-1]+'"))\n'+
                            '{ index<-which(al$V1 %in% paste(g))\n'+
                            '  gene<-subset(x = inp,subset = inp$ALcode %in% g & inp$ann %in% "missense_variant")\n'+
                            '  \n'+
                            '  if (nrow(gene)==0) {high<-0} else \n'+
                            '  {high<-nrow(subset(gene,gene$FstH >= as.numeric(sq)))}\n'+
                            '  mat[index,1]<-high\n'+
                            '}\n'+
                            'write.table(mat,append = F,file="FstHighperGene.Allgenes.' + str(outlier_quantile)+ '.' +pop+'.txt",quote = F, sep = "\t", col.names = T,row.names = T)\n'+
                            'write.table(msum,append = F,file ="FstStatsPerContrast.' + str(outlier_quantile)+ '.' +pop+'.txt",quote = F, sep = "\t", col.names = T,row.names = T)\n')
                rfile.close()
                
                shfile3 = open(recode_dir + 'outlierAASs' +str(outlier_quantile)+pop+'.sh', 'w')
                shfile3.write('#!/bin/bash -e\n' +
                              '#PBS -N outlierAASs' +str(outlier_quantile)+'.'+pop+'\n' +
                              '#PBS -l walltime='+str(time_scratch)+'\n' +
                              '#PBS -l select=1:ncpus='+str(ncpu)+':mem='+str(mem)+'gb:scratch_local='+str(scratch_gb)+'gb\n' +
                              '#PBS -m abe\n' +
                              '#PBS -j oe\n\n' +
                              #'MYR=outlierAASs'+str(outlier_quantile)+pop+'.r\n' +
                              'module add R-3.4.3-gcc\n'+
                              'trap \'clean_scratch\' TERM EXIT\n'+
                              'if [ ! -d "$SCRATCHDIR" ] ; then echo "Scratch not created!" 1>&2; exit 1; fi \n' +
                              'DATADIR="/storage/pruhonice1-ibot/home/holcovam/ScanTools"\n' +
                              'cp $DATADIR/'+ recode_dir + 'outlierAASs'+str(outlier_quantile)+pop+'.r $SCRATCHDIR || exit 1\n'+
                              'cp $DATADIR/'+ alcode_path +' $SCRATCHDIR || exit 1\n'+
                              'cp $DATADIR/'+ recode_dir + pop +'_WS1_MS1_BPM.txt $SCRATCHDIR || exit 1\n'+
                              'cd $SCRATCHDIR || exit 2\n' +
                              'echo data loaded at `date`\n' + 
                              'head -n 1 '+pop+'_WS1_MS1_BPM.txt > '+pop+'synNon_WS1_MS1_BPM.txt \n'+
                              'grep -E "missense_variant|synonymous_variant" '+pop+'_WS1_MS1_BPM.txt >> '+pop+'synNon_WS1_MS1_BPM.txt\n'+                    
                              'Rscript --vanilla outlierAASs'+str(outlier_quantile)+pop+'.r\n'+
                              'rm outlierAASs'+str(outlier_quantile)+pop+'.r\n'+
                              'rm '+ alcode_path.split("/")[-1] +'\n'+
                              'rm *_WS1_MS1_BPM.txt\n'+
                              'cp $SCRATCHDIR/* $DATADIR/'+recode_dir+' || export CLEAN_SCRATCH=false\n'+
                              'printf "\\nFinished\\n\\n"\n')
                shfile3.close()
                if print1 is False:
                    cmd3 = ('qsub ' + recode_dir + 'outlierAASs' +str(outlier_quantile)+pop+'.sh')
                    p3 = subprocess.Popen(cmd3, shell=True)
                    sts3 = os.waitpid(p3.pid, 0)[1]
                else:
                    file3 = open(recode_dir + 'outlierAASs'+str(outlier_quantile)+pop+'.r', 'r')
                    data3 = file3.read()
                    print(data3)
                    file4 = open(recode_dir + 'outlierAASs' +str(outlier_quantile)+pop+'.sh', 'r')
                    data4 = file4.read()
                    print(data4)
                    
                #os.remove(recode_dir + 'outlierAASs'+str(outlier_quantile)+pop+'.r')
                os.remove(recode_dir + 'outlierAASs' +str(outlier_quantile)+pop+'.sh')

    def N_SperGene(self, lineages, recode_dir, print1, time_scratch, ncpu, mem, scratch_gb, alcode_path):

        if recode_dir.endswith("/") is False:
            recode_dir += "/"
        togrep=str('missense_variant|synonymous_variant')
        rfile = open(recode_dir + 'N_SperGene.r', 'w')
        rfile.write('lineages<-c('+lineages+')\n'+
                    'for (lin in  lineages){\n'+
                    #'  system(command = paste("grep -E '+togrep+'",lin,".table.repol.txt > ",lin,".synNon.txt",sep=""))\n'+
                    '  inp<-read.table(paste(lin,".synNon.txt",sep=""),header=F)\n'+
                    '  inp[inp=="-9"]<-NA\n'+
                    '  for (g in readLines("'+alcode_path.split("/")[-1]+'"))\n'+
                    '  { mat <- matrix(nrow = 1, ncol = 6)\n'+
                    '    gene<-subset(x = inp,subset = inp$V7 %in% g)\n'+
                    '    if (nrow(subset(gene,gene$V8 %in% "missense_variant"))==0) {n<-0} else \n'+
                    '    {n<-sum(subset(gene,gene$V8 %in% "missense_variant")[,10:ncol(gene)],na.rm = T)}\n'+
                    '    if (nrow(subset(gene,gene$V8 %in% "synonymous_variant"))==0) {s<-0} else \n'+
                    '    {s<-sum(subset(gene,gene$V8 %in% "synonymous_variant")[,10:ncol(gene)],na.rm = T)}\n'+
                    '    if ((n+s)<5) {check<-NA} else \n'+
                    '    {check<-1}\n'+
                    '    ns<-as.numeric(n)/as.numeric(s)\n'+
                    '    mat[1,1]<-as.character(gene[1,1])\n'+
                    '    mat[1,2]<-g\n'+
                    '    mat[1,3]<-n\n'+
                    '    mat[1,4]<-s\n'+
                    '    mat[1,5]<-ns\n'+
                    '    mat[1,6]<-check\n'+
                    '    write.table(mat,append = T,file = "NSperGene.AllLin.Allgenes.txt",quote = F, sep = "\t",col.names = F,row.names = F)}\n'+
                    '  mat1 <- matrix(nrow = 1, ncol = 6)\n'+
                    '  n<-sum(subset(inp,inp$V8 %in% "missense_variant")[,10:ncol(inp)],na.rm = T)\n'+
                    '  s<-sum(subset(inp,inp$V8 %in% "synonymous_variant")[,10:ncol(inp)],na.rm = T)\n'+
                    '  ns<-as.numeric(n)/as.numeric(s)\n'+
                    '  mat1[1,1]<-as.character(inp[1,1])\n'+
                    '  mat1[1,2]<-"genome"\n'+
                    '  mat1[1,3]<-n\n'+
                    '  mat1[1,4]<-s\n'+
                    '  mat1[1,5]<-ns\n'+
                    '  mat1[1,6]<-"1"\n'+
                    '  write.table(mat1,append = T,file = "NSgenomeWide.AllLin.txt",quote = F, sep = "\t",col.names = F,row.names = F)}\n')
        rfile.close()
        
        shfile3 = open(recode_dir + 'N_SperGene.sh', 'w')
        shfile3.write('#!/bin/bash -e\n' +
                      '#PBS -N N_SperGene\n' +
                      '#PBS -l walltime='+str(time_scratch)+'\n' +
                      '#PBS -l select=1:ncpus='+str(ncpu)+':mem='+str(mem)+'gb:scratch_local='+str(scratch_gb)+'gb\n' +
                      '#PBS -m abe\n' +
                      '#PBS -j oe\n\n' +
                      'module add R-3.4.3-gcc\n'+
                      'trap \'clean_scratch\' TERM EXIT\n'+
                      'if [ ! -d "$SCRATCHDIR" ] ; then echo "Scratch not created!" 1>&2; exit 1; fi \n' +
                      'DATADIR="/storage/pruhonice1-ibot/home/holcovam/ScanTools"\n' +
                      'cp $DATADIR/'+ recode_dir + 'N_SperGene.r $SCRATCHDIR || exit 1\n'+
                      'cp $DATADIR/'+ alcode_path +' $SCRATCHDIR || exit 1\n'+
                      'linbash=('+ lineages.replace(","," ") +')\n'+
                      'for i in "${linbash[@]}"; do\n'+
                      'cp $DATADIR/'+ recode_dir +'$i.table.repol.txt $SCRATCHDIR || exit 1\n'+
                      'done \n'+
                      'cd $SCRATCHDIR || exit 2\n' +
                      'echo data loaded at `date`\n' +
                      #'linbash=('+ lineages.replace(","," ") +')\n'+
                      'for i in "${linbash[@]}"; do\n'+
                      '    grep -E "missense_variant|synonymous_variant" $i.table.repol.txt > $i.synNon.txt\n'+
                      'done \n'+
                      'Rscript --vanilla N_SperGene.r\n'+
                      'rm N_SperGene.r\n'+
                      'rm '+ alcode_path.split("/")[-1] +'\n'+
                      'rm *.table.repol.txt\n'+
                      'cp $SCRATCHDIR/* $DATADIR/'+recode_dir+' || export CLEAN_SCRATCH=false\n'+
                      'printf "\\nFinished\\n\\n"\n')
        shfile3.close()

        if print1 is False:
            cmd3 = ('qsub ' + recode_dir + 'N_SperGene.sh')
            p3 = subprocess.Popen(cmd3, shell=True)
            sts3 = os.waitpid(p3.pid, 0)[1]

        else:
            file3 = open(recode_dir + 'N_SperGene.r', 'r')
            data3 = file3.read()
            print(data3)
            file4 = open(recode_dir + 'N_SperGene.sh', 'r')
            data4 = file4.read()
            print(data4)
            
        #os.remove(recode_dir + 'N_SperGene.r')
        os.remove(recode_dir + 'N_SperGene.sh')



##NOT FINISHED!!!!!
    def calcFnD(self, recode_dir, p1, p2, p3, pO, window_size, min_snps, print1=False, mem=16000, numcores=1, partition="nbi-medium", use_repol=True, keep_intermediates=False, time_scratch="0-12:00", use_scratch=False, scratch_path="/nbi/scratch/monnahap"):
        '''Purpose:  \
           Notes: '''

        if use_repol is True:
            suffix = '.table.repol.txt'
        else:
            suffix = '.table.recode.txt'

        if recode_dir.endswith("/") is False:
            recode_dir += "/"

        if use_scratch is True:
            if scratch_path.endswith("/") is False:
                scratch_path += "/"
            tmpdir = scratch_path + recode_dir.split("/")[-2] + "/"
            if os.path.exists(tmpdir) is False:
                os.mkdir(tmpdir)
        else:
            tmpdir = recode_dir

        output_name = ".".join([p1,p2,p3,pO]) + "_WS" + str(window_size) + "_MS" + str(min_snps)
        if os.path.exists(recode_dir) is True:
            pop_num = 0
            file_string = ""
            for pop in [p1, p2, p3, pO]:
                try:
                    a = open(recode_dir + pop + suffix, 'r')
                    a.close()
                    file_string += recode_dir + pop + suffix + " "
                    pop_num += 1
                except IOError:
                    print("Did not find input file for pop ", pop)
            if pop_num != 4:
                print("Did not find all input files!!  Aborting.")
                os.remove(recode_dir + output_name + '.concat.txt')
            else:
                shfile3 = open(recode_dir + output_name + '.FnD.sh', 'w')

                shfile3.write('#!/bin/bash\n' +
                              '#SBATCH -J ' + output_name + '.FnD.sh' + '\n' +
                              '#SBATCH -e ' + self.oande + output_name + '.FnD.err' + '\n' +
                              '#SBATCH -o ' + self.oande + output_name + '.FnD.out' + '\n' +
                              '#SBATCH -p ' + str(partition) + '\n' +
                              '#SBATCH -n ' + str(numcores) + '\n' +
                              '#SBATCH -t ' + str(time_scratch) + '\n' +
                              '#SBATCH --mem=' + str(mem) + '\n' +
                              'source python-3.5.1\n' +
                              'sort -k3,3 -k4,4n -m ' + file_string + '> ' + tmpdir + output_name + '.concat.txt\n' +
                              'python3 ' + self.code_dir + '/calcFnD.py -i ' + tmpdir + output_name + '.concat.txt' + ' -o ' + recode_dir + ' -prefix ' + output_name + ' -ws ' + str(window_size) + ' -ms ' + str(min_snps) + '\n')
                if keep_intermediates is False:
                    shfile3.write('rm ' + tmpdir + output_name + '.concat.txt')
                shfile3.close()

                if print1 is False:
                    cmd3 = ('sbatch -d singleton ' + recode_dir + output_name + '.FnD.sh')
                    p3 = subprocess.Popen(cmd3, shell=True)
                    sts3 = os.waitpid(p3.pid, 0)[1]

                    self.log_file.write("###  Calculate f_D and D  ###\n" +
                                        "Input Directory: " + recode_dir + "\n" +
                                        "Window Size: " + str(window_size) + "\n" +
                                        "Minimum SNPs in a window: " + str(min_snps) + "\n" +
                                        "Use repolarized data: " + str(use_repol) + "\n" +
                                        "Populations: " + ",".join([p1, p2, p3, pO]) + "\n")

                else:
                    file3 = open(recode_dir + output_name + '.FnD.sh', 'r')
                    data3 = file3.read()
                    print(data3)

                os.remove(recode_dir + output_name + '.FnD.sh')
        else:
            print("Did not find recode_dir.  Must run splitVCFs followed by recode before able to calculate f_D and D")


    def calcAFS(self, recode_dir, data_name, window_size="50000", num_reps="200", sampind="-99", pops="-99", time_scratch='00:30:00', mem="4", use_repol=True, print1=False, allow_one_missing=True,scratch_gb="10", ncpu="2"):
        """Purpose:  calculate allele frequency spectrum for each population in pops list and downsample to number specified in sampind.  If no downsampling is specified the number of indviduals in the population (minus one if allow_one_missing=True) will be used."""
        if recode_dir.endswith("/") is False:
            recode_dir += "/"
        if pops == "-99":
            pops = self.pops

        if use_repol is True:
            suffix = '.table.repol.txt'
        else:
            suffix = '.table.recode.txt'
        for pop in pops:
            if sampind == "-99":
                si = self.samp_nums[pop]
                if allow_one_missing is True:
                    si -= 1
            else:
                si = sampind
            infile = recode_dir + pop + suffix
            shfile3 = open(infile + '.afs.sh', 'w')
            prefix = pop + "_" + data_name
            shfile3.write('#!/bin/bash -e\n' +
                          '#PBS -N '+pop + "_" + data_name+'\n' +
                          '#PBS -l walltime='+time_scratch+'\n' +
                          '#PBS -l select=1:ncpus='+ncpu+':mem='+mem+'gb:scratch_local='+scratch_gb+'gb\n' +
                          '#PBS -m abe\n' +
                          '#PBS -j oe\n\n' +
                          'module add python-3.4.1-gcc\n'+
                          'module add python34-modules-gcc\n'+
                          'trap \'clean_scratch\' TERM EXIT\n'+
                          'if [ ! -d "$SCRATCHDIR" ] ; then echo "Scratch not created!" 1>&2; exit 1; fi \n' +
                          'DATADIR="/storage/pruhonice1-ibot/home/holcovam/ScanTools"\n' +
                          'cp $DATADIR/calcAFS.py $SCRATCHDIR || exit 1\n' +
                          'cp $DATADIR/'+ infile +' $SCRATCHDIR || exit 1\n'+
                          'cd $SCRATCHDIR || exit 2\n' +
                          'echo data loaded at `date`\n\n' +
                          'python3 calcAFS.py -i ' +pop + suffix+ ' -o $SCRATCHDIR/ -p ' + prefix + ' -nrep ' + num_reps + ' -ws ' + window_size + ' -sampind ' + str(si) + '\n\n' +
                          'rm '+pop + suffix+'\n'+
                          'rm calcAFS.py\n'+
                          'cp $SCRATCHDIR/* $DATADIR/'+recode_dir+' || export CLEAN_SCRATCH=false\n'+
                          'printf "\\nFinished\\n\\n"\n')

            shfile3.close()
            if print1 is False:
                cmd3 = ('qsub ' + infile + '.afs.sh')
                p3 = subprocess.Popen(cmd3, shell=True)
                sts3 = os.waitpid(p3.pid, 0)[1]

            else:
                file3 = open(infile + '.afs.sh', 'r')
                data3 = file3.read()
                print(data3)
            os.remove(infile + '.afs.sh')


    def calcFreqs(self, recode_dir, outfile_name, sites_file, pops="-99", time_scratch="0-02:00", partition="nbi-short", mem="8000", use_repol=True, print1=False):
        """Calls calcFreqs_atSites.py.  Takes a list of sites in (sites_file, should be formatted so that each line simply has scaffold and position, with scaffold simply coded as an integer 0-8) and calculates the allele frequency in each population (list_of_populations) at each site"""
        if recode_dir.endswith("/") is False:
            recode_dir += "/"
        if pops == "-99":
            pops = self.pops

        if use_repol is True:
            suffix = '.table.repol.txt'
        else:
            suffix = '.table.recode.txt'

        pop_string = ""
        for pop in pops:
            if os.path.exists(recode_dir + pop + suffix) is True:
                pop_string += pop + ","
        pop_string.strip(",")

        shfile3 = open('CalcFreqs.sh', 'w')
        shfile3.write('#!/bin/bash\n' +
                      '#SBATCH -J CalcFreqs.sh' + '\n' +
                      '#SBATCH -e ' + self.oande + 'CalcFreqs.err' + '\n' +
                      '#SBATCH -o ' + self.oande + 'CalcFreqstest.out' + '\n' +
                      '#SBATCH -p ' + partition + '\n'
                      '#SBATCH -n 1' + '\n' +
                      '#SBATCH -t ' + str(time_scratch) + '\n' +
                      '#SBATCH --mem=' + str(mem) + '\n' +
                      'source python-3.5.1\n' +
                      'python3 ' + self.code_dir + '/calcFreqs_atSites.py -i ' + recode_dir + ' -o ' + recode_dir + ' -of ' + outfile_name + ' -s ' + sites_file + ' -pops ' + pop_string + ' -suffix ' + suffix + '\n')
        shfile3.close()
        if print1 is False:
            cmd3 = ('sbatch CalcFreqs.sh')
            p3 = subprocess.Popen(cmd3, shell=True)
            sts3 = os.waitpid(p3.pid, 0)[1]

        else:
            file3 = open('CalcFreqs.sh', 'r')
            data3 = file3.read()
            print(data3)
        os.remove('CalcFreqs.sh')


    def concatWPM(self, recode_dir, suffix, outname, pops='all'):
        '''Purpose:  Concatenate _WPM.txt files corresponding to populations indicated in pops parameter.'''

        if recode_dir.endswith("/") is False:
            recode_dir += "/"

        if os.path.exists(recode_dir) is True:
            new = open(recode_dir + outname + "_WPM.txt", 'w')
            if pops == 'all':
                pops = self.pops
            head = False
            for i, pop in enumerate(pops):
                try:
                    with open(recode_dir + pop + suffix, 'r') as inf:
                        for j, line in enumerate(inf):
                            if j == 0 and head is False:
                                new.write(line)
                                head = True
                            elif j != 0:
                                new.write(line)
                except FileNotFoundError:
                    print("Did not find _WPM.txt file for population: ", pop)


    def concatBPM(self, recode_dir, suffix, outname, pops='all', strict=False):
        '''Purpose:  Concatenate _WPM.txt files corresponding to populations indicated in pops parameter.'''

        if recode_dir.endswith("/") is False:
            recode_dir += "/"

        if os.path.exists(recode_dir) is True:
            new = open(recode_dir + outname + "_BPM.txt", 'w')
            if pops == 'all':
                pops = self.pops
            head = False
            for file in os.listdir(recode_dir):
                if file.endswith(suffix):
                    if strict is False:
                        if file.split("_")[0][:3] in pops or file.split("_")[0][3:] in pops:
                            with open(recode_dir + file, 'r') as inf:
                                for j, line in enumerate(inf):
                                    if j == 0 and head is False:
                                        new.write(line)
                                        head = True
                                    elif j != 0:
                                        new.write(line)
                    elif strict is True:
                        if file.split("_")[0][:3] in pops and file.split("_")[0][3:] in pops:
                            with open(recode_dir + file, 'r') as inf:
                                for j, line in enumerate(inf):
                                    if j == 0 and head is False:
                                        new.write(line)
                                        head = True
                                    elif j != 0:
                                        new.write(line)


    def findOutliers(self, recode_dir, in_file, column_index_list, percentile, tails='upper'):
        '''Purpose:  Take output from either calcwpm or calcbpm and determine outlier metrics for a given percentile.
           Notes: Output will be two csv files (one containing all sites with outliers indicated by 0 or 1 and another containing just outliers)
                  as well as a bed file to be used in annotateOutliers'''

        if recode_dir.endswith("/") is False:
            recode_dir += "/"

        if percentile > 1.0:
            print("!!percentile parameter should be coded as a proportion, not as a percentage!!")
            return False

        if os.path.exists(recode_dir) is True:
            try:
                print("Finding outliers for " + recode_dir + in_file + "\n")
                data = pandas.read_table(recode_dir + in_file, header=0)
                metrics = []
                for i in column_index_list:
                    try:
                        metrics.append(list(data.columns.values)[i])
                    except IndexError:
                        print("IndexError in " + recode_dir + in_file + "\n")
                data.start = data.start.astype(int)
                data.end = data.end.astype(int)
                for metric in metrics:
                    data[metric + '.out'] = 0
                    if tails == 'both':
                        data[metric + '.out'].loc[(data[metric] > data[metric].quantile(q=percentile))] = 1
                        data[metric + '.out'].loc[(data[metric] < data[metric].quantile(q=1.0 - percentile))] = 1
                    elif tails == 'lower':
                        data[metric + '.out'].loc[(data[metric] < data[metric].quantile(q=1.0 - percentile))] = 1
                    elif tails == 'upper':
                        data[metric + '.out'].loc[(data[metric] > data[metric].quantile(q=percentile))] = 1
                    else:
                        print("Did not specify tails option correctly.  Options are: both, upper, and lower")
                data['num_outliers'] = data.iloc[:, -len(metrics):].sum(1)
                data.to_csv(recode_dir + in_file.replace(".txt", "") + '_' + str(percentile) + 'tile_OutLabelled.csv', index=False)
                # select all windows that are outliers for at least one metric
                df_outlier = data[(data.num_outliers != 0)]
                df_outlier.to_csv(recode_dir + in_file.replace(".txt", "") + '_' + str(percentile) + 'tile_OutOnly.csv', index=False)
                df_outlier.to_csv(recode_dir + in_file.replace(".txt", "") + '_' + str(percentile) + 'tile_OutOnly.bed', index=False, sep='\t', columns=["scaff", "start", "end"], header=False)
            except (UnicodeDecodeError, IndexError, KeyError):
                print('Error with file: ' + recode_dir + in_file + "\n")

    def annotateOutliers(self, recode_dir, in_file, annotation_file="../references/lyrataV2/LyV2.gff", overlap_proportion=0.000001, print1=False, time_scratch="0:10:00", mem_mb=100, ncpu=1,scratch_mb=100):
        '''Purpose: annotate bed file from findOutliers using information in annotation_file
           Notes: The output (suffix ol_genes.gff) only contains the window locations along with annotation info and does not contain
                    the original metric information used to determine outliers.  Use mergeAnnotation to merge original outlier file with annotation info'''

        if recode_dir.endswith("/") is False:
            recode_dir += "/"
        basename = in_file.strip(".bed")
        if os.path.exists(recode_dir) is True:
            shfile1 = open(recode_dir + in_file + 'bedtools_gff.sh', 'w')
            shfile1.write('#!/bin/bash -e\n' +
                          '#PBS -N annotate\n' +
                          '#PBS -l walltime='+str(time_scratch)+'\n' +
                          '#PBS -l select=1:ncpus='+str(ncpu)+':mem='+str(mem_mb)+'mb:scratch_local='+str(scratch_mb)+'mb\n' +
                          '#PBS -m abe\n' +
                          '#PBS -j oe\n\n' +
                          'module add bedtools-2.26.0\n' +
                          'trap \'clean_scratch\' TERM EXIT\n'+
                          'if [ ! -d "$SCRATCHDIR" ] ; then echo "Scratch not created!" 1>&2; exit 1; fi \n' +
                          'DATADIR="/storage/pruhonice1-ibot/home/holcovam/ScanTools"\n' +
                          'cp $DATADIR/'+ annotation_file +' $SCRATCHDIR || exit 1\n'+
                          'cp $DATADIR/'+ recode_dir + in_file +' $SCRATCHDIR || exit 1\n'+
                          'cd $SCRATCHDIR || exit 2\n' +
                          'echo data loaded at `date`\n' +
                          'module add bedtools-2.26.0\n' +
                          'bedtools intersect -a ' + in_file + ' -b ' + annotation_file.split("/")[-1] + ' -f ' + str(overlap_proportion) + ' -wo | grep transcript | grep -v transcription | sort -u |' +
                          """awk '{print $1,$2,$3,$7,$8,$9,$10,$12}'""" +
                          '| ' +
                          """tr ' ' '\t' """
                          '> ' + basename + '_genes.gff\n'+
                          'rm '+annotation_file.split("/")[-1]+'\n'+
                          'rm '+ in_file+'\n'+
                          'cp $SCRATCHDIR/* $DATADIR/'+recode_dir+' || export CLEAN_SCRATCH=false\n'+
                          'printf "\\nFinished\\n\\n"\n')
            shfile1.close()

            if print1 is False:
                cmd3 = ('qsub ' + recode_dir + in_file + 'bedtools_gff.sh')
                p3 = subprocess.Popen(cmd3, shell=True)
                sts3 = os.waitpid(p3.pid, 0)[1]

            else:
                file3 = open(recode_dir + in_file + 'bedtools_gff.sh', 'r')
                data3 = file3.read()
                print(data3)

            os.remove(recode_dir + in_file + 'bedtools_gff.sh')

        else:
            print("recode_dir not found")


    def mergeAnnotation(self, recode_dir, outlier_file):
        '''Purpose: Merge the annotation information with the original outlier file results from findOutliers.'''

        if recode_dir.endswith("/") is False:
            recode_dir += "/"
        annotated_outlier_file = outlier_file.strip('.csv') + '_genes.gff'
        if os.path.exists(recode_dir) is True:
            try:
                outliers = pandas.read_csv(recode_dir + outlier_file, header=0)
                annotation = pandas.read_table(recode_dir + annotated_outlier_file, names=["scaff", "start", "end", "gene_start", "gene_end", "overlap", "strand", "geneName"])
            except IOError:
                print("Did not find either original outlier file or the annotated outlier file")
            merged = outliers.merge(annotation, on=["scaff", "start", "end"],)
            merged.to_csv(recode_dir + outlier_file.replace("_OutOnly.csv", "_OutAnnot.csv"), index=False)
        else:
            print("Did not find recode_dir")
            
    #MAJDA     
    def getfunctionOutliers(self, recode_dir, in_file, annotation_file="../references/lyrataV2/LyV2_TAIR10orth_des_20150927.txt"):
        if recode_dir.endswith("/") is False:
            recode_dir += "/"
        basename = in_file.strip("_OutAnnot.csv")
        print(basename)
        outfile = open(recode_dir + basename + '_GF.txt', 'w')
        #with open(infile) as inf:
            #for j, tline in enumerate(inf):
            #for tline in infile:
        infile = open(recode_dir + in_file, 'r')
        print(infile)
        for infileline in infile:
            tline = infileline.split('=')
            uline = ''.join(tline[-1:])
            line = uline.replace('\n', '')
            print(line)
            genefunction = open(annotation_file, 'r')
            for genefunctionline in genefunction:
                data = genefunctionline.split('\t')
                if line in data[0]:
                    outfile.write(data[0]+'\t'+data[1]+'\t'+data[3]+'\t'+data[4]+'\t'+data[5]+'\t'+data[6]+'\n')
        infile.close()
        outfile.close()

    #MAJDA
    def graphs(self, recode_dir, in_filew, in_filed, in_filegf, pops):
        if recode_dir.endswith("/") is False:
            recode_dir += "/"
        rfile = open(recode_dir + pops+'.r', 'w')
        rfile.write('# load libraries\n'+
                    'library(ggplot2,lib.loc="../programs/Rpackages/")\n'+
                    'library(grid)\n'+
                    'descr <- read.table("'+recode_dir+in_filed+'", sep="\t", header=T)\n'+
                    'window <- read.table("'+recode_dir+in_filew+'", header=T)\n'+
                    'f <- read.table("../references/lyrataV2/Lyrata2genesOriented.out",header=F)\n'+
                    'al <- read.table ("'+recode_dir+in_filegf+'", sep="=", header = F)\n'+
                    'for(i in 1:nrow(al)) {\n'+
                    'locus=as.character(al[i,3])\n'+
                    'scaf=substr(locus, 3, 3)\n'+
                    'winsize = 50000\n'+
                    'locstart <- subset(f[,3] , as.vector(f[,1])==locus)\n'+
                    'locend <- subset(f[,4], as.vector(f[,1])==locus)\n'+
                    'middle=(locend+locstart)/2\n'+
                    'start=middle-winsize\n'+
                    'end=middle+winsize\n'+
                    'genes=subset(f, (f[,2]==as.numeric(scaf) & f[,3]>=start & f[,3]<=end) | (f[,2]==1 & f[,4]>=start & f[,4]<=end))\n'+
                    'genes$geneor <- ifelse(genes[,5]=="+","last","first")\n'+
                    'descr_locus <- descr[descr$start >= start & descr$end <=end,]\n'+
                    'window_locus <- window[window$start >= start & window$end <=end,]\n'+
                    'layout <- theme_bw(base_size = 60, base_family="Helvetica") +\n'+
                    '  theme(axis.title.x = element_blank(),\n'+
                    '        axis.text.x = element_blank(),\n'+
                    '        axis.ticks.x = element_blank(),\n'+
                    '        panel.grid = element_blank())\n'+
                    'layouttix <- theme_bw(base_size=60, base_family="Helvetica") +\n'+
                    '  theme(axis.title.x = element_blank(),\n'+
                    '        panel.grid = element_blank())\n'+
                    'pafds <- ggplot(aes(start, AFD), data=descr_locus) + ylab("AFDsnp") +\n'+
                    '  ggtitle(locus) +\n'+
                    '  scale_y_continuous(limits=c(-1,1.1)) +\n'+
                    '  geom_point(alpha=0.35, size=12) +\n'+
                    '  layouttix\n'+
                    'for(j in 1:nrow(genes)){\n'+
                    '  pafds <- pafds + geom_segment(x=genes[j,3], xend=genes[j,4], y=1.1, yend=1.1, colour="grey30",\n'+
                    '                                arrow=arrow(length=unit(0.03, "npc"), ends=genes[j,"geneor"]),size=4)\n'+
                    ' }\n'+
                    'p.AFDs <- pafds  + geom_segment(x=genes[which(genes[,1]==locus), 3],\n'+
                    '                                xend=genes[which(genes[,1]==locus), 4], y=1.1, yend=1.1, color="red", arrow=arrow(length=unit(0.03, "npc"),\n'+
                    '                                                                                                                  ends=genes[which(genes[,1]==locus),"geneor"]),size=4)\n'+
                    'pdxy <- qplot(start, dxy, data=window_locus, geom="line", ylim=c(0,1)) +\n'+
                    '  geom_rect(aes(xmin=genes[which(genes[,1]==locus), 3], xmax=genes[which(genes[,1]==locus), 4], ymin=-Inf, ymax=Inf), fill="grey90") +\n'+
                    '  geom_line(color=I("green"),size=4) + \n'+
                    '  layout\n'+
                    'p.Dxy <- pdxy  + geom_segment(x=genes[which(genes[,1]==locus), 3], xend=genes[which(genes[,1]==locus), 4],\n'+
                    '                              y=1, yend=1, color="red", arrow=arrow(length=unit(0.03, "npc"), ends=genes[which(genes[,1]==locus), "geneor"]),size=4)\n'+
                    'pfst <- qplot(start, Fst, data=window_locus, geom="line", ylim=c(0, 1)) +\n'+
                    '  geom_rect(aes(xmin=genes[which(genes[,1]==locus), 3], xmax=genes[which(genes[,1]==locus), 4],\n'+
                    '                ymin=-Inf, ymax=Inf), fill="grey90") +\n'+
                    '  geom_line(color=I("blue"),size=4) +\n'+
                    '  layout\n'+
                    'p.Fst <- pfst  + geom_segment(x=genes[which(genes[,1]==locus), 3], xend=genes[which(genes[,1]==locus), 4],\n'+
                    '                              y=1, yend=1, color="red", arrow=arrow(length=unit(0.03, "npc"), ends=genes[which(genes[,1]==locus),"geneor"]),size=4)\n'+
                    'prho <- qplot(start, Rho, data=window_locus, geom="line", ylim=c(-0.2, 0.6)) +\n'+
                    '  geom_rect(aes(xmin=genes[which(genes[,1]==locus), 3], xmax=genes[which(genes[,1]==locus), 4],\n'+
                    '                ymin=-Inf, ymax=Inf), fill="grey90",size=4) +\n'+
                    '  geom_line(color=I("goldenrod"),size=4) +\n'+
                    '  layout\n'+
                    'p.Rho <- prho  + geom_segment(x=genes[which(genes[,1]==locus), 3], xend=genes[which(genes[,1]==locus), 4],\n'+
                    '                              y=0.6, yend=0.6, color="red", arrow=arrow(length=unit(0.03, "npc"), \n'+
                    'ends=genes[which(genes[,1]==locus),"geneor"]),size=4)\n'+
                    'path = c(paste ("'+recode_dir+pops+'_'+'",locus,".png", sep=""))\n'+
                    'png(filename = path, width=1754, height=3508, units = "px", pointsize = 60)\n'+
                    'grid.newpage()\n'+
                    'pushViewport(viewport(layout = grid.layout(5,1)))\n'+
                    'vplayout <- function(x,y)\n'+
                    '  viewport(layout.pos.row=x,layout.pos.col=y)\n'+
                    'print(p.AFDs, vp=vplayout(1:2,1))\n'+
                    'print(p.Dxy, vp=vplayout(3,1))\n'+
                    'print(p.Fst, vp=vplayout(4,1))\n'+
                    'print(p.Rho, vp=vplayout(5,1))\n'+
                    'dev.off()}' )
        rfile.close()
        cmd1 = ('Rscript '+recode_dir+pops+'.r')
        subprocess.call(cmd1, shell=True) 
        os.chmod(recode_dir+pops+'.r', 0o755)




    def Outliers(self, recode_dir, in_file, column_index_list, percentile, tails='upper', annotation_file="../references/lyrataV2/LyV2.gff", overlap_proportion=0.000001):
        """Purpose:  Wraps .findOutliers, .annotateOutliers, and .mergeAnnotation into a single function"""
        print("Be sure that no old versions of gff files are in this directory")
        self.findOutliers(recode_dir, in_file, column_index_list, percentile, tails)
        bed_file = in_file.replace(".txt", "") + '_' + str(percentile) + 'tile_OutOnly.bed'
        self.annotateOutliers(recode_dir, bed_file, annotation_file, overlap_proportion)
        outlier_file = in_file.replace(".txt", "") + '_' + str(percentile) + 'tile_OutOnly.csv'
        annotated_outlier_file = outlier_file.strip('.csv') + '_genes.gff'
        while not os.path.exists(recode_dir + annotated_outlier_file):
            time.sleep(1)

        if os.path.isfile(recode_dir + annotated_outlier_file):
            time.sleep(30)
            self.mergeAnnotation(recode_dir, outlier_file)
        else:
            raise ValueError("File error for mergeAnnotation: %s" % recode_dir + annotated_outlier_file)



    def generateFSC2input(self, recode_dir, pops, output_name, bootstrap_block_size=50000, bootstrap_reps=0, mem=16, time_scratch='4:00:00', ncpu=4, print1=False, use_repol=True, keep_intermediates=False, alphabetical_pop_order='false', scratch_path="$SCRATCHDIR", scratch_gb="10"):
        '''Purpose:  Generate --multiSFS for fastsimcoal2 along with a given number of non-parametric block-bootstrapped replicates
           Notes: Must provide the block size for bootstrapping as well as number of bootstrap replicates
                  As of now, the necessary template files for FSC2 must come from elsewhere.  Beware of running this method with numerous populations'''

        if recode_dir.endswith("/") is False:
            recode_dir += "/"

        outdir = recode_dir + "FSC2input_" + output_name + "/"

        if os.path.exists(outdir) is False:
            os.mkdir(outdir)

        num_pops = len(pops)

        if os.path.exists(recode_dir) is True:
            if use_repol is True:
                suffix1 = '.table.repol.txt'
                suffix2 = 'repol.concat.txt'
            else:
                suffix1 = '.table.recode.txt'
                suffix2 = '.recode.concat.txt'
            concat_name = recode_dir + output_name + suffix2
            missing = []
            if os.path.exists(concat_name) is False:
                print("Concatenating input files")
                concat_file = open(recode_dir + output_name + suffix2, 'w')
                for pop in pops:
                    try:
                        with open(recode_dir + pop + suffix1) as infile:
                            for line in infile:
                                concat_file.write(line)
                    except FileNotFoundError:
                        missing.append(pop)
            if len(missing) != 0:
                print("Did not find input files for the following populations:", missing, ".  Aborting!!")
            else:
                print("Finished preparing input data")

                shfile4 = open(output_name + '.fsc2input.sh', 'w')

                shfile4.write('#!/bin/bash -e\n' +
                              '#PBS -N '+output_name+'\n' +
                              '#PBS -l walltime='+str(time_scratch)+'\n' +
                              '#PBS -l select=1:ncpus='+str(ncpu)+':mem='+str(mem)+'gb:scratch_local='+str(scratch_gb)+'gb\n' +
                              '#PBS -m abe\n' +
                              '#PBS -j oe\n\n' +
                              'module add python-3.6.2-gcc\n'+
                              'module add python36-modules-gcc\n'+
                              'trap \'clean_scratch\' TERM EXIT\n'+
                              'if [ ! -d "$SCRATCHDIR" ] ; then echo "Scratch not created!" 1>&2; exit 1; fi \n' +
                              'DATADIR="/storage/pruhonice1-ibot/home/holcovam/ScanTools/"\n' +
                              'cp $DATADIR/FSC2input.py $SCRATCHDIR || exit 1\n' +
                              'cp $DATADIR/'+ concat_name +' $SCRATCHDIR || exit 1\n'+
                              'cd $SCRATCHDIR || exit 2\n' +
                              'echo data loaded at `date`\n\n' +
                           #   'sort -k3,3 -k4,4g ' + output_name + suffix2 + ' | uniq -u > a.txt\n' +###ADDED 21Nov
                            #  'mv a.txt ' + output_name + suffix2 + '\n' +###ADDED 21Nov
                              'python3 FSC2input.py -i ' + output_name + suffix2 + ' -o ' + "FSC2input_" + output_name + ' -prefix ' + output_name + ' -ws ' + str(bootstrap_block_size) + ' -bs ' + str(bootstrap_reps) + ' -np ' + str(num_pops) + ' -alpha ' + str(alphabetical_pop_order) + '\n' +
                              'rm ' + output_name + suffix2 + '\n'+
                              'rm FSC2input.py\n'+
                              'cp $SCRATCHDIR/* $DATADIR/'+outdir+' || export CLEAN_SCRATCH=false\n'+
                              'printf "\\nFinished\\n\\n"\n' +
                              '#rm $DATADIR/'+concat_name +'\n')
#                if keep_intermediates is False:
#                    shfile4.write('rm ' + concat_name + "\n")
                shfile4.close()

                if print1 is False:
                    cmd1 = ('qsub ' + output_name + '.fsc2input.sh')
                    p1 = subprocess.Popen(cmd1, shell=True)
                    sts1 = os.waitpid(p1.pid, 0)[1]

                    self.log_file.write("###  Generate FSC2 Input  ###\n" +
                                        "Input Directory: " + recode_dir + "\n" +
                                        "Bootstrap Block Size: " + str(bootstrap_block_size) + "\n" +
                                        "Bootstrap Replicate Data Sets: " + str(bootstrap_reps) + "\n" +
                                        "Output Name: " + output_name + "\n" +
                                        "Populations: " + str(pops) + "\n")

                else:
                    file3 = open(output_name + '.fsc2input.sh', 'r')
                    data3 = file3.read()
                    print(data3)

                os.remove(output_name + '.fsc2input.sh')

        else:
            print("!!!Did not find recode_dir!!!!")


    def FSC2(self, input_dir, num_reps=50, min_sims=100000, max_ecm=20, calc_CI=False, time_scratch="01:50:00", scratch_mb=100, mem=100, numcores=1, print1=False, overwrite="None", fsc2_path="/storage/pruhonice1-ibot/home/holcovam/programs/fsc26_linux64/fsc26"):

        """This method parallelises job submission of fastsimcoal2, but requires a very specific set up of input files.  The output of '.generateFSC2input' should be a folder that contains the multi-dimensional SFS.  Place this folder in a new folder that will be the FSC2_Data_Parent_Directory.  This directory should also contain one or more template (.tpl) and estimates (.est) files whose format can be found in the fastsimcoal2 documentation.  For each sub-directory containing input data, this method will re-format and rename the .tpl and .est files to reflect the necessary information in the sub-directory multi-dimensional SFS and then submit these jobs to the cluster.  I've tried to make the code as general as possible, but this is one method that will likely require the user to read and understand the code in order to get things working well for them.  Also, a major potential source of errors is in the correct formatting of the .tpl and .est files, so it is worthwhile to ensure that these are correct (by running FSC2 on a subset of your sub-directories) before launching full-scale"""


        tpl_files = []
        est_files = []

        for path in os.listdir(input_dir):
            if os.path.isdir(input_dir + path) and path.startswith("FSC2input"):
                samp_name = path.split("_")[1]
                if samp_name + "_DSFS.obs" not in os.listdir(input_dir + path):
                    print("Did not find input data file for: ", samp_name)
                if calc_CI is True:
                    num_files = 0
                    for file in os.listdir(input_dir + path):
                        if file.endswith("_DSFS.obs") and file.split("_")[-2].split(".")[-1][0:3] == "rep" and file != samp_name + "_DSFS.obs":
                            num_files += 1
                    if num_files < 1:
                        print("Did not find bootstrap replicates for: ", samp_name)
                    else:
                        print("Found ", num_files, " replicate dsfs files for CI calculation for ", samp_name)
            if path.endswith(".tpl"):
                tpl_files.append(path)
                est_files.append(path.split(".")[0])
        if len(tpl_files) == 0:
            print("Did not find any tpl files!! Aborting!!")
        else:
            if 1 == 0:
                if any(os.path.exists(input_dir + k.split(".tpl")[0] + ".est") is False for k in tpl_files):
                    print("Did not find all est files.  Aborting!!")
            else:
                shfile5 = open("FSC2_submit.sh", 'w')
                shfile5.write('#!/bin/bash -e\n' +
                              'python3 FSC2_submit.py -i ' + input_dir + ' -reps ' + str(num_reps) + ' -minsims ' + str(min_sims) + ' -max_ecm ' + str(max_ecm) + ' -ci ' + str(calc_CI) + ' -nc ' + str(numcores) + ' -mem ' + str(mem) + ' -mb ' + str(scratch_mb) + ' -t ' + time_scratch + ' -print1 ' + str(print1) + ' -Ov ' + str(overwrite) + ' -fsc2path ' + fsc2_path + '\n')

                shfile5.close()
                if print1 is False:
                    cmd1 = ("bash FSC2_submit.sh")
                    p1 = subprocess.Popen(cmd1, shell=True)
                    #here maybe os.waipid?

                else:
                    file3 = open("FSC2_submit.sh", 'r')
                    data3 = file3.read()
                    print(data3)
                    
#                os.remove("FSC2_submit.sh")

    def gatherFSC2output(self, parent_dir):
        """This method collects all information from the '.bestlhood' output files of FSC2 that is buried in the sub-directories and outputs the information into one of two files:  Likelihoods file and parameters file."""
        if parent_dir.endswith("/") is False:
            parent_dir += "/"
        get_header = True
        dirname = parent_dir.strip("/").split("/")[-1]
        Lhoodfile = open(parent_dir + dirname + "_FSC2_Likelihoods.txt", 'w')
        paramsfile = open(parent_dir + dirname + "_FSC2_Params.txt", 'w')
        for root, dirs, files in os.walk(parent_dir):
            for file in files:
                if file.endswith(".bestlhoods"):
                    name = file.split(".best")[0]
                    samp_names = name.split("_")[0]
                    model = ".".join(name.split("_")[1:])
                    if get_header is True:
                        get_header = False
                        with open(os.path.join(root, file), 'r') as ff:
                            for i, line in enumerate(ff):
                                if i == 0:
                                    Lhoodfile.write("Model\tSampleNames\tnum_params\tLhood_est\tLhood_obs\tLhood_diff\tAIC\n")
                                else:
                                    info = line.strip("\n").split("\t")
                                    num_params = len(info) - 2
                                    Lhood_est = float(info[-2])
                                    AIC = (2 * num_params) - (2 * (Lhood_est * math.log(10)))
                                    Lhoodfile.write(model + "\t" + samp_names + "\t" + str(num_params) + "\t" + str(Lhood_est) + "\t" + info[-1] + "\t" + str(Lhood_est - float(info[-1])) + "\t" + str(AIC) + "\n")
                                    paramsfile.write(model + "\t" + samp_names + "\t" + line)
                    else:
                        with open(os.path.join(root, file), 'r') as ff:
                            for i, line in enumerate(ff):
                                if i > 0:
                                    info = line.strip("\n").split("\t")
                                    num_params = len(info) - 2
                                    Lhood_est = float(info[-2])
                                    AIC = (2 * num_params) - (2 * (Lhood_est * math.log(10)))
                                    Lhoodfile.write(model + "\t" + samp_names + "\t" + str(num_params) + "\t" + str(Lhood_est) + "\t" + info[-1] + "\t" + str(Lhood_est - float(info[-1])) + "\t" + str(AIC) + "\n")
                                    paramsfile.write(model + "\t" + samp_names + "\t" + line)
        Lhoodfile.close()
        paramsfile.close()


    def queryFSC2input(self, input_dir, index_list, outname):
        """PROVIDE SUMMARY HERE AND IN README"""
        outfile = open(input_dir + outname + ".txt", 'w')
        outfile.write("Outname\tPops\tProp\n")
        for path in os.listdir(input_dir):
            if os.path.isdir(input_dir + path) and path.startswith("FSC2input"):
                samp_name = path.split("_")[1]
                if samp_name + "_DSFS.obs" in os.listdir(input_dir + path):
                    with open(input_dir + path + '/' + samp_name + "_DSFS.obs") as fsc2input:
                        for i, line in enumerate(fsc2input):
                            if i == 2:
                                line = line.strip("\n").strip("\t").split("\t")
                                tot = sum([int(j) for j in line])
                                sp = 0  # Shared polymorphisms
                                for ix in index_list:
                                    sp += int(line[ix])
                    outfile.write(outname + '\t' + samp_name + '\t' + str(float(sp) / float(tot)) + '\n')
