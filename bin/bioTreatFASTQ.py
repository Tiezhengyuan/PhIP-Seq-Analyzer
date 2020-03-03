# -*- coding: utf-8 -*-
"""
Created on Wed Jan 20 21:47:13 2016

@author: yuan
"""
#standard modules
import getopt
import os
#import re
import sys
import shutil

#personal modules
import myDict
import myGenome
import myIO
import myParallel


def par_command(argv):
    phip_libs = ['human', 'virus', 'allergome', 'provirome', 'toxome', 'mouse', 'PE', 'zika', 'arbo', 'LISH']
    #initiate parameter
    na_str='fq_file,barcode_file,index_file,I1_file,I2_file,dir_raw_data,dir_in,out,dir_result'
    par=dict([(key, 'NA') for key in na_str.split(',')])
    par.update({'ref_libs':phip_libs[:2], 'seq_start':0, 'seq_end':None, 'seq_min':10})
    usage_out = 'Usage:\n' + argv[0] + ' [options] -o <raw data directory> ' + \
                '-f <fastq file> -i <index file> -b <barcode file>\n'
    try:
        opts, args = getopt.getopt(argv[1:],"hf:i:b:o:t:r:l:x:y:m:n:c:",["help",\
            "fastq_file", "index_file", "barcode_file", "dir_raw_data", "trim_5end", 'len_trim',\
            'fixed_end5', 'dir_in', 'out', 'I1_file','I2_file','ref_library'])
    except getopt.GetoptError:
        print(usage_out)
        sys.exit(2)
      
    #get parameters 
    for opt, arg in opts:
        if opt in ('-h', '--help'):
            print(usage_out)
            #common usage
            # python Process_FASTQ.py -f * -i * -b * -o * -y *" 
            print("-h --help\tUsage information of this script.")
            print("-t --trim_len\tTrim sequences from the 5'-end or 3'-end of reads (Optional)")
            print("-f --fastq_file\tFastq file determined by a sequencing analyzer.")
            print("-i --index_file\tIndex file matched with the fastq file.")
            print("-b --barcode_file\tBarcode file matched with the index file.")
            print("-o --raw_data\tDirectory storing demulitplexed *fastq files.")
            print("-y --out\tDirectory storing sample_info.csv and variables.txt.")
            print("-c --ref_library\tReference libraries can be any of {}, default is {}.".format(phip_libs, phip_libs[:2]))
            sys.exit(2)
        elif opt in ("-f", "--fastq_file"):
            par['fq_file'] = os.path.abspath(arg)
        elif opt in ("-i", "--index_file"):
            par['index_file'] = os.path.abspath(arg)
        elif opt in ("-b", "--barcode_file"):
            par['barcode_file'] = os.path.abspath(arg)
        elif opt in ("-o", "--raw_data"):
            par['dir_raw_data'] = myIO.dir_os(os.path.abspath(arg)).create_dir()
        elif opt in ('-x', "--dir_in"):
            par['dir_in'] = os.path.abspath(arg)
            par['fq_files'] = myParallel.samples({}).seek_fq(par['dir_in'])
        elif opt in ('-y', "--out"):
            par['out'] = arg
        elif opt in ("-l", "--min_len"):
            # discard shorter reads due to poor sequencing
            par['seq_min'] = abs(int(arg))
        elif opt in ("-t", "--trim_5end"):
            #trim_end5: length of nt from the 5-end
            par['seq_start'] = abs(int(arg))
        elif opt in ("-r", "--fixed_len"):
            #len_trim: length of nt after trimming 5-end and 3-end
            par['seq_len'] = abs(int(arg))
            par['seq_end'] = par['seq_start'] + par['seq_len']
        elif opt in ("-m", "--I1_file"):
            par['I1_file'] = os.path.abspath(arg)
        elif opt in ("-n", "--I2_file"):
            par['I2_file'] = os.path.abspath(arg)
        elif opt in ("-c", "--ref_library"):
            libs = arg.split(',')
            par['ref_libs'] = [x for x in libs if x in phip_libs]
    #   
    myDict.basic(par).print_dict()
    return par
   
###################
#legalize the parameters
def judge_par(par):   
   #judge illegal parameters: 
   #1: demultiplexing mode: 
   #fastq file
   if not os.path.isfile(par['fq_file']):
      print('Warning: No FASTQ file ', par['fq_file'])
   #barcode file
   if not os.path.isfile(par['barcode_file']):
      print('Warning: No barcode file ', par['barcode_file'])
   #index file
   if not os.path.isfile(par['index_file']):
      print("Warning: NO index files ", par['index_file'])
   if not (os.path.isfile(par['I1_file']) and os.path.isfile(par['I2_file'])):
      print("Warning: NO Index file: {} and {}".format(par['I1_file'], par['I2_file']))

   #trim mode
   if 'fq_files' in par.keys():
       if len(par['fq_files']) == 0:
           print("Error: No fastq files were detected under ", par['dir_in'])
           sys.exit(2)
       if par['seq_start'] == 0 and par['seq_end'] == 0:
           print("Error: Enter the trimmed length !")
           sys.exit(2)
   #3: raw data dir
   if not os.path.isdir(par['dir_raw_data']):
      print("Error: No directory storing seperated fastq files", par['dir_raw_data'])
      sys.exit(2)
       
   #no return
   
##############################
if __name__ == "__main__":
    #home_dir=os.path.expanduser("~")+'/'
    #get parameters from command line
    par = par_command(sys.argv)
    #judge parameters
    judge_par(par)
    
    #combine index files if par['index_file'] or par['I1_file'] and par['I2_file']
    #if os.path.isfile(par['I1_file']) and os.path.isfile(par['I2_file']):
    #    myGenome.genome(par['I1_file']).cbind_fq(par['I2_file'], par['index_file'])

    #demultiplexing: split fastq files based on barcode
    if os.path.isfile(par['index_file']):
        myGenome.genome(par['fq_file']).decompose_fq(par)
    elif os.path.isfile(par['I1_file']) and os.path.isfile(par['I2_file']):
        myGenome.genome(par['fq_file']).decompose_fq2(par)
        
    #trim fastq files
    if 'fq_files' in par.keys():
        for fq in par['fq_files']:
            myGenome.genome(fq).trim_fq(par['dir_raw_data'], par['seq_start'], par['seq_end'])
            
    #generate sample_info file under result dir: 
    if par['out'] != 'NA':
        #current dir
        par['dir_bin'] = os.path.dirname(os.path.realpath(__file__)) + '/'
        par['dir_home'] = os.path.abspath(os.path.join(par['dir_bin'], os.pardir)) + '/'
        print('Home directory of phip pipsline: ', par['dir_home'])
        
        #libraries. default is human and virus
        for lib in par['ref_libs']:
            par['dir_result'] = myIO.dir_os(os.path.abspath(par['out']+'_'+lib)).create_dir()
            if os.path.isdir(par['dir_result']):
                #1: sample_info.csv
                par['file_sample_info'] = par['dir_result'] + 'sample_info.csv'
                print('The sample information file: ', par['file_sample_info'])
                #read sample_info.csv
                myParallel.samples(par).export_sample_info()
                #2: copy template variables.txt into lib folder
                template_file='{}variables_{}.txt'.format(par['dir_bin'], lib)
                var_file = '{}variables.txt'.format(par['dir_result'])
                print('Save {} and then update it.'.format(var_file))
                shutil.copy(template_file, var_file)
                #update parameters of variables.txt
                refresh = {'dir_home':par['dir_home'], 'dir_result':par['dir_result'], \
                           'dir_raw_data':par['dir_raw_data']}
                myIO.file_os(var_file, '=').line_replace(refresh)
            else:
                print('Warning: No results directory (ref_lib:{}): {}.'.format(lib, par['dir_result']))

    #
    print('\n\nThe raw data pre-processing is done!\n\n')
    
#
#test simple mode at local computer
#python ~/phip/bin/Process_FASTQ.py -x ~/raw_data/LISH -o ~/raw_data/LISH/test_raw_data -y ~/raw_data/LISH/LISH2
#normal mode at local computer
#python ~/phip/bin/bioTreatFASTQ.py -f ~/raw_data/LISH/Undetermined_S0_L002_R1_001.fastq -i ~/raw_data/LISH/Undetermined_S0_L002_I1_001.fastq -b ~/raw_data/LISH/LISH_barcode.txt  -o ~/raw_data/LISH/test_raw_data -y ~/raw_data/LISH/test -c human,virus,PE,LISH