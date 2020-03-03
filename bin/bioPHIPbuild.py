# -*- coding: utf-8 -*-
#!/usr/bin/env python
"""
Created on Thu Oct 25 10:06:36 2018

@author: yuan
"""

import os
import sys
import numpy as np
import pandas as pd
import tables #used for hdf5 expport
#personal modules
import myAlign
import myCommon
import myDataframe
import myDict
import myGenome
import myIO
import myParallel
import mySystem
import bioPHIPfunc

sys.path.insert(0, '/home/yuan/phip/bin')
###############################################################################
class phip_annot:
    def __init__(self, par):
        self.par = par

#
    def human(self):
        #export
        self.par['lib']='human'
        self.par['file_hdf5'] = self.par['dir_hdf5'] + self.par['lib'] + '.h5'
        store=pd.HDFStore(self.par['file_hdf5'], 'w')
        #input
        self.par['file_ref_fa'] = self.par['dir_ref_seq'] + 'T7Pep2_human.fa'
        self.par['file_annotation'] = self.par['dir_ref_seq'] + 'T7Pep2_human_annot.txt' 
        self.par['file_NC'] = self.par['dir_ref_seq'] + 'human_pep_RC_BeadsOnly.txt'
                
        #annotation file
        print(self.par['file_annotation'])
        annot_df = pd.read_table(self.par['file_annotation'], index_col=False)
        #print(annot_df)
        store.append('annot_df', annot_df)
        
        
#

#end

###############################################################################
#main
par={'dir_phip':'/home/yuan/phip/'}
par['dir_ref_seq']=par['dir_phip']+'ref_seq/'
par['dir_aligner']=par['dir_phip']+'bowtie1/'
par['dir_hdf5']=par['dir_phip']+'hdf5/'


#human
phip_annot(par).human()




annot_df = pd.read_table('/home/yuan/phip/ref_seq/T7Pep2_human_annot.txt', index_col=False)
store=pd.HDFStore('/home/yuan/phip/hdf5/human.f5', 'w')
store.append('annot_df', annot_df)

