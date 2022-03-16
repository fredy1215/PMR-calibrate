# C:\Users\Meng Lu\Box\Thesis\code\app
import scipy.io as sco
import os
import fnmatch
import numpy as np

def find_files(filename, search_path):
 result = []

 if filename.find('*') != -1:

   for root, dir, files in os.walk(search_path):
       for a in files:
         if fnmatch.fnmatch(a,filename):
            result.append(os.path.join(root, a))
   return result
 else:
   for root, dir, files in os.walk(search_path):
         if filename in files:
             result.append(os.path.join(root, filename))
   return result
File = find_files('*.mat','C:\\Users\\Meng Lu\\Box\\Thesis\\code\\app')
matstruct_contents_1_5 = sco.loadmat(File[0])
matstruct_contents_3 = sco.loadmat(File[1])
TI_1_5 = matstruct_contents_1_5['TI']
TI_3 = matstruct_contents_3['TI']
H = list(range(40,100+1,1))
Nflash = list(range(20,50+1,1))
FA = list(range(10,20+1,1))
Echo = list(np.arange(3,6+0.1,0.1))

