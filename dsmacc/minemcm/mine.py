'''Mining pogram for the mcm. Daniel.Ellis.Research@googlemail.com 2018'''
import re,glob,urllib2
import numpy as np
from mpi4py import MPI


comm = MPI.COMM_WORLD

inchis = re.compile(r'<span class="inchi">(.*)<') 
smiles =  re.compile(r'<span class="smiles">(.*)<') 
synonyms =  re.compile(r'<span>\n\t\s*([\w;\-\s]*);\n')

files = glob.glob('*.txt')

def scour(url,spec,mcm): 
                st= urllib2.urlopen(url).read()
                return [spec,smiles.findall(st)[0],inchis.findall(st)[0],';'.join(synonyms.findall(st)),mcm] 




if comm.rank == 0:
    print " Running on %d cores" % comm.size
    gather = []
   
for f in range(len(files)):
    mcm = re.findall(r'MCM (v[\d\.]+)',str(tuple(open(files[f]))[:20]))[0]

    if comm.rank == 0:    
        species = re.findall(r'\d+[\\t \s]+\b([A-z]\w+)\b \\n',str(tuple(open(files[f]))))
        split = np.array_split(species,comm.size-1)
        
        for i in xrange(comm.size-1):
            comm.send(split[i], dest=i+1, tag=f)

    else:
        data = comm.recv(source=0, tag=f)
        
        data = data[:2]
        print data
        ret  = []
        for spec in data:
            try: 
                ret.append(scour('http://mcm.leeds.ac.uk/MCM/browse.htt?species=%s'%spec,spec,mcm))
            except: 
                try: 
                    ret.append(scour('http://mcm.leeds.ac.uk/MCM/browse.htt?species=%s'%spec,spec,mcm))
                except: 
                    try: 
                        ret.append(scour('http://mcm.leeds.ac.uk/MCMv3.2/browse.htt?species=%s'%spec,spec,mcm))
                      
                        

                    except Exception as e: 
                        print 'Failed on: ' + spec   , e  
                
        comm.send(ret, dest=0, tag=f)
        

    if comm.rank == 0:  
        for i in xrange(comm.size-1):
            gather.extend( comm.recv(source=i+1, tag=f ) )
           

        

   
if comm.rank == 0:     
   
   gather = sorted(gather, key=lambda e: e[0], reverse=False)
   
   import pandas as pd
   
   df = pd.DataFrame(gather,columns= ['name','smiles','inchi','synonyms','v'])

   df.replace('Exception: (1064, "You have an error in your SQL syntax.  Check the manual that corresponds to your MySQL server version for the right syntax to use near \'%s\' at line 1")','',inplace=True)
   
   df['v'][df.duplicated('name')]='v3.2 & v3.3.1'
   df=df.drop_duplicates('name',keep='last')
   print df.head()
   
   df.to_csv('smiles_mined.csv')
       
        
