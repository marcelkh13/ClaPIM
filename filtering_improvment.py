import sys
import ast
from random import randint
from math import ceil

# defual values
#kmer_len= 64
#kmers_inMAT=128
#eth = 4
#num_queries=1000

filename= str(sys.argv[1]) # *.fna file containing the genome
kmer_len = int(sys.argv[2]) # the length of the kmers to extract from the genome
kmers_inMAT = int(sys.argv[3]) # number of kmers to fit in each crossbar array
eth = int(sys.argv[4]) # the filter error threshold
num_queries = int(sys.argv[5]) # how many random queries to check. a larger number provides more precise average


ref=""
with open(filename) as f:
  for l in f:
    if l[0]!=">":
      ref= "".join([ref,l.replace("\n","")])
f.close()


ref_kmers={}
r = len(ref) - kmer_len + 1
kmers_num=0

#extracting the kmers from the reference genome
for i in range(r):
  kmer = ref[i:i+kmer_len]
  if kmer.find('N') == -1:
    if ref_kmers.get(kmer) == None:
      #this is a new kmer, build  its histrogram
      kmers_num += 1
      hist={"A":0 , "T":0, "G":0, "C":0}
      for c in range(kmer_len):
        if kmer[c]=="A" or kmer[c]=="a":
          hist["A"]+=1
        if kmer[c]=="T" or kmer[c]=="t":
          hist["T"]+=1
        if kmer[c]=="G" or kmer[c]=="g":
          hist["G"]+=1
        if kmer[c]=="C" or kmer[c]=="c":
          hist["C"]+=1
      # update hist
      ref_kmers[kmer] = hist



#ref_kmers contains all different kmers in the genome as keys and their histogram as values

hist_list = ref_kmers.values() 
hist_groups = {} 

#for each hist, count the number of kmers with the same hist
for hist in hist_list:
  key = str(hist)
  hist_groups[key]= hist_groups.get(key,0) +1



num_checks_avr=0
num_checks_avr_b=0


checked=0

for i in range(num_queries):
  tot_MATs=0
  checked_MATs=0
  
  #generate random query hist
  Aq=randint(0,kmer_len)
  Tq=randint(0,kmer_len-Aq)
  Gq=randint(0,kmer_len-Aq-Tq)
  Cq=kmer_len-Aq-Tq-Gq
  
  for hist_str,count in hist_groups.items():
    hist = ast.literal_eval(hist_str) 
    #kmer hist
    Ak = int(hist["A"])
    Tk = int(hist["T"])
    Gk = int(hist["G"])
    Ck = int(hist["C"])
    tot_MATs += ( count/kmers_inMAT +1*(count%kmers_inMAT!=0) )#ceil(float(count/kmers_inMAT))
    if ( abs(Aq-Ak) + abs(Tq-Tk) + abs(Gq-Gk) + abs(Cq-Ck) ) < 2*eth+1 :    
      checked_MATs += ( count/kmers_inMAT +1*(count%kmers_inMAT!=0) ) 
      
    
  checked +=  checked_MATs

checked_avr= float(checked)/num_queries


if checked_avr!=0 and tot_MATs!=0:
  print("The power improvment is " , float(tot_MATs)/checked_avr )
else:
  print("ERROR: Try to change the number of queries")

print("%%DONE%%")


















    
    
    
    
    
    
    
  