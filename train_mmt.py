import sys
import os
import subprocess

def mal_string(parnum,parval,correct=True):
    mal="mal:"
    for i in range(6): mal+=("%f%s" % ((parval if parnum==i else 0),("," if i<5 else "")))
    if correct: mal+="::cor"
    return mal

#directories
#this bool: true to make the training; false to use trained correction
train=False
dtag="20151111"
code_dir="/n/atlasfs/atlascode/backedup/stchan/mmtrigger/"
exec_dir=code_dir+"Root/"
ntup_dir=code_dir+"ntuple/dlm_new/"
hist_dir="/n/atlasfs/atlasdata/atlasdata1/stchan/mmtrigger/raw_hist/"+dtag+"/"
scri_dir=code_dir+"scripts/"+dtag+"/"
otag="train"if train else dtag
scr_pat="%s%s_%s.sh" #script directory,output tag, spec tag (the string that's a chain of different property strings)
slurm_cmd=exec_dir+"MMtriggerSim"

#MMtriggerSim command line format: MMtriggerSim [ntuple_file] [output_dir] [output_tag] {property:value_for_par_par()}

#variables we want to scan across
file_pattern="%sNSWPRDValAlgPt%iGeV.digi.ntuple.root" # will be % (ntup_dir,pT_value)
pTs=[10,20,30,50,100,200]
charges=[0,1,2]
CTs = [[2,1],[2,2],[3,2],[3,3],[4,4]]
def_CT=[4,4]
mal_bounds=[[0,5],[0,5],[0,5],[-0.0015,0.0015],[0,0.0015],[0,0.0015]]
#mal_bounds=[[0,5]]#,[-0.0015,0.0015]]
#mal_bounds=[[0.,0.0015],[0.,0.0015]]
pnum=0
ndiv_mal=20
toggle=[not train]
malpts=[200] if train else [20,100]
if not os.path.isdir(code_dir): sys.exit("Code directory %s does not exist. Aborting!" % (code_dir))
if not os.path.isdir(exec_dir): sys.exit("Executable's directory %s does not exist. Aborting!" % (exec_dir))
if not os.path.isdir(ntup_dir): sys.exit("Ntuple directory %s does not exist. Aborting!" % (ntup_dir))
if not os.path.isdir(scri_dir):
    os.makedirs(scri_dir)
    print("Created output directory %s" % scri_dir)
if not os.path.isdir(hist_dir):
    os.makedirs(hist_dir)
    print("Created output directory %s" % hist_dir)
scripts=[]

for mal in mal_bounds:
    increment=(mal[1]-mal[0])/(1.*ndiv_mal)
    for i in range(ndiv_mal+1):
        for pt in malpts:
            for tf in toggle:
                val=mal[0]+i*increment
                snm=scr_pat%(scri_dir,otag,mal_string(pnum,val,tf))
                scripts.append(snm)
                script=open(snm,"w")
                script.write("#!/bin/bash\n")
                if train: script.write(slurm_cmd+" "+file_pattern%(ntup_dir,pt)+" "+hist_dir+" "+otag+" "+mal_string(pnum,val,tf)+" xct:2 uvct:1")
                else: script.write(slurm_cmd+" "+file_pattern%(ntup_dir,pt)+" "+hist_dir+" "+otag+" "+mal_string(pnum,val,tf)+" 3cor:/n/atlasfs/atlasdata/atlasdata1/stchan/mmtrigger/raw_hist/"+dtag+"/,Pt200GeV.digi.ntuple_train")
                script.close()
    pnum+=1

mem="2048" #"16000" 
for scr in scripts:
#    print(["sbatch","--mem-per-cpu","16000","--time","1-0:0:0","-p","pleiades","--error=%sslurmout/train-%%j.out"%exec_dir,"--output=%sslurmout/train-%%j.out"%exec_dir,scr])
    submit=subprocess.Popen(["sbatch","--mem-per-cpu",mem,"--time","1-0:0:0","-p","pleiades",scr])
    submit.wait()


"""    
for pT in pTs:
    for qt in charges:
        for ct in CTs:
            submit=subprocess.Popen(launch_cmd+[file_pattern % (ntup_dir,pT),hist_dir,otag,("qt:%i" % qt),("xct:%i" % ct[0]),("uvct:%i" % ct[1])])
            submit.wait()
"""
