import sys
import os
import subprocess

def mal_string(parnum,parval,correct=True):
    mal="mal:"
    for i in range(6): mal+=("%f%s" % ((parval if parnum==i else 0),("," if i<5 else "")))
    if correct: mal+="::cor"
    return mal

#directories
redo=0
otag="20160422"
code_dir="/n/atlasfs/atlascode/backedup/stchan/mmtrigger/"
exec_dir=code_dir+"Root/"
ntup_dir=code_dir+"ntuple/dlm_new/"
hist_dir="/n/atlasfs/atlasdata/atlasdata1/stchan/mmtrigger/raw_hist/"+otag+"/"
scri_dir=code_dir+"scripts/"+otag+"/"
scr_pat="%s%s_%s.sh" #script directory,output tag, spec tag (the string that's a chain of different property strings)
slurm_cmd=exec_dir+"MMtriggerSim"
h_vmm=0.445*64/7583.5 #slope road size corresponding to one VMM chip (64 x 0.445 mm strips at the first z plane (z=7583.5 mm)
n_per_vmm=[-999,1,2,4]
#vmm_fat=True

#MMtriggerSim command line format: MMtriggerSim [ntuple_file] [output_dir] [output_tag] {property:value_for_par_par()}

#variables we want to scan across
file_pattern="%sNSWPRDValAlgPt%iGeV.digi.ntuple.root" # will be % (ntup_dir,pT_value)
pTs=[10,20,30,50,100,200]
bgs=[0,1]#"qbg:0","qbg:1"]
charges=[0,1,2]
CTs = [[2,1],[2,2],[3,2],[3,3],[4,4]]
def_CT=[4,4]
mal_bounds=[[0,5],[0,5],[0,5],[-0.0015,0.0015],[0,0.0015],[0,0.0015]]
#mal_bounds=[[0,5]]#,[-0.0015,0.0015]]
#mal_bounds=[[0.,0.0015],[0.,0.0015]]
pnum=0
ndiv_mal=20
toggle=[True]

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

for pt in pTs:
    for c in CTs:
        for bg in range(8):
            bgtag="qbg:%i"%bg
            snm=scr_pat%(scri_dir,otag,bgtag+"_pt%i"%pt+"_ctx%i"%c[0]+"_ctuv%i"%c[1])
            scripts.append(snm)
            script=open(snm,"w")
            script.write("#!/bin/bash\n")
            script.write(slurm_cmd+" "+file_pattern%(ntup_dir,pt)+" "+hist_dir+" "+otag+" "+bgtag+(" colsk:%i"%bg)+" xct:%i"%c[0]+" uvct:%i"%c[1]+" redo:%i"%redo+" \n")
            script.close()

mem="8192" #16000
for scr in scripts:
    submit=subprocess.Popen(["sbatch","--mem-per-cpu",mem,"--time","3-0:0:0","-p","pleiades",scr])
#    submit.wait()


"""    
for pT in pTs:
    for qt in charges:
        for ct in CTs:
            submit=subprocess.Popen(launch_cmd+[file_pattern % (ntup_dir,pT),hist_dir,otag,("qt:%i" % qt),("xct:%i" % ct[0]),("uvct:%i" % ct[1])])
            submit.wait()
"""
