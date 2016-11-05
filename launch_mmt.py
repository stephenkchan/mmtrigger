import sys
import os
import subprocess

def mal_string(parnum,parval,correct=True):
    mal="mal:"
    for i in range(6): mal+=("%i%s" % ((parval if parnum==i else 0),("," if i<5 else "")))
    if correct: mal+="::cor"
    return mal

#directories
code_dir="/n/atlascode/backedup/stchan/mmtrigger/"
exec_dir=code_dir+"Root/"
ntup_dir=code_dir+"ntuple/dlm_new/"
hist_dir="/n/atlasdata1/stchan/mmtrigger/raw_hist/20150703/"
otag="15etabins"
launch_cmd=["bsub","-q","pleiades",exec_dir+"MMtriggerSim"]

#MMtriggerSim command line format: MMtriggerSim [ntuple_file] [output_dir] [output_tag] {property:value_for_par_par()}

#variables we want to scan across
file_pattern="%sNSWPRDValAlgPt%iGeV.digi.ntuple.root" # will be % (ntup_dir,pT_value)
pTs=[10,30,50,100,200]
charges=[0,1,2]
CTs = [[2,1],[2,2],[3,2],[3,3],[4,4]]

if not os.path.isdir(code_dir): sys.exit("Code directory %s does not exist. Aborting!" % (code_dir))
if not os.path.isdir(exec_dir): sys.exit("Executable's directory %s does not exist. Aborting!" % (exec_dir))
if not os.path.isdir(ntup_dir): sys.exit("Ntuple directory %s does not exist. Aborting!" % (ntup_dir))
if not os.path.isdir(hist_dir):
    os.makedirs(hist_dir)
    print("Created output directory %s" % hist_dir)

for pT in pTs:
    for qt in charges:
        for ct in CTs:
            submit=subprocess.Popen(launch_cmd+[file_pattern % (ntup_dir,pT),hist_dir,otag,("qt:%i" % qt),("xct:%i" % ct[0]),("uvct:%i" % ct[1])])
            submit.wait()
