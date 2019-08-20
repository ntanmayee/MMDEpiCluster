import argparse
import sys
import os
import csv
import re

def main():    
    class C:
        pass
    opts = C();
    parser = argparse.ArgumentParser(description='Parse SRA Run Table')
    parser.add_argument('runtable', type=str, help="Path to runtable")
    parser.parse_args(args=sys.argv[1:], namespace=opts)
    
    k4me3 ="H3K4me3"
    k27ac ="H3K27ac"
    
    abdict = dict()
    
    with open(opts.runtable,"r") as fd:
        reader = csv.DictReader(fd, delimiter="\t", quotechar='"')
        for row in reader:
            ab1 = re.search("("+k4me3+")|("+k27ac+")",row["chip_antibody1"]).group(0)
            ab2 = re.search("("+k4me3+")|("+k27ac+")",row["chip_antibody2"]).group(0)
                        
            abl = [ab1,ab2]
            #abl.sort() 
            dirname = abl[0]+"_"+abl[1]
            
            if not dirname in abdict:
                abdict[dirname] = dict()
            
            if not row["Experiment"] in abdict[dirname]:
                abdict[dirname][row["Experiment"]] = []
            
            abdict[dirname][row["Experiment"]].append(row["Run"])
    
    for abl, exps in abdict.items():
        if not os.path.exists(abl+"/"):
            os.makedirs(abl)
        
        for srx, srrs in exps.items():
            if not os.path.exists(abl+"/"+srx):
                os.makedirs(abl+"/"+srx)
            
            with open(abl+"/"+srx+"/SRR_Acc_List.txt", mode='wt', encoding='utf-8') as myfile:
                myfile.write('\n'.join(srrs))            

if __name__ == '__main__':
    main()