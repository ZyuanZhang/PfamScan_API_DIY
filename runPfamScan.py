# 调用 PfamScan API 鉴定蛋白的domain/motif
# 注意:每个fasta文件中的序列不超过100 (如果超过100，则会先将fasta分割成多个文件)

import os
from Bio import SeqIO
import time
import sys

def runPfamScan(inFasta, outFile):
    """
        inFasta中序列数目 <= 100
    """
    cmd = 'python pfamscan.py --email "gakkijohn@gmail.com" --database pfam-a --sequence '+inFasta+' --evalue 50 --format txt --outfile '+outFile
    os.system(cmd)


def splitFasta(infile, threshold=100):
    """
       将fasta文件中的序列拆分成多个fasta文件，每个fasta文件中的序列数目<=100
       threshold: 拆分完的fasta文件中蛋白序列的数目的阈值，默认为100
    """
    id_list, seq_list = [], []
    for record in SeqIO.parse(infile, "fasta"):
        id_list.append(str(record.id))
        seq_list.append(str(record.seq))

    if len(id_list) <= threshold:
        return "NoSplit"
    else:
        for i in range(0, len(id_list), threshold):
            if i+threshold <= len(id_list):
                outF = open("/".join(infile.split("/")[0:-1])+"/part_"+infile.split("/")[-1].split(".")[0]+"_"+str(i+1)+"_"+str(i+threshold)+".fasta", "w")
                for j in range(i, i+threshold):
                    outF.write(">%s\n%s\n" % (id_list[j], seq_list[j]))
                outF.close()
            else:
                outF = open("/".join(infile.split("/")[0:-1])+"/part_"+infile.split("/")[-1].split(".")[0]+"_"+str(i+1)+"_"+str(len(id_list))+".fasta", "w")
                for k in range(i, len(id_list)):
                    outF.write(">%s\n%s\n" % (id_list[k], seq_list[k]))
                outF.close()

            
def main(infileFasta):
    res = splitFasta(infile=infileFasta)
    if res == "NoSplit":
        runPfamScan(inFasta = infileFasta,
                    outFile = "./results/"+infileFasta.split("/")[-1].split(".")[0])
    else:
        fastaFile_list = ["./infileFasta/"+f for f in os.listdir("./infileFasta") if f.startswith("part_")]
        for inf in fastaFile_list:
            runPfamScan(inFasta = inf,
                        outFile = "./results/"+inf.split("/")[-1].split(".")[0])


if __name__ == "__main__":
    print("Start: ", time.ctime(), flush=True)
    main(infileFasta = sys.argv[1])
    print("End: ", time.ctime(), flush=True)
