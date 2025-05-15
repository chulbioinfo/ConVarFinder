from gprofiler import gprofiler
import pandas
import sys
import re
import os
import string
import glob
import queue
import threading
import time

def gprofiler_enrichment(gene_list,oNAME):
  enrichment = gprofiler(gene_list, organism='hsapiens')
  try:
    enrichment.to_csv(oNAME, sep='\t')
  except:
    if len(enrichment)==0:
      fpout = open(oNAME,'w')
      fpout.close()
    else:
      print(oNAME)
      print(enrichment)
      sys.exit()



#################################
#################################
# main
#################################

queue = queue.queue()

class ThreadPrank(threading.Thread):
  def __init__(self,queue):
    threading.Thread.__init__(self)
    self.queue = queue
  def run(self):
    while True:
      info = self.queue.get()
      gene_list = info[0]
      oNAME = info[1]
      if len(gene_list)>0:
        gprofiler_enrichment(gene_list,oNAME)
      else:
        fpout = open(oNAME,'w')
        fpout.close()
        self.queue.task_done()

start = time.time()

def main():
  for i in range(6):
    t = ThreadPrank(queue)
    t.setDaemon(True)
    t.start()

  for datatype_PATH in ["iCSAV"]:

      iPATH = "gp_in/" + datatype_PATH + varianttype_PATH
      oPATH = "gp_out/" + datatype_PATH + varianttype_PATH
      flist = glob.glob(iPATH+"*.txt")
      for fNAME in flist:
        oNAME = oPATH + os.path.basename(fNAME)
        fpin = open(fNAME,'r')
        gene_list = []
        for line in fpin:
          gene_list.append(line.strip())
        fpin.close()
        info = [gene_list,oNAME]
        queue.put(info)
  queue.join()
main()
print("Elapsed Time: %s" % (time.time() - start))
