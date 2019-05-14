import os, sys
from optparse import OptionParser

def get_options():
  parser = OptionParser()
  parser.add_option('--input_type', dest='input_type', default='electron', help="Input type [electron,neutrino,pion,photon]" )
  parser.add_option('--clusteringAlgo', dest='clusteringAlgo', default='Histomax_vardr', help="Clustering algorithm" )
  parser.add_option('--geometry', dest='geometry', default='v9', help="HGCal geometry configuration" )
  parser.add_option('--numberOfFiles', dest='numberOfFiles', default=-1, help="Number of files to process" )
  parser.add_option('--queue', dest='queue', default='microcentury', help="HTCondor Queue" )
  return parser.parse_args()

(opt,args) = get_options()

# Define path
path = os.environ['CMSSW_BASE'] + "/src/L1Trigger/egid_analysis/HGCal_L1T_egammaID/run"
f_sub_name = "submit_%s_%s_%s.sub"%(opt.input_type,opt.clusteringAlgo,opt.geometry)
f_out = "%s_%s_$(procID)_%s"%(opt.input_type,opt.clusteringAlgo,opt.geometry)

#Dictionary storing number of files in each directory
f_number = {"electron_v9":400,"photon_v9":400,"pion_v9":394,"neutrino_v9":2599,"electron_v8":270,"neutrino_v8":623}

# If user input number of files = -1, set to f_number
if( opt.numberOfFiles == -1 ): numberOfFiles = f_number["%s_%s"%(opt.input_type,opt.geometry)]
else: numberOfFiles = opt.numberOfFiles

#Create condor submission file
f_sub = open("%s"%f_sub_name,"w+")
f_sub.write("plusone = $(Process) + 1\n")
f_sub.write("procID = $INT(plusone,%d)\n\n")
f_sub.write("executable          = %s/run_cl3d_selection.sh\n"%path)
f_sub.write("arguments           = %s $(procID) %s %s\n"%(opt.input_type,opt.clusteringAlgo,opt.geometry))
f_sub.write("output              = %s/jobs/out/%s.out\n"%(path,f_out))
f_sub.write("error               = %s/jobs/err/%s.err\n"%(path,f_out))
f_sub.write("log                 = %s/jobs/log/%s.log\n"%(path,f_out))
f_sub.write("+JobFlavour         = \"%s\"\n"%opt.queue)
f_sub.write("queue %s\n"%numberOfFiles)
f_sub.close()

#Submit condor file
os.system('condor_submit %s'%f_sub_name)

#Delete submission file
os.system('rm %s'%f_sub_name)
