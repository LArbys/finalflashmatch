import os,sys,commands

cfg = "inter_tool_test.cfg"

data_file = "/home/twongjirad/working/data/larbys/mcc8v6/testdb/bnb5e19/run5925subrun0195"

inter   = "/home/twongjirad/working/data/larbys/final_files/test6/mcc8v6_bnb5e19_test6_inter.root"

ssnet   = data_file+"/ssnetout-larcv-Run005925-SubRun000195.root"
vertex  = data_file+"/vertexout_filter_numu_ana_tree_59250195.root"
shower  = data_file+"/shower_reco_out_59250195.root"
tracker = data_file+"/tracker_reco_59250195.root"

opreco  = data_file+"/opreco-Run005925-SubRun000195.root"
reco2d  = data_file+"/reco2d-Run005925-SubRun000195.root"

args = "%s %s %s %s %s %s %s"%(cfg,inter,ssnet,vertex,shower,tracker,opreco)

os.system("./select_tune_sample %s"%(args))




