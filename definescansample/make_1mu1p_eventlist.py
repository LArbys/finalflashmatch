import os,sys
import ROOT as rt
import pandas as pd

bnbinter_datapath = "~/working/data/larbys/final_files/test6/mcc8v6_bnb5e19_test6_inter.root"
bnbfinal_datapath = "~/working/data/larbys/final_files/test6/mcc8v6_bnb5e19_test6_dllee.root"

tnumu = rt.TChain("numu_ana_tree")
tnumu.Add(bnbfinal_datapath)

nentries = tnumu.GetEntries()

data = {"run":[],"subrun":[],"event":[],"cosmicll":[],"vtxid":[],"goodreco":[]}

for ientry in range(0,nentries):
    tnumu.GetEntry(ientry)

    if ientry%1000==0:
        print "reading entry ",ientry
    
    if tnumu.num_croi>0 and tnumu.num_vertex>0:
        data["run"].append(tnumu.run)
        data["subrun"].append(tnumu.subrun)
        data["event"].append(tnumu.event)
        data["cosmicll"].append(tnumu.CosmicLL)
        data["vtxid"].append(tnumu.vertex_id)
        data["goodreco"].append(int(tnumu.Good3DReco))


        
df = pd.DataFrame.from_dict( data )

df = df.sort_values(by=["cosmicll"],ascending=False)

# dump top 1000 into a list

topnumu = open("topbnb_1mu1p_test6_wgoodrecocut.txt",'w')

irow = 0
for index,row in df.iterrows():
    #print row["run"],row["subrun"],row["cosmicll"]
    if row["goodreco"]==0:
        continue
    print >> topnumu,'\t',int(row["run"]),'\t',int(row["subrun"]),'\t',int(row["event"]),'\t',int(row["vtxid"]),'\t',row["cosmicll"],'\t',int(row["goodreco"])
    irow+=1
    if irow>=1000:
        break

topnumu.close()
