import os
import sys
import ROOT
import PlotUtils
import math
import copy
from array import array
from collections import OrderedDict
from tools.PlotLibrary import PLOT_SETTINGS,VariantPlotsNamingScheme, HistHolder
from config.AnalysisConfig import AnalysisConfig
from config import BackgroundFitConfig
from tools import Utilities,PlotTools

# Get This from Rob. Thanks Rob.
# This helps python and ROOT not fight over deleting something, by stopping ROOT from trying to own the histogram. Thanks, Phil!
# Specifically, w/o this, this script seg faults in the case where I try to instantiate FluxReweighterWithWiggleFit w/ nuE constraint set to False for more than one playlist
ROOT.TH1.AddDirectory(False)

fitPath = '/minerva/data/users/shenry/antinu_e'
nue_file = 'kin_dist_mcFHCScaledSystematics_nx_collab1_fspline.root'.format(fitPath)
nue_datafile = 'kin_dist_dataFHCScaledSystematics_nx_collab1_fspline.root'.format(fitPath)
nue_mcfile = 'kin_dist_mcFHCScaledSystematics_nx_collab1_fspline.root'.format(fitPath)
#nue_file = 'kin_dist_mcFHCFinalScaled_nx_new_fspline.root'.format(fitPath)
#nue_datafile =  'kin_dist_dataFHCFinalScaled_nx_new_fspline.root'.format(fitPath)
#nue_mcfile = 'kin_dist_mcFHCFinalScaled_nx_new_fspline.root'.format(fitPath)

nue_background = ROOT.TFile.Open("{}/{}".format(fitPath,nue_file))
data_background = ROOT.TFile.Open("{}/{}".format(fitPath,nue_datafile))

nuedata_map_path = "{}/{}".format(fitPath,nue_datafile)
nuemc_map_path = "{}/{}".format(fitPath,nue_mcfile)

HISTOGRAMS_TO_NUE_ADD = [
  "Visible Energy vs Q3",
  "Visible Energy vs Lepton Pt",
  "Lepton Energy",
  "Lepton Theta",
  "Q3",
  "Visible Energy",
  "Lepton Pt",
]


def HistAddition(mc_hists): #used incase you want to add to existing hist

    out_h1 = mc_hists.hists["Other"].Clone() 

    out_h2 = mc_hists.hists["NuEElastic"].Clone() #nue events
    out_h3 = mc_hists.hists["NonPhaseSpace"].Clone()
    out_h4 = mc_hists.hists["NonFiducial"].Clone()
    out_h5 = mc_hists.hists["CCPi0"].Clone()
    out_h6 = mc_hists.hists["NCCohPi0"].Clone()
    out_h7 = mc_hists.hists["NCPi0"].Clone()
 
    out_h1.Add(out_h2) #add the nue events to the Other background
    out_h1.Add(out_h3)
    out_h1.Add(out_h4)
    out_h1.Add(out_h5)
    out_h1.Add(out_h6)
    out_h1.Add(out_h7)
    return out_h1

def DataSubtraction(data_hist,nuebar_hist,background_hist,pot_scale):
    out_data = data_hist.Clone()
    out_nuebar = nuebar_hist.Clone()
    out_background = background_hist.Clone()  
    
    out_nuebar.AddMissingErrorBandsAndFillWithCV(nuebar_hist)
    out_background.AddMissingErrorBandsAndFillWithCV(background_hist)
    print(out_nuebar.GetVertErrorBand("Flux"),out_background.GetVertErrorBand("Flux")) 
    out_nuebar.Scale(pot_scale) #FHC data pot normalize
    out_background.Scale(pot_scale) #FHC data pot normalize
    out_data.AddMissingErrorBandsAndFillWithCV(background_hist) 
    out_data.Add(out_nuebar,-1)
    out_data.Add(out_background,-1)
    return out_data

if __name__ == "__main__":
    
    #input knobs
    playlist=AnalysisConfig.playlist
    nue_type_path_map = {} 
    type_path_map = { t:AnalysisConfig.SelectionHistoPath(playlist,t =="data",False) for t in AnalysisConfig.data_types}
    
    nue_type_path_map['data'] = '{}'.format(nuedata_map_path)
    nue_type_path_map['mc'] = '{}'.format(nuemc_map_path)
 
    datafile,mcfile,pot_scale = Utilities.getFilesAndPOTScale(playlist,type_path_map,AnalysisConfig.ntuple_tag)
    #new_playlist_tag = AnalysisConfig.helicity_tag
   
    pots = [None,None]
    nuepots = [None,None]
    for i,t in enumerate(["data","mc"]):
        try:
            path = type_path_map[t]
            nuepath = nue_type_path_map[t]
        except KeyError:
            continue
        pots[i]= Utilities.getPOTFromFile(path)
        nuepots[i] = Utilities.getPOTFromFile(nuepath) 
    nue_pot_scale  = pots[0] / nuepots[1] #RHC_data_pot / FHC_MC_pot
    data_pot_scale = nuepots[0] / nuepots[1] #FHC_data_pot / FHC_MC_pot needed for subtraction
    datanue_pot_scale = pots[0] / nuepots[0] #RHC_data_pot / FHC_data_pot

    regions = ["Signal"] + AnalysisConfig.sidebands
    newplaylist=ROOT.TFile.Open('/minerva/data/users/shenry/antinu_e/kin_dist_mcRHCScaledSystematics_nx_collab1_fspline.root',"UPDATE")
    
    for region in regions:
        for config in HISTOGRAMS_TO_NUE_ADD: 
            nue_hists = HistHolder(config["name"] if "name" in config else config,nue_background,region,True,nue_pot_scale) #this is to grab a different category than what is in the general signal def 
            data_hists = HistHolder(config["name"] if "name" in config else config,data_background,region,True,1.0)
            data = data_hists.hists["Total"]
           
            nue = nue_hists.hists["Nue"] #grab hist from scaled FHC file, make sure binning is the same as in the RHC tuple
            print(nue.GetVertErrorBand("Flux"))
            nuebar = nue_hists.hists["NueBar"] 
            background = nue_hists.hists["Background"] 
            
            nue_prediction = DataSubtraction(data,nuebar,background,data_pot_scale)
            print(nue_prediction.GetVertErrorBand("Flux"))       
            nue_prediction.Scale(1/nue_pot_scale) 
            nue_prediction.Write("{}Other".format(nue.GetName()))
           
            #nue.Scale(1/nue_pot_scale) #scale nue hists, use this for FHC MC prediction   
            #nue.Write("{}OtherNue".format(nue.GetName())) #write as a new hist in the current tuple being used

    newplaylist.Close()
