import ROOT
import os, time, sys, math
import datetime as dt
import PlotUtils
import UnfoldUtils
from multiprocessing import Process

from tools.PlotLibrary import PLOT_SETTINGS,HistHolder
from config.AnalysisConfig import AnalysisConfig
from config.BackgroundFitConfig import CATEGORY_FACTORS
from config.SignalDef import SIGNAL_DEFINATION
from config import PlotConfig
from tools import Utilities,PlotTools
from tools.unfoldingwrapper import DeOverflowWrapper,WithBackgroundWrapper,IdentityWrapper

ROOT.TH1.AddDirectory(False)
baseDir = os.path.dirname(os.path.abspath(__file__))+"/../"

Iterations = "1,2,3,4,5,6,7,8,9,10,15"
#Iterations= "1"
NUniverses =5
MaxChi2 =10000
stepChi2 = 50
remove=""#"163,164,165,166,167,168,169,170,171,172,173,174,175,176,177,178,179,180"
NUniverses_per_job = 100
K=0

def GetPaths(baseDir):
  l = baseDir.split("/")
  index = l.index("cmtuser")
  Release = "/".join(l[0:index+2])
  RelativePath = "/".join(l[index+2:-1])
  return Release,RelativePath

def createTarball( releaseDir,outDir):
  print("I'm inside createTarball()")
  found = os.path.isfile(outDir)
  if found and input("Overwrite tarball? (Y/n)").upper() =="Y":
    os.remove(outDir)
    found = False

  if not found:
    cmd = "tar -czf %s -C %s %s"%(outDir, releaseDir, ".")# remove /Ana/CCNuE from baseDir bc want to tar the release.
    print(cmd)
    os.system(cmd)

  print("I'm done creating the tarballs")

def unpackTarball( mywrapper ,RelativePath):
  # Add lines to wrapper that wil unpack tarball; add additional setup steps here if necessary  
  mywrapper.write("cd $CONDOR_DIR_INPUT\n")
  mywrapper.write("source /cvmfs/minerva.opensciencegrid.org/minerva/setup/setup_minerva_products.sh\n")
  mywrapper.write("tar -xvzf {}\n".format(outdir_tarball.split("/")[-1]))
  # Set up test release
  #mywrapper.write("pushd Tools/ProductionScriptsLite/cmt/\n")#assuming only one test release
  #mywrapper.write("cmt config\n")
  #mywrapper.write(". setup.sh\n")
  #mywrapper.write("popd\n")
  mywrapper.write("pushd Ana/PlotUtils/cmt/\n")
  mywrapper.write("cmt config\n")
  #mywrapper.write("cmt make\n")
  mywrapper.write(". setup.sh\n")
  mywrapper.write("cd ../macros\n")
  mywrapper.write("export LD_LIBRARY_PATH=/cvmfs/minerva.opensciencegrid.org/minerva/software_releases/v22r1p1/lcg/external/GSL/1.10/x86_64-slc7-gcc49-opt/lib:$LD_LIBRARY_PATH\n")
  #mywrapper.write("pushd {}\n".format(RelativePath))
  #mywrapper.write(". setup.sh\n")
  #mywrapper.write("make\n")
  
  #mywrapper.write("popd\n")

def submitJob( tupleName,relative_path,tag):

  # Create wrapper
  wrapper_name = "grid_wrappers/%s/%s_wrapper.sh" % ( processingID , tupleName ) 
  
  my_wrapper = open(wrapper_name,"w")
  my_wrapper.write("#!/bin/sh\n")
  unpackTarball(my_wrapper,relative_path)
    # Write dummy output files to get around dumb copy-back bug
  #makeDummyFile(my_wrapper,'$CONDOR_DIR_LOGS')
  #makeDummyFile(my_wrapper,'$CONDOR_DIR_HISTS')


  # This is the bash line that will be executed on the grid
  my_wrapper.write( "export USER=$(whoami)\n")
  #my_wrapper.write( "export XRD_LOGLEVEL=\"Debug\"\n")

  my_wrapper.write("readlink -f `which TransWarpExtraction` &> $CONDOR_DIR_OUTPUT/log_{tag}_$PROCESS\n".format(tag=tag)) 
  pot_scale=1
  my_wrapper.write("TransWarpExtraction -o $CONDOR_DIR_OUTPUT/transWrap_{tag}_$PROCESS.root -d data -D $CONDOR_DIR_INPUT/{iput} -i data_truth -I $CONDOR_DIR_INPUT/{iput} -m migration -M $CONDOR_DIR_INPUT/{iput} -r reco -R $CONDOR_DIR_INPUT/{iput} -t reco_truth -T $CONDOR_DIR_INPUT/{iput} -n {Iter} -z {dimension} -u {NUniverse} -c {max_chi2} -C {step_chi2} -b -L -f 10.0 &>> $CONDOR_DIR_OUTPUT/log_{tag}_$PROCESS\n".format(tag=tag,iput=ifile,Iter=Iterations,dimension=dimension,NUniverse=NUniverses_per_job,max_chi2=MaxChi2,step_chi2=stepChi2))

  my_wrapper.write("exit $?\n")
  #my_wrapper.write( "python eventSelection.py -p %s --grid --%s-only --ntuple_tag %s --count %d %d  --output $CONDOR_DIR_HISTS %s \n" % (playlist, dataSwitch, gridargs.ntuple_tag, start, count, argstring) )
 
  #addBashLine(my_wrapper,'date')
  my_wrapper.close()
  
  os.system( "chmod 777 %s" % wrapper_name )
  cmd = "jobsub_submit --group=minerva --resource-provides=usage_model=DEDICATED,OPPORTUNISTIC --role=Analysis --memory %dMB -f %s -d OUTPUT %s  -r %s -N %d --expected-lifetime=%dh --cmtconfig=%s -f dropbox://./%s -i /cvmfs/minerva.opensciencegrid.org/minerva/software_releases/%s/ file://%s/%s" % ( memory , outdir_tarball , outdir_hists , os.environ["MINERVA_RELEASE"], njobs, 12, os.environ["CMTCONFIG"], ifile, os.environ["MINERVA_RELEASE"], os.environ["PWD"] , wrapper_name )
  print(cmd)
  os.system(cmd)

def DoFancyTuneToMigration(migration,reco,reco_truth,data_truth,k):
  try:
    scale = reco_truth.Integral(0,-1)/data_truth.Integral(0,-1)
  except:
    scale = reco_truth.Integral(0,-1,0,-1)/data_truth.Integral(0,-1,0,-1)
  for i in range(migration.GetNbinsX()+2):
    for j in range(migration.GetNbinsY()+2):
      reco.AddBinContent(i,-1.0*migration.GetBinContent(i,j))

  for j in range(migration.GetNbinsY()+2):
    target = data_truth.GetBinContent(j)*scale * k+reco_truth.GetBinContent(j)*(1-k)
    correction = target/reco_truth.GetBinContent(j) if reco_truth.GetBinContent(j)!=0 else 1
    reco_truth.SetBinContent(j,correction*reco_truth.GetBinContent(j))
    for i in range(migration.GetNbinsX()+2):
      migration.SetBinContent(i,j,migration.GetBinContent(i,j)*correction)
      migration.SetBinError(i,j,math.sqrt(migration.GetBinContent(i,j)))
      reco.AddBinContent(i,migration.GetBinContent(i,j))


def ConsolidateInputs(Temp_path,mig,reco,reco_truth,data,data_truth):
  f = ROOT.TFile(Temp_path,"RECREATE")
  DoFancyTuneToMigration(mig,reco,reco_truth,data_truth,K)
  for h in [data,data_truth]:
        h = h.Clone("{}_dataclone".format(h.GetName()))
        h.Scale(pot_scale)
        for i in range(h.GetSize()):
          h.SetBinError(i,math.sqrt(h.GetBinContent(i)))
  mig.Write("migration")
  reco_truth.Write("reco_truth")
  reco.Write("reco")
  for h in [data,data_truth]:
    h.Scale(pot_scale)
    for i in range(h.GetSize()):
      h.SetBinError(i,math.sqrt(h.GetBinContent(i)))
  data_truth.Write("data_truth")
  data.Write("data")
  f.Close()


def CalcChi2(h1,h2,cov):
  l1 = h1.GetSize()
  l2 = h2.GetSize()
  
  # if not (l1 == l2 == cov.GetNcols() == cov.GetNrows()):
  #     print (l1,l2,cov.GetNcols(),cov.GetNrows())
  #     raise ValueError("cov matrix/histograms dont have same dimention to calculate correlation.")
  diff = [(i,h1.GetBinContent(i)-h2.GetBinContent(i)) for i in range(l1) if (h1.GetBinContent(i)-h2.GetBinContent(i))!=0]
  v1 = ROOT.TMatrixD(1,len(diff))
  v2 = ROOT.TMatrixD(len(diff),1)
  this_cov = ROOT.TMatrixD(len(diff),len(diff))
  for i,entry in enumerate(diff):
    v1[0][i]=entry[1]
    v2[i][0]=entry[1]
    for j,entry2 in enumerate(diff):
      this_cov[i][j] = cov[entry[0]][entry2[0]]

  SVD = ROOT.TDecompSVD(this_cov)
  #print [SVD.GetV()[i][i] for i in range(len(diff))]
  print(diff[0:10])
  tmp = SVD.Invert()
  print([tmp[i][i] for i in range(len(diff))])
  tmp = tmp*v2
  tmp = v1*tmp
  return tmp[0][0]

def SubtractPoissonHistograms(h,h1):
    h.AddMissingErrorBandsAndFillWithCV (h1)
    errors = []
    for i in range(h.GetSize()):
        errors.append(math.sqrt(h.GetBinError(i)**2 + h1.GetBinError(i)**2))
    h.Add(h1,-1)
    for i in range(h.GetSize()):
        h.SetBinError(i,errors[i])

    return h

def MyUnfoldWrapper(ifile_name):
  ifile = ROOT.TFile.Open(ifile_name)
  response=ifile.Get("migration")
  reco=ifile.Get("reco")
  truth=ifile.Get("reco_truth")
  hlin_measured = ifile.Get("data")
  unfold_old = ifile.Get("data_truth")
  miggg = ROOT.RooUnfoldResponse(reco,truth,response)
  miggg.UseOverflow(1)
  unfolder = ROOT.RooUnfold.New(ROOT.RooUnfold.kBayes,miggg,hlin_measured,2)
  unfolder.IncludeSystematics(0)
  h_unfolded = unfolder.Hreco()
  h_unfolded_cov = unfolder.Ereco()
  #print (h_unfolded_cov[21][21],h_unfolded_cov[21][22])
  #print unfold_old.GetBinContent(181)
  print(CalcChi2(unfold_old,h_unfolded,h_unfolded_cov))
  ifile.Close()

def MnvUnfoldWrapper(ifile_name):
  ifile = ROOT.TFile.Open(ifile_name)
  response=ifile.Get("migration")
  reco=ifile.Get("reco")
  truth=ifile.Get("reco_truth")
  hlin_measured = ifile.Get("data")
  unfold_old = ifile.Get("data_truth")
  mnvunfold=UnfoldUtils.MnvUnfold()
  h_unfolded =truth.Clone()
  covmx = ROOT.TMatrixD(500,500)
  mnvunfold.setUseBetterStatErrorCalc(True)
  mnvunfold.UnfoldHistoWithFakes(h_unfolded,covmx,response,response.ProjectionX("fuck",0,-1),ROOT.nullptr
,ROOT.nullptr
,ROOT.nullptr,1,False,False)
  for i in range(18):
    print ([covmx[i][j] for j in range(18)])
  print ("done")
 
  # for i in range(h_unfolded.GetSize()):
  #   covmx[i][i]=0
  # covmx.ResizeTo(h_unfolded.GetSize(),h_unfolded.GetSize())
  # h_unfolded.PushCovMatrix("unfold",covmx)
  # matrix = h_unfolded.GetTotalErrorMatrix()
  # #print (h_unfolded.GetBinError(20)**2,matrix[20][20])
  # print CalcChi2(unfold_old,h_unfolded,matrix)
  # plotter=PlotUtils.MnvPlotter()
  # plotter.chi2_use_overflow_err=True
  # _ =ROOT.Long()
  # print plotter.Chi2DataMC(h_unfolded,unfold_old,_,1,True,False,False)
  # print plotter.Chi2DataMC(h_unfolded,unfold_old,1,True,False,False)
  ifile.Close()

def RunTransWrapper(ifile,ofile):
  cmd = "TransWarpExtraction -o {output} -d data -D {iput} -i data_truth -I {iput} -m migration -M {iput} -r reco -R {iput} -t reco_truth -T {iput} -n {Iter} -z {dimension} -u {NUniverse} -c {max_chi2} -C {step_chi2} {remove} -b -L ".format(output=ofile,iput=ifile,Iter=Iterations,dimension=dimension,NUniverse=NUniverses,max_chi2=MaxChi2,step_chi2=stepChi2,remove="-x {}".format(remove) if remove else "")
  print(cmd)
  os.system(cmd)


AlternativeModels = {
  #"CV":None,
   #"LowQ2Pi0": ("LowQ2Pi",0),
  # "LowQ2Pi1": ("LowQ2Pi",1),
  # "LowQ2Pi2": ("LowQ2Pi",2),
  # "LowQ2Pi3": ("LowQ2Pi",3),
  # "2p2h0": ("Low_Recoil_2p2h_Tune",0),
  # "2p2h1": ("Low_Recoil_2p2h_Tune",1),
  # "2p2h2": ("Low_Recoil_2p2h_Tune",2),
  # "RPA_highq20" :("RPA_HighQ2",0),
  # "RPA_highq21" :("RPA_HighQ2",1),
  #"RPA_lowq20" :("RPA_LowQ2",0),
  # "RPA_lowq21" :("RPA_LowQ2",1),
  # "MK_Model":("MK_model",0),
  "FSI_Weight0":("fsi_weight",0),
  # "FSI_Weight1":("fsi_weight",1),
  # "FSI_Weight2":("fsi_weight",2),
  # "SuSA2p2h":("SuSA_Valencia_Weight",0),
  # "GenieMaCCQE_UP":("GENIE_MaCCQE",1),
  # "GenieMaCCQE_DOWN":("GENIE_MaCCQE",0),
}

Measurables = [
  "Eavail_q3",
  "Eavail",
  "Eavail_Lepton_Pt"
]

if __name__ == '__main__':
  playlist=AnalysisConfig.playlist
  type_path_map = { t:AnalysisConfig.SelectionHistoPath(playlist,t =="data",False) for t in AnalysisConfig.data_types}
  data_file,mc_file,pot_scale = Utilities.getFilesAndPOTScale(playlist,type_path_map,AnalysisConfig.ntuple_tag)
  #unfolded_file = ROOT.TFile.Open(AnalysisConfig.UnfoldedHistoPath(playlist,False))
  ifile="fin.root"
  release_path,relative_path=GetPaths(baseDir)
  PNFS_switch = "persistent"
  # Automatically generate unique output directory
  #processingID = '%s_%s-%s' % ("CCNUE_selection", dt.date.today() , dt.datetime.today().strftime("%H") )
  processingID = "CCNUE_Warping5"
  
  outdir_hists = "/pnfs/minerva/%s/users/%s/%s_hists" % (PNFS_switch,os.environ["USER"],processingID)
  os.system( "mkdir %s" % outdir_hists )
  outdir_logs = "/pnfs/minerva/%s/users/%s/%s_logs" % (PNFS_switch,os.environ["USER"],processingID)
  os.system( "mkdir %s" % outdir_logs )
  os.system( "mkdir -p grid_wrappers/%s" % processingID )
  outdir_tarball="/pnfs/minerva/resilient/tarballs/hsu-%s.tar.gz" % (processingID)
  createTarball(release_path,outdir_tarball)
  memory =2000
  for prefix in Measurables:
    MnvMigration = Utilities.GetHistogram(mc_file, prefix+ "_migration")
    MnvMCreco = Utilities.GetHistogram(mc_file, prefix)
    MnvMCtruth = Utilities.GetHistogram(mc_file, prefix+ "_migration_truth")
    if not MnvMCtruth:
      MnvMCtruth = MnvMigration.ProjectionY("_truth",0,-1)
    MnvMCSignal = Utilities.GetHistogram(mc_file, prefix+ "_migration_reco")
    if not MnvMCSignal:
      MnvMCSignal = MnvMigration.ProjectionX("_reco",0,-1)
    dimension = MnvMCreco.GetDimension()
    print(dimension)
    if dimension ==1:
      H2D = PlotUtils.MnvH1D
    else:
      H2D = PlotUtils.MnvH2D


    for model,value in AlternativeModels.items():
      data_migration = PlotUtils.MnvH2D(MnvMigration.GetCVHistoWithStatError())
      data_truth = H2D(MnvMCtruth.GetCVHistoWithStatError())
      data_reco = H2D(MnvMCreco.GetCVHistoWithStatError())
      data_signal = H2D(MnvMCSignal.GetCVHistoWithStatError())
      if model=="CV":
        model_migration = PlotUtils.MnvH2D(MnvMigration.GetCVHistoWithStatError())
        model_truth=H2D(MnvMCtruth.GetCVHistoWithStatError())
        model_signal=H2D(MnvMCSignal.GetCVHistoWithStatError())
        model_reco=H2D(MnvMCreco.GetCVHistoWithStatError())
      else:
        try:
          model_migration = PlotUtils.MnvH2D(MnvMigration.GetVertErrorBand(value[0]).GetHist(value[1]))
          model_truth = H2D(MnvMCtruth.GetVertErrorBand(value[0]).GetHist(value[1]))
          model_signal= H2D(MnvMCSignal.GetVertErrorBand(value[0]).GetHist(value[1]))
          model_reco= H2D(MnvMCreco.GetVertErrorBand(value[0]).GetHist(value[1]))
        except TypeError as e:
          print(e)
          continue

      #ConsolidateInputs(ifile,model_migration,model_signal,model_truth,data_signal,data_truth)
      #ConsolidateInputs(ifile,model_migration,data_reco,data_truth,model_reco,model_truth)
      ConsolidateInputs(ifile,model_migration,model_reco,model_truth,data_reco,data_truth)
      #ConsolidateInputs(ifile,model_migration,model_signal,model_truth,data_signal,model_truth)
      #ConsolidateInputs(ifile,model_migration,model_reco,model_truth,data_reco,model_truth)
      #ConsolidateInputs(ifile,model_migration,model_signal,model_truth,model_signal,data_truth)
      
      njobs = int(math.ceil(1.0*NUniverses/NUniverses_per_job))
      cmdString = "CCNuE-%s-%s" % (playlist,"mc")
      #wrapper.setup(migration,cv_reco,cv_bkg)
      #MyUnfoldWrapper(ifile)
      #MnvUnfoldWrapper(ifile)
      #submitJob(cmdString,relative_path,model+prefix)
      RunTransWrapper(ifile,"transWrapT_{}.root".format(model+prefix))
      #break
      # mnvunfold=UnfoldUtils.MnvUnfold()
      # h_test = MnvMCtruth.Clone()
      # cov=ROOT.TMatrixD(1,1)
      # mnvunfold.UnfoldHistoWithFakes(h_test,cov,MnvMigration,data,cv_signal,cv_truth,cv_bkg,5,False,False)

      #print CalcChi2(h_test,data_truth,cov)
