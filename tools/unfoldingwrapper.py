import ROOT
import math
import PlotUtils

class MyResponseMaker(object):
    def setup(self,response, mc_signal = None, mc_background=None,mc_truth = None):
        self.response=response
        self.mc_signal = mc_signal
        self.mc_truth = mc_truth
        self.mc_background = mc_background

    @staticmethod
    def ResponseCoverter(hist,nm,fakes=None):
        if hist.GetNbinsX()+2 == nm or hist.GetNbinsX() == nm:
            name = "new_{}".format(hist.GetName())
            new_hist =PlotUtils.MnvH2D(name,name,nm,0,nm,nm+(fakes is not None),0,nm+(fakes is not None))

            for i in range(nm):
                for j in range(nm):
                    new_hist.SetBinContent(i+1,j+1,hist.GetBinContent(i,j))
                    new_hist.SetBinError(i+1,j+1,hist.GetBinError(i,j))
            if fakes:
                truth_tmp = hist.ProjectionX()
                for i in range(nm):
                    new_hist.SetBinContent(i+1,nm+1,fakes.GetBinContent(i)-truth_tmp.GetBinContent(i))
                    new_hist.SetBinError(i+1,nm+1,math.sqrt(fakes.GetBinError(i)**2+truth_tmp.GetBinError(i)**2))
            return new_hist
        else:
            raise ValueError("index mismatch, res x bins :{}, mes x bins: {}".format(hist.GetNbinsX(),nm))

    @staticmethod
    def Linearlizer(hist):
        nd = hist.GetSize()
        name = "new_{}".format(hist.GetName())
        new_hist =PlotUtils.MnvH1D(name,name,nd,0,nd)

        for i in range(nd):
            new_hist.SetBinContent(i+1,hist.GetBinContent(i))
            new_hist.SetBinError(i+1,hist.GetBinError(i))
        return new_hist,nd

    @staticmethod
    def Delinearlizer(hist_sample,hist_input):
        return_hist=hist_sample.Clone("new_{}".format(hist_sample.GetName()))
        nd = hist_sample.GetSize()
        for i in range(nd):
            return_hist.SetBinContent(i,hist_input.GetBinContent(i+1))
            return_hist.SetBinError(i,hist_input.GetBinError(i+1))
        return return_hist

    @staticmethod
    def SubtractPoissonHistograms(h,h1):
        try:
            h.AddMissingErrorBandsAndFillWithCV (h1)
        except:
            print("Not MnvHXD, dont mind")
        errors = []
        for i in range(h.GetSize()):
            errors.append(math.sqrt(h.GetBinError(i)**2 + h1.GetBinError(i)**2))
        h.Add(h1,-1)
        for i in range(h.GetSize()):
            h.SetBinError(i,errors[i])

        return h

class IdentityWrapper(MyResponseMaker):
    def GetMigration(self):
        return self.response,self.mc_signal,self.mc_truth

    def GetToUnfold(self,data,pot_scale,background = None):
        if background is None:
            background = self.mc_background.Clone()
        background.Scale(pot_scale)
        MyResponseMaker.SubtractPoissonHistograms(data,background)
        return data

    def PostProcessing(self,data):
        return data

    def PrepareTruth(self,data_truth):
        return data_truth


class DeOverflowWrapper(MyResponseMaker):
    def GetMigration(self):
        nm = self.mc_background.GetSize()
        new_migration = MyResponseMaker.ResponseCoverter(self.response,nm)
        return new_migration,new_migration.ProjectionX("reco",0,-1),new_migration.ProjectionY("truth",0,-1)

    def GetToUnfold(self,data,pot_scale,background = None):
        if background is None:
            background = self.mc_background.Clone()
        background.Scale(pot_scale)
        MyResponseMaker.SubtractPoissonHistograms(data,background)
        return MyResponseMaker.Linearlizer(data)[0]

    def PostProcessing(self,data):
        return MyResponseMaker.DeLinearlizer(self.truth,data)

    def PrepareTruth(self,data_truth):
        return MyResponseMaker.Linearlizer(data_truth)[0]


class WithBackgroundWrapper(MyResponseMaker):
    def GetMigration(self):
        mc_reco = self.mc_signal.Clone()
        mc_reco.Add(self.mc_background)
        new_reco,nm = MyResponseMaker.Linearlizer(mc_reco)
        new_migration = MyResponseMaker.ResponseCoverter(self.response,nm)
        return new_migration,new_reco,new_migration.ProjectionY()

    def GetToUnfold(self,data,pot_scale,background=None):
        return MyResponseMaker.Linearlizer(data)[0]

    def PostProcessing(self,data):
        return MyResponseMaker.DeLinearlizer(self.truth,data)

    def PrepareTruth(self,data_truth):
        tmp = MyResponseMaker.Linearlizer(data_truth)[0]
        return tmp






