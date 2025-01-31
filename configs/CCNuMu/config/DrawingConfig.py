import ROOT
from collections import OrderedDict
from config import PlotConfig
from tools.PlotLibrary import HistHolder
from tools import PlotTools
from config.SystematicsConfig import CONSOLIDATED_ERROR_GROUPS,DETAILED_ERROR_GROUPS


Default_Plot_Type="stacked"
Default_Scale = lambda histHolder:histHolder.POTScale(True)
DefaultSlicer = PlotTools.PrepareSlicer


# the order and color of each band in stacked plots defined here.
COLORS=ROOT.MnvColors.GetColors()
Categories = OrderedDict()
Categories["CCQE"]={
    "title" : "QE",
    "color": COLORS[0]
}

Categories["CCDelta"]={
    "title" : "Delta",
    "color": COLORS[1]
}

Categories["CC2p2h"]={
    "title": "2p2h",
    "color": COLORS[3]
}
Categories["CCDIS"]= {
    "title": "DIS",
    "color" : COLORS[2]
}
Categories["CCOther"]= {
    "title": "CC-Other",
    "color" : COLORS[4]
}

# Categories["NuEElastic"] = {
#     "title":"#nu + e elastic",
#     "color" : COLORS[6]
# }
Categories["Other"] = {
    "title":"Others",
    "cate" : {"Other","NonFiducial","NonPhaseSpace","CCNuE","NC"},
    "color" : COLORS[6]
}


SignalDecomposition = {
    "CCQE" :
    {
        "title" : "CC #nu_{mu}-QE",
        "color": COLORS[0]
    },
    "CCDelta" : {
        "title" : "CC #nu_{mu}-Delta",
        "color": COLORS[1]
    },
    "CCDIS" : {
        "title" : "CC #nu_{mu}-DIS",
        "color": COLORS[2]
    },
    "CC2p2h" : {
        "title" : "CC #nu_{mu}-2p2h",
        "color": COLORS[3]
    },
    "Background" : {
        "title" : "Backgrounds",
        "cate": {"Other","NonFiducial","NonPhaseSpace","CCNuE","NC"},
        "color": COLORS[4]
    }
    # "CCNuEAntiNu" : {
    #     "title" : "CC #bar{#nu_{e}}",
    #     "color": COLORS[4]
    # }
}


DefaultPlotters={
    "comp":{"func": PlotTools.PrepareComp},
    "ratio":{"func": PlotTools.PrepareRatio},
    "err": {"func": PlotTools.PrepareErr,
            "args":(True,False,CONSOLIDATED_ERROR_GROUPS)},
    "stacked":{"func":PlotTools.PrepareStack,
               "args": (Categories,)},
    "diff":{"func":PlotTools.PrepareDiff},
    "migration":{"func":PlotTools.PrepareMigration},
    "sigdep":{"func":PlotTools.PrepareSignalDecompose,
              "args": (SignalDecomposition,True,False)},
    "sigdepRatio":{"func":PlotTools.PrepareSignalDecomposeRatio,
                   "args": (SignalDecomposition,True,False)},

}

PLOTS_TO_MAKE = [
    # {"name":"Lepton Energy",
    #  "plot_type":"stacked"},
    # {"name":"Neutrino Energy",
    #  "plot_type":"stacked"},
    # {"name":"Neutrino Energy",
    #   "sideband_group": ("ex_and_pi0",("Excess","Pi0"))},
    # {"name":"Visible Energy",
    #   "sideband_group": ("ex_and_pi0",("Excess","Pi0"))},
    # {"name": "Lepton Theta"},
    # {"name": "Lepton Theta"},
    # {"name": "Psi"},
    # {"name":"Q0"},
    # {"name":"Q3"},
    # {"name":"Visible Energy",
    #   "plot_type" : "comp",},
    # {"name":"Q3",
    #  "plot_type" : "sigdep",},
    # {"name":"Visible Energy",
    #  "plot_type" : "err",},
    # {"name":"Q0",
    #  "plot_type" : "err",},
    # {"name":"Q3",
    #  "plot_type" : "err",},
    # {"name":"Q3 Migration",
    #   "plot_type":"migration",
    #  "slicer": lambda x: [x.Clone()]},
    # {"name":"Q0 Migration",
    #  "plot_type":"migration",
    #  "slicer": lambda x: [x.Clone()]},
     # {"name":"Visible Energy Migration",
     #  "plot_type":"migration",
     #  "slicer": lambda x: [x.Clone()]},
    # {"name":"Visible Energy vs q3"},
    # {"name":"Visible Energy vs q3",
    #  "plot_type":"comp"},
    #  {"name":"Visible Energy vs q3",
    #  "plot_type":"err"},
    #  {"name":"Visible Energy vs q3",
    #  "plot_type":"sigdepRatio"},
    # {"name":"Visible Energy vs q3",
    #  "plot_type":"sigdep"},

    # {"variables":["Q0","True Visible Energy"],
    #  "scale":lambda hist:hist.AreaScale(True),
    #  "slicer":lambda hist: PlotTools.Make2DSlice(hist,True),
    #  "plot_type":"comp",
    #  "args":(True,False,True),
    #  },
    #  {"name":"Visible Energy Migration",
    #   "slicer":lambda hist: PlotTools.Make2DSlice(hist,True),
    #   "scale":lambda hist:hist.AreaScale(True),
    #   "args":(True,False,True),
    #   "plot_type":"comp",
    #   },
    #  {"name":"Q0 Migration",
    #   "slicer":lambda hist: PlotTools.Make2DSlice(hist,True),
    #   "scale":lambda hist:hist.AreaScale(True),
    #   "args":(True,False,True),
    #   "plot_type":"comp",
    #   },
    # {"variables":["Visible Energy","True Q0"],
    #  "scale":lambda hist:hist.AreaScale(True),
    #  "slicer":lambda hist: PlotTools.Make2DSlice(hist,True),
    #  "args":(True,False,True),
    #  "plot_type":"comp",
    # },
    # {"name":"Visible Energy Migration",
    #  "plot_type":"comp",
    #  "slicer":lambda hist: PlotTools.Make2DSlice(hist,True)},
    #  {"name":"Q0 Migration",
    #  "plot_type":"comp",
    #  "slicer":lambda hist: PlotTools.Make2DSlice(hist,True)},
    # {"name":"True Signal Visible Energy vs q3"
    #  "plot_type":"err"},
    {"name":"Visible Energy vs q3",
     "plot_type" : "comp"},
    {"name":"Visible Energy vs q3",
     "plot_type" : "ratio"},
    {"name":"Visible Energy vs q3",
     "plot_type" : "err"},
    {"name":"Visible Energy vs q3",
     "plot_type" : "sigdep",},
    # {"name":"q0 vs q3",
    #  "plot_type" : "sigdep",},
    # {"name":"Visible Energy"},
    # {"variables":["Visible Zoomin","Lepton Energy"],
    #  "slicer": lambda hist: PlotTools.Make2DSlice(hist,interval=5),
    #  "args":(NuEElasticCategory,),
    #  "scale":lambda histHolder:histHolder.POTScale(False)},
    # {"name":"FrontdEdX_position" },
    # {"name":"Front dEdX" },
    # {"name":"Inline Upstream Energy Weighted Position",
    #  "canvasconfig":PlotTools.Logy},
    # {"name":"Inline Upstream Energy",
    #  "canvasconfig":PlotTools.Logx,
    #  "scale" : lambda histHolder:histHolder.POTScale(False)},
    {"variables":["Visible Energy","Lepton Pt"]},
    {"variables":["Visible Energy","Lepton Pt"],
     "plot_type":"comp"},
    {"variables":["Visible Energy","Lepton Pt"],
     "plot_type":"sigdepRatio"},
     {"variables":["Visible Energy","Lepton Pt"],
     "plot_type":"err"},
    # {"variables":["Visible Energy","Lepton Pt"],
    #  "slicer": PlotTools.IdentitySlicer,
    #  "scale":lambda histHolder:histHolder.POTScale(False),
    #  "plot_type":"diff"},
    # {"variables":["Visible Energy","Lepton Pt"],
    #  "sideband_group": ("Excess",("Excess_High_Inline","Excess_Low_Inline")),
    #  "slicer": PlotTools.IdentitySlicer,
    #  "scale":lambda histHolder:histHolder.POTScale(False),
    #  "plot_type":"diff"},
 
    #  {"variables":["Lepton Energy"],
    #  "sideband_group": ("Excess",("Excess_High_Inline","Excess_Low_Inline"))},
    # {"variables":["Sum Visible Energy","Lepton Pt"],
    #  "sideband_group": ("Excess",("Excess_High_Inline","Excess_Low_Inline")),
    #  "slicer":PlotTools.ProjectionX},
    # {"variables":["Sum Visible Energy","Lepton Pt"],
    #  "slicer":PlotTools.ProjectionX},
    # {"variables":["Visible Energy","Lepton Pt"],
    #  "plot_type":"err"},
    # {"variables":["Lepton Pt"]},
    # {"variables":["Lepton Energy"]},

    # {"variables":["Lepton Pt"],
    #  "plot_type":"err"},
    # {"variables":["Visible Energy","Lepton Energy"]},
    # {"variables":["N_EMLike_prongs"]},
    # {"variables":["Proton Angle"]},
    # {"variables":["Proton Angle Residual"],
    #  "plot_type":"sigdep"},
    # {"variables":["Proton Angle Residual"]},
    # {"name":"Lepton Energy Low Inline",
    #  "plot_type" : "comp"},
    # {"name":"Lepton Energy High Inline",
    #  "plot_type" : "comp"},
    #  {"name":"Lepton Energy Low Inline",
    #  "plot_type" : "err"},
    # {"name":"Lepton Energy High Inline",
    #  "plot_type" : "err"},
    # {"name":"Pseudo q3"},
    # {"name":"Psi"},
    # {"variables":["W"]},
    # {"variables":["Qsquare"]},
    #{"variables":["Vertex Multiplicity"]},
    #{"variables":["Start Multiplicity"]},
    #{"variables":["Lepton Visible Energy"],
    # },
    #{"variables":["Sum Visible Energy"],
    # },
    # {"variables":["Exuv"]},
    # {"variables":["Euv"]},
    # {"name":"True Signal Visible Energy vs q3",
    #  "plot_type":"comp"},
    #{"variables":["Visible Energy","Lepton Visible Energy"],
    # "plot_type":"comp"},
    # {"variables":["Sum Visible Energy"],
    #  "plot_type":"stacked",
    #  },
    # {"name":"Sum Visible Energy low inline"},
    # {"name":"Sum Visible Energy high inline"},
    # {"variables":["Front dEdX", "Visible Energy"],
    #  "sideband_group": ("Merged",("Signal","Excess","Pi0"))},
    #  {"variables":["Front dEdX", "Q3"],
    #  "sideband_group": ("Merged",("Signal","Excess","Pi0"))},
    # {"variables":["Front dEdX", "Lepton Pt"],
    #  "sideband_group": ("Merged",("Signal","Excess","Pi0"))},
    # {"variables":["Front dEdX", "Sum Visible Energy"],
    #  "sideband_group": ("Merged",("Signal","Excess","Pi0"))},

    # {"variables":["Lepton Visible Energy","Visible Energy"]},
    # {"variables":["Sum Visible Energy","Visible Energy"]},
    # #  {"variables":["Front dEdX", "Q3"],
    # # "sideband_group": ("Merged",("Excess","Pi0","Signal"))},
    # {"variables":["Vertex Apothem","Visible Energy"]},
    #  #"sideband_group": ("Merged",("Signal","Excess","Pi0"))},
    # {"variables":["Vertex Y","Visible Energy"]},
    # {"variables":["Vertex Z","Visible Energy"]}
    #  {"variables":["NuEFuzz E","Lepton Pt"],
    #  "plot_type" : "sigdep",
    #   "canvasconfig":PlotTools.Logy},
    # {"variables":["NuEFuzz E","Q3"],
    #  "plot_type" : "sigdep",
    #   "canvasconfig":PlotTools.Logy},
    # {"variables":["NuEFuzz Proton","Lepton Pt"],
    #  "plot_type" : "sigdep",
    #  "canvasconfig":PlotTools.Logy},
    # {"variables":["NuEFuzz Muon","Q3"],
    #  "plot_type" : "sigdep",
    #   "canvasconfig":PlotTools.Logy},
    # {"variables":["NuEFuzz Proton","Lepton Pt"],
    #  "plot_type" : "sigdep",
    #   "canvasconfig":PlotTools.Logy},
    # {"variables":["NuEFuzz Proton","Q3"] ,
    #  "plot_type" : "sigdep",
    #   "canvasconfig":PlotTools.Logy},
    # {"variables":["NuEFuzz Meson","Lepton Pt"],
    #  "plot_type" : "sigdep",
    #   "canvasconfig":PlotTools.Logy},
    # {"variables":["NuEFuzz Meson","Q3"],
    #  "plot_type" : "sigdep",
    #   "canvasconfig":PlotTools.Logy},
    #  {"variables":["NuEFuzz Neutron","Lepton Pt"],
    #  "plot_type" : "sigdep",
    #   "canvasconfig":PlotTools.Logy},
    # {"variables":["NuEFuzz Neutron","Q3"],
    #  "plot_type" : "sigdep",
    #   "canvasconfig":PlotTools.Logy},
    #  {"variables":["NuEFuzz EM","Lepton Pt"],
    #  "plot_type" : "sigdep",
    #   "canvasconfig":PlotTools.Logy },
    # {"variables":["NuEFuzz EM","Q3"],
    #  "plot_type" : "sigdep",
    #   "canvasconfig":PlotTools.Logy },
    #  {"variables":["NuEFuzz Xtalk","Lepton Pt"],
    #  "plot_type" : "sigdep",
    #   "canvasconfig":PlotTools.Logy },
    # {"variables":["NuEFuzz Xtalk","Q3"],
    #  "plot_type" : "sigdep",
    #   "canvasconfig":PlotTools.Logy },
    #  {"variables":["NuEFuzz Other","Lepton Pt"],
    #  "plot_type" : "sigdep",
    #   "canvasconfig":PlotTools.Logy},
    # {"variables":["NuEFuzz Other","Q3"],
    #  "plot_type" : "sigdep",
    #   "canvasconfig":PlotTools.Logy},
]

