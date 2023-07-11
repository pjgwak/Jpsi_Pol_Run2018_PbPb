import ROOT
from ROOT import TFile, TCanvas, RooFit, RooRealVar, RooChebychev, RooFit, RooArgList
ROOT.PyConfig.DisableRootLogon = False


#=== Read Inputs ===#
infile = TFile.Open('../skimmedFiles/OniaRooDataSet_isMC1_Psi2S_pp_y0.00_2.40_Effw0_Accw0_PtW0_TnP0_230517.root')
in_dataset = infile.Get('dataset')


#=== Set kinematics ===#
mass_low = 3.3; mass_high = 4.1
pt_low = 6.5; pt_high = 9
y_low = 0; y_high = 1.6
c_low = 0; c_high = 180

# Define cut
kine_cut = f'pt>{pt_low} && pt<{pt_high} && abs(y)>{y_low} && abs(y)<{y_high} && mass>{mass_low} && mass<{mass_high}'
acc_cut = '( ((abs(eta1) <= 1.2) && (pt1 >=3.5)) || ((abs(eta2) <= 1.2) && (pt2 >=3.5)) || ((abs(eta1) > 1.2) && (abs(eta1) <= 2.1) && (pt1 >= 5.47-1.89*(abs(eta1)))) || ((abs(eta2) > 1.2)  && (abs(eta2) <= 2.1) && (pt2 >= 5.47-1.89*(abs(eta2)))) || ((abs(eta1) > 2.1) && (abs(eta1) <= 2.4) && (pt1 >= 1.5)) || ((abs(eta2) > 2.1)  && (abs(eta2) <= 2.4) && (pt2 >= 1.5)) )'
nan_cut = '!TMath::IsNaN(ctau3D) && !TMath::IsNaN(ctau3DRes)'
sign_cut = 'recoQQsign==0'
total_cut = kine_cut + '&&' + acc_cut + '&&' + sign_cut + '&&' + nan_cut

# Applay cuts to dataset
dataset = in_dataset.reduce(total_cut)
#dataset.Print('v')


#=== Set Fit models ===#
# Set paramaters
mass = RooRealVar('mass', 'mass (GeV)',3.3,4.1)
mean = RooRealVar('mean', 'mean of signal', 3.686)
sigma = RooRealVar('sigma', 'width of signal (GeV)', 0.05, 0.01, 0.12)
x_A = RooRealVar("x_A","sigma ratio ", 0.4, 1e-6, 1)
tail_alpha = RooRealVar('tail_alpha', 'tail length of Crystal Ball', 5, 1e-6, 100)
power_n = RooRealVar('power_n', 'power of Crystal Ball (height)', 5, 1e-6,100)
sigma2 = ROOT.RooFormulaVar("sigma2","@0*@1",RooArgList(sigma, x_A) );
tail_alpha2 = ROOT.RooFormulaVar("tail_alpha2","1.0*@0",RooArgList(tail_alpha) );
power_n2 = ROOT.RooFormulaVar("power_n2","1.0*@0",RooArgList(power_n) );
cb_fraction = RooRealVar('cb_fraction', 'fraction of CBs in double CB function', 0.5,0.01,1)

n_signal = RooRealVar('n_signal', 'number of Jpsi', 20000, 0, 350000)

# Build fit model
cb1 = ROOT.RooCBShape('cb1', 'single Crystal Ball', mass, mean, sigma, tail_alpha, power_n)
cb2 = ROOT.RooCBShape('cb2', 'single Crystal Ball', mass, mean, sigma2, tail_alpha2, power_n2)
sig_double_cb =ROOT.RooAddPdf('sig_double_cb', 'Double Crystal Ball', RooArgList(cb1, cb2), RooArgList(cb_fraction))

fit_model = ROOT.RooAddPdf('fit_model', 'dobule CB', RooArgList(sig_double_cb), RooArgList(n_signal))

fit_result = fit_model.fitTo(dataset, Extended=True, Save=True, Minimizer=('Minuit', 'migradimproved '), MaxCalls=10000)


#=== Plotting ===#
c1 = TCanvas('c1', 'c1', 1200, 1200)
c1.cd()
xframe = mass.frame()
dataset.plotOn(xframe)
fit_model.plotOn(xframe)
#dataset.statOn(xframe, Layout=(0.55, 0.99, 0.8))
fit_model.plotOn(xframe, RooFit.Components(cb1), RooFit.LineColor(ROOT.kGreen), RooFit.LineStyle(3)) 
fit_model.plotOn(xframe, RooFit.Components(cb2), RooFit.LineColor(ROOT.kRed), RooFit.LineStyle(3)) 
fit_model.paramOn(xframe)
xframe.Draw()
c1.SaveAs('test.pdf')