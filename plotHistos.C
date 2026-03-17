#include "TClonesArray.h"
#include "TFile.h"
#include "TH1D.h"
#include "TH3D.h"
#include "TH1F.h"
#include "TH2F.h"
#include "THnSparse.h"
#include "TList.h"
#include "TLorentzVector.h"
#include "TMath.h"
#include "TNtuple.h"
#include "TProfile.h"
#include "TRandom3.h"
#include "TString.h"
#include "TTree.h"
#include "TVector3.h"
#include <algorithm>
#include <cstring>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <math.h>
#include <sstream>
#include <stdlib.h>
#include <string>
#include <vector>
#include <cstdio>
#include <ctime>
#include <time.h>
#include <valarray>
#include <chrono>
#include "TCanvas.h"
#include "TPad.h"
#include "TStyle.h"
#include "TLegend.h"

void StyleHisto(TH1F *histo, Float_t Low, Float_t Up, Int_t color, Int_t style, TString TitleX, TString TitleY, TString title, Float_t mSize, Float_t xOffset, Float_t yOffset)
{
  histo->GetYaxis()->SetRangeUser(Low, Up);
  histo->SetLineColor(color);
  histo->SetMarkerColor(color);
  histo->SetMarkerStyle(style);
  histo->SetMarkerSize(mSize);
  histo->GetXaxis()->SetTitle(TitleX);
  histo->GetXaxis()->SetTitleSize(0.05);
  histo->GetXaxis()->SetLabelSize(0.05);
  histo->GetXaxis()->SetTitleOffset(xOffset);
  histo->GetYaxis()->SetTitle(TitleY);
  histo->GetYaxis()->SetTitleSize(0.05);
  histo->GetYaxis()->SetTitleOffset(yOffset);
  histo->GetYaxis()->SetLabelSize(0.05);
  histo->SetTitle(title);
}

void SetFont(TH1F *histo)
{
  histo->GetXaxis()->SetTitleFont(43);
  histo->GetXaxis()->SetLabelFont(43);
  histo->GetYaxis()->SetTitleFont(43);
  histo->GetYaxis()->SetLabelFont(43);
}
void SetTickLength(TH1F *histo, Float_t TickLengthX, Float_t TickLengthY)
{
  histo->GetXaxis()->SetTickLength(TickLengthX);
  histo->GetYaxis()->SetTickLength(TickLengthY);
}

void SetHistoTextSize(TH1F *histo, Float_t XSize, Float_t XLabelSize, Float_t XOffset, Float_t XLabelOffset, Float_t YSize, Float_t YLabelSize, Float_t YOffset, Float_t YLabelOffset)
{
  histo->GetXaxis()->SetTitleSize(XSize);
  histo->GetXaxis()->SetLabelSize(XLabelSize);
  histo->GetXaxis()->SetTitleOffset(XOffset);
  histo->GetXaxis()->SetLabelOffset(XLabelOffset);
  histo->GetYaxis()->SetTitleSize(YSize);
  histo->GetYaxis()->SetLabelSize(YLabelSize);
  histo->GetYaxis()->SetTitleOffset(YOffset);
  histo->GetYaxis()->SetLabelOffset(YLabelOffset);
}

void StyleCanvas(TCanvas *canvas, Float_t TopMargin, Float_t BottomMargin, Float_t LeftMargin, Float_t RightMargin)
{
  canvas->SetFillColor(0);
  canvas->SetTickx(1);
  canvas->SetTicky(1);
  gPad->SetTopMargin(TopMargin);
  gPad->SetLeftMargin(LeftMargin);
  gPad->SetBottomMargin(BottomMargin);
  gPad->SetRightMargin(RightMargin);
  gStyle->SetLegendBorderSize(0);
  gStyle->SetLegendFillColor(0);
  gStyle->SetLegendFont(42);
}

void StylePad(TPad *pad, Float_t LMargin, Float_t RMargin, Float_t TMargin, Float_t BMargin)
{
  pad->SetFillColor(0);
  pad->SetTickx(1);
  pad->SetTicky(1);
  pad->SetLeftMargin(LMargin);
  pad->SetRightMargin(RMargin);
  pad->SetTopMargin(TMargin);
  pad->SetBottomMargin(BMargin);
}

TString TitleRelDistance = "#sqrt{#Delta#phi^{2} + #Delta y^{2}}";
TString TitlePairs = "#it{N}_{#Lambda#bar{#Lambda} pairs}";
TString TitleInvMass = "#it{M}_{#Lambda#bar{#Lambda}} (GeV/#it{c}^{2})";

Float_t xTitle = 30;
Float_t xOffset = 1.3;
Float_t yTitle = 38;
Float_t yOffset = 1.7;

Float_t xLabel = 30;
Float_t yLabel = 30;
Float_t xLabelOffset = 0.015;
Float_t yLabelOffset = 0.01;

Float_t tickX = 0.03;
Float_t tickY = 0.025;

void plotHistos(int nEvents = 1000)
{

  // --- open histo file ---
  TString Sinputfile = Form("../DoubleLambdaHistos_%d.root", nEvents);
  TFile *f = new TFile(Sinputfile, "");
  if (!f || f->IsZombie())
  {
    std::cerr << "Error opening file: " << Sinputfile << std::endl;
    return;
  }
  cout << "Opened file: " << Sinputfile << std::endl;

  THnSparseD *hsparse_same = (THnSparseD *)f->Get("hsparse_same");

  TH1D *hRelativeDistance = hsparse_same->Projection(0);
  TH1D *hPairInvMass = hsparse_same->Projection(1);
  TH1D *hPairIsPrimary = hsparse_same->Projection(2);
  TH1D *hPairIsPair = hsparse_same->Projection(3);

  TCanvas *canvas = new TCanvas("canvas", "canvas", 1200, 800);
  canvas->Divide(2, 2);
  canvas->cd(1);
  gStyle->SetOptStat(0);
  hRelativeDistance->SetTitle("Relative Distance;#sqrt{#Delta#phi^{2} + #Delta y^{2}};#it{N}_{#Lambda#bar{#Lambda} pairs}");
  hRelativeDistance->Draw("HIST");
  canvas->cd(2);
  hPairInvMass->SetTitle("Invariant Mass;M_{#Lambda#bar{#Lambda}} (GeV/#it{c}^{2});#it{N}_{#Lambda#bar{#Lambda} pairs}");
  hPairInvMass->Draw("HIST");
  canvas->cd(3);
  hPairIsPrimary->SetTitle("Pair is Primary;#it{N}_{#Lambda#bar{#Lambda} pairs}");
  hPairIsPrimary->Draw("HIST");
  canvas->cd(4);
  hPairIsPair->SetTitle("Pair Is Pair;#it{N}_{#Lambda#bar{#Lambda} pairs}");
  hPairIsPair->Draw("HIST");

  TH1F *hRelDistanceAll = (TH1F *)hRelativeDistance->Clone("hRelDistanceAll");
  TH1F *hRelDistanceAllNorm = (TH1F *)hRelativeDistance->Clone("hRelDistanceAllNorm");
  hRelDistanceAllNorm->Scale(1. / hRelDistanceAll->Integral());
  TH1F *hInvMassAll = (TH1F *)hPairInvMass->Clone("hInvMassAll");
  TH1F *hInvMassAllNorm = (TH1F *)hPairInvMass->Clone("hInvMassAllNorm");
  hInvMassAllNorm->Scale(1. / hInvMassAll->Integral());
  int Npairs = hRelativeDistance->GetEntries();
  cout << Npairs << " pairs in total" << endl;

  hsparse_same->GetAxis(2)->SetRange(2, 2); // primary pairs
  TH1F *hRelDistancePrimary = (TH1F *)hsparse_same->Projection(0);
  hRelDistancePrimary->SetName("hRelDistancePrimary");
  int Nprim =  hRelDistancePrimary->GetEntries();
  cout << Nprim << " primary pairs in total" << endl;
  TH1F *hRelDistancePrimaryNorm = (TH1F *)hRelDistancePrimary->Clone("hRelDistancePrimaryNorm");
  hRelDistancePrimaryNorm->Rebin(2);
  hRelDistancePrimaryNorm->Scale(1. / 2 / hRelDistancePrimaryNorm->Integral());
  TH1F *hInvMassPrimary = (TH1F *)hsparse_same->Projection(1);
  TH1F *hInvMassPrimaryNorm = (TH1F *)hInvMassPrimary->Clone("hInvMassPrimaryNorm");
  hInvMassPrimaryNorm->Scale(1. / hInvMassPrimary->Integral());
  TH1F *hRatioPrimaryAll = (TH1F *)hRelDistancePrimary->Clone("hRatioPrimaryAll");
  hRatioPrimaryAll->Divide(hRelDistanceAll);
  TH1F *hRatioMassPrimaryAll = (TH1F *)hInvMassPrimary->Clone("hRatioMassPrimaryAll");
  hRatioMassPrimaryAll->Divide(hInvMassAll);

  hsparse_same->GetAxis(3)->SetRange(2, 2); // primary pairs from same mother
  TH1F *hRelDistanceSameMother = (TH1F *)hsparse_same->Projection(0);
  hRelDistanceSameMother->SetName("hRelDistanceSameMother");
  int NprimSM = hRelDistanceSameMother->GetEntries();
  cout << NprimSM << " primary pairs from same mother" << endl;
  TH1F *hRelDistanceSameMotherNorm = (TH1F *)hRelDistanceSameMother->Clone("hRelDistanceSameMotherNorm");
  hRelDistanceSameMotherNorm->Rebin(2);
  hRelDistanceSameMotherNorm->Scale(1. / 2 / hRelDistanceSameMotherNorm->Integral());
  TH1F *hInvMassSameMother = (TH1F *)hsparse_same->Projection(1);
  hInvMassSameMother->SetName("hInvMassSameMother");
  TH1F *hInvMassSameMotherNorm = (TH1F *)hInvMassSameMother->Clone("hInvMassSameMotherNorm");
  hInvMassSameMotherNorm->Scale(1. / hInvMassSameMother->Integral());
  TH1F *hRatioSameMotherPrimary = (TH1F *)hRelDistanceSameMother->Clone("hRatioSameMotherPrimary");
  hRatioSameMotherPrimary->Divide(hRelDistanceAll);
  TH1F *hRatioMassSameMotherPrimary = (TH1F *)hInvMassSameMother->Clone("hRatioMassSameMotherPrimary");
  hRatioMassSameMotherPrimary->Divide(hInvMassAll);

  hsparse_same->GetAxis(3)->SetRange(1, 1); // primary pairs from different mothers
  TH1F *hRelDistanceDiffMother = (TH1F *)hsparse_same->Projection(0);
  hRelDistanceDiffMother->SetName("hRelDistanceDiffMother");
  int NprimDM = hRelDistanceDiffMother->GetEntries();
  cout << NprimDM << " primary pairs from different mothers" << endl;
  TH1F *hRelDistanceDiffMotherNorm = (TH1F *)hRelDistanceDiffMother->Clone("hRelDistanceDiffMotherNorm");
  hRelDistanceDiffMotherNorm->Rebin(2);
  hRelDistanceDiffMotherNorm->Scale(1. / 2 / hRelDistanceDiffMotherNorm->Integral());
  TH1F *hInvMassDiffMother = (TH1F *)hsparse_same->Projection(1);
  hInvMassDiffMother->SetName("hInvMassDiffMother");
  TH1F *hInvMassDiffMotherNorm = (TH1F *)hInvMassDiffMother->Clone("hInvMassDiffMotherNorm");
  hInvMassDiffMotherNorm->Scale(1. / hInvMassDiffMother->Integral());

  hsparse_same->GetAxis(3)->SetRange(0, -1);
  hsparse_same->GetAxis(2)->SetRange(1, 1); // secondary pairs
  TH1F *hRelDistanceSecondary = (TH1F *)hsparse_same->Projection(0);
  hRelDistanceSecondary->SetName("hRelDistanceSecondary");
  TH1F *hRelDistanceSecondaryNorm = (TH1F *)hRelDistanceSecondary->Clone("hRelDistanceSecondaryNorm");
  hRelDistanceSecondaryNorm->Scale(1. / hRelDistanceSecondary->Integral());
  TH1F *hInvMassSecondary = (TH1F *)hsparse_same->Projection(1);
  hInvMassSecondary->SetName("hInvMassSecondary");
  TH1F *hInvMassSecondaryNorm = (TH1F *)hInvMassSecondary->Clone("hInvMassSecondaryNorm");
  hInvMassSecondaryNorm->Scale(1. / hInvMassSecondary->Integral());

  TCanvas *fractionCanvas = new TCanvas("fractionCanvas", "Fraction Canvas", 1200, 800);
  StyleCanvas(fractionCanvas, 0.03, 0.15, 0.18, 0.03);
  TH1F *hDummyFraction = new TH1F("hDummyFraction", "Fraction of pairs", 3, -0.5, 2.5);
  StyleHisto(hDummyFraction, 0, 0.15, 1, 1, "Category", "Fraction", "", 1.15, 1.3, 1.7);
  hDummyFraction->SetBinContent(1, (float)Nprim / Npairs);
  hDummyFraction->SetBinContent(2, (float)NprimSM / Npairs);
  hDummyFraction->SetBinContent(3, (float)NprimDM / Npairs);
  hDummyFraction->GetXaxis()->SetBinLabel(1, "Primary/Total");
  hDummyFraction->GetXaxis()->SetBinLabel(2, "Primary Same Mother/Total");
  hDummyFraction->GetXaxis()->SetBinLabel(3, "Primary Diff Mother/Total");
  hDummyFraction->Draw("HIST");
  fractionCanvas->SaveAs("../FractionPairs.pdf");
  fractionCanvas->SaveAs("../FractionPairs.png");
  fractionCanvas->SaveAs("../FractionPairs.eps");

  TCanvas *relDistanceCanvas = new TCanvas("relDistanceCanvas", "Relative Distance", 800, 600);
  StyleCanvas(relDistanceCanvas, 0.03, 0.15, 0.18, 0.03);
  TH1F *hDummyRelDistance = (TH1F *)hRelativeDistance->Clone("hDummyRelDistance");
  hDummyRelDistance->Reset();
  SetFont(hDummyRelDistance);
  StyleHisto(hDummyRelDistance, 0, 0.03, 1, 1, TitleRelDistance, TitlePairs, "", 1, 1.15, 1.6);
  SetHistoTextSize(hDummyRelDistance, xTitle, xLabel, xOffset, xLabelOffset, yTitle, yLabel, yOffset, yLabelOffset);
  SetTickLength(hDummyRelDistance, tickX, tickY);
  hDummyRelDistance->GetXaxis()->SetRangeUser(0, 3.49);
  StyleHisto(hRelDistanceAllNorm, 0, 1, kRed + 1, 22, TitleRelDistance, TitlePairs, "", 1, 1.15, 1.6);
  StyleHisto(hRelDistancePrimaryNorm, 0, 1, kBlue + 1, 33, TitleRelDistance, TitlePairs, "", 1, 1.15, 1.6);
  StyleHisto(hRelDistanceSameMotherNorm, 0, 1, kGreen + 1, 20, TitleRelDistance, TitlePairs, "", 1, 1.15, 1.6);
  StyleHisto(hRelDistanceDiffMotherNorm, 0, 1, kMagenta + 1, 22, TitleRelDistance, TitlePairs, "", 1, 1.15, 1.6);
  StyleHisto(hRelDistanceSecondaryNorm, 0, 1, kCyan + 1, 20, TitleRelDistance, TitlePairs, "", 1, 1.15, 1.6);
  StyleHisto(hRelDistanceAll, 0, 1, kRed + 1, 22, TitleRelDistance, TitlePairs, "", 1, 1.15, 1.6);
  StyleHisto(hRelDistancePrimary, 0, 1, kBlue + 1, 33, TitleRelDistance, TitlePairs, "", 1, 1.15, 1.6);
  StyleHisto(hRelDistanceSameMother, 0, 1, kGreen + 1, 20, TitleRelDistance, TitlePairs, "", 1, 1.15, 1.6);
  StyleHisto(hRelDistanceDiffMother, 0, 1, kMagenta + 1, 22, TitleRelDistance, TitlePairs, "", 1, 1.15, 1.6);
  StyleHisto(hRelDistanceSecondary, 0, 1, kCyan + 1, 20, TitleRelDistance, TitlePairs, "", 1, 1.15, 1.6);
  relDistanceCanvas->cd();
  hDummyRelDistance->Draw("");
  hRelDistanceAllNorm->Draw("same");
  hRelDistancePrimaryNorm->Draw("same");
  hRelDistanceSameMotherNorm->Draw("same");
  hRelDistanceDiffMotherNorm->Draw("same");
  hRelDistanceSecondaryNorm->Draw("same");
  // hRelDistanceAll->Draw("same");
  // hRelDistancePrimary->Draw("same");
  // hRelDistanceSameMother->Draw("same");
  // hRelDistanceDiffMother->Draw("same");
  // hRelDistanceSecondary->Draw("same");
  TLegend *legendRelDistance = new TLegend(0.37, 0.61, 0.67, 0.91);
  legendRelDistance->SetFillStyle(0);
  legendRelDistance->SetTextAlign(12);
  legendRelDistance->SetTextSize(0.04);
  legendRelDistance->AddEntry(hRelDistanceAllNorm, "All pairs", "pl");
  legendRelDistance->AddEntry(hRelDistancePrimaryNorm, "Primary pairs", "pl");
  legendRelDistance->AddEntry(hRelDistanceSameMotherNorm, "Primary pairs from same mother", "pl");
  legendRelDistance->AddEntry(hRelDistanceDiffMotherNorm, "Primary pairs from different mothers", "pl");
  legendRelDistance->AddEntry(hRelDistanceSecondaryNorm, "Non-primary pairs", "pl");
  legendRelDistance->Draw("");
  relDistanceCanvas->SaveAs("../RelativeDistance.pdf");
  relDistanceCanvas->SaveAs("../RelativeDistance.png");
  relDistanceCanvas->SaveAs("../RelativeDistance.eps");

  TCanvas *invMassCanvas = new TCanvas("invMassCanvas", "Invariant Mass", 800, 600);
  StyleCanvas(invMassCanvas, 0.03, 0.15, 0.18, 0.03);
  TH1F *hDummyInvMass = (TH1F *)hPairInvMass->Clone("hDummyInvMass");
  hDummyInvMass->Reset();
  SetFont(hDummyInvMass);
  StyleHisto(hDummyInvMass, 0, 0.12, 1, 1, TitleInvMass, TitlePairs, "", 1, 1.15, 1.6);
  SetHistoTextSize(hDummyInvMass, xTitle, xLabel, xOffset, xLabelOffset, yTitle, yLabel, yOffset, yLabelOffset);
  SetTickLength(hDummyInvMass, tickX, tickY);
  hDummyInvMass->GetXaxis()->SetRangeUser(2, 8);
  StyleHisto(hInvMassAllNorm, 0, 1, kRed + 1, 22, TitleInvMass, TitlePairs, "", 1, 1.15, 1.6);
  StyleHisto(hInvMassPrimaryNorm, 0, 1, kBlue + 1, 33, TitleInvMass, TitlePairs, "", 1, 1.15, 1.6);
  StyleHisto(hInvMassSameMotherNorm, 0, 1, kGreen + 1, 20, TitleInvMass, TitlePairs, "", 1, 1.15, 1.6);
  StyleHisto(hInvMassDiffMotherNorm, 0, 1, kMagenta + 1, 22, TitleInvMass, TitlePairs, "", 1, 1.15, 1.6);
  StyleHisto(hInvMassSecondaryNorm, 0, 1, kCyan + 1, 20, TitleInvMass, TitlePairs, "", 1, 1.15, 1.6);
  invMassCanvas->cd();
  hDummyInvMass->Draw("");
  hInvMassAllNorm->Draw("same");
  hInvMassPrimaryNorm->Draw("same");
  hInvMassSameMotherNorm->Draw("same");
  hInvMassDiffMotherNorm->Draw("same");
  hInvMassSecondaryNorm->Draw("same");
  legendRelDistance->Draw("");
  invMassCanvas->SaveAs("../InvariantMass.pdf");
  invMassCanvas->SaveAs("../InvariantMass.png");
  invMassCanvas->SaveAs("../InvariantMass.eps");

  TCanvas *fracRelDistanceCanvas = new TCanvas("fracRelDistanceCanvas", "Fraction of Primary Pairs", 800, 600);
  StyleCanvas(fracRelDistanceCanvas, 0.03, 0.15, 0.18, 0.03);
  TH1F *hDummyFracRelDistance = (TH1F *)hRelativeDistance->Clone("hDummyFracRelDistance");
  hDummyFracRelDistance->Reset();
  SetFont(hDummyFracRelDistance);
  StyleHisto(hDummyFracRelDistance, 0, 0.3, 1, 1, TitleRelDistance, "Fraction of all pairs", "", 1, 1.15, 1.6);
  SetHistoTextSize(hDummyFracRelDistance, xTitle, xLabel, xOffset, xLabelOffset, yTitle, yLabel, yOffset, yLabelOffset);
  SetTickLength(hDummyFracRelDistance, tickX, tickY);
  hDummyFracRelDistance->GetXaxis()->SetRangeUser(0, 3.49);
  StyleHisto(hRatioPrimaryAll, 0, 1, kBlue + 1, 33, TitleRelDistance, "", "", 1, 1.15, 1.6);
  StyleHisto(hRatioSameMotherPrimary, 0, 1, kGreen + 1, 20, TitleRelDistance, "", "", 1, 1.15, 1.6);
  fracRelDistanceCanvas->cd();
  hDummyFracRelDistance->Draw("");
  hRatioPrimaryAll->Draw("same");
  hRatioSameMotherPrimary->Draw("same");
  TLegend *legendFracRelDistance = new TLegend(0.23, 0.73, 0.53, 0.92);
  legendFracRelDistance->SetFillStyle(0);
  legendFracRelDistance->SetTextAlign(12);
  legendFracRelDistance->SetTextSize(0.04);
  legendFracRelDistance->AddEntry(hRatioPrimaryAll, "Primary pairs / All pairs", "pl");
  legendFracRelDistance->AddEntry(hRatioSameMotherPrimary, "Primary pairs from same mother / All pairs", "pl");
  legendFracRelDistance->Draw("");
  fracRelDistanceCanvas->SaveAs("../FractionPrimaryPairs_RelDistance.pdf");
  fracRelDistanceCanvas->SaveAs("../FractionPrimaryPairs_RelDistance.png");
  fracRelDistanceCanvas->SaveAs("../FractionPrimaryPairs_RelDistance.eps");

  TCanvas *fracInvMassCanvas = new TCanvas("fracInvMassCanvas", "Fraction of Primary Pairs", 800, 600);
  StyleCanvas(fracInvMassCanvas, 0.03, 0.15, 0.18, 0.03);
  TH1F *hDummyFracInvMass = (TH1F *)hPairInvMass->Clone("hDummyFracInvMass");
  hDummyFracInvMass->Reset();
  SetFont(hDummyFracInvMass);
  StyleHisto(hDummyFracInvMass, 0, 0.3, 1, 1, TitleInvMass, "Fraction of all pairs", "", 1, 1.15, 1.6);
  SetHistoTextSize(hDummyFracInvMass, xTitle, xLabel, xOffset, xLabelOffset, yTitle, yLabel, yOffset, yLabelOffset);
  SetTickLength(hDummyFracInvMass, tickX, tickY);
  hDummyFracInvMass->GetXaxis()->SetRangeUser(2, 6);
  StyleHisto(hRatioMassPrimaryAll, 0, 1, kBlue + 1, 33, TitleInvMass, "", "", 1, 1.15, 1.6);
  StyleHisto(hRatioMassSameMotherPrimary, 0, 1, kGreen + 1, 20, TitleInvMass, "", "", 1, 1.15, 1.6);
  fracInvMassCanvas->cd();
  hDummyFracInvMass->Draw("");
  hRatioMassPrimaryAll->Draw("same");
  hRatioMassSameMotherPrimary->Draw("same");
  TLegend *legendFracInvMass = new TLegend(0.23, 0.73, 0.53, 0.92);
  legendFracInvMass->SetFillStyle(0);
  legendFracInvMass->SetTextAlign(12);
  legendFracInvMass->SetTextSize(0.04);
  legendFracInvMass->AddEntry(hRatioMassPrimaryAll, "Primary pairs / All pairs", "pl");
  legendFracInvMass->AddEntry(hRatioMassSameMotherPrimary, "Primary pairs from same mother / All pairs", "pl");
  legendFracInvMass->Draw("");
  fracInvMassCanvas->SaveAs("../FractionPrimaryPairs_InvMass.pdf");
  fracInvMassCanvas->SaveAs("../FractionPrimaryPairs_InvMass.png");
  fracInvMassCanvas->SaveAs("../FractionPrimaryPairs_InvMass.eps");

  // write to output file
  TString Sfileout = Form("../PlotDoubleLambdaHistos_%d.root", nEvents);
  TFile *fout = new TFile(Sfileout, "recreate");
  fout->cd();

  fout->Close();
  cout << "Histograms saved to file: " << Sfileout << endl;
}