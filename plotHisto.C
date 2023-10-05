#include <TCanvas.h>
#include <TLegend.h>
#include <TH1.h>
#include <TH2.h>
#include <Riostream.h>
#include <TROOT.h>
#include <TVector3.h>
#include <TFile.h>

void plotHisto()
{
    TFile *outputFile = new TFile("myFileProton.root", "READ");

    TCanvas *c2 = new TCanvas("c2", "c2", 2700, 1000);
    c2->Divide(3, 1);
    c2->cd(1);
    dynamic_cast<TH2D *>(outputFile->Get("h2EbVsAngb"))->Draw("colz");
    c2->cd(2);
    dynamic_cast<TH2D *>(outputFile->Get("h2EBVsAngB"))->Draw("colz");
    c2->cd(3);
    dynamic_cast<TH2D *>(outputFile->Get("h2EsVsAngs"))->Draw("colz");
    c2->Update();

    TCanvas *c3 = new TCanvas("c3", "c3", 2700, 1000);
    c3->Divide(3, 1);
    c3->cd(1);
    dynamic_cast<TH2D *>(outputFile->Get("h2AngBVsAngb"))->Draw("colz");
    c3->cd(2);
    dynamic_cast<TH2D *>(outputFile->Get("h2AngBVsAngs"))->Draw("colz");
    c3->cd(3);
    dynamic_cast<TH2D *>(outputFile->Get("h2AngbVsAngs"))->Draw("colz");
    c3->Update();

    TCanvas *c4 = new TCanvas("c4", "c4", 1024, 1024);
    // c4->Divide(2, 1);
    c4->cd(1);
    dynamic_cast<TH1D *>(outputFile->Get("h1ExA"))->Draw();
    // c4->cd(2);
    // dynamic_cast<TH1D *>(outputFile->Get("h1ExA_sub"))->Draw();
    c4->Update();

    TCanvas *c5 = new TCanvas("c5", "c5", 1024, 768);
    c5->Divide(2, 1);
    c5->cd(1);
    dynamic_cast<TH1D *>(outputFile->Get("h1Phib"))->Draw();
    c5->cd(2);
    dynamic_cast<TH1D *>(outputFile->Get("h1Phis"))->Draw();
    c5->Update();
}