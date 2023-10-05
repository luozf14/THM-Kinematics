#include <TCanvas.h>
#include <TLegend.h>
#include <TH1.h>
#include <TH2.h>
#include <TF1.h>
#include <TGraph.h>
#include <Riostream.h>
#include <TROOT.h>
#include <Math/Integrator.h>
#include <Math/SpecFunc.h>
#include <TRandom3.h>
#include <TStopwatch.h>

Double_t h_bar = 197.3; // MeV*fm/c

TGraph *grPhi_r;
TF1 *f1;
Double_t integrator(Double_t *r, Double_t *p);

void GetPsi_p()
{
    ROOT::EnableImplicitMT();

    TStopwatch timer;
    timer.Start();
    ROOT::Math::IntegratorOneDimOptions::SetDefaultIntegrator("Gauss");

    // read bound state wave function phi(r)
    Double_t r, phi_r;
    std::string inputFile = "fort.58";
    std::ifstream dataStream(inputFile.c_str());
    grPhi_r = new TGraph();
    grPhi_r->SetName("grPhi_r");
    grPhi_r->SetTitle("#phi(r)");
    grPhi_r->GetYaxis()->SetTitle("#phi(r)");
    grPhi_r->GetYaxis()->CenterTitle();
    grPhi_r->GetXaxis()->SetTitle("r [fm]");
    grPhi_r->GetXaxis()->CenterTitle();
    Int_t nLines = 0;
    while (!dataStream.eof())
    {
        dataStream >> r >> phi_r;
        grPhi_r->SetPoint(nLines, r, phi_r);
        nLines++;
    }
    Double_t totalPhi_r = grPhi_r->Integral(0., 50.);
    std::cout << "totalPhi_r: " << totalPhi_r << std::endl;

    TGraph *grPhi_rr = new TGraph();
    for (int i = 0; i < 2000; i++)
    {
        Double_t xx = 50. / 2000. * (double)i;
        grPhi_rr->SetPoint(i, xx, grPhi_r->Eval(xx, 0, "S"));
    }

    // Define integration function
    f1 = new TF1("f1", integrator, 0., 50., 1);

    TGraph *grPsi_k = new TGraph();
    grPsi_k->SetName("grPsi_k");
    grPsi_k->SetTitle("#psi(k)");
    grPsi_k->GetYaxis()->SetTitle("#psi(k)");
    grPsi_k->GetYaxis()->CenterTitle();
    grPsi_k->GetXaxis()->SetTitle("k [fm^{-1}]");
    grPsi_k->GetXaxis()->CenterTitle();

    TGraph *grPsi_p2 = new TGraph();
    grPsi_p2->SetName("grPsi_p2");
    grPsi_p2->SetTitle("|#psi(p)|^{2}");
    grPsi_p2->GetYaxis()->SetTitle("|#psi(p)|^{2}");
    grPsi_p2->GetYaxis()->CenterTitle();
    grPsi_p2->GetXaxis()->SetTitle("p [MeV/c]");
    grPsi_p2->GetXaxis()->CenterTitle();

    Int_t nPoints = 500;
    Double_t p[nPoints], k[nPoints], ig[nPoints];
    for (int i = 0; i < nPoints; i++)
    {
        p[i] = 500. / (double)nPoints * (double)i;
        k[i] = p[i] / h_bar;
        f1->SetParameter(0, k[i]);
        ig[i] = f1->Integral(0., 30.);
        grPsi_k->SetPoint(i, k[i], ig[i]);
        grPsi_p2->SetPoint(i, p[i], TMath::Power(ig[i], 2.));
    }

    Double_t totalPsi_p2 = grPsi_p2->Integral(0., 500.);
    std::cout << "totalPsi_p2: " << totalPsi_p2 << std::endl;

    std::ofstream outputFile;
    outputFile.open("Psi_p2_1.txt");
    for (int i = 0; i < nPoints; i++)
    {
        grPsi_p2->SetPoint(i, p[i], TMath::Power(ig[i], 2.) / totalPsi_p2);
        outputFile << p[i] << " " << TMath::Power(ig[i], 2.) / totalPsi_p2 << "\n";
    }
    outputFile.close();

    TCanvas *c1 = new TCanvas("c1", "c1", 1024, 1024);
    c1->cd();
    grPhi_r->SetMarkerSize(2);
    grPhi_r->SetMarkerStyle(21);
    grPhi_r->Draw("ap*");
    grPhi_rr->SetLineColor(kRed);
    grPhi_rr->Draw("same");

    TCanvas *c2 = new TCanvas("c2", "c2", 1024, 1024);
    c2->cd();
    grPsi_k->Draw("alp");

    TCanvas *c3 = new TCanvas("c3", "c3", 1024, 1024);
    c3->cd();
    grPsi_p2->Draw("alp");

    // Sample momentum p using acceptance-rejection sampling method
    // TRandom3 ran(0);
    // TH1D *h1SampleP = new TH1D("h1SampleP", "Sampled p", 100, 0, 500);
    // for (int i = 0; i < 2000; i++)
    // {
    //     Double_t sampleP = ran.Uniform(0, 500);
    //     Double_t temp = ran.Uniform(0, 0.010);
    //     if (temp <= grPsi_p2->Eval(sampleP, 0, "S"))
    //     {
    //         h1SampleP->Fill(sampleP);
    //     }
    // }

    // TCanvas *c4 = new TCanvas("c4", "c4", 1024, 1024);
    // c4->cd();
    // h1SampleP->Draw();
    timer.Stop();
    timer.Print();
}

Double_t integrator(Double_t *x, Double_t *p)
{
    Double_t rr = x[0];
    Double_t kk = p[0];
    return grPhi_r->Eval(rr, 0, "S") * ROOT::Math::sph_bessel(1, rr * kk);
}
