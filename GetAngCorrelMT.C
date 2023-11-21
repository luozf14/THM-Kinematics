//-------------Multithread version of GetAngCorrel-------------//
#include <Math/SpecFunc.h>
#include <ROOT/TProcessExecutor.hxx>
#include <Riostream.h>
#include <TCanvas.h>
#include <TF1.h>
#include <TFile.h>
#include <TGraph.h>
#include <TH1.h>
#include <TH2.h>
#include <TLegend.h>
#include <TLorentzVector.h>
#include <TROOT.h>
#include <TRandom3.h>
#include <TStopwatch.h>
#include <TVector3.h>

// constants
Double_t Ea = 30.;
Double_t EA = 27.6;
Double_t ma = 13. * 931.49 + 3.12500933; // mass of 13C
Double_t mx = 12. * 931.49;              // mass of 12C
Double_t ms = 931.49 + 8.07131806;       // mass of neutron
Double_t mu_sx = ms * mx / (ms + mx);

// 12C+12C
Double_t mA = 12. * 931.49; // mass of 12C
Double_t mb = 4. * 931.49 + 2.42491587;          // mass of 4He
Double_t mB = 20. * 931.49 - 7.04193217; // mass of 20Ne
// Double_t mB = 20. * 931.49 - 7.04193217 + 1.634; // mass of 20Ne_alpha1
// Double_t mb = 931.49 + 7.288971064;      // mass of proton
// Double_t mB = 23. * 931.49 - 9.52985352; // mass of 23Na

// 12C+16O
// Double_t mA = 16. * 931.49 - 4.73700217; // mass of 16O
// Double_t mb = 4. * 931.49 + 2.42491587;  // mass of 4He
// Double_t mB = 24. * 931.49 -13.933578; // mass of 24Mg

// 12C+18O
// Double_t mA = 18. * 931.49 - 0.78281634; // mass of 18O
// Double_t mb = 4. * 931.49 + 2.42491587;  // mass of 4He
// Double_t mB = 26. * 931.49 - 16.214544;  // mass of 26Mg

Double_t mF = mx + mA;
Double_t mu_sF = ms * mF / (ms + mF);
Double_t e_sx = mx + ms - ma;
Double_t ksx = std::sqrt(2. * mu_sx * e_sx);

TGraph *grPsi_p2;

std::vector<TVector3> Kinematics(UInt_t workerID, Double_t &Ecm);
std::vector<TVector3> InverseKinematics(UInt_t workerID, Double_t &Ecm);
Double_t ExA_global;

Double_t MdmAng = 5. / 180. * TMath::Pi();
Double_t SlitDistance = 63.5; //[cm]
Double_t SlitWidth = 2. * std::tan(2. / 180. * TMath::Pi()) * SlitDistance;
TVector3 vecMdm(1, 1, 1);

Double_t SiAng = 110. / 180. * TMath::Pi();
Double_t SiDistance = 14.224; //[cm]
Double_t SiWidth = 5.;        //[cm]
TVector3 vecSi(1, 1, 1);

Double_t flangeOpening = 2. * std::atan(2.5 / 28.67);
Double_t flangeAng = 155. / 180. * TMath::Pi();
TVector3 vecSpectatorDec(1, 1, 1);

const UInt_t nNumbers = 1000000U;
const UInt_t nWorkers = 24U;

void GetAngCorrelMT() {
  TStopwatch timer;
  timer.Start();

  gSystem->Exec("rm myFile.root");
  ROOT::EnableThreadSafety();

  //
  // Read the spectator momentum distribution from txt file
  //
  // Define TGraph of the momentum
  grPsi_p2 = new TGraph();
  grPsi_p2->SetName("grPsi_p2");
  grPsi_p2->SetTitle("|#psi(p)|^{2}");
  grPsi_p2->GetYaxis()->SetTitle("|#psi(p)|^{2}");
  grPsi_p2->GetYaxis()->CenterTitle();
  grPsi_p2->GetXaxis()->SetTitle("p [MeV/c]");
  grPsi_p2->GetXaxis()->CenterTitle();
  // Read file
  std::string inputFile = "Psi_p2.txt";
  std::ifstream dataStream(inputFile.c_str());
  Double_t p, psi_p2;
  Int_t nPoints = 0;
  while (!dataStream.eof()) {
    dataStream >> p >> psi_p2;
    grPsi_p2->SetPoint(nPoints, p, psi_p2);
    nPoints++;
  }
  // Plot momentum distribution
  // TCanvas *c1 = new TCanvas("c1", "c1", 1024, 1024);
  // c1->cd();
  // grPsi_p2->Draw("alp");

  vecMdm.SetMag(SlitDistance);
  vecMdm.SetTheta(MdmAng);
  vecMdm.SetPhi(0.);

  vecSi.SetMag(SiDistance);
  vecSi.SetTheta(SiAng);
  vecSi.SetPhi(TMath::Pi());

  vecSpectatorDec.SetTheta(flangeAng);
  vecSpectatorDec.SetPhi(0.);

  auto workItem = [](UInt_t workerID) {
    TFile f(Form("myFile_%u.root", workerID), "RECREATE");

    TH2D *h2EbVsAngb =
        new TH2D("h2EbVsAngb", "E_{b} vs #theta_{b}", 200, 0, 180, 200, 0, 50);
    h2EbVsAngb->SetXTitle("#theta_{b} [deg]");
    h2EbVsAngb->SetYTitle("E_{b} [MeV]");
    TH2D *h2EBVsAngB =
        new TH2D("h2EBVsAngB", "E_{B} vs #theta_{B}", 200, 0, 35, 200, 0, 2);
    h2EBVsAngB->SetXTitle("#theta_{B} [deg]");
    h2EBVsAngB->SetYTitle("E_{B} [MeV/u]");
    TH2D *h2EsVsAngs =
        new TH2D("h2EsVsAngs", "E_{n} vs #theta_{n}", 200, 0, 180, 200, 0, 5);
    h2EsVsAngs->SetXTitle("#theta_{n} [deg]");
    h2EsVsAngs->SetYTitle("E_{n} [MeV]");

    TH2D *h2AngBVsAngb = new TH2D("h2AngBVsAngb", "#theta_{B} vs #theta_{b}",
                                  200, 0, 180, 200, 0, 35);
    h2AngBVsAngb->SetXTitle("#theta_{b} [deg]");
    h2AngBVsAngb->SetYTitle("#theta_{B} [deg]");
    TH2D *h2AngBVsAngs = new TH2D("h2AngBVsAngs", "#theta_{B} vs #theta_{n}",
                                  200, 0, 180, 200, 0, 35);
    h2AngBVsAngs->SetXTitle("#theta_{n} [deg]");
    h2AngBVsAngs->SetYTitle("#theta_{B} [deg]");
    TH2D *h2AngbVsAngs = new TH2D("h2AngbVsAngs", "#theta_{b} vs #theta_{n}",
                                  200, 0, 180, 200, 0, 180);
    h2AngbVsAngs->SetXTitle("#theta_{n} [deg]");
    h2AngbVsAngs->SetYTitle("#theta_{b} [deg]");

    TH1D *h1Phis = new TH1D("h1Phis", "#phi_{s}", 1000, -180, 180);
    h1Phis->SetXTitle("#phi_{s} [deg]");
    TH1D *h1Phib = new TH1D("h1Phib", "#phi_{b}", 1000, -180, 180);
    h1Phib->SetXTitle("#phi_{b} [deg]");

    TH1D *h1ExA = new TH1D("h1ExA", "ExA", 1000, 0, 20);
    h1ExA->SetXTitle("ExA [MeV]");

    TH1D *h1ExA_sub = new TH1D("h1ExA_sub", "ExA", 1000, 0, 20);
    h1ExA_sub->SetXTitle("ExA [MeV]");

    for (UInt_t i = 0; i < nNumbers; ++i) {
      Double_t Ecm;
      std::vector<TVector3> vecP = Kinematics(workerID, Ecm);
      TVector3 ps_l = vecP[0];
      TVector3 pb_l = vecP[1];
      TVector3 pB_l = vecP[2];
      // if (!(Ecm >= 1. && Ecm <= 2.))
      //     continue;
      // if (pB_l.Mag2() / (2. * mB) / (mB / 931.5) <= 0.75)
      //     continue;

      Double_t theta_b = pb_l.Theta() * 180. / TMath::Pi();
      Double_t theta_B = pB_l.Theta() * 180. / TMath::Pi();
      Double_t theta_s = ps_l.Theta() * 180. / TMath::Pi();

      Double_t phi_b = pb_l.Phi() * 180. / TMath::Pi();
      Double_t phi_B = pB_l.Phi() * 180. / TMath::Pi();
      Double_t phi_s = ps_l.Phi() * 180. / TMath::Pi();

      TVector3 heavyLocal =
          (SlitDistance / std::cos(pB_l.Angle(vecMdm))) * pB_l.Unit() - vecMdm;
      heavyLocal.RotateY(-MdmAng);
      Bool_t heavyDet = std::abs(theta_B - MdmAng * 180. / TMath::Pi()) < 2. &&
                        std::abs(heavyLocal.X()) < SlitWidth / 2. &&
                        std::abs(heavyLocal.Y()) < SlitWidth / 2.;

      TVector3 lightLocal =
          (SiDistance / std::cos(pb_l.Angle(vecSi))) * pb_l.Unit() - vecSi;
      lightLocal.RotateY(SiAng);
      Bool_t lightDet = std::abs(theta_b - SiAng * 180. / TMath::Pi()) < 20. &&
                        std::abs(lightLocal.X()) < SiWidth / 2. &&
                        std::abs(lightLocal.Y()) < SiWidth / 2.;

      Bool_t spectatorDet = ps_l.Angle(vecSpectatorDec) < flangeOpening / 2.;

      if (!(heavyDet && lightDet && spectatorDet))
        continue;
      h1Phib->Fill(phi_b);
      h1Phis->Fill(phi_s);

      h2EbVsAngb->Fill(theta_b, pb_l.Mag2() / (2. * mb));
      h2EBVsAngB->Fill(theta_B, pB_l.Mag2() / (2. * mB) / (mB / 931.5));
      h2EsVsAngs->Fill(theta_s, ps_l.Mag2() / (2. * ms));

      h2AngBVsAngb->Fill(theta_b, theta_B);
      h2AngBVsAngs->Fill(theta_s, theta_B);
      h2AngbVsAngs->Fill(theta_s, theta_b);

      h1ExA->Fill(Ecm);
      if (Ecm > 1. && Ecm < 2.)
        h1ExA_sub->Fill(Ecm);
    }

    h2EbVsAngb->Write();
    h2EBVsAngB->Write();
    h2EsVsAngs->Write();

    h2AngBVsAngb->Write();
    h2AngBVsAngs->Write();
    h2AngbVsAngs->Write();

    h1Phis->Write();
    h1Phib->Write();
    h1ExA->Write();
    h1ExA_sub->Write();

    return 0;
  };

  // Create the collection which will hold the threads, our "pool"
  std::vector<std::thread> workers;

  // Fill the "pool" with workers
  for (auto workerID : ROOT::TSeqI(nWorkers)) {
    workers.emplace_back(workItem, workerID);
  }

  // Now join them
  for (auto &&worker : workers)
    worker.join();

  gSystem->Exec("hadd myFile.root myFile_*.root");
  gSystem->Exec("rm myFile_*.root");

  TFile *outputFile = new TFile("myFile.root", "READ");

  TCanvas *c2 = new TCanvas("c2", "c2", 2700, 1000);
  c2->Divide(3, 1);
  c2->cd(1);
  dynamic_cast<TH2D *>(outputFile->Get("h2EbVsAngb"))->Draw("colz");
  c2->cd(2);
  dynamic_cast<TH2D *>(outputFile->Get("h2EBVsAngB"))->Draw("colz");
  c2->cd(3);
  dynamic_cast<TH2D *>(outputFile->Get("h2EsVsAngs"))->Draw("colz");
  c2->Update();

  TCanvas *c3 = new TCanvas("c3", "c3", 2048, 1024);
  c3->Divide(3, 1);
  c3->cd(1);
  dynamic_cast<TH2D *>(outputFile->Get("h2AngBVsAngb"))->Draw("colz");
  c3->cd(2);
  dynamic_cast<TH2D *>(outputFile->Get("h2AngBVsAngs"))->Draw("colz");
  c3->cd(3);
  dynamic_cast<TH2D *>(outputFile->Get("h2AngbVsAngs"))->Draw("colz");
  c3->Update();

  TCanvas *c4 = new TCanvas("c4", "c4", 1024, 1024);
  c4->Divide(2, 1);
  c4->cd(1);
  dynamic_cast<TH1D *>(outputFile->Get("h1ExA"))->Draw();
  c4->cd(2);
  dynamic_cast<TH1D *>(outputFile->Get("h1ExA_sub"))->Draw();
  c4->Update();

  TCanvas *c5 = new TCanvas("c5", "c5", 1024, 768);
  c5->Divide(2, 1);
  c5->cd(1);
  dynamic_cast<TH1D *>(outputFile->Get("h1Phib"))->Draw();
  c5->cd(2);
  dynamic_cast<TH1D *>(outputFile->Get("h1Phis"))->Draw();
  c5->Update();

  timer.Stop();
  timer.Print();
}

std::vector<TVector3> Kinematics(UInt_t workerID, Double_t &Ecm) {
  //
  // Kinematics for THM a+A->s+F*, F*->b+B
  //
  TRandom3 ran(0);

  // Calculation in lab system
  TVector3 kA(0., 0., std::sqrt(2. * mA * EA)); // momentum vector of A in lab
  // 1. Sample the momentum of s(neutron)
  TVector3 ks(1., 1., 1.); // momentum vector of s in lab
  Double_t ExA = -1000.;
  while (ExA <= 0.) {
    Double_t mag_s = ran.Uniform(0., ksx); // momentum magnitude of s in lab
    Double_t temp = ran.Uniform(0, 0.010);
    while (temp > grPsi_p2->Eval(mag_s, 0, "S")) {
      mag_s = ran.Uniform(0., ksx);
      temp = ran.Uniform(0, 0.010);
    }
    ks.SetMag(mag_s);
    ks.SetTheta(std::acos(ran.Uniform(-1., 1.)));
    ks.SetPhi(ran.Uniform(-TMath::Pi(), TMath::Pi()));
    // 2. Calculate ExA
    ExA = mx / (mx + mA) * EA - ks.Mag2() / (2. * mu_sF) +
          ks.Dot(kA) / (mx + mA) - e_sx;
  }
  // 3. Calculate the momentum of b and B (x+A->b+B)
  TVector3 pF_l = kA - ks; // momentum vector of F(x+A) in lab
  // 3.1 Momentum of b and B in CMS of x+A
  Double_t Q2 = mx + mA - mb - mB; // Q-value of x+A->b+B
  Double_t EbB = ExA + Q2;
  Double_t mag_b = std::sqrt(2. * EbB / (1. / mb + 1. / mB));
  TVector3 kb(1., 1., 1.); // momentum vector of b in CMS of x+A
  kb.SetMag(mag_b);
  kb.SetTheta(std::acos(ran.Uniform(-1., 1.)));
  kb.SetPhi(ran.Uniform(-TMath::Pi(), TMath::Pi()));
  TVector3 kB = -kb; // momentum vector of B in CMS of x+A
  // 3.2 Boost kb and kB to lab system
  TVector3 vF_l = (1. / (mb + mB)) * pF_l; // velocity vector in CMS of b+B
  TVector3 ub = (1. / mb) * kb;            // velocity vector of b in CMS of b+B
  TVector3 uB = (1. / mB) * kB;            // velocity vector of B in CMS of b+B
  TVector3 vb_l = ub + vF_l;               // velocity vector of b in lab
  TVector3 vB_l = uB + vF_l;               // velocity vector of B in lab
  TVector3 pb_l = mb * vb_l;               // momentum vector of b in lab
  TVector3 pB_l = mB * vB_l;               // momentum vector of B in lab

  // TVector3 testP = kA - ks - pb_l - pB_l;
  // std::cout << "momentum excess: " << testP.Mag() << std::endl;
  // std::cout << "energy excess: " << ks.Mag2() / (2. * ms) + pb_l.Mag2() / (2.
  // * mb) + pB_l.Mag2() / (2. * mB) - kA.Mag2() / (2. * mA) - (ma + mA - ms -
  // mb - mB) << std::endl;

  std::vector<TVector3> vec;
  vec.push_back(ks);
  vec.push_back(pb_l);
  vec.push_back(pB_l);
  Ecm = ExA;

  return vec;
}

std::vector<TVector3> InverseKinematics(UInt_t workerID, Double_t &Ecm) {
  TRandom3 ran(0);

  //
  // Kinematics for THM a+A->s+F*, F*->b+B
  //
  // Boost from lab to anti-lab system (Pa=0)
  Double_t mag_a =
      std::sqrt(2. * ma * Ea);      // momentum magnitude of a(13C) in lab
  Double_t va = mag_a / ma;         // velocity of a in lab [c]
  Double_t EA = 0.5 * mA * va * va; // kinetic energy of A(12C) in anti-lab
  TVector3 kA(0., 0.,
              -std::sqrt(2. * mA * EA)); // momentum vector of A in anti-lab
  TVector3 uA = (1. / mA) * kA;          // velocity vector of A in anti-lab

  // Calculation in anti-lab system
  // 1. Sample the momentum of s(neutron)
  TVector3 ks(1., 1., 1.); // momentum vector of s in anti-lab
  Double_t ExA = -1000.;
  while (ExA <= 0.) {
    Double_t mag_s =
        ran.Uniform(0., ksx); // momentum magnitude of s in anti-lab
    Double_t temp = ran.Uniform(0, 0.010);
    while (temp > grPsi_p2->Eval(mag_s, 0, "S")) {
      mag_s = ran.Uniform(0., ksx);
      temp = ran.Uniform(0, 0.010);
    }
    ks.SetMag(mag_s);
    ks.SetTheta(std::acos(ran.Uniform(-1., 1.)));
    ks.SetPhi(ran.Uniform(-TMath::Pi(), TMath::Pi()));
    // 2. Calculate ExA
    ExA = mx / (mx + mA) * EA - ks.Mag2() / (2. * mu_sF) +
          ks.Dot(kA) / (mx + mA) - e_sx;
  }
  // 3. Boost ks to lab system
  TVector3 va_l(0., 0., va);    // momentum vector of a in lab
  TVector3 vs = (1. / ms) * ks; // velocity vector of s in anti-lab
  TVector3 vs_l = vs + va_l;    // velocity vector of s in lab
  TVector3 ps_l = ms * vs_l;    // momentum vector of s in lab
  TVector3 pa_l = ma * va_l;    // momentum vector of a in lab
  TVector3 vA_l = uA + va_l;    // velocity vector of A in lab
  TVector3 pA_l = mA * vA_l;    // momentum vector of A in lab
  TVector3 pF_l = pa_l - ps_l;  // momentum vector of F* in lab
  // 4. Calculate the momentum of b and B (x+A->b+B)
  // 4.1 Momentum of b and B in CMS of x+A
  Double_t Q2 = mx + mA - mb - mB; // Q-value of x+A->b+B
  Double_t EbB = ExA + Q2;
  Double_t mag_b = std::sqrt(2. * EbB / (1. / mb + 1. / mB));
  TVector3 kb(1., 1., 1.); // momentum vector of b in CMS of x+A
  kb.SetMag(mag_b);
  kb.SetTheta(std::acos(ran.Uniform(-1., 1.)));
  kb.SetPhi(ran.Uniform(-TMath::Pi(), TMath::Pi()));
  TVector3 kB = -kb; // momentum vector of B in CMS of x+A
  // 4.2 Boost kb and kB to lab system
  TVector3 vF_l = (1. / (mb + mB)) * pF_l; // velocity vector in CMS of b+B
  TVector3 ub = (1. / mb) * kb;            // velocity vector of b in CMS of b+B
  TVector3 uB = (1. / mB) * kB;            // velocity vector of B in CMS of b+B
  TVector3 vb_l = ub + vF_l;               // velocity vector of b in lab
  TVector3 vB_l = uB + vF_l;               // velocity vector of B in lab
  TVector3 pb_l = mb * vb_l;               // momentum vector of b in lab
  TVector3 pB_l = mB * vB_l;               // momentum vector of B in lab

  // TVector3 testP = pa_l + pA_l - ps_l - pb_l - pB_l;
  // std::cout << "momentum excess: " << testP.Mag() << std::endl;
  // std::cout << "energy excess: " << ps_l.Mag2() / (2. * ms) + pb_l.Mag2() /
  // (2. * mb) + pB_l.Mag2() / (2. * mB) - pa_l.Mag2() / (2. * ma) - pA_l.Mag2()
  // / (2. * mA) - (ma + mA - ms - mb - mB) << std::endl;

  std::vector<TVector3> vec;
  vec.push_back(ps_l);
  vec.push_back(pb_l);
  vec.push_back(pB_l);
  Ecm = ExA;

  return vec;
}
