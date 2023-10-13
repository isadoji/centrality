Double_t NBDFit(Float_t nhits, Float_t meanMul, Float_t k){

 return  (tgamma(nhits+k)/((tgamma(nhits+1))*(tgamma(k))))*(((TMath::Power(meanMul/k,nhits))/(TMath::Power(meanMul/k +1,nhits+k))));

}

Double_t NANCESTORS(Float_t f, Float_t Npart, Float_t Ncoll){

 return  f*Npart + (1-f)*Ncoll;

}

void ExtractGlauber3(Int_t minmul, Double_t fparameter, Double_t kparameter){
    
    TH1F *hnpart = new TH1F("hnpart","hnpart",450,0,450);
    TH1F *hncoll = new TH1F("hncoll","hncoll",1200,0,1200);
    TH1F *hnancestors = new TH1F("hnancestors","hnancestors",1200,0,1200);

    TFile *InputFile = new TFile("/home/pedro/Softwares/Glauber/gmc-BiBi-snn30.7-md0.4-nd-1.0-rc1-smax99.0.root", "READ");
    TTree *tree = (TTree*)InputFile->Get("nt_Bi_Bi;1");

    TFile *f = new TFile("/home/pedro/Storage/centrality/BiBi9GeV/BiBi9GeV/BiBi9GeV.root", "READ");
    TH1F *hmultiplicity = (TH1F*)f->Get("hRefMultSTAR;1");

    TFile *fo = new TFile("Fit.root","RECREATE");

    Float_t maxmul = hmultiplicity->FindLastBinAbove(1)-minmul;

    Float_t fNpart{-1.};
    Float_t fNcoll{-1.};
    Float_t fNhard{-1.};
    Float_t fB{-1.};

    tree->SetBranchAddress("Npart", &fNpart);
    tree->SetBranchAddress("Ncoll", &fNcoll);
    tree->SetBranchAddress("Nhard", &fNhard);
    tree->SetBranchAddress("B", &fB);

    Int_t Nevents = tree->GetEntries();
    Int_t NeventsMul = hmultiplicity->GetEntries();
    cout<<NeventsMul<<endl;

    for(Int_t i=0; i<Nevents; i++){
        tree->GetEntry(i);
        hnpart->Fill(fNpart);
        hncoll->Fill(fNcoll);

        Int_t Nancestors = NANCESTORS(fparameter,fNpart,fNcoll);
        hnancestors->Fill(Nancestors);
    }

    Float_t Nancestors_max = hnancestors->FindLastBinAbove(1);
    Float_t meanMul = maxmul/Nancestors_max;

    cout<<"Max Num Ancestors= "<<Nancestors_max<<" Max Multiplicity= "<<maxmul<<" Mu= "<<meanMul<<endl;

    Int_t nBins = (meanMul+1.)*3 < 10 ? 10 : (meanMul+1.)*3;
    TH1F *hglauberhisto = new TH1F("hglauberhisto","hglauberhisto",nBins,0,nBins);

    for (Int_t i=0; i<maxmul; ++i) {
        Float_t val = NBDFit(i, meanMul, kparameter);
        if (val>1e-10){
	  hglauberhisto->SetBinContent(i+1, val);
	  cout<<i+1<<" "<<val<<endl;
	}
    }

    TH1F *hglauberfit = new TH1F("hglauberfit","hglauberfit",maxmul*1.3,0,maxmul*1.3); // + 50 just to see it completely

    for (Int_t i=0; i<Nevents; i++){
        tree->GetEntry(i);
        Int_t Na = NANCESTORS(fparameter,fNpart,fNcoll);
        Float_t nHits=0;
        for (Int_t j=0; j<Na; j++) nHits += Int_t(hglauberhisto->GetRandom());
        hglauberfit->Fill(nHits);
    }

    Int_t fGlauberFitHistoInt {0}; 
    Int_t fDataHistoInt {0};
    Int_t lowchibin = minmul;
    Int_t highchibin = maxmul;
    for (Int_t i=lowchibin; i<highchibin; i++){

        fGlauberFitHistoInt += hglauberfit->GetBinContent(i+1);
        fDataHistoInt += hmultiplicity->GetBinContent(i+1);
    }
    Float_t ScaleFactor = (Float_t)fDataHistoInt/fGlauberFitHistoInt;
    hglauberfit->Scale(ScaleFactor); 


    cout<<hmultiplicity->Chi2Test(hglauberfit,"UU")<<endl;

    TCanvas *c1 = new TCanvas("c1","Charged particles multiplicity BMD ");
    c1->SetLogy();
    hmultiplicity->SetStats(kFALSE);
    hmultiplicity->GetXaxis()->SetTitle("Multiplicity");
    hmultiplicity->GetYaxis()->SetTitle("Events number");
    hmultiplicity->SetLineColor(4);
    hmultiplicity->GetXaxis()->SetRange(0,350);
    hmultiplicity->Draw();

    hglauberfit->SetLineColor(2);
    hglauberfit->Draw("sames");

    fo->cd();
    hmultiplicity->Write();
    hnancestors->Write();
    hglauberhisto->Write();
    hglauberfit->Write();
    fo->Close();

}


