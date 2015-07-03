#include "TH1.h"
#include "TH2.h"
#include "TROOT.h"
#include "TFile.h"
#include "TF1.h"
#include "TLegend.h"
#include "TCanvas.h"
#include "TColor.h"
#include "TMath.h"
#include "TChain.h"
#include "TProfile.h"
#include "TNtuple.h"
#include <sstream>
#include <iostream>
#include <fstream>

using namespace std;

void set_plot_style()
{
    const Int_t NRGBs = 5;
    const Int_t NCont = 255;

    Double_t stops[NRGBs] = { 0.00, 0.34, 0.61, 0.84, 1.00 };
    Double_t red[NRGBs]   = { 0.00, 0.00, 0.87, 1.00, 0.51 };
    Double_t green[NRGBs] = { 0.00, 0.81, 1.00, 0.20, 0.00 };
    Double_t blue[NRGBs]  = { 0.51, 1.00, 0.12, 0.00, 0.00 };
    TColor::CreateGradientColorTable(NRGBs, stops, red, green, blue, NCont);
    gStyle->SetNumberContours(NCont);
}

void Fit(TH1* h, Double_t& outPos, Double_t& outFWHM, Double_t minpeak, Double_t maxpeak)
{
    // Perform fit.
    
    Char_t tmp[256];

    Double_t x1 = 60;
    Double_t x2 = 300;

    TF1* func = new TF1("func", "gaus(0)+pol2(3)", x1, x2);
    func->SetParameters(h->GetMaximum(), 135, 15, 1, 1, 0.1, 0.1);
    func->SetLineColor(kBlue);
    func->SetParLimits(1, minpeak, maxpeak);
    func->SetParLimits(2, 8, 22);
    h->Fit(func, "RBQO");
    
    // get position and FWHM
    Double_t fPi0Pos = func->GetParameter(1);
    outPos = fPi0Pos;
    outFWHM = 2.35*func->GetParameter(2);

    // indicator line
    TLine line;
    line.SetLineWidth(2);
    line.SetX1(fPi0Pos);
    line.SetX2(fPi0Pos);
    line.SetY1(0);
    line.SetY2(h->GetMaximum());

    TF1* fBG = new TF1("funcBG", "pol2", x1, x2);
    for (Int_t i = 0; i < 3; i++) fBG->SetParameter(i, func->GetParameter(3+i));
    fBG->SetLineColor(kRed);

	delete func;
	delete fBG;
}

// create a new histogram- IM(gg) and Eg as function of a run number
// For each run number: get for each detector element average peak position and average energy. Fill upp ntuples
// When all runs have been looped through, plot as function of detector nr, 
void TAPS_gain()
{
	string path = "/home/adlarson/data2014.07/CaLib/it8newcl/";

	TString pathtohist = "CaLib_TAPS_IM_Neut_TAPS_mult";
	ostringstream s;

	int ifirst_runnr = 4952; // specify first and last run
	int ilast_runnr = 5903;
	int ibinlow = 0;
	int ibinhigh = 438;

	Double_t  peak;
	Double_t  FWHM;
	Double_t* peakrun = new Double_t[(ilast_runnr-ifirst_runnr)*438];

	TH2F** detnr = new TH2F*[438];
	for (int i = 0; i < 438; i++)
	{
		detnr[i] = new TH2F( Form("detnr%03d", i), Form("detnr%03d", i), (ilast_runnr-ifirst_runnr + 20), ifirst_runnr-10, ilast_runnr+10, 400, 100., 170. );	
	}
	TH2F* detnr_allruns	= new TH2F("detnr_allruns", "detnr_allruns", 438, 0, 438, 300, 100., 170.); 
	TH3F* sumAll = new TH3F("sumAll", "sumAll", 500, 0., 1000., 438, 0, 438, 10, 0, 10);

	std::vector<int>	fileNumbers;
	for (int irunnr = ifirst_runnr; irunnr <= ilast_runnr; irunnr++)
	{
		s << path << "Hist_CBTaggTAPS_" << irunnr << ".root";
		TString file1 = s.str();
		f = new TFile( file1 );
		if(f->IsOpen())
		{
			TH3F *rec1 = (TH3F*)f->Get(pathtohist);
			if(rec1)
			sumAll->Add(rec1);
			
			delete rec1;
			fileNumbers.push_back(irunnr);
			f->Close();
		}
		s.str("");
	}
	for( int ibn = ibinlow; ibn < ibinhigh; ibn++)
	{
		TH1F* prpi0av = (TH1F*)sumAll->ProjectionX( "IM(#gamma#gamma) (MeV)", ibn+1, ibn+1, 0, 6)->Clone();
	// insert the average over the full peak projection here
		if( prpi0av->Integral() == 0 ) 
		{
			detnr_allruns->Fill( ibn , 160.00 );
		}
		else
		{
			Double_t help[2];
			Double_t minpeak = 110.0;
			Double_t maxpeak = 147.0;
			Fit(prpi0av, help[0], help[1], minpeak, maxpeak);
			detnr_allruns->Fill( ibn , help[0] );
		}
		delete prpi0av;
	}

	TCanvas* allruns = new TCanvas("allruns", "allruns", 200, 10, 1000, 700);
	allruns->SetGridx();
	allruns->SetGridy();
	detnr_allruns->Draw("BOX");
	allruns->SaveAs("allruns_vs_detel.eps");	// peak position per detector element for all runs specified
	delete allruns;

	ofstream prtofile;
	prtofile.open("TAPS_gain_Jun17_2015.txt");  // name of the gain correction file

	std::vector<int>::iterator it = fileNumbers.begin();
	for( int irunnr = ifirst_runnr; irunnr <= ilast_runnr; irunnr++ )
	{	
		if( irunnr == *(it) &&  it != fileNumbers.end() )
		{
			prtofile << *(it) << "\t" ;
			std::cout << *(it) << std::endl;
			TH3F* summedHist = new TH3F("summedHist", "summedHist", 500, 0, 1000, 438, 0, 438, 10, 0, 10);
			int count = 0;
			for(int i=-5; i<6; i++)   // floating average of 10 runs to get stable peak position for most TAPS elements. 									  
			{							//	  Play around with this to see what works best for you. In the last iteration you may need to increase this for the det elements which still are not stable.
			
				if( (it+i) < fileNumbers.begin() ) continue;
				if( (it+i) >= fileNumbers.end() )	continue;

				count++;
				s << path << "Hist_CBTaggTAPS_" << *(it+i) << ".root";
				TString file1 = s.str();
				f = new TFile( file1 );
				if(f->IsOpen())
				{
					TH3F *rec1 = (TH3F*)f->Get(pathtohist);
					if(rec1)
						summedHist->Add(rec1);
					f->Close();
				}
				s.str("");
			}
			for( int ibn = ibinlow; ibn < ibinhigh; ibn++)
			{

				TH1F* prpi0av = (TH1F*)summedHist->ProjectionX( "IM(#gamma#gamma) (MeV)", ibn+1, ibn+1, 0, 6 )->Clone();  // take only the events where nr of rec clusters in CB + TAPS (TAPS with time cut) < 7. This to avoid the worst of comb background

				if( prpi0av->Integral() == 0 ) 
				{
					detnr[ibn]->Fill( *(it), 160.00 );
					prtofile << std::setprecision(5) << "1.0000" << "\t" ;
					delete prpi0av;
				}

				else
				{	

					prpi0av->Rebin(2); 						// Rebin IM(gg) histogram. This may also need to be changed, depending on statistics. I also tested with Rebin(5)
					prpi0av->Draw();
					Double_t help[2];
					Double_t minpeak = 111.0;				// here I set the min and max value allowed for the peak to be in the fitting procedure. Play around with these values
					Double_t maxpeak = 148.0;


					// if(ibn >= 12 && ibn <= 14)			// it may turn out that for single elements fit will still fail. Perhaps you need to set individual min and max values for these elements for some iteration
					// {									// Not very elegant, but perhaps this also shows that TAPS is a bit more difficult to handle.
					// 	minpeak = 120.0;
					// 	maxpeak = 154.0;
					// }

					Fit(prpi0av, help[0], help[1], minpeak, maxpeak);

					detnr[ibn]->Fill( *(it), help[0] );
					if(prpi0av) delete prpi0av;
				
					if( (help[0] > minpeak) && (help[0] < maxpeak) )
					{	
						detnr[ibn]->Fill( *(it), help[0] );
						prtofile << std::setprecision(5) << (134.9766*134.9766)/(help[0]*help[0]) << "\t" ; 
					}
					else
					{						
						detnr[ibn]->Fill( *(it), 160.0 );
						prtofile << std::setprecision(5) << "1.0000" << "\t" ;
					}
				}
			}
			prtofile << "\n" ;
			if(summedHist) delete summedHist;
			
			++it;
		}
		else
		{
			prtofile << irunnr << "\t" ;
	  		for( int k = 0; k < ibinhigh ; k++ )
	  		{
	  			prtofile << std::setprecision(5) << "1.0000" << "\t" ;
	  		}
	  		prtofile << "\n" ;
	  		if(summedHist) delete summedHist;
	  	}
			
	}
	prtofile.close();

	int nCanvas = 28;
	TCanvas** TAPS_E = new TCanvas*[nCanvas];
	for (int i = 0; i < nCanvas; i++)
	{
		TAPS_E[i] = new TCanvas(Form("TAPS_E%02d", i), Form("TAPS_E%02d", i), 200, 10, 1000, 700);
		TAPS_E[i]->SetFillColor(18);
		TAPS_E[i]->Divide(4,4);
		for( int k = 1; k <= 16; k++)
		{
			TAPS_E[i]->cd(17-k);
			int inr = (438 - (k + 16*i) );
			if(inr < 0) continue;
			gPad->SetGridx();
			gPad->SetGridy();
			detnr[inr]->Draw("BOX");
		}
			TAPS_E[i]->SaveAs(Form("TAPSit4_%02d.eps", i));
	}
	for(int i = 0; i < nCanvas; i++)
		delete TAPS_E[i];
	delete TAPS_E;
	exit(0);
}