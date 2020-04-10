#include "TROOT.h"
#include "TMath.h"
#include "TRandom.h"
#include "TGraph.h"
#include "TF2.h"
#include "TCanvas.h"
#include "TCutG.h"
#include "Math/Minimizer.h"
#include "Math/Factory.h"
#include "Math/Functor.h"
#include <cmath>
#include <chrono>
#include <iostream>
using namespace std;
using namespace TMath;

void ImpactIt()
{
	gSystem->Beep(650, 300);
	using clock = std::chrono::steady_clock;
	clock::time_point start = clock::now();
	gSystem->Beep(500, 300);

	// 1, 2 st
	TFile *f_input = new TFile("/mnt/d/Work/root/builddir/macros/hists_NICA_5000bins.root");
	analzdc_AuAuss11mb_7sect_no_central_module_target_0_24_57440ev.root
	TH2F *hist = (TH2F *)f_input->Get("pEcalcEdep_nt");
	hist->RebinX(2);
	hist->RebinY(2);
	int Bin_Y = hist->GetYaxis()->GetNbins();
	int Bin_X = hist->GetXaxis()->GetNbins();
	cout << Bin_X <<"  "<< Bin_Y << endl;
	//TCutG* graph_cut = new TCutG("mycut", 4);
	TCanvas *canvas = new TCanvas("c_Et_El_a", "c_Et_El_a", 600, 450);
	TCanvas *canvas2 = new TCanvas("2", "2", 600, 450);

	// 3 st
	//double fNevents = 0, gdTEvents, gdTEvents_1 = 0;
	//double edep = 0, impPar_scint = 0;
	TH2F *pElEt2 = new TH2F("pElEt2", "pElEt2", 1200, 0, 30, 1200, 0, 30);
	TH2F *pElEt_1 = new TH2F("pElEt_1", "pElEt_1", 1200, 0, 30, 1200, 0, 30);
	TH2F *pElEt_2 = new TH2F("pElEt_2", "pElEt_2", 1200, 0, 30, 1200, 0, 30);
	TH2F *pElEt_3 = new TH2F("pElEt_3", "pElEt_3", 1200, 0, 30, 1200, 0, 30);
	TH2F *pElEt_4 = new TH2F("pElEt_4", "pElEt_4", 1200, 0, 30, 1200, 0, 30);
	TH2F *pElEt_5 = new TH2F("pElEt_5", "pElEt_5", 1200, 0, 30, 1200, 0, 30);
	TH2F *pElEt_6 = new TH2F("pElEt_6", "pElEt_6", 1200, 0, 30, 1200, 0, 30);
	TH2F *pElEt_7 = new TH2F("pElEt_7", "pElEt_7", 1200, 0, 30, 1200, 0, 30);
	TH2F *pElEt_8 = new TH2F("pElEt_8", "pElEt_8", 1200, 0, 30, 1200, 0, 30);
	TH2F *pElEt_9 = new TH2F("pElEt_9", "pElEt_9", 1200, 0, 30, 1200, 0, 30);
	TH2F *pElEt_10 = new TH2F("pElEt_10", "pElEt_10", 1200, 0, 30, 1200, 0, 30);
	TH2F *pElEt_11 = new TH2F("pElEt_11", "pElEt_11", 1200, 0, 30, 1200, 0, 30);
	TH2F *pElEt_12 = new TH2F("pElEt_12", "pElEt_12", 1200, 0, 30, 1200, 0, 30);
	TH2F *pElEt_13 = new TH2F("pElEt_13", "pElEt_13", 1200, 0, 30, 1200, 0, 30);
	TH2F *pElEt_14 = new TH2F("pElEt_14", "pElEt_14", 1200, 0, 30, 1200, 0, 30);
	TH2F *pElEt_15 = new TH2F("pElEt_15", "pElEt_15", 1200, 0, 30, 1200, 0, 30);
	TH2F *pElEt_16 = new TH2F("pElEt_16", "pElEt_16", 1200, 0, 30, 1200, 0, 30);
	TH2F *pElEt_17 = new TH2F("pElEt_17", "pElEt_17", 1200, 0, 30, 1200, 0, 30);
	TH2F *pElEt_18 = new TH2F("pElEt_18", "pElEt_18", 1200, 0, 30, 1200, 0, 30);
	TH2F *pElEt_19 = new TH2F("pElEt_19", "pElEt_19", 1200, 0, 30, 1200, 0, 30);
	TH2F *pElEt_20 = new TH2F("pElEt_20", "pElEt_20", 1200, 0, 30, 1200, 0, 30);
	TH2F *pElEt_21 = new TH2F("pElEt_21", "pElEt_21", 1200, 0, 30, 1200, 0, 30);
	TH2F *pElEt_22 = new TH2F("pElEt_22", "pElEt_22", 1200, 0, 30, 1200, 0, 30);
	TH2F *pElEt_23 = new TH2F("pElEt_23", "pElEt_23", 1200, 0, 30, 1200, 0, 30);


	//TFile *analzdc = new TFile();
	//analzdc = TFile::Open("analzdc_AuAuss11mb_7sect_no_central_module_target_0_24_57440ev.root");
	//TTree *theTree = (TTree*)analzdc->Get("nt1");
	//gdTEvents=theTree->GetEntries();
	//theTree = SetBranchAddress("impPar",&impPar);

	
	int count = 0.;
	int k = 0;
	int sum = 0;
	int i_cut = 0;

	TMarker *m;

	cout << hist->GetYaxis()->GetNbins() << " " << hist->GetXaxis()->GetNbins() << endl;
	int s = 0;

	//подсчет кол-ва событий

	for (int i = 0; i < hist->GetXaxis()->GetNbins(); ++i) // integral or i = 1 СДЕЛАЙ С 0 !!!!!!!!!!!!!!1
	{
		for (int j = 0; j < hist->GetYaxis()->GetNbins(); ++j)
		{
			count = hist->GetBinContent(i, j);
			//cout << i<< " " << j << " "<< count << endl;
			if (count > 0)
			{
				sum = sum + count;
				//cout << sum << "\r";
			}
		}
	}

	int arr_cuts[25];
	int n_cuts = 10;
	for (int init_cuts_array = n_cuts; init_cuts_array > 0; --init_cuts_array)
	{
		arr_cuts[init_cuts_array] = round(sum/n_cuts*init_cuts_array);
		cout << init_cuts_array << " " << arr_cuts[init_cuts_array] << endl;
	}
	cout << sum << endl;
	//sleep(1000);
	//уравнение полуоси y=0.839x+4.1137

	double xs, ys, xf, yf, koef, xs_d, ys_d, xf_d, yf_d, koef_d, yorth, xorth, ycheck, xs_rot, ys_rot, xf_rot, yf_rot, xs_rot_d, ys_rot_d, xf_rot_d, yf_rot_d;
	double cut_point_d[24][5];
	double delta_x = 0;
	double xorth_prev;
	double a = 26.725;
	double b = 2.96187;
	double th = 0.6981; //наклон эллипса 0,769760013422076672 // 44.1* градуса
	double x0 = 0.9875;
	double y0 = 4.94198;
	/*
	double a = 23.9928;
	double b = 2.64050;
	double th = 0.7698; //наклон эллипса 0,769760013422076672 // 44.1* градуса
	double x0 = 1.195;
	double y0 = 5.0276;
	*/

	bool BinDup[2500][2500];
	//double means[11];
	for (int initBD = 0; initBD < 2500; initBD++)
	{
		for (int initBD2 = 0; initBD2 < 2500; initBD2++)
		{
			BinDup[initBD][initBD2] = false;
		}
	}
	int points_k = 0;
	int points_k2 = 0;

	int binnumX = 0;
	int binnumY = 0;
	bool DupFlag;
	bool flag_v = true;
	int init_cuts_array = n_cuts;

	/////////////////////////////////////
	////                             ////
	////         //                  ////
	////      ////                   ////
	////       //                    ////
	////     ///// st stage          ////
	/////////////////////////////////////
	


	cout << "***************1**************" << endl;
	for (int i = 0; i < 350000; i++) // along x
	{
		// ===========================================
		// под
		xs_d = i * 0.0001; //считaем точки на эллипсе 1
		ys_d = y0 - (b * sqrt(a * a - xs_d * xs_d + 2 * xs_d * x0 - x0 * x0)) / a;
		xf_d = xs_d + 0.00001; //считaем точки на эллипсе 2 (смещение)
		yf_d = y0 - (b * sqrt(a * a - xf_d * xf_d + 2 * xf_d * x0 - x0 * x0)) / a;

		// под
		xs_rot_d = xs_d * cos(th) - ys_d * sin(th); // получили новые коорд в повернутой СК
		ys_rot_d = xs_d * sin(th) + ys_d * cos(th);
		xf_rot_d = xf_d * cos(th) - yf_d * sin(th); // получили новые коорд в повернутой СК (смещенные)
		yf_rot_d = xf_d * sin(th) + yf_d * cos(th);

		koef_d = (yf_rot_d - ys_rot_d) / (xf_rot_d - xs_rot_d); // коэф наклона прямой через 2 точки на эллипсе под

		//под осью
		for (int l = 0; l < 35000; l++) //идем по перпендикуляру
		{
			xorth = l * 0.001;
			yorth = ys_rot_d - (1 / koef_d) * (xorth - xs_rot_d); //получили ур-е перпендикуляра
																  //рисуем перпендикуляры
			//double orth_count = std::fmod(points_k, round(sum/20));
			if (points_k < arr_cuts[init_cuts_array] /*&&  orth_count < 40*/ && init_cuts_array > 0 && xorth > 0 && yorth > 0 && ((yorth - (0.839 * xorth + 4.1137)) < 0)) //0.065 * xorth + 0.256
			{
				init_cuts_array = init_cuts_array - 1;
				//delta_x = xorth - xorth_prev;
				//if (delta_x > 0.05)
				//{
					/*cut_point_d[0][1] = 0;
					cut_point_d[0][2] = 4;
					cut_point_d[0][4] = 0;
					cut_point_d[0][3] = 0;*/
					i_cut++;
					cut_point_d[i_cut][1] = (ys_rot_d + xs_rot_d / koef_d - 4.1137) / (0.839 + 1 / koef_d); //? +- 
					cut_point_d[i_cut][2] = 0.839 * cut_point_d[i_cut][1] + 4.1137;
					cut_point_d[i_cut][4] = 0;
					cut_point_d[i_cut][3] = koef_d * (ys_rot_d - cut_point_d[i_cut][4]) + xs_rot_d;

					cout << " CUTS | " << cut_point_d[i_cut][1] << " " << cut_point_d[i_cut][2]
						 << " " << cut_point_d[i_cut][3] << " " << cut_point_d[i_cut][4] << " | points_k | " << points_k << " | " <<  init_cuts_array+1
						 << " |arr_cuts| " << arr_cuts[init_cuts_array+1] << " | Delta " << points_k-arr_cuts[init_cuts_array+1] << endl;
				//}
				//xorth_prev = xorth;
			}
			
			//cout << points_k << " " << orth_count << endl;

			binnumX = hist->GetXaxis()->FindBin(xorth);
			binnumY = hist->GetYaxis()->FindBin(yorth);
			if (((yorth - (0.839 * xorth + 4.1137)) < 0) && xorth > 0 && yorth > 0) // && binnumX < 4000 && binnumY < 4000)
			{
				binnumX = hist->GetXaxis()->FindBin(xorth);
				binnumY = hist->GetYaxis()->FindBin(yorth);
				count = hist->GetBinContent(binnumX, binnumY);
				if (BinDup[binnumX][binnumY])
				{
					DupFlag = false;
				}
				else
				{
					BinDup[binnumX][binnumY] = true;
					DupFlag = true;
				}

				if (DupFlag && count > 0)
				{
					points_k2 = points_k2 + count;
					points_k = sum - points_k2;
					//cout << points_k2 <<" / " << sum << "\r";
				}
			}
		}
	}

	i_cut++;
	cut_point_d[i_cut][1] = 28; // x
	cut_point_d[i_cut][2] = 0.839*28 + 4.1137; // y=0,065x+0.253 1.03 * xorth - 4.1137
	cut_point_d[i_cut][3] = 28;
	cut_point_d[i_cut][4] = 0;
	
	i_cut++;
	cut_point_d[i_cut][1] = 0; // x
	cut_point_d[i_cut][2] = 20; // y
	cut_point_d[i_cut][3] = 0;
	cut_point_d[i_cut][4] = 5;
	
	/////////////////////////////////////
	////                             ////
	////    ////////                 ////
	////        ///                  ////
	////    ///                      ////
	////   //////// nd stage         ////
	/////////////////////////////////////

	
	points_k = arr_cuts[init_cuts_array];
    cout << "2nd stage " << endl;
	cout << "Points k " << points_k << endl;
	bool flag_v2 = true;
	delta_x = 0;
	xorth_prev = 0;
	for (int i = 0; i<300000; i++) // along x i = 89000
	{
		// ===========================================
		//над осью
		//cout << i << "\r";
		xs = i*0.0001; //считaем точки на эллипсе 1 над осью
		ys = (b*sqrt(a*a-xs*xs+2*xs*x0-x0*x0))/a+y0;
		xf = xs+0.00001; //считaем точки на эллипсе 2 (смещение) над осью
		yf = (b*sqrt(a*a-xf*xf+2*xf*x0-x0*x0))/a+y0;
	
		//над
		xs_rot = xs*cos(th)-ys*sin(th); // получили новые коорд в повернутой СК
		ys_rot = xs*sin(th)+ys*cos(th);
		xf_rot = xf*cos(th)-yf*sin(th); // получили новые коорд в повернутой СК (смещенные)
		yf_rot = xf*sin(th)+yf*cos(th);
		

		koef = (yf_rot-ys_rot)/(xf_rot-xs_rot); // коэф наклона прямой через 2 точки на эллипсе над
		//полуось y=0,065x+0.253

		//над большой полуосью
		for (int l = 0; l<30000; l++) //идем по перпендикуляру l = 8000
		{
			xorth = l*0.001;
			yorth = ys_rot-(1/koef)*(xorth-xs_rot);
			//рисуем перпендикуляры
			//double orth_count = std::fmod(points_k,round(sum/20));
			/*if (points_k < 5700) 
			{
				cout << " %%%%%%%%%%%%%%% 5700 %%%%%%%%%%%%%%%%%" << endl;
			}*/
			if (points_k < arr_cuts[init_cuts_array] && init_cuts_array > 0 && xorth > 0 && yorth > 0 && ((yorth - (0.839 * xorth + 4.1137)) > 0)) 
			{
				init_cuts_array = init_cuts_array - 1;

				//delta_x = xorth - xorth_prev;
				//if (delta_x > 0.05)
				//{
					i_cut++;
					cut_point_d[i_cut][3] = (ys_rot + xs_rot / koef - 4.1137) / (0.839 + 1 / koef);
					cut_point_d[i_cut][4] = 0.839 * cut_point_d[i_cut][3] + 4.1137;
					cut_point_d[i_cut][2] = 30;
					cut_point_d[i_cut][1] = koef * (ys_rot - cut_point_d[i_cut][2]) + xs_rot;

					cout << " CUTS | " << cut_point_d[i_cut][1] << " " << cut_point_d[i_cut][2]
						 << " " << cut_point_d[i_cut][3] << " " << cut_point_d[i_cut][4] << " | points_k | " << points_k << " | " <<  init_cuts_array+1
						 << " |arr_cuts| " << arr_cuts[init_cuts_array+1] << " | Delta " << points_k-arr_cuts[init_cuts_array+1] << endl;

				//}
				xorth_prev = xorth;
				//cout << xorth << " " << yorth << " " << xs_rot_d << " " << ys_rot_d << endl;
			}
			//cout << points_k << " " << orth_count << endl;

			binnumX = hist->GetXaxis()->FindBin(xorth);
			binnumY = hist->GetYaxis()->FindBin(yorth);
			if (((yorth - (0.839*xs_rot+ 4.1137)) > 0) && xorth>0 && yorth >0)// && binnumX<3000 && binnumY<3000)
			{
				binnumX = hist->GetXaxis()->FindBin(xorth);
				binnumY = hist->GetYaxis()->FindBin(yorth);
				count = hist->GetBinContent(binnumX,binnumY);
				if (BinDup[binnumX][binnumY])
				{
					DupFlag = false;
				}
				else 
				{
					BinDup[binnumX][binnumY] = true;
					DupFlag = true;
				}
				
				if (DupFlag && count>0)
				{
					points_k = points_k - count;
					//cout << points_k << "\r";
				}
			}
		}
	}

	i_cut++;
	cut_point_d[i_cut][1] = 30; // x
	cut_point_d[i_cut][2] = 30;// y 
	cut_point_d[i_cut][3] = 30;
	cut_point_d[i_cut][4] = 0.839*30+4.1137; // y=0,065x+0.253;

	for (i_cut = 0; i_cut < 23; i_cut++)
	{
		cout << " CUTS " << i_cut << " | " << cut_point_d[i_cut][1] << " " << cut_point_d[i_cut][2]
						 << " " << cut_point_d[i_cut][3] << " " << cut_point_d[i_cut][4] << endl;
	}
	sleep(5000);

	
	/////////////////////////////////////
	////    ////////                 ////
	////        ///                  ////
	////      ///                    ////
	////        ///                  ////
	////   //////// rd stage (events)////
	/////////////////////////////////////


	//for(Int_t iEventN=0; iEventN<theTree->GetEntries(); iEventN++)// loop on events
	// { theTree->GetEntry(iEventN); }

	//TFile *_file0 = TFile::Open("analzdc_AuAuss11mb_7sect_no_central_module_target_0_24_57440ev.root");
		//TFile *_file0 = TFile::Open("/mnt/d/Work/root/builddir/macros/analzdc_AuAuss11mb_7sect_no_central_module_target_0_24_57440ev.root");
	//TFile *_file0 = TFile::Open("/mnt/d/Work/root/builddir/macros/analzdc_AuAuss11mb_7sect_no_central_module_target_0_24_57440ev.root")
	//TFile *_file0 = TFile::Open("data.root");
		//TTree *nt1 = (TTree *)_file0->Get("nt1");
		//gdTEvents = nt1->GetEntries();
		//cout << "loop 1 -> gdTEvents " << gdTEvents << endl;

	TCutG *graph_cut = new TCutG("mycut", 5);

	TH1F *hImpPar = new TH1F("hImpPar", "hImpPar", 400, 0., 20.);
	TH1F *hImpPar1 = new TH1F("hImpPar1", "hImpPar1", 400, 0., 20.);
	TH1F *hImpPar2 = new TH1F("hImpPar2", "hImpPar1", 400, 0., 20.);
	TH1F *hImpPar3 = new TH1F("hImpPar3", "hImpPar1", 400, 0., 20.);
	TH1F *hImpPar4 = new TH1F("hImpPar4", "hImpPar1", 400, 0., 20.);
	TH1F *hImpPar5 = new TH1F("hImpPar5", "hImpPar1", 400, 0., 20.);
	TH1F *hImpPar6 = new TH1F("hImpPar6", "hImpPar1", 400, 0., 20.);
	TH1F *hImpPar7 = new TH1F("hImpPar7", "hImpPar1", 400, 0., 20.);
	TH1F *hImpPar8 = new TH1F("hImpPar8", "hImpPar1", 400, 0., 20.);
	TH1F *hImpPar9 = new TH1F("hImpPar9", "hImpPar1", 400, 0., 20.);
	TH1F *hImpPar10 = new TH1F("hImpPar10", "hImpPar10", 400, 0., 20.);
	TH1F *hImpPar11 = new TH1F("hImpPar11", "hImpPar11", 400, 0., 20.);
	TH1F *hImpPar12= new TH1F("hImpPar12", "hImpPar12", 400, 0., 20.);
	TH1F *hImpPar13 = new TH1F("hImpPar13", "hImpPar13", 400, 0., 20.);
	TH1F *hImpPar14 = new TH1F("hImpPar14", "hImpPar14", 400, 0., 20.);
	TH1F *hImpPar15 = new TH1F("hImpPar15", "hImpPar15", 400, 0., 20.);
	TH1F *hImpPar16 = new TH1F("hImpPar16", "hImpPar16", 400, 0., 20.);
	TH1F *hImpPar17 = new TH1F("hImpPar17", "hImpPar17", 400, 0., 20.);
	TH1F *hImpPar18 = new TH1F("hImpPar18", "hImpPar18", 400, 0., 20.);
	TH1F *hImpPar19 = new TH1F("hImpPar19", "hImpPar19", 400, 0., 20.);
	TH1F *hImpPar20 = new TH1F("hImpPar1", "hImpPar1", 400, 0., 20.);
	TH1F *hImpPar21 = new TH1F("hImpPar2", "hImpPar1", 400, 0., 20.);
	TH1F *hImpPar22 = new TH1F("hImpPar2", "hImpPar1", 400, 0., 20.);
	TH1F *hImpPar23 = new TH1F("hImpPar2", "hImpPar1", 400, 0., 20.);


///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	
	gStyle->SetStatW(0.3);
  gStyle->SetStatH(0.2);
  gStyle->SetOptStat(0000);
  gStyle->SetOptFit(00110);
//root6 kBird palette
  Double_t stops[9] = { 0.0000, 0.1250, 0.2500, 0.3750, 0.5000, 0.6250, 0.7500, 0.8750, 1.0000};  
  Double_t red[9]   = { 0.2082, 0.0592, 0.0780, 0.0232, 0.1802, 0.5301, 0.8186, 0.9956, 0.9764};
  Double_t green[9] = { 0.1664, 0.3599, 0.5041, 0.6419, 0.7178, 0.7492, 0.7328, 0.7862, 0.9832};
  Double_t blue[9]  = { 0.5293, 0.8684, 0.8385, 0.7914, 0.6425, 0.4662, 0.3499, 0.1968, 0.0539};
  TColor::CreateGradientColorTable(9, stops, red, green, blue, 255);
  Double_t radToDeg = 180./TMath::Pi(); 
  Double_t degToRad = TMath::Pi()/180.; 

  Double_t radTodeg = 180./TMath::Pi();  
  const Int_t nbCBMmods = 90; //NICA 45+45 mods (VETO1 - 1-45, VETO 2 - 46-89)
  const Int_t nbCBMsect = 7; 
   Float_t rNICA[nbCBMmods],xCur, yCur, yCur1,rsumNICA,rsumNICA_1;
  Double_t xxxxNICA[nbCBMmods],yyyyNICA[nbCBMmods];
  Float_t zzzzNICA[nbCBMsect];
  Double_t fphi_mod[nbCBMmods] = {0.}; 
  Double_t fphi_sect[nbCBMmods][nbCBMsect] = {0.}; 
  Int_t modNb=-1, rNICA_int,rsumNICA_int;
  //TString srNICA[nbCBMmods] = {""};
  TString srNICA[nbCBMmods];
  TString srsumNICA = "";
  Double_t fNevents=0, gdTEvents=0;

  Char_t buf[1000],buf1[1000],buff[1000],buf2[1000],buf3[1000]; 
  Double_t impPar,edep_7sect_1,edep_7sect_2,edep_7sect_mod15_31_1,edep_7sect_mod15_31_2,edep_7sect_mod7_39_1,edep_7sect_mod7_39_2;
  Double_t edep;
  Double_t edepMod_1[91],edepMod_2[91],edepMod_3[91],edepMod_4[91],edepMod_5[91],edepMod_6[91],edepMod_7[91];

  Double_t Ebeam_11 = 4.5;
  Double_t Ebeam_5 = 1.5;

  Int_t icd = 0; Int_t ih;
  TPad *disp[5000];
  TCanvas *canv[5000];
  TH2F *hhClone;
    Int_t bbb;

  Double_t errx[100] = {0};

  //NICA 45 mods, 15x15 hole 10cm in mod 23 VETO 1
  xCur = 60; yCur = 45; 
  for(Int_t iii=0; iii<5; iii++) {
    for(Int_t jjj=0; jjj<7; jjj++) {
      modNb = 7*iii+jjj+5;
      xxxxNICA[modNb] = xCur-(jjj+1)*15.;
      yyyyNICA[modNb] = yCur-(iii+1)*15.;
      cout <<"xx " <<iii <<" " <<jjj <<" " <<modNb+1 <<" " <<xxxxNICA[modNb] <<" " <<yyyyNICA[modNb] <<endl;
    }
  }
  xCur = 45;
  for(Int_t jjj=0; jjj<2; jjj++) {
    for(Int_t iii=0; iii<5; iii++) {
      if(jjj==0) {modNb=iii; yyyyNICA[modNb]=yCur;}
      if(jjj==1) {modNb=40+iii; yyyyNICA[modNb]=-yCur;} 
      xxxxNICA[modNb] = xCur-(iii+1)*15.;
      cout <<"xx " <<iii <<" " <<jjj <<" " <<modNb+1 <<" " <<xxxxNICA[modNb] <<" " <<yyyyNICA[modNb] <<endl;
    }
  }
  
  //NICA 45 mods, 15x15 hole 10cm in mod 68 VETO 2
  for(Int_t jjj=0; jjj<nbCBMmods/2; jjj++) {
    xxxxNICA[jjj+45] = -xxxxNICA[jjj]; 
    //xxxxNICA[jjj+45] = xxxxNICA[jjj]; 
    yyyyNICA[jjj+45] = yyyyNICA[jjj];
  }

  Int_t giEventToPrint = 1000;
  Double_t impPar_centr_bins[10] = {0};
  TFile histoFile("hists_new.root","RECREATE");
  //hists
  TH1F *hImpPar_nt = new TH1F("hImpPar_nt","hImpPar_nt",2000,0.,20.);
  TH2F *pImpEdep_nt = new TH2F("pImpEdep_nt","pImpEdep_nt",80,0.,20.,250,0,50);
  TH2F *pXY_nt = new TH2F("pXY_nt","pXY_nt",7,-52.5,52.5,7,-52.5,52.5);
  TH2F *pXY_impLT6_nt = new TH2F("pXY_impLT6_nt","pXY_impLT6_nt",7,-52.5,52.5,7,-52.5,52.5);
  TH2F *pXY_impGE6_nt = new TH2F("pXY_impGE6_nt","pXY_impGE6_nt",7,-52.5,52.5,7,-52.5,52.5);
  TH2F *pAB_nt = new TH2F("pAB_nt","pAB_nt",400,-30.,10., 600,-5,10.);//ax+ball energy  //2det//11 GeV all plots for conf. with this binning !!!!!!!!!!!!!
  TH2F *pBA_nt = new TH2F("pBA_nt","pBA_nt",280,-2,5.,500,-20.,5.);//ax+ball energy  //2det
  TH2F *pImpNspect_nt = new TH2F("pImpNspect_nt","pImpNspect_nt",80,0,20.,300,0,300);
  TH2F *pImpThetaMax_nt = new TH2F("pImpThetaMax_nt","pImpThetaMax_nt",80,0,20.,1000,0,0.5);
  TH2F *pImpA_nt = new TH2F("pImpA_nt","pImpA_nt",80,0,20.,500,-20.,5.);
  TH2F *pImpB_nt = new TH2F("pImpB_nt","pImpB_nt",80,0.,20.,280,-2,5.);
  TH2F *pEdepB_nt = new TH2F("pEdepB_nt","pEdepB_nt",250,0.,50.,280,-2,5.);
  TH2F *pEdepA_nt = new TH2F("pEdepA_nt","pEdepA_nt",250,0.,50.,500,-20.,5.);
  TH2F *pEdepEcalc_nt = new TH2F("pEdepEcalc_nt","pEdepEcalc_nt",250,0.,50.,250,0,50);
  TH2F *pEcalcEdep_nt = new TH2F("pEcalcEdep_nt","pEcalcEdep_nt",250,0.,50.,250,0,50);
  TH2F *pEdepEcalc_nt_norm = new TH2F("pEdepEcalc_nt_norm","pEdepEcalc_nt_norm",250,0.,50.,250,0,1);//norm po max Ecalc (po Y)
  TH2F *pEcalcEdep_nt_norm = new TH2F("pEcalcEdep_nt_norm","pEcalcEdep_nt_norm",250,0.,50.,250,0,1);//norm po max Edep (po Y)
  TH2F *pEdepEcorr_nt = new TH2F("pEdepEcorr_nt","pEdepEcorr_nt",250,0.,50.,250,0,50);
  TH2F *pImpEdepDivA_nt = new TH2F("pImpEdepDivA_nt","pImpEdepDivA_nt",80,0,20.,10,-1000,1000.);
  TH2F *pImpBDivA_nt = new TH2F("pImpBDivA_nt","pImpBDivA_nt",80,0,20.,10,-500,50.);
  TH2F *pImpEdepMinusEcalc_nt = new TH2F("pImpEdepMinusEcalc_nt","pImpEdepMinusEcalc_nt",80,0,20.,500,-50,50.);
  TH2F *pEdepEdepDivA_nt = new TH2F("pEdepEdepDivA_nt","pEdepEdepDivA_nt",250,0,50.,10,-1000,1000.);
  TH2F *pEdepBDivA_nt = new TH2F("pEdepBDivA_nt","pEdepBDivA_nt",250,0,50.,10,-50,50.);
  TH2F *pEdepEdepMinusEcalc_nt = new TH2F("pEdepEdepMinusEcalc_nt","pEdepEdepMinusEcalc_nt",250,0,50.,10,-50,50.);
  TH1F *hh_EdepEcalc[10]; 
  TH2F *pImpThetaMax_nt_centr_bins[10]; 
  TH2F *pAB_nt_centr_bins[10]; 
  TH2F *pImpEdepMinusEcalc_nt_centr_bins[10]; 
  TH2F *pImpA_nt_centr_bins[10]; 
  TH2F *pImpB_nt_centr_bins[10]; 
  TH2F *pEdepEcalc_nt_centr_bins[10]; 
  TH2F *pEdepEcorr_nt_centr_bins[10]; 
  for(Int_t ii=0;ii<10;ii++) { // icut
    sprintf(buf1,"hh_EdepEcalc%i",ii+1);
    sprintf(buf2,"EdepEcalc,bin%i",ii+1);
    hh_EdepEcalc[ii] = new TH1F(buf1,buf2,2000,0,20); 

    sprintf(buf1,"pAB_nt_centr_bins%i",ii+1);
    sprintf(buf2,"pAB_nt_centr_bins,bin%i",ii+1);
    pAB_nt_centr_bins[ii] = new TH2F(buf1,buf2,400,-30.,10., 600,-5,10.); //ax+b all energy (err=1./sqrt())  //2det//11 GeV all plots with this binning
    //pAB_nt_centr_bins[ii] = new TH2F(buf1,buf2,1000,-20.,5.,560,-2,5.); //ax+b all energy (err=1./sqrt())  //2det//5GeV?
 
    Double_t impPar_centr_bins[10] = {0};

    sprintf(buf1,"pImpEdepMinusEcalc_nt_centr_bins%i",ii+1);
    sprintf(buf2,"pImpEdepMinusEcalc_nt_centr_bins,bin%i",ii+1);
    pImpEdepMinusEcalc_nt_centr_bins[ii] = new TH2F(buf1,buf2,80,0.,20.,500,-50,50); 
    
    sprintf(buf1,"pImpThetaMax_nt_centr_bins%i",ii+1);
    sprintf(buf2,"pImpThetaMax_nt_centr_bins,bin%i",ii+1);
    pImpThetaMax_nt_centr_bins[ii] = new TH2F(buf1,buf2,80,0.,20.,1000,0,0.5); 
    sprintf(buf1,"pImpA_nt_centr_bins%i",ii+1);
    sprintf(buf2,"pImpA_nt_centr_bins,bin%i",ii+1);
    pImpA_nt_centr_bins[ii] = new TH2F(buf1,buf2,80,0,20.,500,-20.,5.);

    sprintf(buf1,"pImpB_nt_centr_bins%i",ii+1);
    sprintf(buf2,"pImpB_nt_centr_bins,bin%i",ii+1);
    pImpB_nt_centr_bins[ii] = new TH2F(buf1,buf2,80,0.,20.,280,-2,5.);
    
    sprintf(buf1,"pEdepEcalc_nt_centr_bins%i",ii+1);
    sprintf(buf2,"pEdepEcalc_nt_centr_bins,bin%i",ii+1);
    pEdepEcalc_nt_centr_bins[ii] = new TH2F(buf1,buf2,250,0,50,500,-50.,50.);

    sprintf(buf1,"pEdepEcorr_nt_centr_bins%i",ii+1);
    sprintf(buf2,"pEdepEcorr_nt_centr_bins,bin%i",ii+1);
    pEdepEcorr_nt_centr_bins[ii] = new TH2F(buf1,buf2,250,0,50,500,-50.,50.);
  }


//   NICA rings   45 (44) mods
//   (45 mods - central module has number 23
//    44 mods - no central module)
//   Inner rectangle
//   1.      16, 22, 24 (23), 30 (29)
//   2.      15, 17, 29 (28), 31 (30)
//   Middle rectangle
//   3.        9, 21, 25 (24), 37 (36)
//   4.        8, 10, 14, 18, 28 (27), 32 (31), 36 (35), 38 (37)
//   5.        7, 11, 35 (34), 39 (38)
//   Outer rectangle
//   6.        3, 20, 26 (25), 43 (42)
//   7.        2, 4, 13, 19, 27 (26), 33 (32), 42 (41), 44 (43) 
//   8.        1, 5, 6, 12, 34 (33), 40 (39), 41 (40), 45 (44)


TTree *pTree_zdc=0;
TTree *pTree_wall=0;
TTree *bigtree_en;
TTree *nt1;

//TPad *disp[50];
  Double_t dist = 319.5; //cm
  const Int_t nbRings = 8;// 
  //const Int_t nbModsInRings[nbRings] = {4, 4, 4, 8, 4, 4, 8, 8};
  const Int_t nbModsInRings[nbRings] = {8, 8, 8, 16, 8, 8, 16, 16};//2det

  Int_t modsInRings[8][16];//2det

  Double_t xxx_edep[10] = {0.,0.01,0.03,0.05,0.07,0.09,0.1,0.12,0.14,0.2};//fit with atan(r/l)
  //Double_t xxx_edep[10] = {2.2,2.5,2.7,3.,3.2,3.3,3.5,3.6,3.7,4.};//fit with eta

  Double_t yyy_imp[10] = {0.,10.,20.,30.,40.,50.,60.,70.,80.,90.};
  //Double_t errx[100] = {0};   
  Double_t erry[100] = {0};   

  Double_t energyInRing_sum_event=0;
  Double_t energyInRing_sum[nbRings];
  Double_t energyInRing[nbRings];
  Double_t ring_radius[nbRings]; 
  Double_t theta_ring[nbRings]; 
  Double_t theta_ring_radian[nbRings]; 
  Double_t ratio[nbRings]; 
  Double_t atan_ring_radius_dist[nbRings]; 
  Double_t eta_ring_radius_dist[nbRings]; 

  Int_t iFinish = 0;

  Int_t nbPoints_fit = 0;
  Int_t nbPoints_new = 0;
  Double_t fit_ratio[nbRings];
  Double_t fit_ratio_zap[nbRings];
  Double_t fit_erry[nbRings];
  Double_t fit_theta[nbRings];

  Double_t ccc_lin, ccc_exp, ccc_lin_radius_1, ccc_lin_radius_2, ccc_lin_atan_1, AA, BB, T0, theta_max, CC, CC_zap;

//Fit functions 
  //TGraphErrors *gr0 = new TGraphErrors(10,xxx_edep,yyy_imp,errx,errx);
  TGraphErrors *gr0 = new TGraphErrors();
  TGraphErrors *gr1 = new TGraphErrors();
  TGraphErrors *gr12 = new TGraphErrors();

  TF1 *f2 = new TF1("f2","[0]*x+[1]",0.,0.17); //A*x+B //atan(r/l) //w/o initial pars//last
  TF1 *f12 = new TF1("f12","[0]",0.,0.17);//C-for const funct//w/o initial pars

  modsInRings[0][0] = 16;
  modsInRings[0][1] = 22;
  modsInRings[0][2] = 24;
  modsInRings[0][3] = 30;

  modsInRings[0][4] = 61;
  modsInRings[0][5] = 67;
  modsInRings[0][6] = 69;
  modsInRings[0][7] = 75;

  //15, 17, 29 (28), 31 (30)
  //60, 62, 74, 76 
  modsInRings[1][0] = 15;
  modsInRings[1][1] = 17;
  modsInRings[1][2] = 29;
  modsInRings[1][3] = 31;

  modsInRings[1][4] = 60;
  modsInRings[1][5] = 62;
  modsInRings[1][6] = 74;
  modsInRings[1][7] = 76;

  //9, 21, 25 (24), 37 (36)
  //54, 66, 70, 82 
  modsInRings[2][0] = 9;
  modsInRings[2][1] = 21;
  modsInRings[2][2] = 25;
  modsInRings[2][3] = 37;

  modsInRings[2][4] = 54;
  modsInRings[2][5] = 66;
  modsInRings[2][6] = 70;
  modsInRings[2][7] = 82;

  //8, 10, 14, 18, 28 (27), 32 (31), 36 (35), 38 (37)
  //53, 55, 59, 63, 73,     77,      81,      83 
  modsInRings[3][0] = 8;
  modsInRings[3][1] = 10;
  modsInRings[3][2] = 14;
  modsInRings[3][3] = 18;
  modsInRings[3][4] = 28;
  modsInRings[3][5] = 32;
  modsInRings[3][6] = 36;
  modsInRings[3][7] = 38;

  modsInRings[3][8] = 53;
  modsInRings[3][9] = 55;
  modsInRings[3][10] = 59;
  modsInRings[3][11] = 63;
  modsInRings[3][12] = 73;
  modsInRings[3][13] = 77;
  modsInRings[3][14] = 81;
  modsInRings[3][15] = 83;


  modsInRings[4][0] = 7;
  modsInRings[4][1] = 11;
  modsInRings[4][2] = 35;
  modsInRings[4][3] = 39;

  modsInRings[4][4] = 52;
  modsInRings[4][5] = 56;
  modsInRings[4][6] = 80;
  modsInRings[4][7] = 84;

  modsInRings[5][0] = 3;
  modsInRings[5][1] = 20;
  modsInRings[5][2] = 26;
  modsInRings[5][3] = 43;

  modsInRings[5][4] = 48;
  modsInRings[5][5] = 65;
  modsInRings[5][6] = 71;
  modsInRings[5][7] = 88;

  modsInRings[6][0] = 2;
  modsInRings[6][1] = 4;
  modsInRings[6][2] = 13;
  modsInRings[6][3] = 19;
  modsInRings[6][4] = 27;
  modsInRings[6][5] = 33;
  modsInRings[6][6] = 42;
  modsInRings[6][7] = 44;

  modsInRings[6][8] = 47;
  modsInRings[6][9] = 49;
  modsInRings[6][10] = 58;
  modsInRings[6][11] = 64;
  modsInRings[6][12] = 72;
  modsInRings[6][13] = 78;
  modsInRings[6][14] = 87;
  modsInRings[6][15] = 89;

  modsInRings[7][0] = 1;
  modsInRings[7][1] = 5;
  modsInRings[7][2] = 6;
  modsInRings[7][3] = 12;
  modsInRings[7][4] = 34;
  modsInRings[7][5] = 40;
  modsInRings[7][6] = 41;
  modsInRings[7][7] = 45;

  modsInRings[7][8] = 46;
  modsInRings[7][9] = 50;
  modsInRings[7][10] = 51;
  modsInRings[7][11] = 57;
  modsInRings[7][12] = 79;
  modsInRings[7][13] = 85;
  modsInRings[7][14] = 86;
  modsInRings[7][15] = 90;


  cout <<"modsInRings " <<modsInRings[0][0] <<" " <<modsInRings[0][3] <<endl;
  cout <<"modsInRings " <<modsInRings[7][0] <<" " <<modsInRings[7][7] <<endl;
  //return;

  //rings radius and theta
  for(Int_t im=0;im<nbRings;im++) {
    ring_radius[im] = TMath::Sqrt(xxxxNICA[modsInRings[im][0]-1]*xxxxNICA[modsInRings[im][0]-1]+yyyyNICA[modsInRings[im][0]-1]*yyyyNICA[modsInRings[im][0]-1]);
    theta_ring[im] = TMath::ATan(ring_radius[im]/dist)*radTodeg;//deg.
    theta_ring_radian[im] = TMath::ATan(ring_radius[im]/dist);//rad.
    atan_ring_radius_dist[im] = TMath::ATan(ring_radius[im]/dist);//rad

    eta_ring_radius_dist[im] = -TMath::Log(TMath::Tan(atan_ring_radius_dist[im]/2.));
    //eta_ring_radius_dist[nbRings-1-im] = -TMath::Log(TMath::Tan(atan_ring_radius_dist[im]/2.));

    cout <<"ring_radius " <<im <<" " <<modsInRings[im][0] <<" " <<xxxxNICA[modsInRings[im][0]-1] <<" " <<yyyyNICA[modsInRings[im][0]-1] <<" " <<ring_radius[im] <<" " <<theta_ring[im] <<" " <<atan_ring_radius_dist[im] <<" " <<eta_ring_radius_dist[nbRings-1-im] <<endl;
  }

   ccc_lin = 2.*dist*TMath::Pi()/225.;//225=15*15 //w/o rdr
   ccc_exp = 2.*TMath::Pi()/225.;//225=15*15 //w/o rdr
   ccc_lin_radius_1 = TMath::Pi()/225.;//225=15*15 //w/o rdr
   ccc_lin_radius_2 = TMath::Pi()/(225.*3.);//225=15*15 //with rdr

   ccc_lin_atan_1 = 2.*dist*dist*TMath::Pi()/(225.);//atat(r/l)//with rdr//2. -> 2 det

  	//TFile *_file0 = TFile::Open("/mnt/d/Work/root/builddir/macros/analzdc_QGSM_AuAu_11_mb_500_38000ev.root");  
	TFile *_file0 = TFile::Open("/mnt/d/Work/root/builddir/macros/analzdc_AuAuss11mb_7sect_no_central_module_target_0_24_57440ev.root","UPDATE");  

    nt1=(TTree*)_file0->Get("nt1");
    gdTEvents=nt1->GetEntries();
    cout <<"loop 1 -> gdTEvents " <<gdTEvents <<endl;

    nt1->SetBranchAddress("impPar",&impPar);
    nt1->SetBranchAddress("edep_7sect_1",&edep_7sect_1);
    nt1->SetBranchAddress("edep_7sect_2",&edep_7sect_2);
    nt1->SetBranchAddress("edepMod_7",&edepMod_7);

//goto FIT;

for(Int_t iEventN=0;iEventN<gdTEvents;iEventN++){// loop on events
      nt1->GetEntry(iEventN);
      hImpPar_nt->Fill(impPar);	  
      pImpEdep_nt->Fill(impPar,(edep_7sect_1+edep_7sect_2));	  
      
      for(Int_t im=0;im<nbCBMmods/2;im++) {
	if(edepMod_7[im+1]>0) {
	  pXY_nt->Fill(xxxxNICA[im],yyyyNICA[im],edepMod_7[im+1]);
	  if(impPar<6) pXY_impLT6_nt->Fill(xxxxNICA[im],yyyyNICA[im],edepMod_7[im+1]);
	  if(impPar>=6) pXY_impGE6_nt->Fill(xxxxNICA[im],yyyyNICA[im],edepMod_7[im+1]);
	}
      }      
       
    }
    

   Double_t integral;
   for(Int_t ic=0; ic<10; ic++) {
     integral = 0.1*(ic+1)*hImpPar_nt->Integral();
     cout <<"centr " <<ic+1 <<" " <<integral <<endl;
     for(Int_t i=0; i<1600; i++) {
       bbb = 1*(i+1);
       //if(hImpPar_nt->Integral(0.,bbb)>5700) cout <<"bbb " <<bbb <<" " <<hImpPar_nt->Integral(0.,bbb) <<endl;
       if(TMath::Abs(hImpPar_nt->Integral(0.,bbb)-integral)<=36.) {
	 cout <<"compare " <<bbb <<" " <<integral <<" " <<hImpPar_nt->Integral(0.,bbb) <<" " <<ic+1 <<endl;
	 impPar_centr_bins[ic] = (Double_t)bbb/100.;
	 if(ic==9) impPar_centr_bins[ic] = 17.;
	 break;
       }
     }
     //cout <<"marina " <<endl;
   }
//return;

   cout <<"impPar_centrality " <<endl;
   for(Int_t ic=0; ic<10; ic++) {
     //cout <<impPar_centr_bins[ic] <<" " ;
     cout <<ic+1 <<" " <<impPar_centr_bins[ic] <<" " <<hImpPar_nt->Integral(0.,impPar_centr_bins[ic]*100) <<" " <<100.*hImpPar_nt->Integral(0.,impPar_centr_bins[ic]*100)/hImpPar_nt->Integral() <<endl; 
  }
   //cout <<endl;

   //return;

FIT:
// draw fit results
Double_t Ecalc=0;
Double_t Ecorr=0;
Double_t Edep=0;

//loop on events
   for(Int_t im=0;im<nbRings;im++) {
     energyInRing_sum[im] = 0;
   }
   for(Int_t iEventN=0;iEventN<gdTEvents;iEventN++){// loop on events
 
      nt1->GetEntry(iEventN);
   //fit in each event      
      nbPoints_fit = 0;      
      Ecorr = 0; Ecalc = 0; Edep = 0;
      energyInRing_sum_event = 0;

      if(impPar<20) {//all centralities
      //if(impPar>6) {
	for(Int_t im=0;im<nbRings;im++) {
	  energyInRing[im] = 0;
	  for(Int_t imm=0;imm<nbModsInRings[im];imm++) {
	    energyInRing_sum[im] += edepMod_7[modsInRings[im][imm]];
	    energyInRing[im] += edepMod_7[modsInRings[im][imm]];

	  }//imm
  
	  if(energyInRing[im]>0) {	    
	    energyInRing_sum_event += energyInRing[im];

	    fit_theta[nbPoints_fit] = atan_ring_radius_dist[im];//for fit with lin function//atan(r/l)//last
	    fit_ratio[nbPoints_fit] = energyInRing[im]/nbModsInRings[im];
	    fit_ratio_zap[nbPoints_fit] = energyInRing[im]/nbModsInRings[im];
	    fit_erry[nbPoints_fit] = 1./(energyInRing[im]);//new errors
	    nbPoints_fit++; 
	  }//if(energyInRing[im]>0) 

	}//ring


	if(nbPoints_fit==nbRings) {
	  Ecorr = 0;
	  for(Int_t im=0;im<nbCBMmods;im++) {
	    if(edepMod_7[im+1]>0) {
	      Ecorr += (edepMod_7[im+1] - fit_ratio[nbPoints_fit-1]);  
	      //Ecorr += (edepMod_7[im+1]); //test 
	    }
	  }
	  
	  for(Int_t im=0;im<nbPoints_fit;im++) {
	    //cout <<"before " <<im <<" " <<fit_ratio[im] <<endl; 
	    fit_ratio[im] -= fit_ratio[nbPoints_fit-1];
	    //cout <<"after " <<im <<" " <<fit_ratio[im] <<endl; 
	    //if(im>3) fit_erry[im] /= 10.;
	  }
	  CC = 10.;
	  CC_zap = CC;

	FITFIT:
	  //TGraphErrors *gr1 = new TGraphErrors(nbPoints_fit,fit_theta,fit_ratio,errx,fit_erry); 
	  gr1 = new TGraphErrors(nbPoints_fit,fit_theta,fit_ratio,errx,fit_erry); 
	if(iEventN<gdTEvents-1) gr1->Fit("f2","RQN");//works, not draw fit
	//if(iEventN<41940) gr1->Fit("f2","RQN");//works, not draw fit
	else gr1->Fit("f2","RQ"); //f2 - line or expo (TF1...)
	  AA = f2->GetParameter(0);
	  BB = f2->GetParameter(1);
	  cout <<"iEventN,AA,BB before " <<iEventN <<" " <<AA <<" " <<BB <<endl;
	  //if(CC<=0.02) goto finish; //ev 25

	//from fit with lin funct
	AA = f2->GetParameter(0);//comment if use procedure with bg
	BB = f2->GetParameter(1);//comment if use procedure with bg
    ////////////////////
    ////////////////////
  	cout <<"marina 1 " <<endl;	  
	  T0 = TMath::Tan(-BB/AA);
	  theta_max = T0;
	  if(AA!=0) {
	    if((-AA)>=0) {
	      Ecalc = ccc_lin_atan_1*(AA*((T0*T0+1)*TMath::ATan(T0) - T0) + BB*T0*T0);
	    }//if((-AA)>=0)
	  }
	  else Ecalc = 0;
    Edep = edep_7sect_1+edep_7sect_2;
	//pElEt_1->Fill(Ecalc, Edep);
	hImpPar->Fill(impPar);
	for (int iii = 0; iii < 12; iii++)
		{ //проверяем, к какому сектору относится точка
			graph_cut->SetPoint(0, cut_point_d[iii][1], cut_point_d[iii][2]);
			graph_cut->SetPoint(1, cut_point_d[iii][3], cut_point_d[iii][4]);
			graph_cut->SetPoint(2, cut_point_d[iii + 1][3], cut_point_d[iii + 1][4]);
			graph_cut->SetPoint(3, cut_point_d[iii + 1][1], cut_point_d[iii + 1][2]);
			graph_cut->SetPoint(4, cut_point_d[iii][1], cut_point_d[iii][2]);

			if (graph_cut->IsInside(Ecalc, Edep) && iii == 0)
			{
				hImpPar1->Fill(impPar);
				pElEt_1->Fill(Ecalc, Edep);
				cout << "++1+++++++++++++++" << endl;
			}
			else if (graph_cut->IsInside(Ecalc, Edep) && iii == 1)
			{
				hImpPar2->Fill(impPar);
				pElEt_2->Fill(Ecalc, Edep);
				cout << "+++++2++++++++++++" << endl;
			}
			else if (graph_cut->IsInside(Ecalc, Edep) && iii == 2)
			{
				hImpPar3->Fill(impPar);
				pElEt_3->Fill(Ecalc, Edep);
				cout << "++++++++3+++++++++" << endl;
			}
			else if (graph_cut->IsInside(Ecalc, Edep) && iii == 3)
			{
				hImpPar4->Fill(impPar);
				pElEt_4->Fill(Ecalc, Edep);
				cout << "+++++++++++4++++++" << endl;
			}
			else if (graph_cut->IsInside(Ecalc, Edep) && iii == 4)
			{
				hImpPar5->Fill(impPar);
				pElEt_5->Fill(Ecalc, Edep);
				cout << "+++++++++++++5++++" << endl;
			}
			else if (graph_cut->IsInside(Ecalc, Edep) && iii == 5)
			{
				hImpPar6->Fill(impPar);
				pElEt_6->Fill(Ecalc, Edep);
				cout << "+++++++++++++++6++" << endl;
			}
			else if (graph_cut->IsInside(Ecalc, Edep) && iii == 6)
			{
				hImpPar7->Fill(impPar);
				pElEt_7->Fill(Ecalc, Edep);
				cout << "++++++++++++++++7+" << endl;
			}
			/*else if (graph_cut->IsInside(Ecalc, Edep) && iii == 7)
			{
				hImpPar8->Fill(impPar);
				pElEt_8->Fill(Ecalc, Edep);
				cout << "++++++++++++++++8+" << endl;
			}*/
			else if (graph_cut->IsInside(Ecalc, Edep) && iii == 8)
			{
				hImpPar9->Fill(impPar);
				pElEt_9->Fill(Ecalc, Edep);
				cout << "+++++++++++9++++++" << endl;
			}
			else if (graph_cut->IsInside(Ecalc, Edep) && iii == 9)
			{
				hImpPar10->Fill(impPar);
				pElEt_10->Fill(Ecalc, Edep);
				cout << "+++++10++++++++++++" << endl;
			}
			else if (graph_cut->IsInside(Ecalc, Edep) && iii == 10)
			{
				hImpPar11->Fill(impPar);
				pElEt_11->Fill(Ecalc, Edep);
				cout << "++11+++++++++++++++" << endl;
			}
			else if (graph_cut->IsInside(Ecalc, Edep) && iii == 11)
			{
				hImpPar12->Fill(impPar);
				pElEt_12->Fill(Ecalc, Edep);
				cout << "++11+++++++++++++++" << endl;
			}
			else if (graph_cut->IsInside(Ecalc, Edep) && iii == 12)
			{
				hImpPar13->Fill(impPar);
				pElEt_13->Fill(Ecalc, Edep);
				cout << "++11+++++++++++++++" << endl;
			}
			/*else if (graph_cut->IsInside(Ecalc, Edep) && iii == 13)
			{
				hImpPar14->Fill(impPar);
				pElEt_14->Fill(Ecalc, Edep);
				cout << "++11+++++++++++++++" << endl;
			}*/
			else if (graph_cut->IsInside(Ecalc, Edep) && iii == 14)
			{
				hImpPar15->Fill(impPar);
				pElEt_15->Fill(Ecalc, Edep);
				cout << "++11+++++++++++++++" << endl;
			}
			else if (graph_cut->IsInside(Ecalc, Edep) && iii == 15)
			{
				hImpPar16->Fill(impPar);
				pElEt_16->Fill(Ecalc, Edep);
				cout << "++11+++++++++++++++" << endl;
			}
			else if (graph_cut->IsInside(Ecalc, Edep) && iii == 16)
			{
				hImpPar17->Fill(impPar);
				pElEt_17->Fill(Ecalc, Edep);
				cout << "++11+++++++++++++++" << endl;
			}
			else if (graph_cut->IsInside(Ecalc, Edep) && iii == 17)
			{
				hImpPar18->Fill(impPar);
				pElEt_18->Fill(Ecalc, Edep);
				cout << "++11+++++++++++++++" << endl;
			}
			else if (graph_cut->IsInside(Ecalc, Edep) && iii == 18)
			{
				hImpPar19->Fill(impPar);
				pElEt_19->Fill(Ecalc, Edep);
				cout << "++11+++++++++++++++" << endl;
			}
			else if (graph_cut->IsInside(Ecalc, Edep) && iii == 19)
			{
				hImpPar20->Fill(impPar);
				pElEt_20->Fill(Ecalc, Edep);
				cout << "++11+++++++++++++++" << endl;
			}
			else if (graph_cut->IsInside(Ecalc, Edep) && iii == 20)
			{
				hImpPar21->Fill(impPar);
				pElEt_21->Fill(Ecalc, Edep);
				cout << "++11+++++++++++++++" << endl;
			}
			else if (graph_cut->IsInside(Ecalc, Edep) && iii == 21)
			{
				hImpPar22->Fill(impPar);
				pElEt_22->Fill(Ecalc, Edep);
				cout << "++11+++++++++++++++" << endl;
			}
			/*else if (graph_cut->IsInside(Ecalc, Edep) && iii == 22)
			{
				hImpPar23->Fill(impPar);
				pElEt_23->Fill(Ecalc, Edep);
				cout << "++11+++++++++++++++" << endl;
			}*/
		}

	  //pEdepEcalc_nt->Fill((edep_7sect_1+edep_7sect_2),Ecalc);
	}//if(nbPoints_fit>1)
}//if(impPar...)

}//for(Int_t iEventN=0;iEventN<gdTEvents;iEventN++)
///////////////////////////////////////////////////////////////////////////////////////////////////////////
		
		
		//}

	//cout << Et << " " << El << endl;
	//pElEt2->Fill(El, Et);
	//hImpPar->Write();

	//hImpPar1->Write();
		
	hImpPar1->Fit("gaus");
	hImpPar2->Fit("gaus");
	hImpPar3->Fit("gaus");
	hImpPar4->Fit("gaus");
	hImpPar5->Fit("gaus");
	hImpPar6->Fit("gaus");
	hImpPar7->Fit("gaus"); //,"","",5,11);
	//hImpPar8->Fit("gaus"); //,"","",4,10);
	hImpPar9->Fit("gaus");
	hImpPar10->Fit("gaus"); //,"","",0,5);
	hImpPar11->Fit("gaus");
	hImpPar12->Fit("gaus");
	/*hImpPar13->Fit("gaus");
	hImpPar14->Fit("gaus");
	hImpPar15->Fit("gaus");
	hImpPar16->Fit("gaus");
	hImpPar17->Fit("gaus");
	hImpPar18->Fit("gaus"); 
	hImpPar19->Fit("gaus");
	hImpPar20->Fit("gaus");
	hImpPar21->Fit("gaus");
	hImpPar22->Fit("gaus");*/
	//hImpPar23->Fit("gaus");

		hImpPar->SetLineColor(1);
		hImpPar->Draw();
	hImpPar1->SetLineColor(1);
	hImpPar1->Draw("same");
	hImpPar2->SetLineColor(2);
	hImpPar2->Draw("same");
	hImpPar3->SetLineColor(3);
	hImpPar3->Draw("same");
	hImpPar4->SetLineColor(4);
	hImpPar4->Draw("same");
	hImpPar5->SetLineColor(5);
	hImpPar5->Draw("same");
	hImpPar6->SetLineColor(6);
	hImpPar6->Draw("same");
	hImpPar7->SetLineColor(14);
	hImpPar7->Draw("same");	
	//hImpPar8->SetLineColor(8);
	//hImpPar8->Draw("same");
	hImpPar9->SetLineColor(1);
	hImpPar9->Draw("same");
	hImpPar10->SetLineColor(2);
	hImpPar10->Draw("same");
	hImpPar11->SetLineColor(3);
	hImpPar11->Draw("same");
	hImpPar12->SetLineColor(8);
	hImpPar12->Draw("same");
	/*hImpPar13->SetLineColor(5);
	hImpPar13->Draw("same");
	hImpPar14->SetLineColor(6);
	hImpPar14->Draw("same");
	hImpPar15->SetLineColor(7);
	hImpPar15->Draw("same");
	hImpPar16->SetLineColor(8);
	hImpPar16->Draw("same");*/

	//hImpPar4->Draw("same");
	//hImpPar5->Draw("same");
	//pElEt2->Draw();
	/*pElEt_1->SetMarkerColor(1);
	pElEt_1->Draw();
	//TFile histoFileFull("hists_NICA.root","RECREATE");
	//pElEt_1->Write();
	pElEt_2->SetMarkerColor(2);
	pElEt_2->Draw("same");
	pElEt_3->SetMarkerColor(3);
	pElEt_3->Draw("same");
	pElEt_4->SetMarkerColor(4);
	pElEt_4->Draw("same");
	pElEt_5->SetMarkerColor(5);
	pElEt_5->Draw("same");
	pElEt_6->SetMarkerColor(6);
	pElEt_6->Draw("same");
	pElEt_7->SetMarkerColor(7);
	pElEt_7->Draw("same");
	pElEt_8->SetMarkerColor(8);
	pElEt_8->Draw("same");
	pElEt_9->SetMarkerColor(9);
	pElEt_9->Draw("same");
	pElEt_10->SetMarkerColor(1);
	pElEt_10->Draw("same");
	pElEt_11->SetMarkerColor(2);
	pElEt_11->Draw("same");
	pElEt_12->SetMarkerColor(3);
	pElEt_12->Draw("same");
	pElEt_13->SetMarkerColor(1);
	pElEt_13->Draw("same");
	//pElEt_14->SetMarkerColor(5);
	//pElEt_14->Draw("same");
	pElEt_15->SetMarkerColor(6);
	pElEt_15->Draw("same");
	pElEt_16->SetMarkerColor(8);
	pElEt_16->Draw("same");
	pElEt_17->SetMarkerColor(9);
	pElEt_17->Draw("same");
	pElEt_18->SetMarkerColor(1);
	pElEt_18->Draw("same");
	pElEt_19->SetMarkerColor(2);
	pElEt_19->Draw("same");
	pElEt_20->SetMarkerColor(3);
	pElEt_20->Draw("same");
	pElEt_21->SetMarkerColor(4);
	pElEt_21->Draw("same");
	pElEt_22->SetMarkerColor(5);
	pElEt_22->Draw("same");
*/
	//canvas->SetCansvasSize(600,600);
	//canvas->SetWindowSize(500, 500);

	//g1->Draw("AP");
	//g2->Draw("AP");

	//canvas->SaveAs("test.png");
	clock::time_point end = clock::now();
	clock::duration execution_time = end - start;
	gSystem->Beep(432, 1000);
	cout << "Complete in " << (execution_time.count() / (1000000000)) << " s" << endl;
}
