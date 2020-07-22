//Macro to read geant4 simulated data
//Author: Luis Manzanillas

void macro_read_data(string stringDataName = "output/Bi207-17_04_2020-14-10-21", int nCores = 1){
       //declare tree to read data
       TTree *t_geant4_data = new TTree("t_geant4_data","data from ascii file");	
       for(int i = 0; i < nCores; i++){
		TString totalName = "";
		totalName +=stringDataName;
		totalName +="/";
		totalName +="PenAttenuationSetup_nt_PEN_t";
		totalName +=i;
		totalName +=".csv";
 		cout<<" name "<<totalName<<endl;
		if(i == 0){
			t_geant4_data->ReadFile(totalName,"event:EdepTrigger:EdepPEN1:EdepPEN2:EdepPEN3:EdepPEN4:EdepOther:NTrigger:NPMT1:NPMT2:NPMT3:NPMT4:NPMT5");
		}else{
			t_geant4_data->ReadFile(totalName);
		}	
       }
       float event, EdepTrigger, EdepPEN1, EdepPEN2, EdepPEN3, EdepPEN4, EdepOther, NTrigger, NPMT1, NPMT2, NPMT3, NPMT4, NPMT5;
       t_geant4_data -> SetBranchAddress("event",&event);
       t_geant4_data -> SetBranchAddress("EdepTrigger",&EdepTrigger);
       t_geant4_data -> SetBranchAddress("EdepPEN1",&EdepPEN1);
       t_geant4_data -> SetBranchAddress("EdepPEN2",&EdepPEN2);
       t_geant4_data -> SetBranchAddress("EdepPEN3",&EdepPEN3);
       t_geant4_data -> SetBranchAddress("EdepPEN4",&EdepPEN4);
       t_geant4_data -> SetBranchAddress("EdepOther",&EdepOther);
       t_geant4_data -> SetBranchAddress("NTrigger",&NTrigger);
       t_geant4_data -> SetBranchAddress("NPMT1",&NPMT1);
       t_geant4_data -> SetBranchAddress("NPMT2",&NPMT2);
       t_geant4_data -> SetBranchAddress("NPMT3",&NPMT3);
       t_geant4_data -> SetBranchAddress("NPMT4",&NPMT4);
       t_geant4_data -> SetBranchAddress("NPMT5",&NPMT5);

       //set min and max for histograms, units are in keV
       double min_E_x = 10.;       
       double max_E_x = 2010.;
       double bin_E_width = 5.;
       int n_E_bins = (int)((max_E_x-min_E_x)/bin_E_width);

              
       TH1D* h_edep_trigger = new TH1D("h_edep_trigger","h_edep_PEN_trigger;Edep [keV]; entries / keV ",200,1.,201.);
       TH1D* h_edep_PEN1 = new TH1D("h_edep_PEN1","Edep Sample 1;Edep [keV]; entries / 5 keV ",n_E_bins,min_E_x,max_E_x);
       TH1D* h_edep_PEN2 = new TH1D("h_edep_PEN2","Edep Sample 2;Edep [keV]; entries / 5 keV ",n_E_bins,min_E_x,max_E_x);
       TH1D* h_edep_PEN3 = new TH1D("h_edep_PEN3","Edep Sample 3;Edep [keV]; entries / 5 keV ",n_E_bins,min_E_x,max_E_x);
       TH1D* h_edep_PEN4 = new TH1D("h_edep_PEN4","Edep Sample 4;Edep [keV]; entries / 5 keV ",n_E_bins,min_E_x,max_E_x);
       TH1D* h_edep_PEN_3_4 = new TH1D("h_edep_PEN_3_4","Edep Samples 4+3;Edep [keV]; entries / 1.0 keV ",n_E_bins,min_E_x,max_E_x);
       TH1D* h_edep_all_PEN = new TH1D("h_edep_all_PEN","Edep all samples;Edep [keV]; entries / 1.0 keV ",n_E_bins,min_E_x,max_E_x);
       TH1D* h_edep_all_PEN_res_ten = new TH1D("h_edep_all_PEN_res_ten","Edep all samples 10% Eres;Edep [keV]; entries / 1.0 keV ",n_E_bins,min_E_x,max_E_x);
       TH1D* h_edep_all_PEN_res_five = new TH1D("h_edep_all_PEN_res_five","Edep all samples 5% Eres;Edep [keV]; entries / 1.0 keV ",n_E_bins,min_E_x,max_E_x);
       TH1D* h_edep_all_materials = new TH1D("h_edep_all_materials","Edep all materials;Edep [keV]; entries / 1.0 keV ",n_E_bins,min_E_x,max_E_x);

       double min_Charge_x = 5.;
       double max_Charge_x = 1105.;
       double bin_Charge_width = 5.;
       int n_Charge_bins = (int)((max_Charge_x-min_Charge_x)/bin_Charge_width);
       double max_Charge_x_all = 3000.;
       int n_Charge_bins_all = (int)((3000.-min_Charge_x)/bin_Charge_width);

       TH1D* h_charge_pmt1 = new TH1D("h_total_pmt1","PMT 1; charge [# photons];entries / 5 ph ",n_Charge_bins,min_Charge_x,max_Charge_x);
       TH1D* h_charge_pmt2 = new TH1D("h_total_pmt2","PMT 2; charge [# photons];entries / 5 ph ",n_Charge_bins,min_Charge_x,max_Charge_x);
       TH1D* h_charge_pmt3 = new TH1D("h_total_pmt3","PMT 3; charge [# photons];entries / 5 ph ",n_Charge_bins,min_Charge_x,max_Charge_x);
       TH1D* h_charge_pmt4 = new TH1D("h_total_pmt4","PMT 4; charge [# photons];entries / 5 ph ",n_Charge_bins,min_Charge_x,max_Charge_x);
       TH1D* h_charge_pmt5 = new TH1D("h_total_pmt5","PMT 5; charge [# photons];entries / 5 ph ",n_Charge_bins,min_Charge_x,max_Charge_x);
       TH1D* h_total_charge = new TH1D("h_total_charge","All PMTs; charge [# photons];entries / 5 ph ",n_Charge_bins_all,min_Charge_x,max_Charge_x_all);
       TH1D* h_total_charge_sidePMTs = new TH1D("h_total_charge_sidePMTs","Lateral PMTs; charge [# photons];entries / 5 ph ",n_Charge_bins_all,min_Charge_x,max_Charge_x_all);
       TH1D* h_ratio = new TH1D("h_ratio","h_ratio; Q1/(Q1+Q2);entries / bin ",150,0.,1.);
       TH2D* h_Edep_TotalCharge = new TH2D("h_Edep_TotalCharge","h_Edep_TotalCharge;Edep PEN4[MeV]; # detected photons",n_E_bins,min_E_x,max_E_x,n_Charge_bins,min_Charge_x,max_Charge_x_all);

       int nEntries = t_geant4_data->GetEntries();
       TRandom3 *myRandoms = new TRandom3 ();
       double pmt_resolution = 0.7;
       double scale_factor = 2.0;
       double E_res_ten = 0.1;
       double E_res_five = 0.05;
       for(int i = 0; i < nEntries; ++i){
	       t_geant4_data->GetEntry(i);
	       if(EdepTrigger>15. && EdepTrigger<100.){
		       h_edep_PEN1->Fill(EdepPEN1);
		       h_edep_PEN2->Fill(EdepPEN2);
		       h_edep_PEN3->Fill(EdepPEN3);
		       h_edep_PEN4->Fill(EdepPEN4);
		       h_edep_PEN_3_4->Fill(EdepPEN3+EdepPEN4);
		       double Edep_all_pen = (EdepPEN1+EdepPEN2+EdepPEN3+EdepPEN4)/1000.;
		       h_edep_all_PEN->Fill(Edep_all_pen*1000);
		       h_edep_all_materials->Fill(Edep_all_pen*1000+EdepOther);
		       double Edep_with_res = 1000.*(myRandoms -> Gaus(Edep_all_pen, E_res_ten * sqrt (Edep_all_pen)));
		       h_edep_all_PEN_res_ten->Fill(Edep_with_res);
		       Edep_with_res = 1000.*(myRandoms -> Gaus(Edep_all_pen, E_res_five * sqrt (Edep_all_pen)));
		       h_edep_all_PEN_res_five->Fill(Edep_with_res);
		       //Apply spe resolution
		       double n_ph_pmt1 = myRandoms -> Gaus(NPMT1, pmt_resolution * sqrt (NPMT1))/scale_factor;
		       //double n_ph_pmt1 = myRandoms -> Poisson(NPMT1);
		       double n_ph_pmt2 = myRandoms -> Gaus(NPMT2, pmt_resolution * sqrt (NPMT2))/scale_factor;
		       //double n_ph_pmt2 = myRandoms -> Poisson(NPMT2);
		       double n_ph_pmt3 = myRandoms -> Gaus(NPMT3, pmt_resolution * sqrt (NPMT3))/scale_factor;
		       //double n_ph_pmt3 = myRandoms -> Poisson(NPMT3);
		       double n_ph_pmt4 = myRandoms -> Gaus(NPMT4, pmt_resolution * sqrt (NPMT4))/scale_factor;
		       //double n_ph_pmt4 = myRandoms -> Poisson(NPMT4);
		       double n_ph_pmt5 = myRandoms -> Gaus(NPMT5, pmt_resolution * sqrt (NPMT5))/scale_factor;
		       //double n_ph_pmt5 = myRandoms -> Poisson(NPMT5);
		       h_charge_pmt1->Fill(n_ph_pmt1);
		       h_charge_pmt2->Fill(n_ph_pmt2);
		       h_charge_pmt3->Fill(n_ph_pmt3);
		       h_charge_pmt4->Fill(n_ph_pmt4);
		       h_charge_pmt5->Fill(n_ph_pmt5);
		       h_total_charge->Fill(n_ph_pmt1+n_ph_pmt2+n_ph_pmt3+n_ph_pmt4+n_ph_pmt5);
		       h_total_charge_sidePMTs->Fill(n_ph_pmt1+n_ph_pmt2+n_ph_pmt5+n_ph_pmt4);

	       }
       }
	t_geant4_data->Draw("NPMT1+NPMT2:EdepPEN4>>h_Edep_TotalCharge","","colz");
	TCanvas* c_edeps = new TCanvas();
	gStyle->SetOptStat(0);
	h_edep_PEN2->Draw("HIST PLC");
	h_edep_PEN3->Draw("SAME HIST PLC");
	h_edep_PEN4->Draw("SAME HIST PLC");
	h_edep_PEN1->Draw("SAME HIST PLC");
        TCanvas* c_all = new TCanvas();
        h_edep_all_materials->Draw("HIST PLC");
	h_edep_all_PEN->Draw("SAME HIST PLC");
	h_edep_PEN_3_4->Draw("SAME HIST PLC");
      
        TCanvas* c_eres = new TCanvas();	
	h_edep_all_PEN_res_ten->SetLineColor(kRed);
	h_edep_all_PEN_res_five->Draw("");
	h_edep_all_PEN_res_ten->Draw("HISTSAMES");
     

        TCanvas* c_trigger = new TCanvas();
        t_geant4_data->Draw("EdepTrigger>>h_edep_trigger");
       
        TCanvas* c_charges = new TCanvas();
	// Number of PADS
        const Int_t Nx = 3;
        const Int_t Ny = 2;
	// Margins
        Float_t lMargin = 0.02;
        Float_t rMargin = 0.02;
        // Canvas setup
        c_charges->Divide(Nx,Ny,lMargin,rMargin);
        c_charges->cd(1);
        h_charge_pmt1->Draw(""); 	
	h_charge_pmt1->GetXaxis()->SetRangeUser(5.,800.);
        c_charges->cd(2);
        h_charge_pmt2->Draw(""); 	
	h_charge_pmt2->GetXaxis()->SetRangeUser(5.,800.);
        c_charges->cd(3);
        h_charge_pmt3->Draw(""); 	
	h_charge_pmt3->GetXaxis()->SetRangeUser(5.,800.);
        c_charges->cd(4);
        h_charge_pmt4->Draw(""); 	
	h_charge_pmt4->GetXaxis()->SetRangeUser(5.,800.);
        c_charges->cd(5);
        h_charge_pmt5->Draw(""); 	
	h_charge_pmt5->GetXaxis()->SetRangeUser(5.,800.);
        c_charges->cd(6);
	h_total_charge_sidePMTs->Draw("HIST PLC");
        h_total_charge->Draw("HIST SAME PLC");

	TCanvas* c_lateral = new TCanvas();
        //h_total_charge->Draw("HIST PLC");
	h_total_charge_sidePMTs->Draw("HIST PLC");
        h_total_charge->Draw("HIST SAME PLC");

	std::string str=stringDataName;
	std::size_t pos_start = str.find("_x_");
	std::size_t pos_ends = str.find("_mm_");
	std::string pos_x = str.substr(pos_start+3,pos_ends-pos_start-3);
	cout<<"pos_x: "<<pos_x<<endl;

	TF1* f_q1_qtot = new TF1("f_q1_qtot","gaus",0,400);
        Int_t binmax = h_ratio -> GetMaximumBin();
        Double_t x_max = h_ratio -> GetXaxis()->GetBinCenter(binmax);
	Double_t x_min_fit = 0.;
	Double_t x_max_fit = 0.;
	if(x_max > 0.9 || x_max < 0.1){
        	x_min_fit = x_max - 0.05;
        	x_max_fit = x_max + 0.05;
	}else{
        	x_min_fit = x_max - 0.15;
        	x_max_fit = x_max + 0.15;
	}
        h_ratio -> Fit(f_q1_qtot,"","",x_min_fit,x_max_fit);	

 
	std::size_t pos_abs_start = str.find("_ABS_");
	std::size_t pos_abs_ends = str.find("_SigAlpha_");
	std::string abs_file = str.substr(pos_abs_start+5,pos_abs_ends-pos_abs_start-5);

	std::size_t pos_ly_start = str.find("_LY_");
	std::size_t pos_ly_ends = str.find("ph_MeV");
	std::string ly_file = str.substr(pos_ly_start+4,pos_ly_ends-pos_ly_start-4);

	std::size_t pos_SigmaAlpha_start = str.find("_SigAlpha_");
	std::string SigmaAlpha_file = str.substr(pos_SigmaAlpha_start+10,4);

	TString nameFileLYtxt ="other_values_abs_";
	nameFileLYtxt +=abs_file;
	nameFileLYtxt +="_ly_";
	nameFileLYtxt +=ly_file;
	nameFileLYtxt +="_SigmaAlpha_";
	nameFileLYtxt +=SigmaAlpha_file;
	nameFileLYtxt +=".txt";
	ofstream myfileSaveQAData(nameFileLYtxt,std::ios_base::app);
        if (myfileSaveQAData.is_open())
        {
             myfileSaveQAData<<pos_x<<" "<<1.<<" "<<f_q1_qtot -> GetParameter(1)<<" "<<f_q1_qtot -> GetParError(1)<<" "<<f_q1_qtot -> GetParameter(2)<<" "<<f_q1_qtot -> GetParError(2)<<"\n";
             myfileSaveQAData.close();
        }
        else cout << "Unable to open file";

}
