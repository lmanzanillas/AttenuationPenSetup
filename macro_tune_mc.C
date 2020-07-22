//Macro to read geant4 simulated data
//Author: Luis Manzanillas

void macro_tune_mc(string stringDataName = "output/Bi207-17_04_2020-14-10-21", int nCores = 1){
       //declare trees to read data
       TTree *t_data_pmt1 = new TTree("t_data_pmt1","data txt pmt1");
       TTree *t_data_pmt2 = new TTree("t_data_pmt2","data txt pmt2");
       TTree *t_data_pmt3 = new TTree("t_data_pmt3","data txt pmt3");
       TTree *t_data_pmt4 = new TTree("t_data_pmt4","data txt pmt4");
       TTree *t_data_pmt5 = new TTree("t_data_pmt5","data txt pmt5");
       TTree *t_data_pmt6 = new TTree("t_data_pmt6","data txt pmt6");


       //TString Base_name = "/remote/ceph/group/gedet/data/pen/2020/2020-01-28_8182cf65_lt_6pmt_first_data/data_luis_settings2/text_data_for_fit/Spectrum_hv_optimal_setup_1_EJ200_30_30_10mm3_settings2_5PMTs-center-with-reflector-grease2-measurement-1_pmt_";
       TString Base_name = "/remote/ceph/group/gedet/data/pen/2020/2020-01-28_8182cf65_lt_6pmt_first_data/data_luis_settings2/text_data_for_fit/Spectrum_hv_optimal_setup_1_EJ200_30_30_10mm3_settings2_5PMTs-center-collimator-with-reflector-grease2-measurement-4_pmt_";
       TString name_pmt1 = Base_name;
       name_pmt1 +="1.txt";
       t_data_pmt1->ReadFile(name_pmt1,"HV:amplitude:charge");         
       TString name_pmt2 = Base_name;
       name_pmt2 +="2.txt";
       t_data_pmt2->ReadFile(name_pmt2,"HV:amplitude:charge");         
       TString name_pmt3 = Base_name;
       name_pmt3 +="3.txt";
       t_data_pmt3->ReadFile(name_pmt3,"HV:amplitude:charge");         
       TString name_pmt4 = Base_name;
       name_pmt4 +="4.txt";
       t_data_pmt4->ReadFile(name_pmt4,"HV:amplitude:charge");         
       TString name_pmt5 = Base_name;
       name_pmt5 +="5.txt";
       t_data_pmt5->ReadFile(name_pmt5,"HV:amplitude:charge");         
       TString name_pmt6 = Base_name;
       name_pmt6 +="6.txt";
       t_data_pmt6->ReadFile(name_pmt6,"HV:amplitude:charge");        
   
       float HV, amplitude, charge;
       t_data_pmt1 -> SetBranchAddress("charge",&charge);
       t_data_pmt2 -> SetBranchAddress("charge",&charge);
       t_data_pmt3 -> SetBranchAddress("charge",&charge);
       t_data_pmt4 -> SetBranchAddress("charge",&charge);
       t_data_pmt5 -> SetBranchAddress("charge",&charge);
       t_data_pmt6 -> SetBranchAddress("charge",&charge);

       double min_Charge_x = 75.;
       double max_Charge_x = 1105.;
       double bin_Charge_width = 2.;
       int n_Charge_bins = (int)((max_Charge_x-min_Charge_x)/bin_Charge_width);
       double min_Charge_x_all = 100.;
       double max_Charge_x_all = 4000.;
       double bin_Charge_width_all = 5.;
       int n_Charge_bins_all = (int)((max_Charge_x_all-min_Charge_x_all)/bin_Charge_width_all);

       TH1D* h_charge_data_pmt1 = new TH1D("h_charge_data_pmt1","PMT 1; charge [# photons];entries / 5 ph ",n_Charge_bins,min_Charge_x,max_Charge_x);
       h_charge_data_pmt1 -> SetLineColor(kRed);
       TH1D* h_charge_data_pmt2 = new TH1D("h_charge_data_pmt2","PMT 2; charge [# photons];entries / 5 ph ",n_Charge_bins,min_Charge_x,max_Charge_x);
       h_charge_data_pmt2 -> SetLineColor(kRed);
       TH1D* h_charge_data_pmt3 = new TH1D("h_charge_data_pmt3","PMT 3; charge [# photons];entries / 5 ph ",n_Charge_bins,min_Charge_x,max_Charge_x);
       h_charge_data_pmt3 -> SetLineColor(kRed);
       TH1D* h_charge_data_pmt4 = new TH1D("h_charge_data_pmt4","PMT 4; charge [# photons];entries / 5 ph ",n_Charge_bins,min_Charge_x,max_Charge_x);
       h_charge_data_pmt4 -> SetLineColor(kRed);
       TH1D* h_charge_data_pmt5 = new TH1D("h_charge_data_pmt5","PMT 5; charge [# photons];entries / 5 ph ",n_Charge_bins,min_Charge_x,max_Charge_x);
       h_charge_data_pmt5 -> SetLineColor(kRed);
       TH1D* h_charge_data_pmt6 = new TH1D("h_charge_data_pmt6","PMT 6; charge [# photons];entries / 5 ph ",n_Charge_bins,min_Charge_x,max_Charge_x);
       h_charge_data_pmt6 -> SetLineColor(kRed);

       TH1D* h_total_charge_data = new TH1D("h_total_charge_data","All PMTs data; charge [# photons];entries / 5 ph ",n_Charge_bins_all,min_Charge_x_all,max_Charge_x_all);
       h_total_charge_data -> SetLineColor(kRed);
       TH1D* h_total_charge_data_sidePMTs = new TH1D("h_total_charge_data_sidePMTs","Lateral PMTs data; charge [# photons];entries / 5 ph ",n_Charge_bins_all,min_Charge_x_all,max_Charge_x_all);
       h_total_charge_data_sidePMTs -> SetLineColor(kRed-5);

       int nEvents = t_data_pmt1->GetEntries();
       double charge_pmt1, charge_pmt2, charge_pmt3, charge_pmt4, charge_pmt5, charge_pmt6;
       for(int i = 0; i < nEvents; i++){

          t_data_pmt1->GetEntry(i);
          charge_pmt1 = charge;

          t_data_pmt2->GetEntry(i);
          charge_pmt2 = charge;

          t_data_pmt3->GetEntry(i);
          charge_pmt3 = charge;

          t_data_pmt4->GetEntry(i);
          charge_pmt4 = charge;

          t_data_pmt5->GetEntry(i);
          charge_pmt5 = charge;

          t_data_pmt6->GetEntry(i);
          charge_pmt6 = charge;
          
          h_charge_data_pmt1->Fill(charge_pmt1);
          h_charge_data_pmt2->Fill(charge_pmt2);
          h_charge_data_pmt3->Fill(charge_pmt3);
          h_charge_data_pmt4->Fill(charge_pmt4);
          h_charge_data_pmt5->Fill(charge_pmt5);
          h_charge_data_pmt6->Fill(charge_pmt6);
          
          double charge_all_pmts = charge_pmt1 + charge_pmt2 + charge_pmt3 + charge_pmt4 + charge_pmt6;
          double charge_side_pmts = charge_pmt2 + charge_pmt3 + charge_pmt4 + charge_pmt6;
        
          if(charge_all_pmts > min_Charge_x_all && charge_pmt5 > 1.){  
          	h_total_charge_data ->Fill(charge_all_pmts);
          }
          if(charge_side_pmts > min_Charge_x_all && charge_pmt5 > 1.){  
	        h_total_charge_data_sidePMTs ->Fill(charge_side_pmts);
          }
       } 
      
       //Fit distribution to estimate the scale factor in hte simulations
       TF1* f_gaus_lateral = new TF1("f_gaus_lateral","gaus",900.,1200.);
       h_total_charge_data_sidePMTs->Fit(f_gaus_lateral,"R","",750.,1100.);
       double peak_1MeV_lateral_data = f_gaus_lateral->GetParameter(1);
       
       TF1* f_gaus_all = new TF1("f_gaus_all","gaus",1400.,1700.);
       h_total_charge_data -> Fit(f_gaus_all,"R","",1300.,1700.);
       double peak_1MeV_all_data = f_gaus_all->GetParameter(1);
       
        TCanvas* c_charges_data = new TCanvas();
	// Number of PADS
        const Int_t Nx = 3;
        const Int_t Ny = 2;
	// Margins
        Float_t lMargin = 0.02;
        Float_t rMargin = 0.02;
        // Canvas setup
        c_charges_data->Divide(Nx,Ny,lMargin,rMargin);
        c_charges_data->cd(1);
        h_charge_data_pmt1->Draw(""); 	
	h_charge_data_pmt1->GetXaxis()->SetRangeUser(min_Charge_x,800.);
        c_charges_data->cd(2);
        h_charge_data_pmt2->Draw(""); 	
	h_charge_data_pmt2->GetXaxis()->SetRangeUser(min_Charge_x,800.);
        c_charges_data->cd(3);
        h_charge_data_pmt3->Draw(""); 	
	h_charge_data_pmt3->GetXaxis()->SetRangeUser(min_Charge_x,800.);
        c_charges_data->cd(4);
        h_charge_data_pmt4->Draw(""); 	
	h_charge_data_pmt4->GetXaxis()->SetRangeUser(min_Charge_x,800.);
        c_charges_data->cd(5);
        h_charge_data_pmt6->Draw(""); 	
	h_charge_data_pmt6->GetXaxis()->SetRangeUser(min_Charge_x,800.);
        c_charges_data->cd(6);
	h_total_charge_data_sidePMTs->Draw("");
        h_total_charge_data_sidePMTs->GetXaxis()->SetRangeUser(50.,2000.);
        h_total_charge_data->Draw("HIST SAME");




       //declare tree to read MC
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

       TH1D* h_charge_mc_pmt1 = new TH1D("h_total_pmt1","PMT mc 1; charge [# photons];entries / 5 ph ",n_Charge_bins,min_Charge_x,max_Charge_x);
       TH1D* h_charge_mc_pmt2 = new TH1D("h_total_pmt2","PMT mc 2; charge [# photons];entries / 5 ph ",n_Charge_bins,min_Charge_x,max_Charge_x);
       TH1D* h_charge_mc_pmt3 = new TH1D("h_total_pmt3","PMT mc 3; charge [# photons];entries / 5 ph ",n_Charge_bins,min_Charge_x,max_Charge_x);
       TH1D* h_charge_mc_pmt4 = new TH1D("h_total_pmt4","PMT mc 4; charge [# photons];entries / 5 ph ",n_Charge_bins,min_Charge_x,max_Charge_x);
       TH1D* h_charge_mc_pmt5 = new TH1D("h_total_pmt5","PMT mc 5; charge [# photons];entries / 5 ph ",n_Charge_bins,min_Charge_x,max_Charge_x);
       TH1D* h_total_charge_mc_bruto = new TH1D("h_total_charge_mc_bruto","All PMTs; charge [# photons];entries / 5 ph ",n_Charge_bins_all,min_Charge_x_all,max_Charge_x_all);
       TH1D* h_total_charge_mc = new TH1D("h_total_charge_mc","All PMTs mc; charge [# photons];entries / 5 ph ",n_Charge_bins_all,min_Charge_x_all,max_Charge_x_all);
       TH1D* h_total_charge_mc_sidePMTs = new TH1D("h_total_charge_mc_sidePMTs","Lateral PMTs mc; charge [# photons];entries / 5 ph ",n_Charge_bins_all,min_Charge_x_all,max_Charge_x_all);
       TH1D* h_total_charge_mc_sidePMTs_bruto = new TH1D("h_total_charge_mc_sidePMTs_bruto","Lateral PMTs; charge [# photons];entries / 5 ph ",n_Charge_bins_all,min_Charge_x_all,max_Charge_x_all);
       h_total_charge_mc_sidePMTs -> SetLineColor(kBlue-5);

       TH1D* h_ratio = new TH1D("h_ratio","h_ratio; Q1/(Q1+Q2);entries / bin ",150,0.,1.);
       TH2D* h_Edep_TotalCharge = new TH2D("h_Edep_TotalCharge","h_Edep_TotalCharge;Edep PEN4[MeV]; # detected photons",n_E_bins,min_E_x,max_E_x,n_Charge_bins,min_Charge_x_all,max_Charge_x_all);

       int nEntries = t_geant4_data->GetEntries();
       TRandom3 *myRandoms = new TRandom3 ();
       double pmt_resolution1 = 0.5;
       double pmt_resolution2 = 0.5;
       double pmt_resolution3 = 0.5;
       double pmt_resolution4 = 0.3;
       double pmt_resolution5 = 0.5;
       double pmt_resolution6 = 0.5;
       double pmt_resolution_bottom = 0.9;
       double E_res_ten = 0.1;
       double E_res_five = 0.05;
       
       for(int i = 0; i < nEntries; ++i){
	       t_geant4_data->GetEntry(i);
	       if(EdepTrigger>15.){
                    h_total_charge_mc_bruto->Fill(NPMT1+NPMT2+NPMT3+NPMT4+NPMT5);
                    h_total_charge_mc_sidePMTs_bruto->Fill(NPMT1+NPMT2+NPMT4+NPMT5);
               }
       }
       //Fit the histograms
       TCanvas* c_mc_all = new TCanvas();
       int bin_max_mc_total = h_total_charge_mc_bruto->GetMaximumBin();
       double x_bin_max_mc_total = h_total_charge_mc_bruto->GetXaxis()->GetBinCenter(bin_max_mc_total);
       TF1* f_gaus_total_mc = new TF1("f_gaus_total_mc","gaus",x_bin_max_mc_total-200.,x_bin_max_mc_total+200.);
       h_total_charge_mc_bruto->Fit(f_gaus_total_mc,"R","",x_bin_max_mc_total-200.,x_bin_max_mc_total+200.);
       //h_total_charge_mc_bruto->Fit(f_gaus_total_mc,"R","",3100,3400.);
       double peak_1MeV_total_mc = f_gaus_total_mc->GetParameter(1);

       TCanvas* c_mc_side = new TCanvas();
       int bin_max_mc_lateral = h_total_charge_mc_sidePMTs_bruto->GetMaximumBin();
       double x_bin_max_mc_lateral = h_total_charge_mc_sidePMTs_bruto->GetXaxis()->GetBinCenter(bin_max_mc_lateral);
       //TF1* f_gaus_lateral_mc = new TF1("f_gaus_lateral_mc","gaus",x_bin_max_mc_lateral-200.,x_bin_max_mc_lateral+200.);
       TF1* f_gaus_lateral_mc = new TF1("f_gaus_lateral_mc","gaus",1500.,2000.);
       h_total_charge_mc_sidePMTs_bruto->Fit(f_gaus_lateral_mc,"R","",x_bin_max_mc_lateral-200.,x_bin_max_mc_lateral+200.);
       //h_total_charge_mc_sidePMTs_bruto->Fit(f_gaus_lateral_mc,"R","",1800,2000.);
       double peak_1MeV_lateral_mc = f_gaus_lateral_mc->GetParameter(1);

       cout<<" data peak lateral "<<peak_1MeV_lateral_data<<" data peak all "<<peak_1MeV_all_data<<" data all/lateral "<<peak_1MeV_all_data/peak_1MeV_lateral_data<<endl;
       cout<<" mc peak lateral "<<peak_1MeV_lateral_mc<<" mc peak all "<<peak_1MeV_total_mc<<" mc all/lateral "<<peak_1MeV_total_mc/peak_1MeV_lateral_mc<<endl;       

       double scale_factor_lateral = peak_1MeV_lateral_mc/peak_1MeV_lateral_data;
       if(scale_factor_lateral < 0.2){ scale_factor_lateral = 1.0;}
       cout<<" scale_factor_lateral "<<scale_factor_lateral<<endl;
       double scale_factor_bottom = (peak_1MeV_total_mc-peak_1MeV_lateral_mc)/(peak_1MeV_all_data-peak_1MeV_lateral_data);
       if(scale_factor_bottom < 0.1) {scale_factor_bottom = 1.0;}
       cout<<" scale_factor_bottom "<<scale_factor_bottom<<endl;
      

 
       for(int i = 0; i < nEntries; ++i){
	       t_geant4_data->GetEntry(i);
	       //if(EdepTrigger>20. && EdepTrigger<100.){
	       if(EdepTrigger>15.){
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
		       double n_ph_pmt1 = myRandoms -> Gaus(NPMT1, pmt_resolution1 * sqrt (NPMT1))/scale_factor_lateral;
		       //double n_ph_pmt1 = myRandoms -> Poisson(NPMT1);
		       double n_ph_pmt2 = myRandoms -> Gaus(NPMT2, pmt_resolution2 * sqrt (NPMT2))/scale_factor_lateral;
		       //double n_ph_pmt2 = myRandoms -> Poisson(NPMT2);
		       double n_ph_pmt3 = myRandoms -> Gaus(NPMT3, pmt_resolution_bottom * sqrt (NPMT3))/scale_factor_bottom;
		       //double n_ph_pmt3 = myRandoms -> Poisson(NPMT3);
		       double n_ph_pmt4 = myRandoms -> Gaus(NPMT4, pmt_resolution4 * sqrt (NPMT4))/scale_factor_lateral;
		       //double n_ph_pmt4 = myRandoms -> Poisson(NPMT4);
		       double n_ph_pmt5 = myRandoms -> Gaus(NPMT5, pmt_resolution5 * sqrt (NPMT5))/scale_factor_lateral;
		       //double n_ph_pmt5 = myRandoms -> Poisson(NPMT5);
		       h_charge_mc_pmt1->Fill(n_ph_pmt1);
		       h_charge_mc_pmt2->Fill(n_ph_pmt2);
		       h_charge_mc_pmt3->Fill(n_ph_pmt3);
		       h_charge_mc_pmt4->Fill(n_ph_pmt4);
		       h_charge_mc_pmt5->Fill(n_ph_pmt5);
		       h_total_charge_mc->Fill(n_ph_pmt1+n_ph_pmt2+n_ph_pmt3+n_ph_pmt4+n_ph_pmt5);
		       h_total_charge_mc_sidePMTs->Fill(n_ph_pmt1+n_ph_pmt2+n_ph_pmt5+n_ph_pmt4);

	       }
       }
	TCanvas* c_edeps = new TCanvas();
	gStyle->SetOptStat(0);
	h_edep_PEN1->Draw("HIST PLC");
	h_edep_PEN3->Draw("SAME HIST PLC");
	h_edep_PEN4->Draw("SAME HIST PLC");
	h_edep_PEN2->Draw("SAME HIST PLC");
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
      
        TCanvas* c_mc_only = new TCanvas();
        c_mc_only->Divide(Nx,Ny,lMargin,rMargin);
        c_mc_only->cd(3);
        h_charge_mc_pmt1->DrawNormalized("");
        c_mc_only->cd(2);
        h_charge_mc_pmt2->DrawNormalized("");
        c_mc_only->cd(1);
        h_charge_mc_pmt3->DrawNormalized("");
        c_mc_only->cd(4);
        h_charge_mc_pmt4->DrawNormalized("");
        c_mc_only->cd(5);
        h_charge_mc_pmt5->DrawNormalized("");
        c_mc_only->cd(6);
	h_total_charge_mc_sidePMTs->DrawNormalized("");
        h_total_charge_mc->DrawNormalized("HIST SAME ");

 
        TCanvas* c_charges_mc = new TCanvas();
        // Canvas setup
        c_charges_mc->Divide(Nx,Ny,lMargin,rMargin);
        c_charges_mc->cd(3);
        h_charge_mc_pmt1->DrawNormalized(""); 	
        h_charge_data_pmt3->DrawNormalized("HISTSAME");
	//h_charge_mc_pmt1->GetXaxis()->SetRangeUser(5.,800.);
        c_charges_mc->cd(2);
        h_charge_mc_pmt2->DrawNormalized(""); 	
        h_charge_data_pmt2->DrawNormalized("HISTSAME");
	//h_charge_mc_pmt2->GetXaxis()->SetRangeUser(5.,800.);
        c_charges_mc->cd(1);
        h_charge_mc_pmt3->DrawNormalized(""); 	
        h_charge_data_pmt1->DrawNormalized("HISTSAME");
	//h_charge_mc_pmt3->GetXaxis()->SetRangeUser(5.,800.);
        c_charges_mc->cd(4);
        h_charge_mc_pmt4->DrawNormalized(""); 	
        h_charge_data_pmt4->DrawNormalized("HISTSAME");
	//h_charge_mc_pmt4->GetXaxis()->SetRangeUser(5.,800.);
        c_charges_mc->cd(5);
        h_charge_mc_pmt5->DrawNormalized(""); 	
	//h_charge_mc_pmt5->GetXaxis()->SetRangeUser(5.,800.);
        h_charge_data_pmt6->DrawNormalized("HISTSAME");

        c_charges_mc->cd(6);
	h_total_charge_mc_sidePMTs->DrawNormalized("");
	h_total_charge_data_sidePMTs->DrawNormalized("HIST SAME");
        h_total_charge_mc->DrawNormalized("HIST SAME ");
        h_total_charge_data->DrawNormalized("HIST SAME ");
        h_total_charge_data->GetXaxis()->SetRangeUser(50.,1800.);


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
