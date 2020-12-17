//Macro to read geant4 simulated data
//Author: Luis Manzanillas

void macro_tune_mc(string stringDataName = "output/Bi207-17_04_2020-14-10-21", int nCores = 1, double x_min_lat = 0., double x_max_lat = 0., double x_min_all = 0., double x_max_all = 0.){
       //declare trees to read data
       TTree *t_data_pmt1 = new TTree("t_data_pmt1","data txt pmt1");
       TTree *t_data_pmt2 = new TTree("t_data_pmt2","data txt pmt2");
       TTree *t_data_pmt3 = new TTree("t_data_pmt3","data txt pmt3");
       TTree *t_data_pmt4 = new TTree("t_data_pmt4","data txt pmt4");
       TTree *t_data_pmt5 = new TTree("t_data_pmt5","data txt pmt5");
       TTree *t_data_pmt6 = new TTree("t_data_pmt6","data txt pmt6");


       TString Base_name = "/remote/ceph/group/gedet/data/pen/2020/2020-01-28_8182cf65_lt_6pmt_first_data/data_luis_settings2/text_data_for_fit/Spectrum_hv_optimal_setup_2_EJ200_SoLid_30_30_1mm3_settings2_5PMTs-center-collimator-without-reflector-tape-grease1-measurement-2_pmt_";
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

       double min_Charge_x = 5.;
       double max_Charge_x = 995.;
       double bin_Charge_width = 5.;
       int n_Charge_bins = (int)((max_Charge_x-min_Charge_x)/bin_Charge_width);
       double min_Charge_x_all = 20.;
       double max_Charge_x_all = 3000.;
       double bin_Charge_width_all = 10.;
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
      
       //Fit distribution to estimate the scale factor in the simulations
       TF1* f_gaus_lateral = new TF1("f_gaus_lateral","gaus",900.,1200.);
       int bin_max_data_side = h_total_charge_data_sidePMTs->GetMaximumBin();
       double x_bin_max_data_side = h_total_charge_data_sidePMTs->GetXaxis()->GetBinCenter(bin_max_data_side);
       h_total_charge_data_sidePMTs->Fit(f_gaus_lateral,"R","",x_bin_max_data_side-70.,x_bin_max_data_side+100.);
       double peak_1MeV_lateral_data = f_gaus_lateral->GetParameter(1);
       
       TF1* f_gaus_all = new TF1("f_gaus_all","gaus",1400.,1700.);
       int bin_max_data_all = h_total_charge_data->GetMaximumBin();
       double x_bin_max_data_all = h_total_charge_data->GetXaxis()->GetBinCenter(bin_max_data_all);
       h_total_charge_data -> Fit(f_gaus_all,"R","",x_bin_max_data_all-70.,x_bin_max_data_all+100.);
       double peak_1MeV_all_data = f_gaus_all->GetParameter(1);
      
       double xminFitData = peak_1MeV_lateral_data/4. - 60;
       double xmaxFitData = peak_1MeV_lateral_data/4. + 50;
       TF1* f_gaus_data = new TF1("f_gaus_data","gaus",0.,1500.);
       double mean_data_pmt1, mean_data_pmt2,mean_data_pmt3,mean_data_pmt4,mean_data_pmt6;
       h_charge_data_pmt1 -> Fit(f_gaus_data,"R","",peak_1MeV_all_data-peak_1MeV_lateral_data-50,peak_1MeV_all_data-peak_1MeV_lateral_data+50); 
       mean_data_pmt1 = f_gaus_data->GetParameter(1);
       h_charge_data_pmt2 -> Fit(f_gaus_data,"R","",xminFitData,xmaxFitData); 
       mean_data_pmt2 = f_gaus_data->GetParameter(1);
       h_charge_data_pmt3 -> Fit(f_gaus_data,"R","",xminFitData,xmaxFitData); 
       mean_data_pmt3 = f_gaus_data->GetParameter(1);
       h_charge_data_pmt4 -> Fit(f_gaus_data,"R","",xminFitData,xmaxFitData); 
       mean_data_pmt4 = f_gaus_data->GetParameter(1);
       h_charge_data_pmt6 -> Fit(f_gaus_data,"R","",xminFitData,xmaxFitData); 
       mean_data_pmt6 = f_gaus_data->GetParameter(1);

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
       h_charge_data_pmt1->GetXaxis()->SetRangeUser(min_Charge_x,990.);
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
       std::string str=stringDataName;
       std::size_t pos_base_name_ends = str.find("0.csv");
       std::string base_name = str.substr(0,pos_base_name_ends);
       cout<<" base_name "<<base_name<<endl;

       TTree *t_geant4_data = new TTree("t_geant4_data","data from ascii file");	
       for(int i = 0; i < nCores; i++){
		TString totalName = "";
		totalName +=base_name;
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
       double bin_E_width = 15.;
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

       TH1D* h_charge_mc_bruto_pmt1 = new TH1D("h_total_bruto_pmt1","PMT mc 1; charge [# photons];entries / 5 ph ",n_Charge_bins,min_Charge_x,max_Charge_x);
       TH1D* h_charge_mc_bruto_pmt2 = new TH1D("h_total_bruto_pmt2","PMT mc 2; charge [# photons];entries / 5 ph ",n_Charge_bins,min_Charge_x,max_Charge_x);
       TH1D* h_charge_mc_bruto_pmt3 = new TH1D("h_total_bruto_pmt3","PMT mc 3; charge [# photons];entries / 5 ph ",n_Charge_bins,min_Charge_x,max_Charge_x);
       TH1D* h_charge_mc_bruto_pmt4 = new TH1D("h_total_bruto_pmt4","PMT mc 4; charge [# photons];entries / 5 ph ",n_Charge_bins,min_Charge_x,max_Charge_x);
       TH1D* h_charge_mc_bruto_pmt5 = new TH1D("h_total_bruto_pmt5","PMT mc 5; charge [# photons];entries / 5 ph ",n_Charge_bins,min_Charge_x,max_Charge_x);
       TH1D* h_charge_mc_bruto_pmt6 = new TH1D("h_total_bruto_pmt6","PMT mc 6; charge [# photons];entries / 5 ph ",n_Charge_bins,min_Charge_x,max_Charge_x);

       int nScales = 10;
       TH1D* h_charge_mc_scale_pmt1[nScales];
       TH1D* h_charge_mc_scale_pmt2[nScales];
       TH1D* h_charge_mc_scale_pmt3[nScales];
       TH1D* h_charge_mc_scale_pmt4[nScales];
       TH1D* h_charge_mc_scale_pmt5[nScales];
       
       TString name_h_scale("");
       for(int scale = 0; scale < nScales; scale++){
             name_h_scale = "h_scale_";
             name_h_scale +=scale;
             name_h_scale +="_pmt1";
             h_charge_mc_scale_pmt1[scale] = new TH1D(name_h_scale,"PMT mc; charge [# photons];entries / 5 ph ",n_Charge_bins,min_Charge_x,max_Charge_x);
             name_h_scale = "h_scale_";
             name_h_scale +=scale;
             name_h_scale +="_pmt2";
             h_charge_mc_scale_pmt2[scale] = new TH1D(name_h_scale,"PMT mc; charge [# photons];entries / 5 ph ",n_Charge_bins,min_Charge_x,max_Charge_x);
             name_h_scale = "h_scale_";
             name_h_scale +=scale;
             name_h_scale +="_pmt3";
             h_charge_mc_scale_pmt3[scale] = new TH1D(name_h_scale,"PMT mc; charge [# photons];entries / 5 ph ",n_Charge_bins,min_Charge_x,max_Charge_x);
             name_h_scale = "h_scale_";
             name_h_scale +=scale;
             name_h_scale +="_pmt4";
             h_charge_mc_scale_pmt4[scale] = new TH1D(name_h_scale,"PMT mc; charge [# photons];entries / 5 ph ",n_Charge_bins,min_Charge_x,max_Charge_x);
             name_h_scale = "h_scale_";
             name_h_scale +=scale;
             name_h_scale +="_pmt5";
             h_charge_mc_scale_pmt5[scale] = new TH1D(name_h_scale,"PMT mc; charge [# photons];entries / 5 ph ",n_Charge_bins,min_Charge_x,max_Charge_x);
       }

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
       double pmt_resolution1 = 0.6;
       double pmt_resolution2 = 0.6;
       double pmt_resolution3 = 0.6;
       double pmt_resolution4 = 0.6;
       double pmt_resolution5 = 0.6;
       double pmt_resolution6 = 0.6;
       double pmt_resolution_bottom = 0.6;
       double E_res_ten = 0.1;
       double E_res_five = 0.05;
       
       for(int i = 0; i < nEntries; ++i){
	       t_geant4_data->GetEntry(i);
	       if(EdepTrigger>15.){
                    h_total_charge_mc_bruto->Fill(NPMT1+NPMT2+NPMT3+NPMT4+NPMT5);
                    h_total_charge_mc_sidePMTs_bruto->Fill(NPMT1+NPMT2+NPMT4+NPMT5);
		    h_charge_mc_bruto_pmt1->Fill(NPMT1);
		    h_charge_mc_bruto_pmt2->Fill(NPMT2);
		    h_charge_mc_bruto_pmt3->Fill(NPMT3);
		    h_charge_mc_bruto_pmt4->Fill(NPMT4);
		    h_charge_mc_bruto_pmt5->Fill(NPMT5);
               }
       }
       //Fit the histograms
       TCanvas* c_mc_all = new TCanvas();
       int bin_max_mc_total = h_total_charge_mc_bruto->GetMaximumBin();
       double x_bin_max_mc_total = h_total_charge_mc_bruto->GetXaxis()->GetBinCenter(bin_max_mc_total);
       TF1* f_gaus_total_mc = new TF1("f_gaus_total_mc","gaus",x_bin_max_mc_total-100.,x_bin_max_mc_total+200.);
       h_total_charge_mc_bruto->Fit(f_gaus_total_mc,"R","",x_bin_max_mc_total-100.,x_bin_max_mc_total+200.);
       //h_total_charge_mc_bruto->Fit(f_gaus_total_mc,"R","",x_min_all,x_max_all);
       double peak_1MeV_total_mc = f_gaus_total_mc->GetParameter(1);

       TCanvas* c_mc_side = new TCanvas();
       int bin_max_mc_lateral = h_total_charge_mc_sidePMTs_bruto->GetMaximumBin();
       double x_bin_max_mc_lateral = h_total_charge_mc_sidePMTs_bruto->GetXaxis()->GetBinCenter(bin_max_mc_lateral);
       //TF1* f_gaus_lateral_mc = new TF1("f_gaus_lateral_mc","gaus",x_bin_max_mc_lateral-200.,x_bin_max_mc_lateral+200.);
       TF1* f_gaus_lateral_mc = new TF1("f_gaus_lateral_mc","gaus",2500.,3000.);
       h_total_charge_mc_sidePMTs_bruto->Fit(f_gaus_lateral_mc,"R","",x_bin_max_mc_lateral-100.,x_bin_max_mc_lateral+100.);
       //h_total_charge_mc_sidePMTs_bruto->Fit(f_gaus_lateral_mc,"R","",x_min_lat,x_max_lat);
       double peak_1MeV_lateral_mc = f_gaus_lateral_mc->GetParameter(1);
       TCanvas* c_fits = new TCanvas();
       double xminFit = peak_1MeV_lateral_mc/4. -70.;
       double xmaxFit = peak_1MeV_lateral_mc/4. +50.;
       double mean_pmt1, mean_pmt2, mean_pmt3, mean_pmt4, mean_pmt5;
       TF1* f_gaus = new TF1("f_gaus","gaus",0.,1500.);
       h_charge_mc_bruto_pmt1->Fit(f_gaus,"R","",xminFit,xmaxFit);
       mean_pmt1=f_gaus->GetParameter(1);

       h_charge_mc_bruto_pmt2->Fit(f_gaus,"R","",xminFit,xmaxFit);
       mean_pmt2=f_gaus->GetParameter(1);

       h_charge_mc_bruto_pmt3->Fit(f_gaus,"R","",peak_1MeV_total_mc-peak_1MeV_lateral_mc - 50,peak_1MeV_total_mc-peak_1MeV_lateral_mc + 50);
       mean_pmt3=f_gaus->GetParameter(1);

       h_charge_mc_bruto_pmt4->Fit(f_gaus,"R","",xminFit,xmaxFit);
       mean_pmt4=f_gaus->GetParameter(1);

       h_charge_mc_bruto_pmt5->Fit(f_gaus,"R","",xminFit,xmaxFit);
       mean_pmt5=f_gaus->GetParameter(1);

       cout<<" data peak lateral "<<peak_1MeV_lateral_data<<" data peak all "<<peak_1MeV_all_data<<" data all/lateral "<<peak_1MeV_all_data/peak_1MeV_lateral_data<<endl;
       cout<<" mc peak lateral "<<peak_1MeV_lateral_mc<<" mc peak all "<<peak_1MeV_total_mc<<" mc all/lateral "<<peak_1MeV_total_mc/peak_1MeV_lateral_mc<<endl;       

       double scale_factor_lateral = peak_1MeV_lateral_mc/peak_1MeV_lateral_data;
       if(scale_factor_lateral < 0.5){ scale_factor_lateral = 1.0;}
       cout<<" scale_factor_lateral "<<scale_factor_lateral<<endl;
       double scale_factor_bottom = (peak_1MeV_total_mc-peak_1MeV_lateral_mc)/(peak_1MeV_all_data-peak_1MeV_lateral_data);
       if(scale_factor_bottom < 0.5) {scale_factor_bottom = 1.0;}
       cout<<" scale_factor_bottom "<<scale_factor_bottom<<endl;
       
       double scale_factor_pmt1 = mean_pmt1/mean_data_pmt3; 
       double scale_factor_pmt2 = mean_pmt2/mean_data_pmt2; 
       double scale_factor_pmt3 = mean_pmt3/mean_data_pmt1; 
       double scale_factor_pmt4 = mean_pmt4/mean_data_pmt4; 
       double scale_factor_pmt5 = mean_pmt5/mean_data_pmt6; 

 
       for(int i = 0; i < nEntries; ++i){
	       t_geant4_data->GetEntry(i);
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
		       double n_ph_pmt1 = myRandoms -> Gaus(NPMT1, pmt_resolution1 * sqrt (NPMT1))/scale_factor_pmt1;
		       double n_ph_pmt2 = myRandoms -> Gaus(NPMT2, pmt_resolution2 * sqrt (NPMT2))/scale_factor_pmt2;
		       double n_ph_pmt3 = myRandoms -> Gaus(NPMT3, pmt_resolution_bottom * sqrt (NPMT3))/scale_factor_pmt3;
		       double n_ph_pmt4 = myRandoms -> Gaus(NPMT4, pmt_resolution4 * sqrt (NPMT4))/scale_factor_pmt4;
		       double n_ph_pmt5 = myRandoms -> Gaus(NPMT5, pmt_resolution5 * sqrt (NPMT5))/scale_factor_pmt5;
		       double n_ph_pmt1_all = myRandoms -> Gaus(NPMT1, pmt_resolution1 * sqrt (NPMT1))/scale_factor_lateral;
		       double n_ph_pmt2_all = myRandoms -> Gaus(NPMT2, pmt_resolution2 * sqrt (NPMT2))/scale_factor_lateral;
		       double n_ph_pmt3_all = myRandoms -> Gaus(NPMT3, pmt_resolution_bottom * sqrt (NPMT3))/scale_factor_bottom;
		       double n_ph_pmt4_all = myRandoms -> Gaus(NPMT4, pmt_resolution4 * sqrt (NPMT4))/scale_factor_lateral;
		       double n_ph_pmt5_all = myRandoms -> Gaus(NPMT5, pmt_resolution5 * sqrt (NPMT5))/(scale_factor_lateral);
		       h_charge_mc_pmt1->Fill(n_ph_pmt1_all);
		       h_charge_mc_pmt2->Fill(n_ph_pmt2_all);
		       h_charge_mc_pmt3->Fill(n_ph_pmt3_all);
		       h_charge_mc_pmt4->Fill(n_ph_pmt4_all);
		       h_charge_mc_pmt5->Fill(n_ph_pmt5_all);
		       h_total_charge_mc->Fill(n_ph_pmt1_all + n_ph_pmt2_all + n_ph_pmt3_all + n_ph_pmt4_all + n_ph_pmt5_all);
		       h_total_charge_mc_sidePMTs->Fill(n_ph_pmt1_all + n_ph_pmt2_all + n_ph_pmt5_all + n_ph_pmt4_all);
	       }
       }
	TCanvas* c_edeps = new TCanvas();
	gStyle->SetOptStat(0);
	h_edep_PEN1->Draw("HIST PLC");
	//h_edep_PEN3->Draw("SAME HIST PLC");
	//h_edep_PEN4->Draw("SAME HIST PLC");
	h_edep_PEN2->Draw("SAME HIST PLC");
        TCanvas* c_all = new TCanvas();
        //h_edep_all_materials->Draw("HIST PLC");
	h_edep_all_PEN->Draw("HIST PLC");
	h_edep_PEN1->Draw("SAME HIST PLC");
	h_edep_PEN2->Draw("SAME HIST PLC");
	//h_edep_PEN_3_4->Draw("SAME HIST PLC");
      
        TCanvas* c_eres = new TCanvas();	
	h_edep_all_PEN_res_ten->SetLineColor(kRed);
	h_edep_all_PEN_res_five->Draw("");
	//h_edep_all_PEN_res_ten->Draw("HISTSAMES");
     

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
	h_total_charge_mc_sidePMTs->DrawNormalized("HISTE");
        h_total_charge_mc->DrawNormalized("HIST SAME ");

 
        TCanvas* c_charges_mc = new TCanvas();
	h_charge_mc_pmt1->Scale(1./h_charge_mc_pmt1->Integral());
	h_charge_mc_pmt2->Scale(1./h_charge_mc_pmt2->Integral());
	h_charge_mc_pmt3->Scale(1./h_charge_mc_pmt3->Integral());
	h_charge_mc_pmt4->Scale(1./h_charge_mc_pmt4->Integral());
	h_charge_mc_pmt5->Scale(1./h_charge_mc_pmt5->Integral());
	h_charge_data_pmt3->Scale(1./h_charge_data_pmt3->Integral());
	h_charge_data_pmt1->Scale(1./h_charge_data_pmt1->Integral());
	h_charge_data_pmt2->Scale(1./h_charge_data_pmt2->Integral());
	h_charge_data_pmt4->Scale(1./h_charge_data_pmt4->Integral());
	h_charge_data_pmt6->Scale(1./h_charge_data_pmt6->Integral());
	h_total_charge_mc_sidePMTs->Scale(1./h_total_charge_mc_sidePMTs->Integral());
	h_total_charge_data_sidePMTs->Scale(1./h_total_charge_data_sidePMTs->Integral());
        h_total_charge_mc->Scale(1./h_total_charge_mc->Integral());
        h_total_charge_data->Scale(1./h_total_charge_data->Integral());

        // Canvas setup
        c_charges_mc->Divide(Nx,Ny,lMargin,rMargin);
        c_charges_mc->cd(3);
        h_charge_mc_pmt1->Draw("HISTE"); 	
        h_charge_data_pmt3->Draw("HISTESAME");
	//h_charge_mc_pmt1->GetXaxis()->SetRangeUser(5.,800.);
        c_charges_mc->cd(2);
        h_charge_mc_pmt2->Draw("HISTE"); 	
        h_charge_data_pmt2->Draw("HISTESAME");
	//h_charge_mc_pmt2->GetXaxis()->SetRangeUser(5.,800.);
        c_charges_mc->cd(1);
        h_charge_mc_pmt3->Draw("HISTE"); 	
        h_charge_data_pmt1->Draw("HISTESAME");
	//h_charge_mc_pmt3->GetXaxis()->SetRangeUser(5.,800.);
        c_charges_mc->cd(4);
        h_charge_mc_pmt4->Draw("HISTE"); 	
        h_charge_data_pmt4->Draw("HISTESAME");
	//h_charge_mc_pmt4->GetXaxis()->SetRangeUser(5.,800.);
        c_charges_mc->cd(5);
        h_charge_mc_pmt5->Draw("HISTE"); 	
	//h_charge_mc_pmt5->GetXaxis()->SetRangeUser(5.,800.);
        h_charge_data_pmt6->DrawNormalized("HISTESAME");

        c_charges_mc->cd(6);
	h_total_charge_mc_sidePMTs->Draw("HISTE");
	h_total_charge_data_sidePMTs->Draw("HISTE SAME");
        h_total_charge_mc->Draw("HISTE SAME ");
        h_total_charge_data->Draw("HISTE SAME ");
        h_total_charge_data->GetXaxis()->SetRangeUser(50.,1800.);

        //Perform a chi2 test to validate best set of paramters
        double chi2_pmt1 = 0;
        double chi2_pmt2 = 0;
        double chi2_pmt3 = 0;
        double chi2_pmt4 = 0;
        double chi2_pmt5 = 0;
        double chi2_lateral = 0;
        double chi2_all = 0;
        int dof_pmt1 = 0;
        int dof_pmt2 = 0;
        int dof_pmt3 = 0;
        int dof_pmt4 = 0;
        int dof_pmt5 = 0;
        int dof_lateral = 0;
        int dof_total = 0;
        double n_Events_data_pmt1 = (double)h_charge_data_pmt1->Integral();
        double n_Events_data_pmt2 = (double)h_charge_data_pmt2->Integral();
        double n_Events_data_pmt3 = (double)h_charge_data_pmt3->Integral();
        double n_Events_data_pmt4 = (double)h_charge_data_pmt4->Integral();
        double n_Events_data_pmt6 = (double)h_charge_data_pmt6->Integral();
        double n_Events_mc_pmt1 = (double)h_charge_mc_pmt1->Integral();
        double n_Events_mc_pmt2 = (double)h_charge_mc_pmt2->Integral();
        double n_Events_mc_pmt3 = (double)h_charge_mc_pmt3->Integral();
        double n_Events_mc_pmt4 = (double)h_charge_mc_pmt4->Integral();
        double n_Events_mc_pmt5 = (double)h_charge_mc_pmt5->Integral();
        
        double n_Events_mc_lateral = (double)h_total_charge_mc_sidePMTs->Integral();
        double n_Events_data_lateral = (double)h_total_charge_data_sidePMTs->Integral();
        double n_Events_mc_all = (double)h_total_charge_mc->Integral();
        double n_Events_data_all = (double)h_total_charge_data->Integral();

	double ks_lateral = h_total_charge_mc_sidePMTs -> KolmogorovTest(h_total_charge_data_sidePMTs);
	double ks_all = h_total_charge_mc -> KolmogorovTest(h_total_charge_data);
        
        double ks_pmt1 = h_charge_data_pmt1 -> KolmogorovTest(h_charge_mc_pmt3);
        double ks_pmt2 = h_charge_data_pmt2 -> KolmogorovTest(h_charge_mc_pmt2);
        double ks_pmt3 = h_charge_data_pmt3 -> KolmogorovTest(h_charge_mc_pmt1);
        double ks_pmt4 = h_charge_data_pmt4 -> KolmogorovTest(h_charge_mc_pmt4);
        double ks_pmt6 = h_charge_data_pmt6 -> KolmogorovTest(h_charge_mc_pmt5);

        for(int i = 0; i < n_Charge_bins; i++ ){
		    int bin_i = i + 1;
                    //bottom pmt
		    if(h_charge_mc_pmt3->GetBinContent(bin_i) > 0. || h_charge_data_pmt1->GetBinContent(bin_i) > 0.){chi2_pmt1 += pow(h_charge_data_pmt1->GetBinContent(bin_i)/n_Events_data_pmt1 - h_charge_mc_pmt3->GetBinContent(bin_i)/n_Events_mc_pmt3,2)/(h_charge_mc_pmt3->GetBinContent(bin_i)/pow(n_Events_mc_pmt3,2) + h_charge_data_pmt1->GetBinContent(bin_i)/pow(n_Events_data_pmt1,2)); dof_pmt1++;}
		    if(h_charge_mc_pmt2->GetBinContent(bin_i) > 0. || h_charge_data_pmt2->GetBinContent(bin_i) > 0.){chi2_pmt2 += pow(h_charge_data_pmt2->GetBinContent(bin_i)/n_Events_data_pmt2 - h_charge_mc_pmt2->GetBinContent(bin_i)/n_Events_mc_pmt2,2)/(h_charge_mc_pmt2->GetBinContent(bin_i)/pow(n_Events_mc_pmt2,2) + h_charge_data_pmt2->GetBinContent(bin_i)/pow(n_Events_data_pmt2,2)); dof_pmt2++;}
		    if(h_charge_mc_pmt1->GetBinContent(bin_i) > 0. || h_charge_data_pmt3->GetBinContent(bin_i) > 0.){chi2_pmt3 += pow(h_charge_data_pmt3->GetBinContent(bin_i)/n_Events_data_pmt3 - h_charge_mc_pmt1->GetBinContent(bin_i)/n_Events_mc_pmt1,2)/(h_charge_mc_pmt1->GetBinContent(bin_i)/pow(n_Events_mc_pmt1,2) + h_charge_data_pmt3->GetBinContent(bin_i)/pow(n_Events_data_pmt3,2)); dof_pmt3++;}
		    if(h_charge_mc_pmt4->GetBinContent(bin_i) > 0. || h_charge_data_pmt4->GetBinContent(bin_i) > 0.){chi2_pmt4 += pow(h_charge_data_pmt4->GetBinContent(bin_i)/n_Events_data_pmt4 - h_charge_mc_pmt4->GetBinContent(bin_i)/n_Events_mc_pmt4,2)/(h_charge_mc_pmt4->GetBinContent(bin_i)/pow(n_Events_mc_pmt4,2) + h_charge_data_pmt4->GetBinContent(bin_i)/pow(n_Events_data_pmt4,2)); dof_pmt4++;}
		    if(h_charge_mc_pmt5->GetBinContent(bin_i) > 0. || h_charge_data_pmt6->GetBinContent(bin_i) > 0.){chi2_pmt5 += pow(h_charge_data_pmt6->GetBinContent(bin_i)/n_Events_data_pmt6 - h_charge_mc_pmt5->GetBinContent(bin_i)/n_Events_mc_pmt5,2)/(h_charge_mc_pmt5->GetBinContent(bin_i)/pow(n_Events_mc_pmt5,2) + h_charge_data_pmt6->GetBinContent(bin_i)/pow(n_Events_data_pmt6,2)); dof_pmt5++;}
                    //cout<<" chi2_pmt1 "<<chi2_pmt1<<endl;
	}
        for(int i = 0; i < n_Charge_bins_all; i++){
		    int bin_i = i + 1;
                    if( h_total_charge_mc_sidePMTs ->GetBinContent(bin_i) > 0. || h_total_charge_data_sidePMTs ->GetBinContent(bin_i) > 0.){chi2_lateral += pow(h_total_charge_data_sidePMTs ->GetBinContent(bin_i) - h_total_charge_mc_sidePMTs ->GetBinContent(bin_i),2)/(h_total_charge_mc_sidePMTs ->GetBinContent(bin_i) + h_total_charge_data_sidePMTs ->GetBinContent(bin_i)); dof_lateral++;}
                    if( h_total_charge_mc ->GetBinContent(bin_i) > 0. || h_total_charge_data ->GetBinContent(bin_i) > 0.){chi2_all += pow(h_total_charge_data ->GetBinContent(bin_i) - h_total_charge_mc ->GetBinContent(bin_i),2)/(h_total_charge_mc ->GetBinContent(bin_i) + h_total_charge_data ->GetBinContent(bin_i)); dof_total++;}
        }
	cout<<" pow(3.,2) = "<<pow(3.,2)<<" chi lat "<<chi2_lateral<<" chi all "<<chi2_all<<" h_total_charge_data->Integral() "<<h_total_charge_data->Integral()<<endl;
        	
	std::size_t pos_start = str.find("_x_");
	std::size_t pos_ends = str.find("_mm_");
	std::string pos_x = str.substr(pos_start+3,pos_ends-pos_start-3);
	cout<<"pos_x: "<<pos_x<<endl;


	//Bi207_25_09_2020-01-16-16_x_0_mm_LY_9920ph_MeV_ABS_1.50_Det_4_2_samples_1.75mm_PVT_structure/PenAttenuationSetup_SigAlphaSides_0.15_SigAlphaBottom_0.15_PmtReflecSide_0.17_PmtReflecBottom_0.07_nt_PEN_t24 
	std::size_t pos_abs_start = str.find("_ABS_");
	std::size_t pos_abs_ends = str.find("_Det_");
	std::string abs_file = str.substr(pos_abs_start+5,pos_abs_ends-pos_abs_start-5);
        cout<<" pos_abs_start "<<pos_abs_start<<" pos_abs_ends "<<pos_abs_ends<<" abs_file "<<abs_file<<endl;

	std::size_t pos_ly_start = str.find("_LY_");
	std::size_t pos_ly_ends = str.find("ph_MeV");
	std::string ly_file = str.substr(pos_ly_start+4,pos_ly_ends-pos_ly_start-4);

	std::size_t pos_SigmaAlphaSides_start = str.find("_SigAlphaBottom_");
	std::string SigmaAlpha_Sides = str.substr(pos_SigmaAlphaSides_start-4,4);
        cout<<" SigmaAlpha_Sides "<<SigmaAlpha_Sides<<endl;

	std::size_t pos_SigmaAlphaBottom_start = str.find("_PmtReflecSide_");
	std::string SigmaAlpha_Bottom = str.substr(pos_SigmaAlphaBottom_start-4,4);
        cout<<" SigmaAlpha_Bottom "<<SigmaAlpha_Bottom<<endl;

	std::size_t pos_PMTSides_start = str.find("_PmtReflecBottom_");
	std::string str_PMT_sides = str.substr(pos_PMTSides_start-4,4);
        cout<<" str_PMT_sides "<<str_PMT_sides<<endl;

	std::size_t pos_PMTBottom_start = str.find("_nt_");
	std::string str_PMT_Bottom = str.substr(pos_PMTBottom_start-4,4);
        cout<<" str_PMT_Bottom "<<str_PMT_Bottom<<endl;

	TString nameFileLYtxt ="mc_Tuning_results_pvt_2_samples_17_mm.txt";
	ofstream myfileSaveQAData(nameFileLYtxt,std::ios_base::app);
	if (myfileSaveQAData.is_open())
	{
		myfileSaveQAData<<abs_file<<" "<<ly_file<<" "<<SigmaAlpha_Sides<<" "<<SigmaAlpha_Bottom<<" "<<str_PMT_sides<<" "<<str_PMT_Bottom<<" "<<chi2_pmt1<<" "<<dof_pmt1<<" "<<chi2_pmt2<<" "<<dof_pmt2<<" "<<chi2_pmt3<<" "<<dof_pmt3<<" "<<chi2_pmt4<<" "<<dof_pmt4<<" "<<chi2_pmt5<<" "<<dof_pmt5<<" "<<chi2_lateral<<" "<<dof_lateral<<" "<<chi2_all<<" "<<dof_total<<" "<<scale_factor_lateral<<" "<<scale_factor_bottom<<" "<<ks_lateral<<" "<<ks_all<<" "<<ks_pmt1<<" "<<ks_pmt2<<" "<<ks_pmt3<<" "<<ks_pmt4<<" "<<ks_pmt6<<"\n";
	        myfileSaveQAData.close();
	}
	else cout << "Unable to open file";

}
