//Macro to read geant4 simulated data
//Author: Luis Manzanillas

void macro_read_data_det3(string stringDataName = "output/Bi207-17_04_2020-14-10-21", int nCores = 1){
       //declare tree to read data
       TTree *t_geant4_data = new TTree("t_geant4_data","data from ascii file");	
       for(int i = 0; i < nCores; i++){
		TString totalName = "";
		totalName +=stringDataName;
		totalName +=i;
		totalName +=".csv";
 		cout<<" name "<<totalName<<endl;
		if(i == 0){
			t_geant4_data->ReadFile(totalName,"wl_photon:distance:flag");
		}else{
			t_geant4_data->ReadFile(totalName);
		}	
       }
       float wl_photon, distance, flag;
       t_geant4_data -> SetBranchAddress("wl_photon",&wl_photon);
       t_geant4_data -> SetBranchAddress("distance",&distance);
       t_geant4_data -> SetBranchAddress("flag",&flag);

       double min_wl = 380.;       
       double max_wl = 680.;
       double bin_E_width = 1.;
       int n_E_bins = (int)((max_wl-min_wl)/bin_E_width);

              
       TH1D* h_emitted_photons = new TH1D("h_emitted_photons","Emitted photons; photon wl [nm]; entries / 1 nm ",n_E_bins,min_wl,max_wl);
       TH1D* h_detected_photons = new TH1D("h_detected_photons","Detected photons; photon wl [nm]; entries / 1 nm ",n_E_bins,min_wl,max_wl);
       TH1D* h_absorbed_photons = new TH1D("h_absorbed_photons","Absorbed photons; photon wl [nm]; entries / 1 nm ",n_E_bins,min_wl,max_wl);

       int nEntries = t_geant4_data->GetEntries();
       for(int i = 0; i < nEntries; ++i){
	       t_geant4_data->GetEntry(i);
               float wl_nm = wl_photon/1000000.;   
	       if(int(flag) ==0 ){ h_emitted_photons->Fill(wl_nm);}
	       if(int(flag) ==1 ){ h_detected_photons->Fill(wl_nm);}
	       if(int(flag) ==2 ){ h_absorbed_photons->Fill(wl_nm);}
       }
     
         
        TCanvas* c_stacked = new TCanvas();	
	h_emitted_photons->SetLineColor(kRed);
	h_detected_photons->SetLineColor(kViolet);
	h_emitted_photons->Draw("");
	h_detected_photons->Draw("HISTSAMES");
	h_absorbed_photons->Draw("HISTSAMES");

}
