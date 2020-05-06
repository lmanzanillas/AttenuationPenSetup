//Macro to read geant4 simulated data
//Author: Luis Manzanillas

void macro_read_data_det2(TString stringDataName = "output/Bi207-17_04_2020-14-10-21", int nCores = 1){ 
	TNtuple calls("calls","calls","event:EdepTrigger:EdepPEN1:EdepPEN2:EdepPEN3:EdepPEN4:EdepOther:NTrigger:NPMT1:NPMT2:NPMT3:NPMT4:NPMT5");
	//calls.ReadFile("./output/Bi207-17_04_2020-10-12-58/PenAttenuationSetup_position_0.000000mm_nt_PEN_t0.csv");
        for(int i = 0; i < nCores; i++){
		TString totalName = "";
		totalName +=stringDataName;
		totalName +="/";
		totalName +="PenAttenuationSetup_nt_PEN_t";
		totalName +=i;
		totalName +=".csv";
 		cout<<" name "<<totalName<<endl;
		calls.ReadFile(totalName);
	
	}
	
	TH1D* h_edep_trigger = new TH1D("h_edep_trigger","h_edep_PEN_trigger;Edep [keV]; entries / keV ",200,1.,201.);
	TH1D* h_edep_PEN4 = new TH1D("h_edep_PEN4","h_edep_PEN4;Edep [keV]; entries / 5 keV ",400,10.,2010.);
	TH1D* h_edep_PEN_all = new TH1D("h_edep_PEN_all","h_edep_PEN_all;Edep [keV]; entries / 2.5 keV ",600,10.,2010.);
        TH1D* h_total_charge = new TH1D("h_total_charge","h_total_charge; charge [number of photons];entries / photon ",600,5.,3605.);
        TH1D* h_ratio = new TH1D("h_ratio","h_ratio; Q1/(Q1+Q2);entries / bin ",100,0.,1.);
	calls.Draw("EdepPEN4>>h_edep_PEN4");
        TCanvas* c_all = new TCanvas();
	calls.Draw("(EdepPEN1+EdepPEN2+EdepPEN3)>>h_edep_PEN_all","(NPMT1+NPMT2)>1.&&NPMT1>1&&NPMT2>1&&EdepTrigger>20&&EdepTrigger<100");
        TCanvas* c_trigger = new TCanvas();
        calls.Draw("EdepTrigger>>h_edep_trigger");
        TCanvas* c_charge = new TCanvas();
        calls.Draw("(NPMT1+NPMT2+NPMT3+NPMT4+NPMT5)>>h_total_charge","(NPMT1+NPMT2+NPMT3+NPMT4+NPMT5)>1.&&NPMT1>1&&NPMT2>1&&EdepTrigger>20&&EdepTrigger<100");
	TCanvas* c_ratio = new TCanvas();
        calls.Draw("NPMT1/(NPMT1+NPMT2)>>h_ratio","(NPMT1+NPMT2)>50.&&NPMT1>1&&NPMT2>1&&EdepTrigger>20&&EdepTrigger<100");
}
