{
    
    TF1 * C1 = new TF1("C1","[0]*exp(-x/[1])",0,20000);
    C1->SetParameters(191.1,22.);

    TF1 * C2 = new TF1("C2","[0]*exp(-x/[1])",0,20000);
    C2->SetParameters(229.5,74.);

    TF1 * C3 = new TF1("C3","[0]*exp(-x/[1])",0,20000);
    C3->SetParameters(88.3,208.);

    TF1 * C4 = new TF1("C4","[0]*exp(-x/[1])",0,20000);
    C4->SetParameters(497.,881.);

    TF1 * C5 = new TF1("C5","[0]*exp(-x/[1])",0,20000);
    C5->SetParameters(25.1,4300.);

    TF1 * C6 = new TF1("C6","[0]*exp(-x/[1])",0,20000);
    C6->SetParameters(5.9,18100.);

    TF1 * C7 = new TF1("C7","[0]*exp(-x/[1])",0,20000);
    C7->SetParameters(1.2,87700);

    TF1 * Ctot_div = new TF1("Ctot_div","(C1+C2+C3+C4+C5+C6+C7)",0,20000);

    TF1 * Ctot = new TF1("Ctot","(C1+C2+C3+C4+C5+C6+C7)/[0]",0,20000);
    Ctot->SetParameter(0,double(Ctot_div->Eval(0)));
    
    Ctot->Draw();
    
    TFile * f1 = new TFile("decay.root","RECREATE");
    Ctot->Write();
    f1->Close();
    
}