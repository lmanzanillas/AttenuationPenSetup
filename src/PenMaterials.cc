/**
 * @file PenMaterials.cc
 * @brief Define PenMaterials
 * @author: (modified by) Luis Manzanillas
 * @date 2020 Pen - Max Planck Institut fur Physik
 */


#include "PenMaterials.hh"

#include "G4Material.hh"
#include "G4PhysicalConstants.hh"
#include "G4SystemOfUnits.hh"
#include "G4LogicalSkinSurface.hh"
#include "G4NistManager.hh"
#include "G4MTRunManager.hh"

using namespace std;
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
PenMaterials::PenMaterials()
{
    lightYieldAntracene=20000/MeV; //Anthracene
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
PenMaterials::~PenMaterials()
{
}
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void PenMaterials::Construct()
{
    G4NistManager* nistManager = G4NistManager::Instance();

    // ------------------------------------------------------------------------
    // Elements
    // ------------------------------------------------------------------------
//    G4Element* Ge = new G4Element("Germanium", "Ge", 32., 72.64*g/mole);
    G4Element* Al = new G4Element("Aluminium", "Al", 13., 26.98*g/mole);
    G4Element* Pb = new G4Element("Lead", "Pb", 82., 207.2*g/mole);
    G4Element* O = new G4Element("Oxygen", "O", 8., 16.00*g/mole);
    G4Element* N = new G4Element("Nitrogen", "N", 7., 14.00674*g/mole);
    G4Element* Si = new G4Element("Silicon",  "Si", 14., 28.09*g/mole);
    G4Element* C = new G4Element("Carbon", "C", 6., 12.011*g/mole);
    G4Element* C_graphite = new G4Element("TS_C_of_Graphite", "C_graphite", 6., 12.011*g/mole);
    G4Element* H = new G4Element("Hydrogen","H", 1., 1.00794*g/mole);
    G4Element* H_water = new G4Element("TS_H_of_Water","H_water", 1., 1.00794*g/mole);
    G4Element* H_polyethylene = new G4Element("TS_H_of_Polyethylene", "H_poly", 1., 1.00794*g/mole);
//    G4Element* Cl = new G4Element("Chlorine", "Cl", 17., 35.45*g/mole);
    G4Element* F = new G4Element("Fluorine", "F", 9., 18.998*g/mole);
    G4Element* Zn = new G4Element("Zinc","Zn", 30., 65.38*g/mole);
    G4Element* S = new G4Element("Sulphur","S", 16., 32.06*g/mole);
    G4Element* Fe = new G4Element("Iron", "Fe", 26., 55.847*g/mole);

    G4Element* P = new G4Element("Phosphorus", "P", 15., 30.974*g/mole);
    G4Element* Ti = new G4Element("Titanium", "Ti", 22., 47.867*g/mole);
    G4Element* Cr = new G4Element("Chromium", "Cr", 24., 51.996*g/mole);
    G4Element* Cu = new G4Element("Copper", "Cu", 29., 63.546*g/mole);
//    G4Element* Ni = new G4Element("Nickel", "Ni", 28., 58.693*g/mole);
//    G4Element* Mo = new G4Element("Molybdenum", "Mo", 42., 95.94*g/mole);

    G4Element* Mg = new G4Element("Magnesium", "Mg", 12., 24.305*g/mole);
    G4Element* Mn = new G4Element("Manganese", "Mn", 25., 54.938*g/mole);

    G4Element* Sb = new G4Element("Antimony", "Sb", 51., 121.760*g/mole);
    G4Element* Cs = new G4Element("Cesium", "Cs", 55., 132.905*g/mole);
    G4Element* K = new G4Element("Potassium", "K", 19., 39.098*g/mole);
//    G4Element* B = new G4Element("Boron", "B", 5., 10.811*g/mole);


    // ------------------------------------------------------------------------
    // Isotopes
    // ------------------------------------------------------------------------
    G4Isotope* Li6 = new G4Isotope("Li6", 3, 6, 6.*g/mole);
    G4Isotope* Li7 = new G4Isotope("Li7", 3, 7, 7.*g/mole);

    G4Isotope* B10 = new G4Isotope("B10", 5, 10, 10.*g/mole);
    G4Isotope* B11 = new G4Isotope("B11", 5, 11, 11.*g/mole);
    

    // ------------------------------------------------------------------------
    // Enriched Li6 - 95% by mass
    // ------------------------------------------------------------------------
    G4Element *Li6el = new G4Element("Li6el","Li", 2);
    Li6el->AddIsotope(Li6, 95.*perCent);
    Li6el->AddIsotope(Li7, 5.*perCent);
    
    // ------------------------------------------------------------------------
    // Pure Li6
    // ------------------------------------------------------------------------
    G4Element *Li6pure = new G4Element("Li6pure","Li", 1);
    Li6pure->AddIsotope(Li6, 100.0*perCent);

    // ------------------------------------------------------------------------
    // Enriched B10 - 95% by mass
    // ------------------------------------------------------------------------
//    G4Element *B10el = new G4Element("B10el","B", 2);
//    B10el->AddIsotope(B10, 95.*perCent);
//    B10el->AddIsotope(B11, 5.*perCent);

    // ------------------------------------------------------------------------
    // Enriched B10 - 95% by mass
    // ------------------------------------------------------------------------
    G4Element *B10el = new G4Element("B10el","B", 2);
    B10el->AddIsotope(B10, 95.*perCent);
    B10el->AddIsotope(B11, 5.*perCent);

    
    // ------------------------------------------------------------------------
    // PenMaterials
    // ------------------------------------------------------------------------
    G4double density;
    G4int nel;
    
    // temperature of experimental hall is controlled at 20 degree.
    const G4double expTemp = STP_Temperature+20.*kelvin;
    
    density = universe_mean_density;
    
    // ------------------------------------------------------------------------
    // Vacuum *************************************************************
    // ------------------------------------------------------------------------
    
    G4Material* Vacuum = new G4Material ("Vacuum",
                                         1.0,
                                         1.01*g/mole,
                                         universe_mean_density,
                                         kStateGas,
                                         3.e-18*pascal,
                                         2.73*kelvin);
    G4int vacEntries=0;
    G4double vacEmit[500];
    G4double vacIndex[500];
    G4double vacAbsorb[500];
    G4double vacEnergy[500];
    G4double vacAbsorbconst=100*m;
    G4double pWavelength;
    
    std::ifstream ReadVac;

    G4String Vac="../properties/PVTEmission.dat";
    ReadVac.open(Vac);
    if(ReadVac.is_open()){
        while(!ReadVac.eof()){
            G4String filler;
            ReadVac>>pWavelength>>filler>>vacEmit[vacEntries];
            vacEnergy[vacEntries]=(1240/pWavelength)*eV; //convert wavelength to eV
            vacIndex[vacEntries]=1.0;
            vacAbsorb[vacEntries]=vacAbsorbconst;
            vacEntries++;
        }
    }
    else
        G4cout<<"Error opening file: "<<Vac<<G4endl;
    ReadVac.close();
    
    G4MaterialPropertiesTable* vacMPT=new G4MaterialPropertiesTable();
    vacMPT->AddProperty("RINDEX",vacEnergy,vacIndex,vacEntries);
    vacMPT->AddProperty("ABSLENGTH",vacEnergy,vacAbsorb,vacEntries);
    Vacuum->SetMaterialPropertiesTable(vacMPT);
    
    // ------------------------------------------------------------------------
    // Air
    // ------------------------------------------------------------------------
    
    density = 1.2929e-03 *g/cm3;  // at 20 degree
    
    G4Material* Air = new G4Material("Air", density, 2, kStateGas, expTemp);
    
    
    Air-> AddElement(N,  0.8);
    Air-> AddElement(O,  0.2);
    
    G4int airEntries = 18;
    G4int air_ref_index_Entries = 0;
    G4double air_ref_index_Energy[18];
    G4double air_ref_index_value[18];
    G4double airbulkAbsorb[18];
    
    
    std::ifstream Read_air_ref_index;
    G4String air_ref_index = "../properties/air_ref_index.dat";
    Read_air_ref_index.open(air_ref_index);
    if(Read_air_ref_index.is_open()){
        while(!Read_air_ref_index.eof()){
            G4String filler;
            Read_air_ref_index >> pWavelength >> filler >> air_ref_index_value[air_ref_index_Entries];
            air_ref_index_Energy[air_ref_index_Entries] = (1240/pWavelength)*eV;
	    airbulkAbsorb[air_ref_index_Entries] = 1.0e6*mm;
            air_ref_index_Entries++;
	    if(air_ref_index_Entries > (airEntries-1)){break;}
        }
    }
    else
        G4cout << "Error opening file: " << air_ref_index << G4endl;
    Read_air_ref_index.close();
    

    
    G4MaterialPropertiesTable* airMPT = new G4MaterialPropertiesTable();
    airMPT->AddProperty("RINDEX", air_ref_index_Energy, air_ref_index_value, airEntries);
    airMPT->AddProperty("ABSLENGTH", air_ref_index_Energy, airbulkAbsorb, airEntries);
    Air->SetMaterialPropertiesTable(airMPT);

//
// DEFINE WATER IN MULTIPLE WAYS
//
    G4int ncomponents, natoms;

    // NIST Water
    nistManager->FindOrBuildMaterial("G4_WATER");
   
   	// Water mixture
    density = 1.000 *g/cm3;
    G4Material* water_mixture =
    		new G4Material("Water_mixture", density, ncomponents=2,
    						kStateLiquid, expTemp,1*atmosphere);
    water_mixture->AddElement(H, natoms=2);
    water_mixture->AddElement(O, natoms=1);
    water_mixture->GetIonisation()->SetMeanExcitationEnergy(78.0*eV); // ICRU 73, revised

   	// Water structure
    G4Material* water_structure =
    		new G4Material("Water_structure", density, ncomponents=2,
    						kStateLiquid, expTemp,1*atmosphere);
    water_structure->AddElement(H_water, natoms=2);
    water_structure->AddElement(O, natoms=1);
    water_structure->GetIonisation()->SetMeanExcitationEnergy(78.0*eV); // ICRU 73, revised


    // Boric Acid
    density = 1.435 *g/cm3;
    G4Material* boric_acid =
            new G4Material("Boric_Acid", density, ncomponents=3,
                            kStateLiquid, expTemp,1*atmosphere);
    boric_acid->AddElement(H, natoms=3);
    boric_acid->AddElement(O, natoms=3);
//    boric_acid->AddElement(B, natoms=1);
    boric_acid->AddElement(B10el, natoms=1);
    boric_acid->GetIonisation()->SetMeanExcitationEnergy(87.5*eV); // ESTAR calculated mean excitation energy

    // Borated Water
    G4Material* water_borated =
	    new G4Material("water_borated", 1.02 *g/cm3, ncomponents=2,
    						kStateLiquid, expTemp,1*atmosphere);
    water_borated->AddMaterial(water_mixture, 97.*perCent);
    water_borated->AddMaterial (boric_acid , 3.*perCent);
//    water_borated->AddElement (B10el , 3.*perCent);
    water_borated->GetIonisation()->SetMeanExcitationEnergy(75.6*eV); // ESTAR calculated mean excitation energy


    //
    // DEFINE GRAPHITE IN MULTIPLE WAYS
    //

    //NIST Graphite
    nistManager->FindOrBuildMaterial("G4_GRAPHITE");

    // Graphite mixture
    density = 2.27*g/cm3;
    G4Material* graphite_mixture =
    		new G4Material("Graphite_mixture", density, ncomponents=1,
    						kStateSolid, expTemp, 1*atmosphere);
    graphite_mixture->AddElement(C, natoms=1);
    graphite_mixture->GetIonisation()->SetMeanExcitationEnergy(78.0*eV); // NIST material database

    // Graphite structure
    G4Material* graphite_structure =
    		new G4Material("Graphite_structure", density, ncomponents=1,
    						kStateSolid, expTemp, 1*atmosphere);
    graphite_structure->AddElement(C_graphite, natoms=1);
    graphite_structure->GetIonisation()->SetMeanExcitationEnergy(78.0*eV); // NIST material database

    
    // ------------------------------------------------------------------------
    // Epoxy (binders)
    // ------------------------------------------------------------------------
    
    
    G4Material * Epoxy = new G4Material("Epoxy",1*g/cm3,4,kStateSolid);
    Epoxy->AddElement(C,2);
    Epoxy->AddElement(H,6);
    Epoxy->AddElement(O,1);
    Epoxy->AddElement(Si,1);
    
    G4double ref_index_Energy[500];
    G4double ref_index_value[500];
    ;
    
    for (int i = 0; i < 500; i++){
        ref_index_Energy[i] = 0;
        ref_index_value[i] = 0;
    }

    
    // Read scintillator refractive index
    G4int ref_index_Entries = 0;
    std::ifstream Read_ref_index;
    //  Read_ref_index;
    G4String ref_index_emit = "../properties/epoxy_ref_index.dat";
    Read_ref_index.open(ref_index_emit);
    if(Read_ref_index.is_open()){
        while(!Read_ref_index.eof()){
            G4String filler;
            Read_ref_index >> pWavelength >> filler >> ref_index_value[ref_index_Entries];
            ref_index_Energy[ref_index_Entries] = (1240/pWavelength)*eV;
            ref_index_Entries++;
        }
    }
    else
        G4cout << "Error opening file: "<< ref_index_emit << G4endl;
    Read_ref_index.close();
    
//    G4MaterialPropertiesTable * EpoxyMPT = new G4MaterialPropertiesTable();
//    EpoxyMPT->AddProperty("RINDEX",ref_index_Energy,ref_index_value,ref_index_Entries);
//    Epoxy->SetMaterialPropertiesTable(EpoxyMPT);


    
    //lead for shielding
    G4Material* lead = new G4Material("Lead", 11.340*g/cm3, nel=1);
    lead->AddElement(Pb,1);
    
    //Aluminium for Structure
    G4Material* aluminium = new G4Material("Aluminium", 2.70*g/cm3, 1);
    aluminium->AddElement(Al,1);
    
    //Carbon Fibre between layers
//    G4Material* Cf = new G4Material("CarbonFiber", 0.145*g/cm3, nel=1);
//    Cf->AddElement(C,1);


    
    //HDPE for screen and reflector bars
    G4Material* HDPE = new G4Material("HDPE",0.97*g/cm3,nel=2);
    HDPE->AddElement(C,1);
    HDPE->AddElement(H,2);

    // Borated HDPE
    G4Material* HDPE_borated = new G4Material("HDPE_borated", 1.01*g/cm3, 2);
    HDPE_borated->AddMaterial(HDPE, 95.*perCent);
    HDPE_borated->AddElement(B10el, 5.*perCent);

    //PPL for shielding
    G4Material* PPL = new G4Material("PPL",0.855*g/cm3,nel=2);
    PPL->AddElement(C,3);
    PPL->AddElement(H,6);
    
    //Iron for Table
    G4Material* iron = new G4Material("Iron", 8.05*g/cm3, 1);
    iron->AddElement(Fe,1);

    //Graphite
//    G4Material* Graphite = new G4Material("Graphite", 2.14*g/cm3, 1);
//    Graphite->AddElement(C, 1);
    
    // ------------------------------------------------------------------------
    // Lithium Fluoride: ZnS - 50% Binder
    // ------------------------------------------------------------------------
    
    //Litium Fluoride
    G4Material* LiF = new G4Material("LiF", 2.64*g/cm3, 2);
    LiF->AddElement(Li6el,1);
    LiF->AddElement(F,  1);
    
   
    
    //Zinc Sulphide
    G4Material* ZnS = new G4Material("ZnS", 4.09*g/cm3, 2);
    ZnS->AddElement(Zn, 1);
    ZnS->AddElement(S,  1);

    G4Material* Li6ZnS = new G4Material("Li6ZnS", 2.229*g/cm3, 3);
    Li6ZnS->AddMaterial(LiF,26.14*perCent);
    Li6ZnS->AddMaterial(ZnS,52.28*perCent);
    Li6ZnS->AddMaterial(Epoxy,21.58*perCent);

    // G4cout << *(G4Material::GetMaterialTable()); // print the list of materials

    // Pure Li-6
    G4Material* Li6ZnS_Li6pure = new G4Material("Li6ZnS_Li6pure", 0.133*g/cm3, 1);
    Li6ZnS_Li6pure->AddElement(Li6pure,100.0*perCent);
    
    G4int absorbEntries = 0;
    G4double varabsorblength;
    G4double absorbEnergy[500];
    G4double Absorb[500];
    G4double scintEnergyZnS[500];
    G4double scintEmitZnS[500];
    G4double scintReflect[500];
    
    // G4double wlsabsorblength;
    for (int i = 0; i < 500; i++){
        scintEnergyZnS[i] = 0;
        scintEmitZnS[i] = 0;
        ref_index_Energy[i] = 0;
        ref_index_value[i] = 0;
    }
    
    // Read primary emission spectrum
    G4int scintEntries = 0;
    G4String Scint_file="../properties/lif_zns_emission_spk2.txt";
    std::ifstream ReadScint;
    ReadScint.open(Scint_file);
    if(ReadScint.is_open()){
        while(!ReadScint.eof()){
            G4String filler;
            ReadScint>>pWavelength>>filler>>scintEmitZnS[scintEntries];
            scintEnergyZnS[scintEntries] = (1240/pWavelength)*eV;         //convert wavelength to eV
            // scintSlow[scintEntries] = 0.0;  //arbitrary test value
            scintEntries++;
        }
    }
    else
        G4cout << "Error opening file: " << Scint_file << G4endl;
    ReadScint.close();

    // Read primary bulk absorption
    
    absorbEntries = 0;
    G4String Readabsorblength = "../properties/LiFZnS2.cfg";
    
    std::ifstream Readabsorb;
    Readabsorb.open(Readabsorblength);
    if (Readabsorb.is_open()){
        while(!Readabsorb.eof()){
            G4String filler;
            Readabsorb >> pWavelength >> filler >> varabsorblength;
            absorbEnergy[absorbEntries] = (1240/pWavelength)*eV;
            Absorb[absorbEntries] = varabsorblength * m;
            absorbEntries++;

        }
    }
    else
        G4cout << "----- Error opening file: "<< Readabsorblength << G4endl;
    Readabsorb.close();

    
    // Read scintillator refractive index
    
    ref_index_Entries = 0;
    
    //  Read_ref_index;
    G4String ref_index_emit2 = "../properties/epoxy_ref_index.dat";
      
    // if (!Read_ref_index2)
     std::ifstream  Read_ref_index2;
    
    Read_ref_index2.open("../properties/epoxy_ref_index.dat");
    if(Read_ref_index2.is_open()){
     while(!Read_ref_index2.eof()){
       G4String filler;
       Read_ref_index2 >> pWavelength  >> filler >> ref_index_value[ref_index_Entries];
	 ref_index_Energy[ref_index_Entries] = (1240/pWavelength)*eV;
	 scintReflect[ref_index_Entries] = 1.0;

	 ref_index_Entries++;
	   }
    }
        else
        G4cout << "Error opening file: "<< ref_index_emit << G4endl;
     Read_ref_index2.close();
    

    
    /// Now apply the properties table
    G4MaterialPropertiesTable* scintMPT2 = new G4MaterialPropertiesTable();
    scintMPT2->AddProperty("RINDEX",ref_index_Energy,ref_index_value,ref_index_Entries);
    scintMPT2->AddProperty("ABSLENGTH",absorbEnergy,Absorb,absorbEntries);
    scintMPT2->AddProperty("FASTCOMPONENT",scintEnergyZnS,scintEmitZnS,scintEntries);
    scintMPT2->AddProperty("SLOWCOMPONENT",scintEnergyZnS,scintEmitZnS,scintEntries);
    //G4double efficiency = 1.0;
    //scintMPT2->AddConstProperty("EFFICIENCY",efficiency);
    scintMPT2->AddConstProperty("SCINTILLATIONYIELD",75000/MeV);
    G4double scintRes = 1;
    scintMPT2->AddConstProperty("RESOLUTIONSCALE",scintRes);
    G4double scintSlowconst=1*ms;
    scintMPT2->AddConstProperty("FASTTIMECONSTANT",scintSlowconst);
    scintMPT2->AddConstProperty("SLOWTIMECONSTANT",scintSlowconst);
    scintMPT2->AddConstProperty("YIELDRATIO",1.0);
    scintMPT2->AddProperty("REFLECTIVITY",ref_index_Energy,scintReflect,ref_index_Entries);

    
    Li6ZnS->SetMaterialPropertiesTable(scintMPT2);
    
    // Set the Birks Constant for the neutron scintillator
    //quenching effect for heavy particle
    Li6ZnS->GetIonisation()->SetBirksConstant(0.008*mm/MeV);
    

    // ------------------------------------------------------------------------
    // LiF:ZnS backing => MELINEX 339, aka Polyethylene terephthalate
    // http://www.tekra.com/products/films/polyester-films/polyester-pet/melinex-339
    // https://en.wikipedia.org/wiki/Polyethylene_terephthalate
    // ------------------------------------------------------------------------
    G4Material* MELINEX = new G4Material("MELINEX", 1.38*g/cm3, nel=3);
    MELINEX->AddElement(C,10);
    MELINEX->AddElement(H,8);
    MELINEX->AddElement(O,4);



    // ------------------------------------------------------------------------
    // PVT Scintillator
    // ------------------------------------------------------------------------
    

    // NIST PVT
    G4Material* pvt_nist = nistManager->FindOrBuildMaterial("G4_PLASTIC_SC_VINYLTOLUENE");
    
    // PVT mixture
    density=1.023*g/cm3; // From EJ-200 data sheet
    G4Material* pvt_mixture = new G4Material("PVT_mixture",density,ncomponents=2,
    								kStateSolid, expTemp,1*atmosphere);
    pvt_mixture->AddElement(C,natoms=27);
    pvt_mixture->AddElement(H,natoms=30);
    pvt_mixture->GetIonisation()->SetMeanExcitationEnergy(64.7*eV); // NIST G4_PLASTIC_SC_VINYLTOLUENE
    
    // PVT structure
    G4Material* pvt_structure = new G4Material("PVT_structure",density,ncomponents=2,
    								kStateSolid, expTemp,1*atmosphere);
    pvt_structure->AddElement(C,natoms=27);
    pvt_structure->AddElement(H_polyethylene,natoms=30);
    pvt_structure->GetIonisation()->SetMeanExcitationEnergy(64.7*eV); // NIST G4_PLASTIC_SC_VINYLTOLUENE


    G4double scintEnergyPVT[115];
    G4double scintEmitPVT[115];
    G4int scintEntries_pvt = 115;
    // Read primary emission spectrum
    Scint_file ="../properties/PVTEmission.dat";
    std::ifstream ReadPVTScint;
    ReadPVTScint.open(Scint_file);
    scintEntries = 0;
    if(ReadPVTScint.is_open()){
        while(!ReadPVTScint.eof()){
	 
            G4String filler;
            ReadPVTScint>>pWavelength>>filler>>scintEmitPVT[scintEntries];
            scintEnergyPVT[scintEntries] = (1240./pWavelength)*eV;         //convert wavelength to eV
            scintEntries++;
	    if (scintEntries > (scintEntries_pvt-1)){ G4cout << " ERROR < scint entries  out of range > " << G4endl; break;}
        }
    }
    else
        G4cout << "Error opening file: " << Scint_file << G4endl;
    ReadPVTScint.close();
     
    // Read primary bulk absorption
    G4int abs_entries_pvt = 500;
    G4double absorbEnergy_pvt[500];
    G4double Absorb_pvt[500];
    //G4double Rayleigh_pvt[500];
    Readabsorblength = "../properties/PVTAbsorption.dat";
    
    Readabsorb.open(Readabsorblength);
    absorbEntries = 0;
    if (Readabsorb.is_open()){
        while(!Readabsorb.eof()){
            G4String filler;
            Readabsorb >> pWavelength >> filler >> varabsorblength;
            absorbEnergy_pvt[absorbEntries] = (1240./pWavelength)*eV;
            Absorb_pvt[absorbEntries] = varabsorblength * m;
            //Rayleigh_pvt[absorbEntries] = varabsorblength/5. * m;
            absorbEntries++;
	    if(absorbEntries > (abs_entries_pvt-1)){G4cout << " ERROR < entries abs  out of range > " << G4endl; break;}
        }
    }
    else
        G4cout<<"Error opening file: "<<Readabsorblength<<G4endl;
    Readabsorb.close();

    // Read scintillator refractive index
    G4int entries_pvt_rindex = 11;
    G4double ref_index_Energy_pvt[11];
    G4double ref_index_value_pvt[11];

    std::ifstream  Read_ref_index_pvt;
    ref_index_emit = "../properties/PVTRefIndex.dat";
    Read_ref_index_pvt.open(ref_index_emit);
    ref_index_Entries = 0;
    if(Read_ref_index_pvt.is_open()){
        while(!Read_ref_index_pvt.eof()){
            G4String filler;
            Read_ref_index_pvt >> pWavelength >> filler >> ref_index_value_pvt[ref_index_Entries];
            ref_index_Energy_pvt[ref_index_Entries] = (1240./pWavelength)*eV;
	    //G4cout<<ref_index_Entries<<" rindex "<<ref_index_value_pvt[ref_index_Entries]<<" energy "<<ref_index_Energy_pvt[ref_index_Entries]<<G4endl;
            ref_index_Entries++;
	    if(ref_index_Entries > (entries_pvt_rindex-1)){G4cout << " ERROR < entries ref abs  out of range > " << G4endl;break;}
        }
    }
    else
        G4cout << "Error opening file: "<< ref_index_emit << G4endl;
    Read_ref_index_pvt.close();
    
    // Now apply the properties table
    G4MaterialPropertiesTable* scintMPT = new G4MaterialPropertiesTable();
    scintMPT->AddProperty("RINDEX", ref_index_Energy_pvt, ref_index_value_pvt, entries_pvt_rindex);
    scintMPT->AddProperty("ABSLENGTH", absorbEnergy_pvt, Absorb_pvt, abs_entries_pvt);
    //scintMPT->AddProperty("RAYLEIGH", absorbEnergy_pvt, Rayleigh_pvt, abs_entries_pvt);
    scintMPT->AddProperty("FASTCOMPONENT", scintEnergyPVT, scintEmitPVT, scintEntries_pvt);
    scintMPT->AddProperty("SLOWCOMPONENT", scintEnergyPVT, scintEmitPVT, scintEntries_pvt);
    //efficiency = 1.0;
    //scintMPT->AddConstProperty("EFFICIENCY",efficiency);
    G4double LY_PVT = lightYieldAntracene*0.64;
    scintMPT->AddConstProperty("SCINTILLATIONYIELD",LY_PVT);
    scintRes=1.0;
    scintMPT->AddConstProperty("RESOLUTIONSCALE",scintRes);
    G4double scintFastconst=2.0*ns;  //fluorescence
    scintMPT->AddConstProperty("FASTTIMECONSTANT",scintFastconst);
    scintSlowconst=16.8*ns;  //phosphorescence
    scintMPT->AddConstProperty("SLOWTIMECONSTANT",scintSlowconst);
    scintMPT->AddConstProperty("YIELDRATIO",1.0); //was 1.0
   
    G4cout<<" pvt table ok "<<G4endl; 
    pvt_nist->SetMaterialPropertiesTable(scintMPT);
    pvt_mixture->SetMaterialPropertiesTable(scintMPT);
    pvt_structure->SetMaterialPropertiesTable(scintMPT);
    
    // Set the Birks Constant for the Polystyrene scintillator
    //quenching effect for heavy particle
    pvt_nist->GetIonisation()->SetBirksConstant(0.15*mm/MeV);
    pvt_mixture->GetIonisation()->SetBirksConstant(0.15*mm/MeV);
    pvt_structure->GetIonisation()->SetBirksConstant(0.15*mm/MeV);
    

    // ------------------------------------------------------------------------
    // Tyvek Wrapping
    // ------------------------------------------------------------------------
    
    G4Material* Tyvek = new G4Material("Tyvek", 0.38*g/cm3, 2);
    Tyvek->AddElement(C,2);
    Tyvek->AddElement(H,4);
    
    G4int teflon_entries = 0;
    G4double teflon_energy[500];
    G4double teflon_reflect[500];
    G4double zero[500];
    
    std::ifstream Read_teflon;
    G4String teflon_file = "../properties/teflon.dat";
    
    Read_teflon.open(teflon_file);
    if (Read_teflon.is_open()){
        while(!Read_teflon.eof()){
            G4String filler;
            G4double wavelength;
            G4double teflon_ref_coeff;
            Read_teflon >> wavelength >> filler >> teflon_ref_coeff;
            teflon_energy[teflon_entries] = (1240/wavelength)*eV;
            teflon_reflect[teflon_entries] = teflon_ref_coeff;
            zero[teflon_entries] = 1e-6;
            teflon_entries++;
        }
    }
    else
        G4cout << "Error opening file: " << teflon_file << G4endl;
    Read_teflon.close();
    
    G4MaterialPropertiesTable *mptTyvek;
    mptTyvek = new G4MaterialPropertiesTable();
    
    mptTyvek->AddProperty("REFLECTIVITY",teflon_energy,teflon_reflect,teflon_entries);
    mptTyvek->AddProperty("SPECULARLOBECONSTANT",teflon_energy,zero,teflon_entries);
    mptTyvek->AddProperty("SPECULARSPIKECONSTANT",teflon_energy,zero,teflon_entries);
    mptTyvek->AddProperty("BACKSCATTERCONSTANT",teflon_energy,zero,teflon_entries);
    
    Tyvek->SetMaterialPropertiesTable(mptTyvek);
    

    // ------------------------------------------------------------------------
    // WLS Fibre
    // ------------------------------------------------------------------------
    
    //Polystyrene WLS fibre core
    G4Material* Polystyrene = new G4Material("Polystyrene",1.05*g/cm3,2,kStateSolid,273.15*kelvin,1.0*atmosphere );
    Polystyrene->AddElement( H, 1 ); // Was 0.498
    Polystyrene->AddElement( C, 1 ); // Was 0.502
    
    //PMMA  WLS fibre cladding: PMMA (polymethylmethacrylate, C5H8O2) 
    G4Material* PMMA = new G4Material("PMMA",1.19*g/cm3,3,kStateSolid,273.15*kelvin,1.0*atmosphere );
    PMMA->AddElement( H, 8 ); // Was 0.532
    PMMA->AddElement( C, 5 ); // Was 0.336
    PMMA->AddElement( O, 2 ); // was 0.0132
    
    //Fiber WLS Absorption  ************************************
    
    G4double coreIndex[500];
    G4double cladIndex1[500];
    // G4double cladIndex2[500];
    G4double coreIndexconst = 1.59;
    G4double cladIndexconst1 = 1.49;
    // G4double cladIndexconst2 = 1.42;
    G4double wls_fiber_Energy[500];
    G4double wls_fiber_Absorb[500];
    G4int wls_fiber_Entries = 0;
    
    G4double wls_fiber_absorblength;
    
    
    std::ifstream ReadWLSfiberabsorb;
    G4String fname="../properties/BCF-91AAbsLength.dat";
    ReadWLSfiberabsorb.open(fname);
    if (ReadWLSfiberabsorb.is_open()){
        while(!ReadWLSfiberabsorb.eof()){
            G4String filler;
            ReadWLSfiberabsorb >> pWavelength >> filler >> wls_fiber_absorblength;
            wls_fiber_Energy[wls_fiber_Entries] = (1240/pWavelength)*eV;
            wls_fiber_Absorb[wls_fiber_Entries] = 1.0*wls_fiber_absorblength*m;
            
            coreIndex[wls_fiber_Entries] = coreIndexconst;
            cladIndex1[wls_fiber_Entries] = cladIndexconst1;
            // cladIndex2[wls_fiber_Entries] = cladIndexconst2;
            
            wls_fiber_Entries++;
        }
    }
    
    else
        G4cout << "Error opening file: " <<"../properties/BCF-91AAbsLength.dat" << G4endl;
    ReadWLSfiberabsorb.close();
    
    //Fiber WLS Emission  ****************************************
    G4double wls_fiber_emit_Energy[500];
    G4double wls_fiber_Emit[500];
    G4int wls_fiber_emit_Entries = 0;
    
    std::ifstream ReadWLSfibemit;
    //  G4String WLSfibemit = "../properties/WLSemit.cfg";
    
    fname="../properties/BCF-91AEmit.dat";
    ReadWLSfibemit.open(fname);
    if(ReadWLSfibemit.is_open()){
        while(!ReadWLSfibemit.eof()){
            G4String filler;
            ;
            ReadWLSfibemit >> pWavelength >> filler >> wls_fiber_Emit[wls_fiber_emit_Entries];
            wls_fiber_emit_Energy[wls_fiber_emit_Entries] = (1240/pWavelength)*eV;
            wls_fiber_emit_Entries++;
        }
    }
    else
        G4cout << "Error opening file: " << "../properties/BCF-91AEmit.dat" << G4endl;
    ReadWLSfibemit.close();
    
    
    //Fiber Bulk Absorption ****************************************
    G4int bulk_fiber_Entries = 0;
    G4double bulk_fiber_absorblength;
    G4double bulk_fiber_Absorb[500];
    G4double bulk_fiber_Energy[500];
    std::ifstream ReadFiberBulk;
    
    //  G4String FiberBulk = "../properties/fiberPSTabsorb.dat";
    
    fname="../properties/fiberPSTabsorb.dat";
    ReadFiberBulk.open(fname);
    if(ReadFiberBulk.is_open()){
        while(!ReadFiberBulk.eof()){
            G4String filler;
            ;
            ReadFiberBulk >> pWavelength >> filler >> bulk_fiber_absorblength;
            bulk_fiber_Energy[bulk_fiber_Entries] = (1240/pWavelength)*eV;
            bulk_fiber_Absorb[bulk_fiber_Entries] = 1.0*bulk_fiber_absorblength*m;
            bulk_fiber_Entries++;
        }
    }
    else
        G4cout << "Error opening file: " << "../properties/fiberPSTabsorb.dat" << G4endl;
    ReadFiberBulk.close();
    
    
    //PMMA Fiber Bulk Absorption  ******************************************
    G4double pmma_fiber_bulkAbsorb[500];
    G4double pmma_fiber_Energy[500];
    G4double pmma_fiber_absorblength;
    G4int pmma_fiber_Entries = 0;
    
    std::ifstream Read_pmma_fib_Bulk;
    //  G4String pmma_Bulk = "../properties/PMMABulkAbsorb.dat";
    
    fname="../properties/PMMABulkAbsorb.dat";
    Read_pmma_fib_Bulk.open(fname);
    if(Read_pmma_fib_Bulk.is_open()){
        while(!Read_pmma_fib_Bulk.eof()){
            G4String filler;
            ;
            Read_pmma_fib_Bulk >> pWavelength >> filler >> pmma_fiber_absorblength;
            pmma_fiber_Energy[pmma_fiber_Entries] = (1240/pWavelength)*eV;
            pmma_fiber_bulkAbsorb[pmma_fiber_Entries] = 1.0*pmma_fiber_absorblength*m;
            pmma_fiber_Entries++;
        }
    }
    else
        G4cout << "Error opening file: " << "../properties/PMMABulkAbsorb.dat" << G4endl;
    Read_pmma_fib_Bulk.close();
    
    
    G4MaterialPropertiesTable* fiberwlsMPT = new G4MaterialPropertiesTable();
    fiberwlsMPT->AddProperty("ABSLENGTH",bulk_fiber_Energy,bulk_fiber_Absorb,bulk_fiber_Entries);
    fiberwlsMPT->AddProperty("RINDEX",wls_fiber_Energy,coreIndex,wls_fiber_Entries);
    fiberwlsMPT->AddProperty("WLSABSLENGTH",wls_fiber_Energy,wls_fiber_Absorb,wls_fiber_Entries);
    fiberwlsMPT->AddProperty("WLSCOMPONENT",wls_fiber_emit_Energy,wls_fiber_Emit,wls_fiber_emit_Entries);
    fiberwlsMPT->AddConstProperty("WLSTIMECONSTANT",12*ns);
    
    G4MaterialPropertiesTable* innercladMPT = new G4MaterialPropertiesTable();
    innercladMPT->AddProperty("ABSLENGTH",pmma_fiber_Energy,pmma_fiber_bulkAbsorb,pmma_fiber_Entries);
    innercladMPT->AddProperty("RINDEX",wls_fiber_Energy,cladIndex1,wls_fiber_Entries);
    
    Polystyrene->SetMaterialPropertiesTable(fiberwlsMPT);
    PMMA->SetMaterialPropertiesTable(innercladMPT);
    
    
    // ------------------------------------------------------------------------
    // MPPC
    // ------------------------------------------------------------------------
    
    G4Material* MPPCFilm = new G4Material("MPPCFilm",1.05*g/cm3,1,kStateSolid,273.15*kelvin,1.0*atmosphere);
    MPPCFilm->AddElement(Si,1);
    
   
    // ------------------------------------------------------------------------
    // Mirror
    // ------------------------------------------------------------------------
    
    G4Material* Mirror = new G4Material("Mirror", 1.00*g/cm3,1);
    Mirror->AddElement(Al,1);
    // top mirror optical properties
    // NOTE: mirror_surface is attached to the surface not the material!
    G4double mirror_PP[2] = { 2.0*eV, 4.0*eV };
    //G4double mirror_REFL[2]  = { 1.0,     1.0 };
    G4double mirror_RIND[2]  = { 100,     100 };
    G4MaterialPropertiesTable *mirror_surface = new G4MaterialPropertiesTable();
    //mirror_surface->AddProperty("REFLECTIVITY", mirror_PP, mirror_REFL, 2);
    mirror_surface->AddProperty("RINDEX", mirror_PP, mirror_RIND, 2);
    Mirror->SetMaterialPropertiesTable(mirror_surface);
    
    // ------------------------------------------------------------------------
    // Grease
    // ------------------------------------------------------------------------
    G4Material* Grease = new G4Material("Grease",
                            1.05*g/cm3,
                            2,
                            kStateSolid,
                            273.15*kelvin,
                            1.0*atmosphere );
    
    Grease->AddElement( H, 0.498 );
    Grease->AddElement( C, 0.502 );
    
    G4int GreaseEntries=0;
    G4double GreaseEnergy[500];
    G4double Greaseabsorblength;
    G4double GreasebulkAbsorb[500];
    G4double GreaseIndex[500];
    G4double GreaseIndexconst=1.59;
    
    for (int i = 0; i < 500; i++){
        GreaseEnergy[i]=0;
        GreasebulkAbsorb[i]=0;
        GreaseIndex[i]=0;
    }
    
    GreaseEntries=0;
    std::ifstream ReadGreaseBulk;
    G4String GreaseBulk="../properties/GreaseBulkAbsorb.cfg";
    ReadGreaseBulk.open(GreaseBulk);
    if(ReadGreaseBulk.is_open()){
        while(!ReadGreaseBulk.eof()){
            G4String filler;
            ;
            ReadGreaseBulk>>pWavelength>>filler>>Greaseabsorblength;
            GreaseEnergy[GreaseEntries]=(1240/pWavelength)*eV;
            GreasebulkAbsorb[GreaseEntries]=Greaseabsorblength*m;
            GreaseIndex[GreaseEntries]=GreaseIndexconst;
            GreaseEntries++;
        }
    }
    else
        G4cout<<"Error opening file: "<<GreaseBulk<<G4endl;
    
    ReadGreaseBulk.close();
    G4MaterialPropertiesTable *GreaseMPT=new G4MaterialPropertiesTable();
    GreaseMPT->AddProperty("RINDEX",GreaseEnergy,GreaseIndex,GreaseEntries);
    GreaseMPT->AddProperty("ABSLENGTH",GreaseEnergy,GreasebulkAbsorb,GreaseEntries);
    Grease->SetMaterialPropertiesTable(GreaseMPT);


    //
    // DEFINE Phase1 Container MATERIALS
    //

    //Steel 080M40 (EN8)
    density = 7.87*g/cm3;
    G4Material* steel_EN8 = new G4Material("Steel_EN8", density, ncomponents=6,
                                              kStateSolid, expTemp, 1*atmosphere);
    steel_EN8->AddElement(C,   0.400*perCent);
    steel_EN8->AddElement(Si,  0.200*perCent);
    steel_EN8->AddElement(P,   0.035*perCent);
    steel_EN8->AddElement(S,   0.035*perCent);
    steel_EN8->AddElement(Mn,  0.800*perCent);
    steel_EN8->AddElement(Fe, 98.530*perCent);
    steel_EN8->GetIonisation()->SetMeanExcitationEnergy(320.7*eV); // ESTAR calculated mean excitation energy

    //Aluminium 6061 (EN 573-3)
    density = 2.72*g/cm3;
    G4Material* aluminium_6061 = new G4Material("Aluminium_6061", density, ncomponents=9,
                                              kStateSolid, expTemp, 1*atmosphere);
    aluminium_6061->AddElement(Mg,  1.000*perCent);
    aluminium_6061->AddElement(Al, 97.305*perCent);
    aluminium_6061->AddElement(Si,  0.600*perCent);
    aluminium_6061->AddElement(Ti,  0.075*perCent);
    aluminium_6061->AddElement(Cr,  0.195*perCent);
    aluminium_6061->AddElement(Mn,  0.075*perCent);
    aluminium_6061->AddElement(Fe,  0.350*perCent);
    aluminium_6061->AddElement(Cu,  0.275*perCent);
    aluminium_6061->AddElement(Zn,  0.125*perCent);
    aluminium_6061->GetIonisation()->SetMeanExcitationEnergy(188.6*eV); // ESTAR calculated mean excitation energy

    //CELOTEX foam (Polyisocyanurate) -->  C3H3N3O3
    density = 0.03*g/cm3;
    G4Material* celotex = new G4Material("CELOTEX", density, ncomponents=4,
                                              kStateSolid, expTemp, 1*atmosphere);
    celotex->AddElement(H, natoms=3);
    celotex->AddElement(C, natoms=3);
    celotex->AddElement(N, natoms=3);
    celotex->AddElement(O, natoms=3);
    celotex->GetIonisation()->SetMeanExcitationEnergy(46.1*eV); // ESTAR calculated mean excitation energy



    //
    // DEFINE BR2 MATERIALS
    //


    //NIST Concrete, for BR2 floor and ceil
    nistManager->FindOrBuildMaterial("G4_CONCRETE");

    //NIST STAINLESS-STEEL
    nistManager->FindOrBuildMaterial("G4_STAINLESS-STEEL");

    //NIST PARAFFIN
    nistManager->FindOrBuildMaterial("G4_PARAFFIN");


    //NIST Mylar
    nistManager->FindOrBuildMaterial("G4_MYLAR");

    //NIST Cadmium
    nistManager->FindOrBuildMaterial("G4_Cd");


    //
    // Define CROSS materials
    //

    //NIST Carbon
    nistManager->FindOrBuildMaterial("G4_C");

    //NIST Carbon
    nistManager->FindOrBuildMaterial("G4_NYLON-6/6");

    //NIST POM, for CROSS
    nistManager->FindOrBuildMaterial("G4_POLYOXYMETHYLENE");

    //NIST PE, for CROSS source capsule
    nistManager->FindOrBuildMaterial("G4_POLYETHYLENE");


    ////////////////////////////////////////////////////////////////////////////////////////////////
    //                  PEN staff copied here
    ////////////////////////////////////////////////////////////////////////////////////////////////
    G4int nelements;

    G4Material* fPOM = new G4Material("POM",density=1.41*g/cm3,nelements=3);
    fPOM->AddElement(O,1);
    fPOM->AddElement(C,1);
    fPOM->AddElement(H,2);
  
    G4Material* fABS = new G4Material("ABS",density=1.07*g/cm3,nelements=3);
    fABS->AddElement(C,15);
    fABS->AddElement(H,17);
    fABS->AddElement(N,1);
  
    // Scintillators
    G4int number_of_atoms;
    G4Material* PenMaterial = new G4Material("PEN", density= 1.3*g/cm3, nelements=3,kStateSolid, expTemp,1*atmosphere);
    //G4Material* PenMaterial = new G4Material("PEN", density= 1.023*g/cm3, nelements=3,kStateSolid, expTemp,1*atmosphere);
    PenMaterial->AddElement(O, number_of_atoms=4);
    PenMaterial->AddElement(H, number_of_atoms=10);
    PenMaterial->AddElement(C, number_of_atoms=14);
    PenMaterial->GetIonisation()->SetMeanExcitationEnergy(64.7*eV);
    PenMaterial->GetIonisation()->SetBirksConstant(0.15*mm/MeV);

    //Reflector 
    G4Material* materialVikuiti = new G4Material("Vikuiti", density = 1.29*g/cm3, nelements=3,kStateSolid, expTemp,1*atmosphere);
    materialVikuiti->AddElement(O, number_of_atoms=4);
    materialVikuiti->AddElement(H, number_of_atoms=8);
    materialVikuiti->AddElement(C, number_of_atoms=10);

    G4int vikuiti_entries = 0;
    G4double vikuiti_energy[202];
    G4double vikuiti_reflect[202];
    G4double zero_vikuiti[202];

    std::ifstream Read_vikuiti;
    G4String vikuiti_file = "../properties/vikuiti.dat";

    Read_vikuiti.open(vikuiti_file);
    if (Read_vikuiti.is_open()){
        while(!Read_vikuiti.eof()){
            G4String filler;
            G4double wavelength;
            G4double vikuiti_ref_coeff;
            Read_vikuiti >> wavelength >> filler >> vikuiti_ref_coeff;
            vikuiti_energy[vikuiti_entries] = (1240./wavelength)*eV;
            vikuiti_reflect[vikuiti_entries] = vikuiti_ref_coeff;
            zero_vikuiti[vikuiti_entries] = 1e-6;
            vikuiti_entries++;
            if(vikuiti_entries > 201){break ;}
        }
    }
    else
        G4cout << "Error opening file: " << vikuiti_file << G4endl;
    Read_vikuiti.close();
    vikuiti_entries = 202;

    G4MaterialPropertiesTable *mptVikuiti;
    mptVikuiti = new G4MaterialPropertiesTable();

    mptVikuiti->AddProperty("REFLECTIVITY",vikuiti_energy,vikuiti_reflect,vikuiti_entries);
    mptVikuiti->AddProperty("SPECULARLOBECONSTANT",vikuiti_energy,zero_vikuiti,vikuiti_entries);
    mptVikuiti->AddProperty("SPECULARSPIKECONSTANT",vikuiti_energy,zero_vikuiti,vikuiti_entries);
    mptVikuiti->AddProperty("BACKSCATTERCONSTANT",vikuiti_energy,zero_vikuiti,vikuiti_entries);

    
    materialVikuiti->SetMaterialPropertiesTable(mptVikuiti);

    //Glass PMT
    nistManager->FindOrBuildMaterial("G4_Pyrex_Glass");
    G4Material* silica_SiO2 = nistManager->FindOrBuildMaterial("G4_SILICON_DIOXIDE");
    G4Material* boronOxide_B2O3 = nistManager->FindOrBuildMaterial("G4_BORON_OXIDE");
    G4Material* sodiumOxide_Na2O = nistManager->FindOrBuildMaterial("G4_SODIUM_MONOXIDE");
    G4Material* potasiumOxide_K2O = nistManager->FindOrBuildMaterial("G4_POTASSIUM_OXIDE");
    G4Material* aluminumOxide_Al2O3 = nistManager->FindOrBuildMaterial("G4_ALUMINUM_OXIDE");

    
    density = 2.23*g/cm3;
    G4Material* matererialGlassPMT = new G4Material("BorosilicateGlass", density, ncomponents=5,
		    kStateSolid, expTemp, 1*atmosphere);
    matererialGlassPMT->AddMaterial(silica_SiO2, 81.*perCent);
    matererialGlassPMT->AddMaterial(boronOxide_B2O3, 13.*perCent);
    matererialGlassPMT->AddMaterial(sodiumOxide_Na2O, 2.*perCent);
    matererialGlassPMT->AddMaterial(potasiumOxide_K2O, 2.*perCent);
    matererialGlassPMT->AddMaterial(aluminumOxide_Al2O3, 2.*perCent);


    //Glass PMT
    //G4Material* matererialGlassPMT =  nistManager->FindOrBuildMaterial("G4_Pyrex_Glass");
    G4MaterialPropertiesTable *MPT_GlassPMT = new G4MaterialPropertiesTable();
    G4int glassPMTEntries = 400;
    G4double glassPMTRindex[400];
    G4double glassPMTEnergy[400];
    G4double glassPMTAbs[400];
    std::ifstream Read_glassPMT;
    G4String glassPMT_file = "../properties/Rindex_real_glass_pmt.txt";

    G4int counter_glass = 0;
    Read_glassPMT.open(glassPMT_file);
    if(Read_glassPMT.is_open()){
        while(!Read_glassPMT.eof()){
                G4double wavelength, glass_rindex;
                Read_glassPMT >> wavelength >> glass_rindex;
                glassPMTEnergy[counter_glass] = (1240./wavelength)*eV;
                glassPMTRindex[counter_glass] = glass_rindex;
                glassPMTAbs[counter_glass] = 1*m;
                //G4cout<<" "<<counter_glass<<" glassPMTRindex "<<glassPMTRindex[counter_glass]<<G4endl;
                counter_glass++;
                if(counter_glass > (glassPMTEntries - 1)){break;}
        }
    }else{G4cout << "Error opening file: " << glassPMT_file << G4endl;}
    Read_glassPMT.close();
    MPT_GlassPMT->AddProperty("RINDEX",glassPMTEnergy,glassPMTRindex,glassPMTEntries);
    MPT_GlassPMT->AddProperty("ABSLENGTH",glassPMTEnergy,glassPMTAbs,glassPMTEntries);
    matererialGlassPMT -> SetMaterialPropertiesTable(MPT_GlassPMT);


    nistManager->FindOrBuildMaterial("G4_Si");
  
    nistManager->FindOrBuildMaterial("G4_TEFLON");
  
  
    G4Material* materialTriggerFoilEJ212 = new G4Material("EJ212", density= 1.023*g/cm3, 2);
    materialTriggerFoilEJ212->AddElement(C, 0.475);
    materialTriggerFoilEJ212->AddElement(H, 0.525);

    //A complete guess of the photocathode material ... shouldn't impact the simulations
    G4Material* materialBialkali = new G4Material("Bialkali", density = 1.4 *g/cm3, 3);
    materialBialkali->AddElement(Sb, 0.32);
    materialBialkali->AddElement(Cs, 0.32);
    materialBialkali->AddElement(K, 0.36);
  
    //Titanium foil for Bi source
    G4Material* materialTitanium = new G4Material("titanium", 4.54*g/cm3, 2);
    materialTitanium->AddElement(Ti,0.99);
    materialTitanium->AddElement(O,0.01);

}
