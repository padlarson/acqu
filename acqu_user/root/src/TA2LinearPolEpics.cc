//--Author	Ken Livingston   17th March 2004
//--Rev 	Ken Livingston... sometime 1st production version
//--Update	Ken Livingston...19th Feb 2004...User code
//--Description
//                *** Acqu++ <-> Root ***
// Online/Offline Analysis of Sub-Atomic Physics Experimental Data
//
// TA2LinearPolEpics
//
// User coded version of Crystal Ball apparatus

#include "TA2LinearPolEpics.h"
#include "TF1.h"
#include <iostream>

enum { ELpInitLevel0,ELpInitLevel1,ELpInitLevel2 };  	// allow multiple calls to PostInitialise();
enum { ELpMiscApp, // options for ParseMisc
       ELpMiscPVNames,ELpMiscRads,ELpMiscPolLookupTable, ELpMiscPolCuts, ELpMiscRunRanges };
enum { ELpPolLookupPhotonEnergy,ELpPolLookupEdge, 	//parsing pol lookup table
       ELpPolLookupPol,ELpPolLookupPlane};
enum { ELpPolFromParams, ELpPolFromFile };



static Map_t kValidDetectors[] = {
 {NULL, 		-1}
};



static const Map_t kValidMiscParams[] = {
  {"app",      	ELpMiscApp},
  {"pvnames", 	ELpMiscPVNames},
  {"rads", 	ELpMiscRads},
  {"poltable",  ELpMiscPolLookupTable},
  {"polcuts",   ELpMiscPolCuts},
  {"runs",      ELpMiscRunRanges},
  {NULL,      	-1}
};

static const Map_t kValidPolLookupTableParams[] = {
  {"PhotonEnergy:",      	ELpPolLookupPhotonEnergy},
  {"Edge:",    			ELpPolLookupEdge},
  {"Pol:", 			ELpPolLookupPol}, 
  {"Plane:", 			ELpPolLookupPlane}, 
  {NULL,      	-1}
};

const Double_t k = 26.5601; 
const Int_t VECTORS[]={2,4,6,8,10};    //list of the vectors to be included (022,044);

//fit of a gaussian on a baseline, for gausFit,
Double_t GausOnBase(Double_t *x, Double_t *par) {
  Double_t arg = 0;
  if (par[2] != 0) arg = (x[0] - par[1])/par[2];
   Double_t fitval = par[3] + par[0]*TMath::Exp(-0.5*arg*arg);
   return fitval;
}


//-----------------------------------------------------------------------------
TA2LinearPolEpics::TA2LinearPolEpics( const char* name, TA2System* fAnalysis  )
  :TA2Apparatus( name, fAnalysis, kValidDetectors )
{
  fInitLevel=EInitLevel0;	//Flag that no init done yet
  fPolPlane = (UShort_t)ENullHit;
  fPlane = -1;
  fLastEdge = -1.0;
  fEdgeSetting = -1.0;
  fEdgeRange   = -1.0;
  sprintf(fPlaneString,"Epics"); //Assume plane information from data stream

  sprintf(fCurrentRunFileName,"none"); 	//set to some dummy value
  fIsNewFile = ETrue;			//to begin with it's a new file
  
  fRunNo=0;
  
  fDInc = EBufferEnd;		//make all these Buffer ends, 
  fDCoh = EBufferEnd;		//since they don't ever get filled	
  fDEnh = EBufferEnd;		//by the system, but in the Reconstruct loop directly.
  fDCohPara = EBufferEnd;
  fDEnhPara = EBufferEnd;
  fDCohPerp = EBufferEnd;
  fDEnhPerp = EBufferEnd;
  
  fDEdge = EBufferEnd;
  fDEdgePerp = EBufferEnd;
  fDEdgePara = EBufferEnd;
  fDEdgeDistPara = EBufferEnd;
  fDEdgeDistPerp = EBufferEnd;
  fDEdgeEpics = EBufferEnd;
  fDEdgePerpEpics = EBufferEnd;
  fDEdgeParaEpics = EBufferEnd;
  fDEdgeDistParaEpics = EBufferEnd;
  fDEdgeDistPerpEpics = EBufferEnd;
  fDPolTableEnh = EBufferEnd;
  fDPolTablePol = EBufferEnd;
  fDPolMean =  EBufferEnd;
  fDPolCount = EBufferEnd;

  fDoingScalers = EFalse;
  fHaveIncScaler = EFalse;	//Flag that Inc ref scaler is not loaded
  fHaveTaggerApp = EFalse;	//Flag that tagger apparatus is not loaded
  fHaveLadderDet = EFalse;	//Flag that ladder detector is not loaded
  fScalerCount = 0;	       	//Running counter for scaler reads

  fNPolLookupEdges[ETablePara] = 0;  //init no of edge positions for which there are tables
  fNPolLookupEdges[ETablePerp] = 0;  //init no of edge positions for which there are tables
  fNPolLookupEnergies =0;	//init no of photon energies in each table
  fPolLookupEdges[ETablePara] = NULL;      	//no pol lookup edge array yet
  fPolLookupEdges[ETablePerp] = NULL;      	//no pol lookup edge array yet
  fPolLookupEnergies = NULL;	//no pol lookup energy table ceated yet
  fPolLookupPolarisations[ETablePara] = NULL; //no pol  lookup pol tables ceated yet
  fPolLookupPolarisations[ETablePerp] = NULL; //no pol  lookup pol tables ceated yet
    
  fEpics           = NULL;      // init epics stuff
  fEpicsInit       = 0;
  fNEpics          = 0;
  fEpicsBuff       = NULL;
  fEpicsIndex      = NULL;
  fEpicsCount      = 0; 


  fIsEpics         = kTRUE;   //assume epics unless overridden by setup file
  fIsDiamond       = kFALSE;  //assume no diamond to start

  fHaveApp         = kFALSE; //init flags for options in the setup file
  fHaveScaler      = kFALSE; 
  fHaveCbrem       = kFALSE; 
  fHavePVNames     = kFALSE; 
  fHaveRads        = kFALSE; 
  fHaveFit         = kFALSE; 
  fHavePolTable    = kFALSE; 

  fPolTableNRanges[ETablePara] = 0;
  fPolTableNRanges[ETablePerp] = 0;
  

  fNEpicsVariables = 8;
  fLinPolEpicsVariables[0]=&fA2LinPolCohEdgeEpics;  //store the addresses of the epics
  fLinPolEpicsVariables[1]=&fA2LinPolCohRadiatorEpics;  //values in array
  fLinPolEpicsVariables[2]=&fA2LinPolCohPlaneEpics;
  fLinPolEpicsVariables[3]=&fA2GoniAxis1Epics;
  fLinPolEpicsVariables[4]=&fA2GoniAxis2Epics;
  fLinPolEpicsVariables[5]=&fA2GoniAxis3Epics;
  fLinPolEpicsVariables[6]=&fA2GoniAxis4Epics;
  fLinPolEpicsVariables[7]=&fA2GoniAxis5Epics;

  //default values for all in setup file
  //app
  fNormEnergy      = 700.0;  // defaults for Ee=855, peak ~ 350;
  fEdgeMin         = 280.0;
  fEdgeMax         = 420.0;
  sprintf(fTaggerName,"TAGG"); //default to TAGG, FPD
  sprintf(fTaggerName,"FPD");

  //scaler
  sprintf(fRunRefFiles[0],"dats/scaler.dmp");

  //cbrem
  sprintf(fEdgeString,"epics"); //default to using epics for plane and coh edge
  fForceFit=kFALSE;
  sprintf(fPlaneString,"epics");
  fForcePlane=kFALSE;

  //pvnames relating to variables above in  fLinPolEpicsVariables array
  sprintf(fLinPolEpicsPVNames[0],"A2LinPolCohEdge");
  sprintf(fLinPolEpicsPVNames[1],"A2LinPolCohRadiator");
  sprintf(fLinPolEpicsPVNames[2],"A2LinPolCohPlane");
  sprintf(fLinPolEpicsPVNames[3],"A2GoniAxis1");
  sprintf(fLinPolEpicsPVNames[4],"A2GoniAxis2");
  sprintf(fLinPolEpicsPVNames[5],"A2GoniAxis3");
  sprintf(fLinPolEpicsPVNames[6],"A2GoniAxis4");
  sprintf(fLinPolEpicsPVNames[7],"A2GoniAxis5");

  //rads
  fDiamondID       = 4;       //default rads from recent beamtime, override in setup file
  fAmoID           = 2;


  //
  fBeforeEdge      = 150.0;
  fAfterEdge       =  20.0;
  fPolMin          =   0.1;

}


//-----------------------------------------------------------------------------
TA2LinearPolEpics::~TA2LinearPolEpics()
{
  // Free up allocated memory
  // (Note: It is safe to delete NULL pointers)
  delete [] fPolLookupEnergies;

  for(int i=0; i<2; ++i) {
      delete fPolLookupEdges[i];
  }

  for(int m=0;m<2;m++){
    if(fPolLookupPolarisations[m]){	
      for(int n=0;n<fNPolLookupEdges[m];n++){
        delete fPolLookupPolarisations[m][n];
      }
      delete [] fPolLookupPolarisations[m];
    }
  }
}

void TA2LinearPolEpics::LoadVariable( )
{
  // Input name - variable pointer associations for any subsequent
  // cut or histogram setup.
  // LoadVariable( "name", pointer-to-variable, type-spec );
  // NB scaler variable pointers need the preceeding &
  //    array variable pointers do not.
  // type-spec ED prefix for a Double_t variable
  //           EI prefix for an Int_t variable
  // type-spec SingleX for a single-valued variable
  //           MultiX  for a multi-valued variable

  //                            name        		pointer      	type-spec
  TA2DataManager::LoadVariable("Incoherent", 		&fDInc,    	EDSingleX);
  TA2DataManager::LoadVariable("Coherent", 		&fDCoh,       	EDSingleX);
  TA2DataManager::LoadVariable("Enhancement", 		&fDEnh,      	EDSingleX);
  TA2DataManager::LoadVariable("CoherentPara", 		&fDCohPara,    	EDSingleX);
  TA2DataManager::LoadVariable("EnhancementPara", 	&fDEnhPara,    	EDSingleX);
  TA2DataManager::LoadVariable("CoherentPerp", 		&fDCohPerp,    	EDSingleX);
  TA2DataManager::LoadVariable("EnhancementPerp", 	&fDEnhPerp,    	EDSingleX);
  TA2DataManager::LoadVariable("CohEdge",        	&fDEdge,  EDSingleX);
  TA2DataManager::LoadVariable("CohEdgePara",     &fDEdgePara,  	EDSingleX);
  TA2DataManager::LoadVariable("CohEdgePerp",    &fDEdgePerp,  	EDSingleX);
  TA2DataManager::LoadVariable("CohEdgeDistPara",      	&fDEdgeDistPara,       	EDSingleX);
  TA2DataManager::LoadVariable("CohEdgeDistPerp",      	&fDEdgeDistPerp,       	EDSingleX);
  TA2DataManager::LoadVariable("CohEdgeEpics",   	&fDEdgeEpics,  EDSingleX);
  TA2DataManager::LoadVariable("CohEdgeParaEpics", &fDEdgeParaEpics,  	EDSingleX);
  TA2DataManager::LoadVariable("CohEdgePerpEpics", &fDEdgePerpEpics,  	EDSingleX);
  TA2DataManager::LoadVariable("CohEdgeDistParaEpics",  &fDEdgeDistParaEpics,       	EDSingleX);
  TA2DataManager::LoadVariable("CohEdgeDistPerpEpics",  &fDEdgeDistPerpEpics,       	EDSingleX);
  TA2DataManager::LoadVariable("PolPlane",        	&fPolPlane,  	EDSingleX);
  TA2DataManager::LoadVariable("PolTableEnh",        	&fDPolTableEnh,  	EDSingleX);
  TA2DataManager::LoadVariable("PolTablePol",        	&fDPolTablePol,  	EDSingleX);
  TA2DataManager::LoadVariable("PolMean",        	&fDPolMean,  	EDSingleX);
  TA2DataManager::LoadVariable("PolCount",        	&fDPolCount,  	EDSingleX);
  //
}
//-----------------------------------------------------------------------------
void TA2LinearPolEpics::PostInitialise( )
{
	
  //we need more than one pass at this (ie an Initialise: must come before and after the Display lines in
  //the setup file for this if we're to get the histograms filled.
  Char_t histName[100];

  const Double_t *ladderECal;
  Double_t energy;

  switch(fInitLevel){
  case ELpInitLevel0:	//1st pass before Display: to do most of the init
    

    //look for the Tagger
    if((fTagger=(TA2Tagger *)fParent->GetChild(fTaggerName))!=NULL){
      fHaveTaggerApp=ETrue;
    }
    else{
      fprintf(stderr,"Warning: Couldn't find TA2Tagger *%s\n",fTaggerName);
      break;
    }
    if((fLadder=(TA2Ladder *)fTagger->GetChild(fLadderName,"TA2Detector"))!=NULL){
      fHaveLadderDet=ETrue;
      fPolArray = new Double_t[fLadder->GetNelem()];
      fPolArray[0] =  (Double_t)ENullHit;
      fLadderHits = fLadder->GetHits();

      //This version will have the pol for the multihits tagged on at the end
      fPolArrayM = new Double_t[fLadder->GetNelem()];
      fPolArrayM[0] =  (Double_t)ENullHit;
    }
    else{
      fprintf(stderr,"Warning: Couldn't find TA2Ladder *%s\n",fLadderName);
      break;
    }
    
    fBeamEnergy = fTagger->GetBeamEnergy();	//get the beam energy
    fTaggerChannels = fLadder->GetNelem();	//get the no of elements in the Ladder
    fCurrentPolTable = new Double_t[fLadder->GetNelem()]; //get the no of elements in the Ladder
    fCurrentPolTable_TC = new Double_t[fLadder->GetNelem()]; //get the no of elements in the Ladder
    fCurrentEnhTable = new Double_t[fLadder->GetNelem()]; //get the no of elements in the Ladder
    fCurrentEnhTable_TC = new Double_t[fLadder->GetNelem()]; //get the no of elements in the Ladder
    for(UInt_t n=0;n<fLadder->GetNelem();n++){
      fCurrentPolTable[n]	=-1.0;
      fCurrentPolTable_TC[n]=-1.0;
      fCurrentEnhTable[n]	= 0.0;      
      fCurrentEnhTable_TC[n]= 0.0;
    }
    // to allow for summing consecutive scaler buffers to get better stats.
    fAccScaler=new Double_t*[fNScalerBuffers];
    for(int n=0;n<fNScalerBuffers;n++){
      fAccScaler[n]=new Double_t[fLadder->GetNelem()];
      for(UInt_t m=0;m<fLadder->GetNelem();m++) fAccScaler[n][m]=0.0;
    }
    fScalerEvent=0;
    
    fNormChannel=fTaggerChannels/2;			// default norm chan to .5 of range

    if((ladderECal = fLadder->GetECalibration())){	//get the ladder energy calibration
      fEnergyCalib = new Double_t[fTaggerChannels];
      fLadderPhotonEnergy = new Double_t[fTaggerChannels];
      fEnergyBins = new Double_t[fTaggerChannels+1];
      for(int n=0;n<fTaggerChannels;n++){		// fill the photon energy calib
	energy=fBeamEnergy-ladderECal[n];
	fLadderPhotonEnergy[n]=energy;
	fEnergyCalib[fTaggerChannels-1-n]=energy;
      }
      //Make array of bin lower limits for hist axes
      //first bin low edge
      fEnergyBins[0]=fEnergyCalib[0]-((fEnergyCalib[1]-fEnergyCalib[0])/2.0);
      //last bin high edge
      fEnergyBins[fTaggerChannels]=fEnergyCalib[fTaggerChannels-1]+
	((fEnergyCalib[fTaggerChannels-1]-fEnergyCalib[fTaggerChannels-2])/2.0);
      for(int n=1;n<fTaggerChannels;n++){		// fill the photon energy calib
	fEnergyBins[n]=0.5*(fEnergyCalib[n-1]+fEnergyCalib[n]);
      }
    }

    if((fDoingScalers)&&(fLadder->IsScaler())){		// If using scalers
      fIncSpectrum = new Double_t[fTaggerChannels];	// Array to hold Inc ref. spectrum
      fCohSpectrum = new Double_t[fTaggerChannels];	// Array to hold Coh spectrum
      fEnhSpectrum = new Double_t[fTaggerChannels];	// Array to hold Enhancement spectrum
      fBadScalerChan = new Bool_t[fTaggerChannels];	// Array to hold bad scaler channels 
      
      fScalerCurr=fLadder->GetScalerCurr();

      LoadAmoRef(fRunRefFiles[0]);					//Load reference amorphous scaler scectrum 
    }

    // Only do this if we have an event
    //    if(fIsPolLookupTable){
    // if(LoadPolLookupTable() != 0){
    //	fprintf(stderr,"Warning failed to load PolLookupTable, GetPol() will always return 0\n");
    //  }
    // }
      
    TA2DataManager::PostInit();	//call the Init fron the base class
    break;

  case ELpInitLevel1:	//2nd pass after display to get the ptrs to the histograms
    if((!fHaveTaggerApp)||(!fHaveLadderDet)) break;
    if(!f1Dhist) break;
    sprintf(histName,"%s_CohEdge",fName.Data());
    fHEdge=(TH1F *)f1Dhist->FindObject(histName);
    if(fHEdge!=NULL){
      fHEdge->SetStats(kFALSE);
      fHEdge->GetXaxis()->SetTitle("Scaler Reads");
      fHEdge->GetYaxis()->SetTitle("Edge Position (MeV)");
      fHEdge->SetMinimum(fEdgeMin);
      fHEdge->SetMaximum(fEdgeMax);
    }
    sprintf(histName,"%s_CohEdgePerp",fName.Data());
    fHEdgePerp=(TH1F *)f1Dhist->FindObject(histName);
    if(fHEdgePerp!=NULL){
      fHEdgePerp->SetStats(kFALSE);
      fHEdgePerp->GetXaxis()->SetTitle("Scaler Reads");
      fHEdgePerp->GetYaxis()->SetTitle("Edge Position (MeV)");
      fHEdgePerp->SetMinimum(fEdgeMin);
      fHEdgePerp->SetMaximum(fEdgeMax);
      fHEdgePerp->SetLineColor(4);
    }
    sprintf(histName,"%s_CohEdgePara",fName.Data());
    fHEdgePara=(TH1F *)f1Dhist->FindObject(histName);
    if(fHEdgePara!=NULL){
      fHEdgePara->SetStats(kFALSE);
      fHEdgePara->GetXaxis()->SetTitle("Scaler Reads");
      fHEdgePara->GetYaxis()->SetTitle("Edge Position (MeV)");
      fHEdgePara->SetMinimum(fEdgeMin);
      fHEdgePara->SetMaximum(fEdgeMax);
      fHEdgePara->SetLineColor(2);
    }
    sprintf(histName,"%s_CohEdgeDistPerp",fName.Data());
    fHEdgeDistPerp=(TH1F *)f1Dhist->FindObject(histName);
    if(fHEdgeDistPerp!=NULL){
      fHEdgeDistPerp->SetStats(kFALSE);
      fHEdgeDistPerp->GetXaxis()->SetTitle("Coherent Edge Position (MeV)");
      fHEdgeDistPerp->GetYaxis()->SetTitle("Scaler Reads");
      fHEdgeDistPerp->GetXaxis()->Set(fHEdgeDistPerp->GetNbinsX(),fEdgeMin,fEdgeMax);
      fHEdgeDistPerp->SetLineColor(4);
    }
    sprintf(histName,"%s_CohEdgeDistPara",fName.Data());
    fHEdgeDistPara=(TH1F *)f1Dhist->FindObject(histName);
    if(fHEdgeDistPara!=NULL){
      fHEdgeDistPara->SetStats(kFALSE);
      fHEdgeDistPara->GetXaxis()->SetTitle("Coherent Edge Position (MeV)");
      fHEdgeDistPara->GetYaxis()->SetTitle("Scaler Reads");
      fHEdgeDistPara->GetXaxis()->Set(fHEdgeDistPara->GetNbinsX(),fEdgeMin,fEdgeMax);
      fHEdgeDistPara->SetLineColor(2);
    }

    sprintf(histName,"%s_CohEdgeEpics",fName.Data());
    fHEdgeEpics=(TH1F *)f1Dhist->FindObject(histName);
    if(fHEdgeEpics!=NULL){
      fHEdgeEpics->SetStats(kFALSE);
      fHEdgeEpics->GetXaxis()->SetTitle("Epics Reads");
      fHEdgeEpics->GetYaxis()->SetTitle("Edge Position (MeV)");
      fHEdgeEpics->SetMinimum(fEdgeMin);
      fHEdgeEpics->SetMaximum(fEdgeMax);
    }
    sprintf(histName,"%s_CohEdgePerpEpics",fName.Data());
    fHEdgePerpEpics=(TH1F *)f1Dhist->FindObject(histName);
    if(fHEdgePerpEpics!=NULL){
      fHEdgePerpEpics->SetStats(kFALSE);
      fHEdgePerpEpics->GetXaxis()->SetTitle("Epics Reads");
      fHEdgePerpEpics->GetYaxis()->SetTitle("Edge Position (MeV)");
      fHEdgePerpEpics->SetMinimum(fEdgeMin);
      fHEdgePerpEpics->SetMaximum(fEdgeMax);
      fHEdgePerpEpics->SetLineColor(4);
    }
    sprintf(histName,"%s_CohEdgeParaEpics",fName.Data());
    fHEdgeParaEpics=(TH1F *)f1Dhist->FindObject(histName);
    if(fHEdgeParaEpics!=NULL){
      fHEdgeParaEpics->SetStats(kFALSE);
      fHEdgeParaEpics->GetXaxis()->SetTitle("Epics Reads");
      fHEdgeParaEpics->GetYaxis()->SetTitle("Edge Position (MeV)");
      fHEdgeParaEpics->SetMinimum(fEdgeMin);
      fHEdgeParaEpics->SetMaximum(fEdgeMax);
      fHEdgeParaEpics->SetLineColor(2);
    }
    sprintf(histName,"%s_CohEdgeDistPerpEpics",fName.Data());
    fHEdgeDistPerpEpics=(TH1F *)f1Dhist->FindObject(histName);
    if(fHEdgeDistPerpEpics!=NULL){
      fHEdgeDistPerpEpics->SetStats(kFALSE);
      fHEdgeDistPerpEpics->GetXaxis()->SetTitle("Coherent Edge Position (MeV)");
      fHEdgeDistPerpEpics->GetYaxis()->SetTitle("Epics Reads");
      fHEdgeDistPerpEpics->GetXaxis()->Set(fHEdgeDistPerpEpics->GetNbinsX(),fEdgeMin,fEdgeMax);
      fHEdgeDistPerpEpics->SetLineColor(4);
    }
    sprintf(histName,"%s_CohEdgeDistParaEpics",fName.Data());
    fHEdgeDistParaEpics=(TH1F *)f1Dhist->FindObject(histName);
    if(fHEdgeDistParaEpics!=NULL){
      fHEdgeDistParaEpics->SetStats(kFALSE);
      fHEdgeDistParaEpics->GetXaxis()->SetTitle("Coherent Edge Position (MeV)");
      fHEdgeDistParaEpics->GetYaxis()->SetTitle("Epics Reads");
      fHEdgeDistParaEpics->GetXaxis()->Set(fHEdgeDistParaEpics->GetNbinsX(),fEdgeMin,fEdgeMax);
      fHEdgeDistParaEpics->SetLineColor(2);
    }

    if((!fHaveTaggerApp)||(!fHaveLadderDet)) break;
    sprintf(histName,"%s_Incoherent",fName.Data());
    fHInc=(TH1F *)f1Dhist->FindObject(histName);
    if(fHInc!=NULL){
      fHInc->GetXaxis()->Set(fTaggerChannels,fEnergyBins);
      fHInc->GetXaxis()->SetTitle("Photon Energy (MeV)");
    }
    sprintf(histName,"%s_Coherent",fName.Data());
    fHCoh=(TH1F *)f1Dhist->FindObject(histName);
    if(fHCoh!=NULL){
      fHCoh->GetXaxis()->Set(fTaggerChannels,fEnergyBins);
      fHCoh->GetXaxis()->SetTitle("Photon Energy (MeV)");
    }
    sprintf(histName,"%s_Enhancement",fName.Data());
    fHEnh=(TH1F *)f1Dhist->FindObject(histName);
    if(fHEnh!=NULL){
      fHEnh->SetStats(kFALSE);
      fHEnh->SetMinimum(80.0);
      fHEnh->SetMaximum(300.0);
      fHEnh->GetXaxis()->Set(fTaggerChannels,fEnergyBins);
      fHEnh->GetXaxis()->SetTitle("Photon Energy (MeV)");
      fHEnh->SetLineColor(1);
    }
    sprintf(histName,"%s_CoherentPara",fName.Data());
    fHCohPara=(TH1F *)f1Dhist->FindObject(histName);
    if(fHCohPara!=NULL){
      fHCohPara->GetXaxis()->Set(fTaggerChannels,fEnergyBins);
      fHCohPara->GetXaxis()->SetTitle("Photon Energy (MeV)");
      fHCohPara->SetLineColor(2);
    }
    sprintf(histName,"%s_EnhancementPara",fName.Data());
    fHEnhPara=(TH1F *)f1Dhist->FindObject(histName);
    if(fHEnhPara!=NULL){
      fHEnhPara->SetStats(kFALSE);
      fHEnhPara->SetMinimum(80.0);
      fHEnhPara->SetMaximum(300.0);
      fHEnhPara->GetXaxis()->Set(fTaggerChannels,fEnergyBins);
      fHEnhPara->GetXaxis()->SetTitle("Photon Energy (MeV)");
      fHEnhPara->SetLineColor(2);
   }

    sprintf(histName,"%s_CoherentPerp",fName.Data());
    fHCohPerp=(TH1F *)f1Dhist->FindObject(histName);
    if(fHCohPerp!=NULL){
      fHCohPerp->GetXaxis()->Set(fTaggerChannels,fEnergyBins);
      fHCohPerp->GetXaxis()->SetTitle("Photon Energy (MeV)");
      fHCohPerp->SetLineColor(4);
    }
    sprintf(histName,"%s_EnhancementPerp",fName.Data());
    fHEnhPerp=(TH1F *)f1Dhist->FindObject(histName);
    if(fHEnhPerp!=NULL){
      fHEnhPerp->SetStats(kFALSE);
      fHEnhPerp->SetMinimum(80.0);
      fHEnhPerp->SetMaximum(300.0);
      fHEnhPerp->GetXaxis()->Set(fTaggerChannels,fEnergyBins);
      fHEnhPerp->GetXaxis()->SetTitle("Photon Energy (MeV)");
      fHEnhPerp->SetLineColor(4);
    }
    
    if(fHaveIncScaler){	// if we have a reference
      if(fHInc!=NULL){	// and we have setup a histogram for it
	for(int n=0;n<fTaggerChannels;n++){
	  fHInc->Fill(fEnergyBins[n],fIncSpectrum[n]);	//fill the histogram
	}
      }
    }
    
    if(fHavePolTable){     //If we have pol lookup tables, make the relevant histograms
      fHistE        = new TH1F("PolTableEnhancement", "Pol Table Enhancement", fTaggerChannels,fEnergyBins );
      fHistP        = new TH1F("PolTablePol", "Pol Table Polarization",fTaggerChannels,fEnergyBins);
      fHistE->SetMinimum(0);
      fHistP->SetMinimum(0);
      fHistP->SetMaximum(1);
      fWeightHist   = new TH1F("weightHist",  "weightHist", THETASTEPS+1, 0, THETASTEPS+1 );
      fThetaPol     = new TH2F("thetaPol",    "thetaPol",   fHistE->GetNbinsX(), fHistE->GetXaxis()->GetXbins()->GetArray(), THETASTEPS+1,0, THETASTEPS+1);
      fThetaItot    = new TH2F("thetaItot",   "thetaItot",  fHistE->GetNbinsX(), fHistE->GetXaxis()->GetXbins()->GetArray(), THETASTEPS+1,0, THETASTEPS+1);
      
    }
    
    sprintf(histName,"%s_PolTableEnh",fName.Data());
    fHPolTableEnh=(TH1F *)f1Dhist->FindObject(histName);
    if(fHPolTableEnh!=NULL){
      fHPolTableEnh->SetStats(kFALSE);
      fHPolTableEnh->SetMinimum(0.0);
      fHPolTableEnh->SetMaximum(6.0);
      fHPolTableEnh->GetXaxis()->Set(fTaggerChannels,fEnergyBins);
      fHPolTableEnh->GetXaxis()->SetTitle("Photon Energy (MeV)");
      fHPolTableEnh->SetLineColor(4);
    }
    sprintf(histName,"%s_PolTablePol",fName.Data());
    fHPolTablePol=(TH1F *)f1Dhist->FindObject(histName);
    if(fHPolTablePol!=NULL){
      fHPolTablePol->SetStats(kFALSE);
      fHPolTablePol->SetMinimum(0.0);
      fHPolTablePol->SetMaximum(1.0);
      fHPolTablePol->GetXaxis()->Set(fTaggerChannels,fEnergyBins);
      fHPolTablePol->GetXaxis()->SetTitle("Photon Energy (MeV)");
      fHPolTablePol->SetLineColor(4);
    }
    sprintf(histName,"%s_PolMean",fName.Data());
    fHPolMean=(TH1F *)f1Dhist->FindObject(histName);
    if(fHPolMean!=NULL){
      fHPolMean->SetStats(kFALSE);
      fHPolMean->SetMinimum(0.0);
      fHPolMean->SetMaximum(1.0);
      fHPolMean->GetXaxis()->Set(fTaggerChannels,fEnergyBins);
      fHPolMean->GetXaxis()->SetTitle("Photon Energy (MeV)");
    }
    sprintf(histName,"%s_PolCount",fName.Data());
    fHPolCount=(TH1F *)f1Dhist->FindObject(histName);
    if(fHPolCount!=NULL){
      fHPolCount->SetStats(kFALSE);
      fHPolCount->SetMinimum(0.0);
      fHPolCount->GetXaxis()->Set(fTaggerChannels,fEnergyBins);
      fHPolCount->GetXaxis()->SetTitle("Photon Energy (MeV)");
    }
	 
  default:
    break;
  }
  fInitLevel++;
  fIsInit=ETrue;  
}

//-----------------------------------------------------------------------------
void TA2LinearPolEpics::Reconstruct( ){
  
  Double_t normValue=0.0;
  int ncount=0;
  TF1 *edgeFit=NULL;
  Double_t fitedge;
  Double_t maxgrad;
  Double_t coh_sum=0.0;
  int binx=0;
  Double_t xmax=0.0; 
  Double_t ymax=0.0;
  
  fLastPolRangeIndex[ETablePara]=-1;
  fLastPolRangeIndex[ETablePerp]=-1;
  fLastRunRangeIndex=-1;
  
  
  if(gAR->IsEpicsRead()||gAR->IsScalerRead()){
  
    if(strcmp(gAR->GetFileName(),fCurrentRunFileName)!=0){    //see if it's a new file
      strcpy(fCurrentRunFileName,gAR->GetFileName()); 
      fIsNewFile=ETrue;	//this stays flagged until we get a scaler or epics read
      
      fPolRangeIndex[ETablePara]=-1;
      fPolRangeIndex[ETablePerp]=-1;
      fRunRangeIndex=-1;
      
      sscanf(strrchr(fCurrentRunFileName,'_'),"_%d.dat",&fRunNo);  //get the run no from filename    
      if(fHavePolTable){
	for(int index=0;index<fPolTableNRanges[ETablePara];index++){ //find the pol table run range for this run
	  if((fRunNo>=fPolTableRangesMin[index][ETablePara])&&(fRunNo<=fPolTableRangesMax[index][ETablePara])){
	    fPolRangeIndex[ETablePara]=index;
	    break;
	  }
	}
	if(fPolRangeIndex[ETablePara]<0){
	fprintf(stderr,"Fatal Error: No para pol Table run range for run %d",fRunNo);
	// not a good way to exit (it will segfault)
	exit(1);
	}
	for(int index=0;index<fPolTableNRanges[ETablePerp];index++){ //find the pol table run range for this run
	  if((fRunNo>=fPolTableRangesMin[index][ETablePerp])&&(fRunNo<=fPolTableRangesMax[index][ETablePerp])){
	    fPolRangeIndex[ETablePerp]=index;
	    break;
	  }
	}
	if(fPolRangeIndex[ETablePerp]<0){
	fprintf(stderr,"Fatal Error: No perp Pol Table run range for run %d",fRunNo);
	// not a good way to exit (it will segfault)
	exit(1);
	}
      }
      for(int index=0;index<fNRunRanges;index++){ //find the run table run range for this run
	if((fRunNo>=fRunRangeMin[index])&&(fRunNo<=fRunRangeMax[index])){
	  fRunRangeIndex=index;
	  break;
	}
      }
      if(fRunRangeIndex>-1) {
	LoadAmoRef(fRunRefFiles[fRunRangeIndex]);  //Load new reference amorphous scaler scectrum 
	if(fForcePlane){
	  fPlane=fRunPlanes[fRunRangeIndex];
	}
      }
      else{
		fprintf(stderr,"Fatal Error: No plane info range for run %d",fRunNo);
		// not a good way to exit (it will segfault)
		exit(1);
      }
    }
    
    if(fIsEpics){   //if wer're using epics buffers to get coherent brem info
      ////////////////////////////////////////////////////////////////////////////
      // Deal with EPICS event
      //
      ////////////////////////////////////////////////////////////////////////////
      
      if (gAR->IsEpicsRead()){                   // check if epics buffer
	if(fEpicsInit == kFALSE){                // initialise if not done
	  fNEpics = gAR->GetNEpics();
	  fEpicsBuff = gAR->GetEpicsBuffer();
	  fEpicsIndex = gAR->GetEpicsIndex();
	  fEpics = new TEPICSmodule();
	  fEpicsInit = kTRUE;
	}
	
	// the rest is specific to what epics buffer and channels you need to read.
	if(fEpicsIndex[fA2LinPolEpicsIndex]){  // if the epics buffer has the index required for lin pol data
	  
	  //loop over all the variables 
	  for(Int_t nepics=0;nepics<fNEpicsVariables;nepics++){
	    //returns pointer to the data for the channel if you want to do something smart
	    fEpicsChannelBuffer = fEpics->GetChannel((Char_t*)fLinPolEpicsPVNames[nepics],     //pv name
						     &fEpicsType,                     //pv type
						     fEpicsBuff[fA2LinPolEpicsIndex], //start of epics buffer
						     fLinPolEpicsVariables[nepics],   //address of the variable to be filled
						     &fEpicsNElem);                   //no of elements in the channel (usually singles)
	  }
	  //fEpics->DumpBuffer(fEpicsBuff[fA2LinPolEpicsIndex]);
	}
	
	if(!fForcePlane){                              // if not forcing plane in setup file 
	  fPlane=(Int_t)fA2LinPolCohPlaneEpics;        //set plane and rad from epics
	  fRadiator=(Int_t)fA2LinPolCohRadiatorEpics;  
	  if(fRadiator==fDiamondID){
	    fIsDiamond=kTRUE;
	  }
	  else{
	    fIsDiamond=kFALSE;
	  }
	}
	if(!fForceFit){                              // if not fitting enhancements 
	  fEdge=fA2LinPolCohEdgeEpics;               // get edge from epics
	}
      }
      
      
      if(fHEdgeEpics){
	if(fEpicsCount%(fHEdgeEpics->GetNbinsX()-1)==0){
	  fHEdgeEpics->Reset("ICE");
	}
	fHEdgeEpics->Fill(fEpicsCount%(fHEdgeEpics->GetNbinsX()-1),fA2LinPolCohEdgeEpics);
      }
      if(fHEdgeParaEpics){
	if(fEpicsCount%(fHEdgeEpics->GetNbinsX()-1)==0){
	  fHEdgeParaEpics->Reset("ICE");
	}
	if(fPlane==ETablePara){
	  fHEdgeParaEpics->Fill(fEpicsCount%(fHEdgeParaEpics->GetNbinsX()-1),fA2LinPolCohEdgeEpics);
	}
	else{
	  fHEdgeParaEpics->Fill(fEpicsCount%(fHEdgeParaEpics->GetNbinsX()-1),0);
	}
	
      }
      if(fHEdgePerpEpics){
	if(fEpicsCount%(fHEdgePerpEpics->GetNbinsX()-1)==0){
	  fHEdgePerpEpics->Reset("ICE");
	}
	if(fPlane==ETablePerp){
	  fHEdgePerpEpics->Fill(fEpicsCount%(fHEdgePerpEpics->GetNbinsX()-1),fA2LinPolCohEdgeEpics);
	}
	else{
	  fHEdgePerpEpics->Fill(fEpicsCount%(fHEdgePerpEpics->GetNbinsX()-1),0);
	}
      }
      if((fPlane==ETablePara)&&(fHEdgeDistParaEpics)){
	fHEdgeDistParaEpics->Fill(fA2LinPolCohEdgeEpics);
      }
      else if((fPlane==ETablePerp)&&(fHEdgeDistPerpEpics)){
	fHEdgeDistPerpEpics->Fill(fA2LinPolCohEdgeEpics);
      }
      
      
      fIsNewFile=EFalse;                              //no longer a new file
      fEpicsCount++;                                  //increment the no of epics events
    } 
    // end of epics section
    
    if(gAR->IsScalerRead()&&(!fIsEpics)){            //if scaler read, and no epics events
      fIsNewFile=EFalse;	                           //no longer a new file
    }

    if((fDoingScalers)&&(gAR->IsScalerRead())){ 	   // only do this bit for scaler events
      
      //Fill the various arrays in ascending E_g order
      for(int n=0;n<fTaggerChannels;n++){	           //fill various hists (in ascending E_g order)
	fAccScaler[fScalerEvent%fNScalerBuffers][fTaggerChannels-1-n]=fScalerCurr[n];
      }
      
      for(int n=0;n<fTaggerChannels;n++){
	fCohSpectrum[n]=0.0;
      }
      for(int b=0;b<fNScalerBuffers;b++){
	for(int n=0;n<fTaggerChannels;n++){
	  fCohSpectrum[n]+=fAccScaler[b][n];
	  coh_sum+=fAccScaler[b][n];
	}
      }
      fScalerEvent++;
      //Fill the various arrays in ascending E_g order
      //      for(int n=0;n<fTaggerChannels;n++){	           //fill various hists (in ascending E_g order)
      //	fCohSpectrum[fTaggerChannels-1-n]=fScalerCurr[n];
      //	coh_sum+=fScalerCurr[n];
      //}
      
      for(int v=fNormChannel-10;v<fNormChannel+10;v++){
	if((fCohSpectrum[v]>1.0)&&(fIncSpectrum[v]>1.0)){
	  normValue+=100.0/(fCohSpectrum[v]/fIncSpectrum[v]);
	  ncount++;
	}
      }
      normValue/=ncount;
      
      if(fHCoh!=NULL)fHCoh->Reset("ICE");	//These have to be reset (ie zeroed) first
      if(fHEnh!=NULL)fHEnh->Reset("ICE");
      if((fHCohPara!=NULL)&&(fPlane==ETablePara))fHCohPara->Reset("ICE");
      if((fHEnhPara!=NULL)&&(fPlane==ETablePara))fHEnhPara->Reset("ICE");
      if((fHCohPerp!=NULL)&&(fPlane==ETablePerp))fHCohPerp->Reset("ICE");
      if((fHEnhPerp!=NULL)&&(fPlane==ETablePerp))fHEnhPerp->Reset("ICE");
      
      for(int n=0;n<fTaggerChannels;n++){		    		//fill various hists
	if(fHCoh!=NULL){
	  fHCoh->Fill(fEnergyBins[n],fCohSpectrum[n]);		//main raw one
	}
	if((fHCohPara!=NULL)&&(fPlane==ETablePara)){		//para raw one
	  fHCohPara->Fill(fEnergyBins[n],fCohSpectrum[n]);	
	}
	if((fHCohPerp!=NULL)&&(fPlane==ETablePerp)){		//perp raw one
	  fHCohPerp->Fill(fEnergyBins[n],fCohSpectrum[n]);	
	}
	
	if(fHaveIncScaler){		//Now normalise if we have reference
	  if((fCohSpectrum[n]<1.0)||(fIncSpectrum[n]<1.0)||(fBadScalerChan[n])){
	    if(n==0)fEnhSpectrum[n]=0;
	    else fEnhSpectrum[n]=  fEnhSpectrum[n-1];
	  }
	  else{
	    fEnhSpectrum[n]=normValue*fCohSpectrum[n]/fIncSpectrum[n];
	  }
	  
	  //	if(fHEnh!=NULL)	fHEnh->Fill(fEnergyBins[n],fEnhSpectrum[n]);	
	  if(fHEnh!=NULL)	fHEnh->Fill(fEnergyBins[n],fEnhSpectrum[n]);	
	  
	  
	  switch(fPlane){	//now fill hists etc according to mode
	  case ETablePara:
	    if(fHEnhPara!=NULL){	//para enhancement one
	      fHEnhPara->Fill(fEnergyBins[n],fEnhSpectrum[n]);	
	    }
	    break;
	  case ETablePerp:
	    if(fHEnhPerp!=NULL){	//para enhancement one
	      fHEnhPerp->Fill(fEnergyBins[n],fEnhSpectrum[n]);	
	    }
	    break;
	  default:
	    break;
	  }
	}
      }
      
      if((fHaveIncScaler)&&(fHEnh)){			//find the coherent edge
		 
		  if(fPlane == ETableAmo) fEdge = 0;
		  else
		  {
			   
	binx=0;
	ymax=0;
	for(int b=fHEnh->FindBin(fEdgeMin);b<fHEnh->FindBin(fEdgeMax);b++){
	  if(fHEnh->GetBinContent(b)>ymax){
	    ymax=fHEnh->GetBinContent(b);
	    binx=b;
	  }
	}
	xmax=fHEnh->GetBinCenter(binx);
	ymax=fHEnh->GetBinContent(binx);
	
	
	if(!edgeFit)edgeFit=new TF1("edgeFit",GausOnBase,0,100,4);
	
	edgeFit->SetRange(xmax,xmax+40.0);
	edgeFit->SetParameter(1,xmax);
	edgeFit->SetParameter(2,10.0);
	edgeFit->SetParameter(3,100.0);
	fHEnh->Fit(edgeFit,"QNR");
	
	maxgrad=0;
	for(fitedge=xmax+1;fitedge<xmax+39;
	    fitedge+=(60.0)/500.0){
	  if(fabs(edgeFit->Derivative(fitedge))>maxgrad){
	    maxgrad=fabs(edgeFit->Derivative(fitedge));
	    fEdge = fitedge;
	  }
	}
	}
	if(fHEdge){
	  if(fScalerCount%(fHEdge->GetNbinsX()-1)==0){
	    fHEdge->Reset("ICE");
	  }
	  fHEdge->Fill(fScalerCount%(fHEdge->GetNbinsX()-1),fEdge);
	}
	if(fHEdgePara){
	  if(fScalerCount%(fHEdgePara->GetNbinsX()-1)==0){
	    fHEdgePara->Reset("ICE");
	  }
	  if(fPlane==ETablePara){
	    fHEdgePara->Fill(fScalerCount%(fHEdgePara->GetNbinsX()-1),fEdge);
	  }
	  else{
	    fHEdgePara->Fill(fScalerCount%(fHEdgePara->GetNbinsX()-1),0);
	  }	  
	}
	if(fHEdgePerp){
	  if(fScalerCount%(fHEdgePerp->GetNbinsX()-1)==0){
	    fHEdgePerp->Reset("ICE");
	  }
	  if(fPlane==ETablePerp){
	    fHEdgePerp->Fill(fScalerCount%(fHEdgePerp->GetNbinsX()-1),fEdge);
	  }
	  else{
	    fHEdgePerp->Fill(fScalerCount%(fHEdgePerp->GetNbinsX()-1),0);
	  }
	}
	if((fPlane==ETablePara)&&(fHEdgeDistPara)){
	  fHEdgeDistPara->Fill(fEdge);
	}
	else if((fPlane==ETablePerp)&&(fHEdgeDistPerp)){
	  fHEdgeDistPerp->Fill(fEdge);
	}
	if(!fForceFit){                              // if not fitting enhancements 
	  fEdge=fA2LinPolCohEdgeEpics;               // get edge from epics
	}
	fScalerCount++;
	GetPolDegree(200.0);                         //force it to set up a table
	
      }
    }
  }

  else{
    if(fIsNewFile==ETrue) return;
    FillPolArray();
  }
}


void TA2LinearPolEpics::ParseMisc(char *line){	// read parameters in the setup file

  Char_t miscType[10];
  Char_t tempString[80];
  Char_t planeString[10];
  Int_t type=-1;
  Int_t plane=-1;
  Int_t rangeMin=-1;
  Int_t rangeMax=-1;

  // Use "wildcard" input line to read params for linear pol monitoring
  if(sscanf(line, "%s", miscType) != 1) PrintError( line, "Linear Pol  parameters" );

  type=Map2Key(miscType,kValidMiscParams);
  switch(type){
  case ELpMiscApp:
    if(sscanf(line, "%*s%lf%lf%lf%s%s%s%d%s%s%lf",&fNormEnergy,&fEdgeMin,&fEdgeMax,fTaggerName,fLadderName,
	      fRunRefFiles[0],&fNScalerBuffers,fEdgeString,fPlaneString,&fDeadband)!=10){
      PrintError( line, "Linear Pol app  parameters" );
    }
    fHaveApp=kTRUE;
    fDoingScalers = kTRUE;
    fHaveScaler   = kTRUE;
    if(strcmp(fEdgeString,"fit")==0){	// get edge from fit
      fForceFit=kTRUE;
    }
    else{
      fForceFit=kFALSE;
    }
    if(strcmp(fPlaneString,"file")==0){	// get plane from setup file
      fForcePlane=kTRUE;
    }
    else{
      fForcePlane=kFALSE;
    }
    if((fForceFit)&&(fForcePlane)){    //  in this case we don't use any epics
      fIsEpics=kFALSE;
    }
    fHaveCbrem = kTRUE;
    break;
    
  case ELpMiscRunRanges:
    if(sscanf(line, "%*s%d%d%s%s%lf%lf",&fRunRangeMin[fNRunRanges],&fRunRangeMax[fNRunRanges],
	      fRunRefFiles[fNRunRanges],planeString,&fEdgeSetting, &fEdgeRange)!=6){
      PrintError( line, "Linear Pol Run Range parameters" );
    }
    if(strcmp(planeString,"para")==0){
      fRunPlanes[fNRunRanges]=ETablePara;
    }
    else if(strcmp(planeString,"perp")==0){
      fRunPlanes[fNRunRanges]=ETablePerp;
    }
    else{
      fRunPlanes[fNRunRanges]=ETableAmo;
    }
    fNRunRanges++;
    
    break;
    
  case ELpMiscPolCuts:
    if(sscanf(line, "%*s%lf%lf%lf%lf",&fBeforeEdge,&fAfterEdge,&fPolMin,&fDeadband)!=4){
      PrintError( line, "Linear Pol PolLookupTable parameters" );
    }
    break;
    
  case ELpMiscPolLookupTable:
    //like this:
    //#Misc:   poltable  params       12345     12346 PARA 0.001790 0.000223 6.339087 0.008871 1557.000000 2.000000 6.045563 1.341427 0.000000 0.000000 4.271005
    //#Misc:   poltable  filename.dat 12345     12346
    
    if(sscanf(line, "%*s%d%d%s",&rangeMin,&rangeMax,tempString)!=3){
      PrintError( line, "Linear Pol PolLookupTable parameters" );
    }
    
    if(strcmp(tempString,"params")==0){
      sscanf(line,"%*s%*d%*d%*s%s",tempString);
      if(strcmp(tempString,"para")==0){
	plane=ETablePara;
      }
      else{
	plane=ETablePerp;
      }
      
      fPolTableTypes[fPolTableNRanges[plane]][plane]=ELpPolFromParams;
      fPolTableRangesMin[fPolTableNRanges[plane]][plane]=rangeMin;
      fPolTableRangesMax[fPolTableNRanges[plane]][plane]=rangeMax;
      
      sscanf(line,"%*s%*d%*d%*s%*s%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf",
	     &fPolTableParams[fPolTableNRanges[plane]][plane][0],
	     &fPolTableParams[fPolTableNRanges[plane]][plane][1],
	     &fPolTableParams[fPolTableNRanges[plane]][plane][2],
	     &fPolTableParams[fPolTableNRanges[plane]][plane][3],
	     &fPolTableParams[fPolTableNRanges[plane]][plane][4],
	     &fPolTableParams[fPolTableNRanges[plane]][plane][5],
	     &fPolTableParams[fPolTableNRanges[plane]][plane][6],
	     &fPolTableParams[fPolTableNRanges[plane]][plane][7],
	     &fPolTableParams[fPolTableNRanges[plane]][plane][8],
	     &fPolTableParams[fPolTableNRanges[plane]][plane][9],
	     &fPolTableParams[fPolTableNRanges[plane]][plane][10]);
    }
    else{
      fPolTableTypes[fPolTableNRanges[plane]][plane]=ELpPolFromFile;
      //do file stuff
    }
    
    fPolTableNRanges[plane]++;  //increment the no of ranges for this plane
    fHavePolTable = kTRUE;
    
    break;
  case ELpMiscPVNames:
    if(sscanf(line,"%*s%s%s%s%s%s%s%s%s",       //get read into these variables
	      fLinPolEpicsPVNames[0],           //&fA2LinPolCohEdgeEpics
	      fLinPolEpicsPVNames[1],           //&fA2LinPolCohRadiatorEpics
	      fLinPolEpicsPVNames[2],           //&fA2LinPolCohPlaneEpics
	      fLinPolEpicsPVNames[3],           //&fA2GoniAxis1Epics1
	      fLinPolEpicsPVNames[4],           //&fA2GoniAxis1Epics2
	      fLinPolEpicsPVNames[5],           //&fA2GoniAxis1Epics3
	      fLinPolEpicsPVNames[6],           //&fA2GoniAxis1Epics4
	      fLinPolEpicsPVNames[7])!=8){      //&fA2GoniAxis1Epics5
      PrintError( line, "Linear Pol pvname parameters" );
    }
    fHavePVNames = kTRUE;    
    break;
  case ELpMiscRads:
    if(sscanf(line,"%*s%d%d",&fDiamondID,&fAmoID)!=2){
      PrintError( line, "Linear Pol rads parameters" );
    }
    fHaveRads = kTRUE;
    break;
  default:
    break;
  }
}


int TA2LinearPolEpics::LoadAmoRef(Char_t *refFileName){
  FILE *fp;
  Char_t line[82];
  Int_t channel;
  Double_t content;

  if((fp=fopen(refFileName,"r"))!=NULL){
    channel=0;
    for(;;){                            		// scan lines from scaler dump file
      if(fgets(line,80,fp) == NULL) break;		// check got a line from file
      if((line[0] == '*')||(line[0] == '#')) continue;
      sscanf(line,"%*d%lf",&content);         // get chan and content
      
      if(channel>=fTaggerChannels) break; 		// break if past last chan
      fIncSpectrum[fTaggerChannels-1-channel]=content; //fill array
      channel++;
    } 						 //(increasing E_g order)
    fHaveIncScaler=ETrue;	// flag that we have the Inc. scaler reference
    fclose(fp);		// close the scaler dump file
    //run though the tagger scalers to find bad channels
    for(int n=0;n<fTaggerChannels-4;n++){ //flag if less that 1/9 of neighbours	
      if((18*fIncSpectrum[n])<(fIncSpectrum[n+2]+fIncSpectrum[n+3]+fIncSpectrum[n+4])){
	fBadScalerChan[n]=kTRUE;
      }
      else{
	fBadScalerChan[n]=kFALSE;
      }
      if(fIncSpectrum[n]<10.0)fBadScalerChan[n]=ETrue; // or less than 10.0
    }
    // work out the normalisation channel based on the norm energy
    
    for(int n=0;n<fTaggerChannels;n++){
      if(fEnergyCalib[n]>=fNormEnergy){
	fNormChannel=n;
	break;
      }
    }
    while(fIncSpectrum[fNormChannel--]<10.0); 
    fprintf(stderr,"Using %lf as the fNormEnergy\n",fEnergyCalib[fNormChannel]);
  }
  
  else {
    fprintf(stderr,"Warning: Couldn't open %s \n",refFileName);
    return -1;
  }
      
  return 0;
}

Double_t TA2LinearPolEpics::GetPolDegree(Double_t energy){
  Double_t pol;
  if(fIsNewFile) return -1;	//no meaningful pol info until 1st scaler read,
  if((fPlane!=ETablePara)&&(fPlane!=ETablePerp)) return -1; //no meanungful plane 
  // if((fEdge<(fLastEdge-fDeadband))||(fEdge>(fLastEdge+fDeadband))){
  //    fLastEdge=fEdge;
  //    enhFromParams();
  //}
  //return -1.0 if too far from edge.
  if((energy<fEdge-fBeforeEdge)||(energy>fEdge+fAfterEdge)) return -1;  
  pol = fHistP->GetBinContent(fHistP->FindBin(energy));
  //return -1 if below threshold
  if(pol<fPolMin) return -1.0;
  return pol;
}


Double_t *TA2LinearPolEpics::FillPolArray(){
  Int_t n=0;
  Int_t bin=-1;
  Int_t counter=-1;
  Double_t oldmean=-1.0;
  Double_t pol=-1.0;
  
  if(fIsNewFile) return NULL;	//no meaningful pol info until 1st scaler read,

  //if no valid plane or edge
  if(((fPlane!=ETablePara)&&(fPlane!=ETablePerp))|| ((fEdge<(fEdgeSetting-fEdgeRange))||(fEdge>(fEdgeSetting+fEdgeRange)))){
    while(fLadderHits[n]!=(Int_t)ENullHit){
      fPolArray[n]=-1.0;
      fPolArrayM[n++]=-1.0;
    }
    fPolArray[n]=(Double_t)ENullHit;
    return NULL;
  }
  
  if((fEdge<(fLastEdge-fDeadband))||(fEdge>(fLastEdge+fDeadband))){
    fLastEdge=fEdge;
    fBeforeEdgeBin=fHistP->FindBin(fEdge-fBeforeEdge); //find channels
    fAfterEdgeBin=fHistP->FindBin(fEdge+fAfterEdge);
    enhFromParams();
  }
  while(fLadderHits[n]!=(Int_t)ENullHit){
    //ladder bins go from 0->351 in reverse ord from hist bins, going from 1-352 
    bin=fHistP->GetNbinsX()-fLadderHits[n];
    if((bin<fBeforeEdgeBin)||(bin>fAfterEdgeBin)){
      fPolArray[n]=-1.0;
    } 
    else{
      pol=fHistP->GetBinContent(bin);
      if(pol<fPolMin){
	fPolArray[n]=-1.0;
      }
      else{
	fPolArray[n]=fHistP->GetBinContent(bin); //fill
	if((fHPolCount)&&(fHPolMean)){           //if the hists are there, fill the meam and count
	  counter = fHPolCount->GetBinContent(bin);
	  oldmean = fHPolMean->GetBinContent(bin);
	  fHPolMean->SetBinContent(bin, oldmean + ((fHistP->GetBinContent(bin) - oldmean)/(counter+1)));
	  fHPolCount->SetBinContent(bin,counter+1);
	}
      }
    }
    n++;
  }
  fPolArray[n]=(Double_t)ENullHit;
  
  return fPolArray; 
}

void  TA2LinearPolEpics::enhFromParams(){
  //make an enhancement and corresponding polarization from some the parameters as defined in the CLAS note.
  //this function is can be called stand alone, but will also ba called many times from the fitting function

  Double_t *par;
  Double_t xd[10];
  Double_t xc[10];
  Double_t Q[10];
  Double_t cohContrib;
  Double_t cohTotal;
  Double_t phiTotal;
  Double_t etotal;
  Double_t ptotal;
  Double_t x=0.0;
  Int_t    g=0;
  Double_t weight=0.0;
  Double_t weightSum=0.0;
  Double_t polSum=0.0;
  Double_t phi,chi,cd;
  Double_t amo;
  Int_t jbin=0;
  Double_t E0=fBeamEnergy;
  Double_t Eg=fEdge;

  par=fPolTableParams[fPolRangeIndex[fPlane]][fPlane];
  //
  //for(int p=0;p<11;p++){
  //  cout << p << ": " << par[p] << ", ";
  //  par[p]=fPar[p];
  // }
  //cout << endl;

  //get theta for the current edge position
  par[THETA]  = k/(2.0*E0*E0*((1/Eg)-(1/E0)));                                  //theta from edge and beam energy


  //reset them all for fresh filling
  fHistE->Reset("ICE");
  fHistP->Reset("ICE");
  fThetaPol->Reset("ICE");
  fThetaItot->Reset("ICE");
  fWeightHist->Reset("ICE");


  for(Double_t j=par[THETA]-3.0*par[SIGMA];j<=par[THETA]+3.001*par[SIGMA];j+=(6.0*par[SIGMA])/THETASTEPS){
    
    weight=TMath::Gaus(j,par[THETA],par[SIGMA]);   //get the weight from the gaussian
    weightSum+=weight;                             //add to sum      
    
    //find the discontinuity for each vector
    for(int v=0;v<par[NVEC];v++){
      g=VECTORS[v];
      xd[v]=1.0/((k/(g*par[E0MEV]*j))+1.0);
      Q[v]=(1.0-xd[v])/xd[v];
      xc[v]=xd[v]/(1+((par[THETAR]*par[THETAR])*(1-xd[v])));
    }

    //loop over all bins in the histogram
    for(int bin=1;bin<=fHistE->GetNbinsX();bin++){
      x=fHistE->GetBinCenter(bin)/par[E0MEV];            //find the value of the bin
      amo=1/x;                                    //assume amo = inc = 1/x over regio of interest
      
      cohTotal=0.0;
      phiTotal=0.0;
      
      //loop over all the vectors
      for(int v=0;v<par[NVEC];v++){
	if(x>xd[v]) continue;           //only do up to x_dg
	 
	//work out chi and phi
	phi=(2*Q[v]*Q[v]*x*x)/((1-x)*(1+((1-x)*(1-x))-((4*Q[v]*Q[v]*x*x/(1-x))*(((1-x)/(Q[v]*x))-1))));
	chi=((Q[v]*Q[v]*x)/(1-x))*(1+((1-x)*(1-x))-((4*Q[v]*Q[v]*x*x/(1-x))*(((1-x)/(Q[v]*x))-1)));
	//	cout << j  << "  " << chi << endl;
	cd=0.5*(1+TMath::Erf((x-xc[v])/(TMath::Sqrt(2)*par[SIGMAR])));

	//get coherent contrib for the vector
	cohContrib=cd*par[IVEC+v]*chi;

	//add to the total and update the phi total
	cohTotal+=cohContrib;
	phiTotal+=cohContrib*phi;

      }
      if(cohTotal>0.0) {
	phiTotal/=cohTotal;   //divide by the cohTotal to get the weighted dmean phi
	//cout << x << " " << phiTotal << " " << cohTotal << " " << weight << endl;	 
      }

      //enhancement = coherent total + inc (or amo).
      etotal=(amo+cohTotal)/amo;
      //and pol like this
      //      ptotal=phiTotal*cohTotal/(cohTotal + amo);
      ptotal=phiTotal*cohTotal;

      //add the weighted contribution to the enhancement
      fHistE->Fill(x*par[E0MEV],weight*etotal);

      //keep the pol for this x,theta coord
      fThetaPol->Fill(x*par[E0MEV],jbin,ptotal);

      //keep the total intensity for this x,theta coord
      fThetaItot->Fill(x*par[E0MEV],jbin,cohTotal+amo);
    }
    
    //save the weight for this theta point
    fWeightHist->Fill(jbin,weight);
    jbin++;

  }
  //normalize the sum of the weighted enhancements
  fHistE->Scale(1.0/weightSum);
  
  
  //loop over each x bin, adding the weighted contribs from each theta pos
  for(int bin=1; bin<=fHistP->GetNbinsX(); bin++){
    weightSum=0.0;
    polSum=0.0;
    
    for(int jb=1;jb<=fWeightHist->GetNbinsX();jb++){
      weight=fWeightHist->GetBinContent(jb);

      //      polSum+=fThetaPol->GetBinContent(bin,jb)*fThetaItot->GetBinContent(bin,jb)*weight;
      polSum+=fThetaPol->GetBinContent(bin,jb)*weight;
      weightSum+=fThetaItot->GetBinContent(bin,jb)*weight;
      //polSum+=fThetaPol->GetBinContent(bin,jb)*weight;
      //weightSum+=weight;
    }
    polSum/=weightSum;
    fHistP->Fill(fHistP->GetBinCenter(bin),polSum);
  } 
  //
    for(int n=0;n<=fTaggerChannels;n++){
      if(fHPolTableEnh) fHPolTableEnh->SetBinContent(n+1,fHistE->GetBinContent(n+1));
      if(fHPolTablePol)fHPolTablePol->SetBinContent(n+1,fHistP->GetBinContent(n+1));
      fCurrentPolTable[n]	=fHistP->GetBinContent(n+1);
      fCurrentEnhTable[n]	=fHistE->GetBinContent(n+1);      

      fCurrentPolTable_TC[fTaggerChannels-n]=fHistP->GetBinContent(n+1);
      fCurrentEnhTable_TC[fTaggerChannels-n]=fHistE->GetBinContent(n+1);
    }
}

ClassImp(TA2LinearPolEpics)
