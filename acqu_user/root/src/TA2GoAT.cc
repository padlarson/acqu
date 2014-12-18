#include "TA2GoAT.h"


TA2GoAT::TA2GoAT(const char* Name, TA2Analysis* Analysis) : TA2AccessSQL(Name, Analysis),
                                                                    file(0),
                                                                    treeTracks(0),
                                                                    treeTagger(0),
                                                                    treeLinPol(0),
                                                                    treeTrigger(0),
                                                                    treeDetectorHits(0),
                                                                    treeScaler(0),
                                                                    nParticles(0),
                                                                    clusterEnergy(0),
                                                                    theta(0),
                                                                    phi(0),
                                                                    time(0),
                                                                    clusterSize(0),
                                                                    centralCrystal(0),
                                                                    centralVeto(0),
                                                                    apparatus(0),
                                                                    vetoEnergy(0),
                                                                    MWPC0Energy(0),
                                                                    MWPC1Energy(0),
                                                                    nTagged(0),
                                                                    taggedEnergy(0),
                                                                    taggedChannel(0),
                                                                    taggedTime(0),
                                                                    plane(0),
                                                                    edge(0),
                                                                    edgeSetting(0),
                                                                    nNaIHits(0),
                                                                    NaIHits(0),
                                                                    NaICluster(0),
                                                                    nPIDHits(0),
                                                                    PIDHits(0),
                                                                    nMWPCHits(0),
                                                                    MWPCHits(0),
                                                                    nBaF2PbWO4Hits(0),
                                                                    BaF2PbWO4Hits(0),
                                                                    BaF2PbWO4Cluster(0),
                                                                    nVetoHits(0),
                                                                    VetoHits(0),
                                                                    energySum(0),
                                                                    multiplicity(0),
                                                                    nTriggerPattern(0),
                                                                    triggerPattern(0),
                                                                    nHelicityBits(0),
                                                                    helicity(0),
                                                                    helicityInverted(0),
                                                                    helicityADC(0),
                                                                    nErrors(0),
                                                                    errorModuleID(0),
                                                                    errorModuleIndex(0),
                                                                    errorCode(0),
                                                                    eventNumber(0),
                                                                    eventID(0)														
{
    	strcpy(outputFolder,"");
    	strcpy(inputName,"");
    	strcpy(fileName,"Acqu");

    	AddCmdList(RootTreeConfigKeys);
}


TA2GoAT::~TA2GoAT()
{
    if(treeTracks)
        delete treeTracks;
	if(treeTagger)
		delete treeTagger;
	if(treeLinPol)
        delete treeLinPol;
	if(treeTrigger)
		delete treeTrigger;
	if(treeDetectorHits)
		delete treeDetectorHits;
	if(treeScaler)
		delete treeScaler;
    if(file)
		delete file;
}

void    TA2GoAT::LoadVariable()
{
	// Including histogram output for testing purposes (quick check of variables)
   	TA2AccessSQL::LoadVariable();

    TA2DataManager::LoadVariable("nParticles", 	  &nParticles,   EISingleX);
    TA2DataManager::LoadVariable("clusterEnergy", clusterEnergy, EDMultiX);
    TA2DataManager::LoadVariable("theta", 		  theta,         EDMultiX);
    TA2DataManager::LoadVariable("phi", 		  phi,           EDMultiX);
    TA2DataManager::LoadVariable("time", 		  time,          EDMultiX);

    TA2DataManager::LoadVariable("nTagged", 	  &nTagged,	     EISingleX);
    TA2DataManager::LoadVariable("taggedEnergy",  taggedEnergy,  EDMultiX);
    TA2DataManager::LoadVariable("taggedChannel", taggedChannel, EIMultiX);
    TA2DataManager::LoadVariable("taggedTime", 	  taggedTime,    EDMultiX);

    TA2DataManager::LoadVariable("vetoEnergy",    vetoEnergy,    EDMultiX);
    TA2DataManager::LoadVariable("MWPC0Energy",   MWPC0Energy,   EDMultiX);
    TA2DataManager::LoadVariable("MWPC1Energy",   MWPC1Energy,   EDMultiX);

    TA2DataManager::LoadVariable("energySum",     &energySum,    EDSingleX);

    return;
    
}

void    TA2GoAT::SetConfig(Char_t* line, Int_t key)
{
    	switch(key)
    	{
    	case EG_OUTPUT_FOLDER:
     	   strcpy(outputFolder,line);
     	   while(outputFolder[strlen(outputFolder)-1]=='\n' 
						|| outputFolder[strlen(outputFolder)-1]=='\r')
			outputFolder[strlen(outputFolder)-1]='\0';
        return;
    	case EG_INPUT_NAME:
        	strcpy(inputName,line);
        	while(inputName[strlen(inputName)-1]=='\n' 
						|| inputName[strlen(inputName)-1]=='\r')
			inputName[strlen(inputName)-1]='\0';
        return;
    	case EG_FILE_NAME:
        	strcpy(fileName,line);
        	while(fileName[strlen(fileName)-1]=='\n' 
						|| fileName[strlen(fileName)-1]=='\r')
			fileName[strlen(fileName)-1]='\0';
        return;
    	case EG_BEAM_HELICITY:
                nHelicityBits = sscanf(line, "%i%s%s%s%s%s%s%s%s", &helicityADC, helicityBits[0], helicityBits[1], helicityBits[2], helicityBits[3], helicityBits[4], helicityBits[5], helicityBits[6], helicityBits[7]);
                nHelicityBits--;
                if(nHelicityBits < 2) Error("SetConfig", "Not enough information to construct beam helicity!");
    	    	else
    	    	{
			printf("Helicity");
            for(Int_t i=0; i<nHelicityBits; i++)
			{
                    helicityInhibit[i] = false;
                    if(!strcmp(helicityBits[i],"I") || !strcmp(helicityBits[i],"i")) helicityInhibit[i] = true;
                    else if(!strcmp(helicityBits[i],"L") || !strcmp(helicityBits[i],"l")) helicityBeam[i] = false;
                    else if(!strcmp(helicityBits[i],"H") || !strcmp(helicityBits[i],"h")) helicityBeam[i] = true;
                printf(" - %s %i %i",helicityBits[i], helicityInhibit[i], helicityBeam[i]);
			}
			printf("\n");
    	    	}
	return;
    	default:
        	TA2AccessSQL::SetConfig(line, key);
    	}
}

void    TA2GoAT::PostInit()
{
    clusterEnergy    = new Double_t[TA2GoAT_MAX_PARTICLE];
    theta            = new Double_t[TA2GoAT_MAX_PARTICLE];
    phi              = new Double_t[TA2GoAT_MAX_PARTICLE];
    time             = new Double_t[TA2GoAT_MAX_PARTICLE];
    clusterSize      = new Int_t[TA2GoAT_MAX_PARTICLE];
    centralCrystal   = new Int_t[TA2GoAT_MAX_PARTICLE];
    centralVeto      = new Int_t[TA2GoAT_MAX_PARTICLE];
	
    taggedEnergy     = new Double_t[TA2GoAT_MAX_TAGGER];
    taggedChannel    = new Int_t[TA2GoAT_MAX_TAGGER];
    taggedTime       = new Double_t[TA2GoAT_MAX_TAGGER];
    
    apparatus        = new Int_t[TA2GoAT_MAX_PARTICLE];
    vetoEnergy       = new Double_t[TA2GoAT_MAX_PARTICLE];
    MWPC0Energy      = new Double_t[TA2GoAT_MAX_PARTICLE];
    MWPC1Energy      = new Double_t[TA2GoAT_MAX_PARTICLE];
        
    NaIHits	         = new Int_t[TA2GoAT_MAX_HITS];
    NaICluster       = new Int_t[TA2GoAT_MAX_HITS];
    PIDHits	         = new Int_t[TA2GoAT_MAX_HITS];
    MWPCHits		 = new Int_t[TA2GoAT_MAX_HITS];
    BaF2PbWO4Hits	 = new Int_t[TA2GoAT_MAX_HITS];
    BaF2PbWO4Cluster = new Int_t[TA2GoAT_MAX_HITS];
    VetoHits         = new Int_t[TA2GoAT_MAX_HITS];
    
    triggerPattern   = new Int_t[32];

    errorModuleID 	 = new Int_t[TA2GoAT_MAX_ERROR];
    errorModuleIndex = new Int_t[TA2GoAT_MAX_ERROR];
    errorCode        = new Int_t[TA2GoAT_MAX_ERROR];

	// Default SQL-physics initialisation
        TA2AccessSQL::PostInit();

   	printf("---------\n");
   	printf("Init Tree\n");
   	printf("---------\n");
    
   	// Append input filename to output tree name.
   	TString fullName;
   	if(gAR->GetProcessType() == EMCProcess) fullName = gAR->GetTreeFileList(0);        
   	else  fullName = gAR->GetFileName();
		
	while(fullName.Contains("/")) fullName.Remove(0,1+fullName.Index("/"));
	fullName.ReplaceAll(".dat",".root");
	if(strlen(inputName) && fullName.BeginsWith(inputName)) fullName.Remove(0,strlen(inputName));
	else fullName.Prepend("_");
	fullName.Prepend(fileName);
	if(!strcmp(outputFolder,"") && (gAR->GetTreeDir())) strcpy(outputFolder, gAR->GetTreeDir());
	if((strlen(outputFolder)>0) && strcmp(outputFolder+strlen(outputFolder)-1,"/")) strcat(outputFolder, "/");
	fullName.Prepend(outputFolder);
   	printf("Root file saved to %s\n", fullName.Data());  

    file		     = new TFile(fullName.Data(),"RECREATE");
    treeTracks       = new TTree("tracks",       "tracks");
    treeTagger       = new TTree("tagger",       "tagger");
    treeTrigger	     = new TTree("trigger",      "trigger");
    treeDetectorHits = new TTree("detectorHits", "detectorHits");
	
    treeTracks->Branch("nTracks", &nParticles, "nTracks/I");
    treeTracks->Branch("clusterEnergy", clusterEnergy, "clusterEnergy[nTracks]/D");
    treeTracks->Branch("theta", theta, "theta[nTracks]/D");
    treeTracks->Branch("phi", phi, "phi[nTracks]/D");
    treeTracks->Branch("time", time, "time[nTracks]/D");
    treeTracks->Branch("clusterSize", clusterSize, "clusterSize[nTracks]/I");
    treeTracks->Branch("centralCrystal", centralCrystal, "centralCrystal[nTracks]/I");
    treeTracks->Branch("centralVeto", centralVeto, "centralVeto[nTracks]/I");
    treeTracks->Branch("apparatus", apparatus, "apparatus[nTracks]/I");
    treeTracks->Branch("vetoEnergy", vetoEnergy, "vetoEnergy[nTracks]/D");
    treeTracks->Branch("MWPC0Energy", MWPC0Energy, "MWPC0Energy[nTracks]/D");
    treeTracks->Branch("MWPC1Energy", MWPC1Energy, "MWPC1Energy[nTracks]/D");
	
	treeTagger->Branch("nTagged", &nTagged,"nTagged/I");
    treeTagger->Branch("taggedEnergy", taggedEnergy, "taggedEnergy[nTagged]/D");
    treeTagger->Branch("taggedChannel", taggedChannel, "taggedChannel[nTagged]/I");
    treeTagger->Branch("taggedTime", taggedTime, "taggedTime[nTagged]/D");

    treeTrigger->Branch("energySum", &energySum, "energySum/D");
    treeTrigger->Branch("multiplicity", &multiplicity, "multiplicity/I");
    if(nHelicityBits > 1) treeTrigger->Branch("helicity", &helicity, "helicity/O");
    treeTrigger->Branch("nErrors", &nErrors, "nErrors/I");
    treeTrigger->Branch("errorModuleID", errorModuleID, "errorModuleID[nErrors]/I");
    treeTrigger->Branch("errorModuleIndex", errorModuleIndex, "errorModuleIndex[nErrors]/I");
    treeTrigger->Branch("errorCode", errorCode, "errorCode[nErrors]/I");
	treeTrigger->Branch("nTriggerPattern", &nTriggerPattern, "nTriggerPattern/I");
    treeTrigger->Branch("triggerPattern", triggerPattern, "triggerPattern[nTriggerPattern]/I");
	
    treeDetectorHits->Branch("nNaIHits", &nNaIHits, "nNaIHits/I");
    treeDetectorHits->Branch("NaIHits", NaIHits, "NaIHits[nNaIHits]/I");
    treeDetectorHits->Branch("NaICluster", NaICluster, "NaICluster[nNaIHits]/I");
    treeDetectorHits->Branch("nPIDHits", &nPIDHits, "nPIDHits/I");
    treeDetectorHits->Branch("PIDHits", PIDHits, "PIDHits[nPIDHits]/I");
    treeDetectorHits->Branch("nMWPCHits", &nMWPCHits, "nMWPCHits/I");
    treeDetectorHits->Branch("MWPCHits", MWPCHits, "MWPCHits[nMWPCHits]/I");
    treeDetectorHits->Branch("nBaF2PbWO4Hits", &nBaF2PbWO4Hits, "nBaF2PbWO4Hits/I");
    treeDetectorHits->Branch("BaF2PbWO4Hits", BaF2PbWO4Hits, "BaF2PbWO4Hits[nBaF2PbWO4Hits]/I");
    treeDetectorHits->Branch("BaF2PbWO4Cluster", BaF2PbWO4Cluster, "BaF2PbWO4Cluster[nBaF2PbWO4Hits]/I");
    treeDetectorHits->Branch("nVetoHits", &nVetoHits, "nVetoHits/I");
    treeDetectorHits->Branch("VetoHits", VetoHits, "VetoHits[nVetoHits]/I");

	// Store Scalers for non-MC process
	if (gAR->GetProcessType() != EMCProcess) 
	{
		treeScaler = new TTree("scaler", "scaler");	
		treeScaler->Branch("eventNumber", &eventNumber, "eventNumber/I");
		treeScaler->Branch("eventID", &eventID, "eventID/I");
		printf("GetMaxScaler: %d\n", GetMaxScaler());
	        Char_t str[256];
        sprintf(str, "scalers[%d]/i", GetMaxScaler());
        treeScaler->Branch("scalers", fScaler, str);

		// Store Lin Pol if class is active
		if(fLinPol)
		{
			treeLinPol = new TTree("linPol", "linPol");		
			treeLinPol->Branch("plane", &plane, "plane/I");
			treeLinPol->Branch("edge", &edge, "edge/D");
			treeLinPol->Branch("edgeSetting", &edgeSetting, "edgeSetting/D");
            treeLinPol->Branch("polarizationTable", fLinPol->GetPolTable_TC(), "polarizationTable[352]/D");
            treeLinPol->Branch("enhancementTable", fLinPol->GetEnhTable_TC(), "enhancementTable[352]/D");
		}
	}

	// Define Histograms which will be saved to root tree
	DefineHistograms();

	gROOT->cd();
	
	eventNumber	= 0;

   	printf("---------\n");
   	printf("Running\n");
   	printf("---------\n");	

}

void    TA2GoAT::Reconstruct()
{
	// Fill standard data check histograms
	DataCheckHistograms();

    DataCheckHistogramsAdlarson();

	// Output scaler info on scaler read events
	if((gAR->IsScalerRead()) && (gAR->GetProcessType() != EMCProcess))
	{
		eventID	= gAN->GetNDAQEvent();
		if(treeScaler) treeScaler->Fill();		
		
		if(fLinPol)
		{
			plane 	= fLinPol->GetPolPlane();
			edge 	= fLinPol->GetEdge();
			edgeSetting = fLinPol->GetEdgeSetting();
			if(treeLinPol) treeLinPol->Fill();
		}
	}

	nTagged = 0;
	if(fTagger && fLadder)
	{
                // Get the conversion of tagger channel to photon energy
	        Double_t electron_E = fTagger->GetBeamEnergy();
                const Double_t* ChToE = fLadder->GetECalibration();
	        Int_t fNmult = 1;
	        if ( gAR->IsOnline() ) fNmult = fLadder->GetNMultihit();

		if ( fNmult <= 1 )
		{
        		// Collect Tagger Hits without Multihits
        		nTagged	= fLadder->GetNhits();
                for(Int_t i=0; i<nTagged; i++)
        		{
                    taggedChannel[i]	= fLadder->GetHits(i);
                    taggedTime[i]	= (fLadder->GetTimeOR())[i];
                    taggedEnergy[i] = electron_E - ChToE[taggedChannel[i]];
	        	}
		}
		
		else
		{
        		// Collect Tagger Hits with Multihits
        		for(UInt_t m=0; m<(UInt_t)fNmult; m++)
        		{
        			for(UInt_t i=0; i<fLadder->GetNhitsM(m); i++)
        			{
                        taggedChannel[nTagged+i] 	= (fLadder->GetHitsM(m))[i];
                        taggedTime[nTagged+i]	= (fLadder->GetTimeORM(m))[i];
                        taggedEnergy[nTagged+i] = electron_E - ChToE[taggedChannel[nTagged+i]];
        			}
        			nTagged	+= fLadder->GetNhitsM(m);
	        	}
		}
	}
	
	// Gather particle information
	nParticles = 0;
	if(fCB)
	{
	// Collect CB Hits
    	nParticles	= fCB->GetNParticle();      
        for(Int_t i=0; i<nParticles; i++)
		{
			TA2Particle part = fCB->GetParticles(i);
			
			part.SetParticleID(kRootino); // Set mass to 0 (rootino)
			part.SetMass(0.0);

			// Reset a bunch of inconsistant "no-value" markers
            if(TMath::Abs(part.GetT()) >= TA2GoAT_NULL) clusterEnergy[i] = 0.0;
            else clusterEnergy[i] = part.GetT();

			if(TMath::Abs(part.GetTime()) >= TA2GoAT_NULL) time[i] = 0.0;
			else time[i] = part.GetTime();
			
            if(TMath::Abs(part.GetVetoEnergy()) >= TA2GoAT_NULL) vetoEnergy[i] = 0.0;
            else vetoEnergy[i]	= part.GetVetoEnergy();
			
            if(TMath::Abs(part.GetEnergyMwpc0()) >= TA2GoAT_NULL) MWPC0Energy[i] = 0.0;
            else MWPC0Energy[i] = part.GetEnergyMwpc0();
			
            if(TMath::Abs(part.GetEnergyMwpc1()) >= TA2GoAT_NULL) MWPC1Energy[i] = 0.0;
            else MWPC1Energy[i] = part.GetEnergyMwpc1();
						
            if(part.GetCentralIndex() == ENullHit) centralCrystal[i] = -1;
            else centralCrystal[i] = part.GetCentralIndex();
						
			if(part.GetVetoIndex() == ENullHit) centralVeto[i] = -1;
			else centralVeto[i]	= part.GetVetoIndex();		
			
			// Store other values which don't have this "no-value" option
            apparatus[i]	= (Int_t)EAppCB;
            theta[i]		= part.GetThetaDg();
            phi[i]			= part.GetPhiDg();
            clusterSize[i]  = part.GetClusterSize();

		}
	}

	if(fTAPS)
	{
		// Collect TAPS Hits
        for(Int_t i=0; i<fTAPS->GetNParticle(); i++)
		{
			TA2Particle part = fTAPS->GetParticles(i);
			
			part.SetParticleID(kRootino); // Set mass to 0 (rootino)
			part.SetMass(0.0);				

			// Reset a bunch of inconsistant "no-value" markers
            if(TMath::Abs(part.GetT()) >= TA2GoAT_NULL) clusterEnergy[nParticles+i] = 0.0;
            else clusterEnergy[nParticles+i] = part.GetT();

			if(TMath::Abs(part.GetTime()) >= TA2GoAT_NULL) time[nParticles+i] = 0.0;
			else time[nParticles+i] = part.GetTime();
			
            if(TMath::Abs(part.GetVetoEnergy()) >= TA2GoAT_NULL) vetoEnergy[nParticles+i] = 0.0;
            else vetoEnergy[nParticles+i]	= part.GetVetoEnergy();
		
            if(part.GetCentralIndex() == ENullHit) centralCrystal[nParticles+i] = -1;
            else centralCrystal[nParticles+i]	= part.GetCentralIndex();
			
			if(part.GetVetoIndex() == ENullHit) centralVeto[nParticles+i] = -1;
			else centralVeto[nParticles+i]	= part.GetVetoIndex();		
			
			// Set WC values to NULL
            MWPC0Energy[nParticles+i] = 0.0;
            MWPC1Energy[nParticles+i] = 0.0;
			
			// Store other values which don't have this "no-value" option
            apparatus[nParticles+i]		= (Int_t)EAppTAPS;
            theta[nParticles+i]			= part.GetThetaDg();
            phi[nParticles+i]			= part.GetPhiDg();
			time[nParticles+i]			= part.GetTime();	
            clusterSize[nParticles+i]  	= part.GetClusterSize();

		}
		nParticles += fTAPS->GetNParticle(); // update number of particles
	}

	UInt_t *clhits;
	HitCluster_t *cl;
	UInt_t *hits;
	Int_t clindex[720];

	// Get Detector Hits
	if(fNaI)
	{
        for(Int_t i=0; i<720; i++)
		{
		        clindex[i] = -1;
		}
                clhits = fNaI->GetClustHit();
        for(UInt_t i=0; i<fNaI->GetNCluster(); i++)
		{
		        cl = fNaI->GetCluster(clhits[i]);
			hits = cl->GetHits();
            for(UInt_t j=0; j<(cl->GetNhits()); j++)
			{
			        clindex[hits[j]] = i;
			}
		}
			
        nNaIHits = fNaI->GetNhits();
        for(Int_t i=0; i<nNaIHits; i++)
		{
            NaIHits[i] = fNaI->GetHits(i);
            NaICluster[i] = clindex[NaIHits[i]];
		}
	}

	if(fPID)
	{
        nPIDHits = fPID->GetNhits();
        for(Int_t i=0; i<nPIDHits; i++)
        { PIDHits[i] = fPID->GetHits(i); }
	}

	if(fMWPC)
	{
        nMWPCHits = fMWPC->GetNhits();
        for(Int_t i=0; i<nMWPCHits; i++)
        { MWPCHits[i] = fMWPC->GetHits(i); }
	}

	if(fBaF2PWO)
	{
            for(Int_t i=0; i<720; i++)
		{
		        clindex[i] = -1;
		}
                clhits = fBaF2PWO->GetClustHit();
        for(UInt_t i=0; i<fBaF2PWO->GetNCluster(); i++)
		{
		        cl = fBaF2PWO->GetCluster(clhits[i]);
			hits = cl->GetHits();
            for(UInt_t j=0; j<(cl->GetNhits()); j++)
			{
			        clindex[hits[j]] = i;
			}
		}
			
        nBaF2PbWO4Hits = fBaF2PWO->GetNhits();
        for(Int_t i=0; i<nBaF2PbWO4Hits; i++)
		{
            BaF2PbWO4Hits[i] = fBaF2PWO->GetHits(i);
            BaF2PbWO4Cluster[i] = clindex[BaF2PbWO4Hits[i]];
		}
	}

	if(fVeto)
	{
        nVetoHits = fVeto->GetNhits();
        for(Int_t i=0; i<nVetoHits; i++)
            { VetoHits[i] = fVeto->GetHits(i);}
	}
	
	// Get Trigger information
	TriggerReconstruction();

    nErrors = gAR->GetHardError();
	ReadErrorMk2_t *ErrorBlock = gAR->GetHardwareError();
	ReadErrorMk2_t *Error;
    for(Int_t i=0; i<nErrors; i++)
	{
		Error = ErrorBlock + i;
        errorModuleID[i] = Error->fModID;
        errorModuleIndex[i] = Error->fModIndex;
        errorCode[i] = Error->fErrCode;
	}

    if(nHelicityBits > 1)
	{
        Bool_t helicityBit;
        helicity = true;
        helicityInverted = true;
        for(Int_t i=0; i<nHelicityBits; i++)
		{
            helicityBit = (fADC[helicityADC] & 1<<i);
            if(helicityInhibit[i] && helicityBit)
			{
                errorCode[nErrors] = 9;
                nErrors++;
				break;
			}
            else if(helicityInhibit[i]) continue;
            helicity = (helicity && (helicityBeam[i] == helicityBit));
            helicityInverted = (helicityInverted && (helicityBeam[i] != helicityBit));
            if(helicity == helicityInverted)
			{
                errorCode[nErrors] = 10;
                nErrors++;
				break;
			}
		}
	} 

        if(fMulti[400])
	{
        Int_t eventIDcheck = fMulti[400]->GetHit(0);
        for(Int_t i=1; i<23; i++)
		{
			if(eventIDcheck != fMulti[400]->GetHit(i))
			{
                errorCode[nErrors] = 11;
                nErrors++;
				break;
			}
		}
	}

	//Apply EndBuffer
    clusterEnergy[nParticles] = EBufferEnd;
    theta[nParticles]         = EBufferEnd;
    phi[nParticles]           = EBufferEnd;
    time[nParticles]          = EBufferEnd;
    MWPC0Energy[nParticles]   = EBufferEnd;
    MWPC1Energy[nParticles]   = EBufferEnd;
    vetoEnergy[nParticles] 	  = EBufferEnd;
    taggedChannel[nTagged] 	  = EBufferEnd;
    taggedTime[nTagged] 	  = EBufferEnd;
	
	//Fill Trees
    if(treeTracks) 	treeTracks->Fill();
	if(treeTagger)			treeTagger->Fill();
	if(treeTrigger)  		treeTrigger->Fill();
	if(treeDetectorHits)	treeDetectorHits->Fill();

	//increment event number
	eventNumber++;	
}

void	TA2GoAT::DefineHistograms()
{
	// Define new data check histograms
	Check_CBdE_E		= new TH2F("Check_CBdE_E", "dE_E (all CB clusters compared to PID hits)", 	400, 0, 400, 100, 0, 10);
	Check_CBPhiCorr 	= new TH2F("Check_CBPhiCorr","PID-NaI phi correlation (all CB clusters compared to PID hits)", 24,  0,  24, 180, -180, 180);

	Check_CBdE_E_1PID 	 = new TH2F("Check_CBdE_E_1PID", "dE_E (all CB clusters compared to single PID hits)", 	400, 0, 400, 100, 0, 10);
	Check_CBPhiCorr_1PID = new TH2F("Check_CBPhiCorr_1PID","PID-NaI phi correlation (all CB clusters compared to single PID hits)", 24,  0,  24, 180, -180, 180);

	Check_CBdE_E_pi0 	= new TH2F("Check_CBdE_E_pi0", 		"dE_E after pi0 identification", 	400, 0, 400, 100, 0, 10);
	Check_CBPhiCorr_pi0 = new TH2F("Check_CBPhiCorr_pi0",	"PID-NaI phi correlation after pi0 identification", 24,  0,  24, 180, -180, 180);

	Check_TAPSdE_E		= new TH2F("Check_TAPSdE_E", "dE_E (all TAPS clusters compared to Veto hits)", 	400, 0, 400, 100, 0, 10);
	Check_TAPSPhiCorr 	= new TH2F("Check_TAPSPhiCorr","TAPS-Veto phi correlation (all TAPS clusters compared to Veto hits)", 438,  0,  438, 180, -180, 180);

	Check_TAPSdE_E_1Veto	= new TH2F("Check_TAPSdE_E_1Veto", "dE_E (all TAPS clusters compared to single Veto hits)", 	400, 0, 400, 100, 0, 10);
	Check_TAPSPhiCorr_1Veto = new TH2F("Check_TAPSPhiCorr_1Veto","TAPS-Veto phi correlation (all TAPS clusters compared to single Veto hits)", 438,  0,  438, 180, -180, 180);

	Check_CBHits 		= new TH2F("Check_CBHits", 		"CB Hits by event number", 		10000,  0,  10000000, 720, 0, 720);	
	Check_CBADCHits 	= new TH2F("Check_CBADCHits", 	"CB ADC Hits by event number", 	10000,  0,  10000000, 720, 0, 720);
	Check_CBTDCHits 	= new TH2F("Check_CBTDCHits",	"CB TDC Hits by event number", 	10000,  0,  10000000, 720, 0, 720);

	Check_PIDHits 		= new TH2F("Check_PIDHits", 	"PID Hits by event number", 	10000,  0,  10000000, 24, 0, 24);	
	Check_PIDADCHits 	= new TH2F("Check_PIDADCHits", 	"PID ADC Hits by event number", 10000,  0,  10000000, 24, 0, 24);
	Check_PIDTDCHits 	= new TH2F("Check_PIDTDCHits",	"PID TDC Hits by event number", 10000,  0,  10000000, 24, 0, 24);

	Check_TAPSHits 		= new TH2F("Check_TAPSHits", 	"TAPS Hits by event number", 	10000,  0,  10000000, 438, 0, 438);	
	Check_TAPSADCHits 	= new TH2F("Check_TAPSADCHits", "TAPS ADC Hits by event number",10000,  0,  10000000, 438, 0, 438);
	Check_TAPSTDCHits 	= new TH2F("Check_TAPSTDCHits", "TAPS TDC Hits by event number",10000,  0,  10000000, 438, 0, 438);

	Check_VetoHits 		= new TH2F("Check_VetoHits", 	"Veto Hits by event number", 	10000,  0,  10000000, 438, 0, 438);	
	Check_VetoADCHits 	= new TH2F("Check_VetoADCHits", "Veto ADC Hits by event number",10000,  0,  10000000, 438, 0, 438);
	Check_VetoTDCHits 	= new TH2F("Check_VetoTDCHits", "Veto TDC Hits by event number",10000,  0,  10000000, 438, 0, 438);


    // PA define histograms
    // NaI
    IMgg_vs_CBnr        = new TH2F("IMgg_vs_CB", "IM(gg) vs CB det nr", 1000, 0, 1000, 720, 0, 720);
    IMgg_vs_Eg          = new TH2F("IMgg_vs_Eg", "IM(gg) vs Eg", 1000, 0, 1000, 1000, 0, 1000);
    Eg_vs_CBnr          = new TH2F("Eg_vs_CB", "Eg vs CB det nr", 1000, 0, 1000, 720, 0, 720);
    Th_vs_CBnr          = new TH2F("Th_vs_CBnr", "Theta versus det nr", 720, 0, 180, 720, 0, 720);
    IMgg_etapi0         = new TH1F("IMgg_etapi0", "IMgg for eta and pi0, p eta pi0", 1000, 0, 1000);
    IMetapi0cand        = new TH1F("IMetapi0cand", "IM(4g) candidate p eta pi0", 200, 200, 1200);
    MM4g                = new TH1F("MM4g", "Missing mass calc for 4g",200, 0, 2000);

    // TAPS
    TAPSn               = new TH1F("TAPSn", "Nr of neutrals TAPS", 10, 0, 10);
    IMgg_vs_TAPS1n      = new TH2F("IMgg_vs_TAPS1n", "IM(gg) vs TAPS det nr", 1000, 0, 1000, 500, 0, 500);
    EdEprTAPS           = new TH2F("EdEprTAPS", "EdE proton candidate", 100, 0, 1000, 50, 0, 10);
    EdEnTAPS            = new TH2F("EdEnTAPS", "EdE gamma candidate", 100, 0, 1000, 50, 0, 10);



	
}

void	TA2GoAT::WriteHistograms()
{
	// Write check histograms
	file->mkdir("CheckCB");
	file->cd("CheckCB");
	Check_CBdE_E->Write();
	Check_CBPhiCorr->Write();	
	
	Check_CBdE_E_1PID->Write();
	Check_CBPhiCorr_1PID->Write();

	Check_CBdE_E_pi0->Write();
	Check_CBPhiCorr_pi0->Write();

	file->cd();
	file->mkdir("CheckTAPS");
	file->cd("CheckTAPS");
	Check_TAPSdE_E->Write();
	Check_TAPSPhiCorr->Write();	
	
	Check_TAPSdE_E_1Veto->Write();
	Check_TAPSPhiCorr_1Veto->Write();
	
	file->cd();
	file->mkdir("CheckHitPatterns");
	file->cd("CheckHitPatterns");	
	// Zoom in on all the plots with eventNumber on x-Axis
	Check_CBHits->GetXaxis()->SetRangeUser(0,eventNumber);
	Check_CBHits->Write();
	
	Check_CBADCHits->GetXaxis()->SetRangeUser(0,eventNumber);
	Check_CBADCHits->Write();
	
	Check_CBTDCHits->GetXaxis()->SetRangeUser(0,eventNumber);
	Check_CBTDCHits->Write();

	Check_PIDHits->GetXaxis()->SetRangeUser(0,eventNumber);
	Check_PIDHits->Write();
	
	Check_PIDADCHits->GetXaxis()->SetRangeUser(0,eventNumber);
	Check_PIDADCHits->Write();
	
	Check_PIDTDCHits->GetXaxis()->SetRangeUser(0,eventNumber);
	Check_PIDTDCHits->Write();

	Check_TAPSHits->GetXaxis()->SetRangeUser(0,eventNumber);
	Check_TAPSHits->Write();
	
	Check_TAPSADCHits->GetXaxis()->SetRangeUser(0,eventNumber);
	Check_TAPSADCHits->Write();
	
	Check_TAPSTDCHits->GetXaxis()->SetRangeUser(0,eventNumber);
	Check_TAPSTDCHits->Write();	
	
	Check_VetoHits->GetXaxis()->SetRangeUser(0,eventNumber);
	Check_VetoHits->Write();
	
	Check_VetoADCHits->GetXaxis()->SetRangeUser(0,eventNumber);
	Check_VetoADCHits->Write();
	
	Check_VetoTDCHits->GetXaxis()->SetRangeUser(0,eventNumber);
    Check_VetoTDCHits->Write();


    file->cd();
    file->mkdir("AdlarsonChecks");
    file->cd("AdlarsonChecks");

    IMgg_vs_CBnr->Write();
    Eg_vs_CBnr->Write();
    IMgg_vs_Eg->Write();
    IMgg_etapi0->Write();
    IMetapi0cand->Write();
    MM4g->Write();
    TAPSn->Write();
    EdEprTAPS->Write();
    EdEnTAPS->Write();


	file->cd();

}

void 	TA2GoAT::DataCheckHistograms()
{
	
	if((fCB) && (fPID))
	{
		// Find the charged particle (if possible) and fill dE_E for this
		// Fill the charged particle phi against the PID hit
        for(Int_t i=0; i<fCB->GetNParticle(); i++)
		{
			TA2Particle part = fCB->GetParticles(i);

			if (part.GetVetoEnergy() < 0)  continue;
			if (part.GetVetoEnergy() > 10) continue;			

			for(UInt_t j=0; j<fPID->GetNhits(); j++)
			{
				Check_CBdE_E->Fill(part.GetT(),part.GetVetoEnergy());
				Check_CBPhiCorr->Fill(fPID->GetHits(j),part.GetPhiDg());
			}

			if(fPID->GetNhits() == 1)
			{	
				Check_CBdE_E_1PID->Fill(part.GetT(),part.GetVetoEnergy());
				Check_CBPhiCorr_1PID->Fill(fPID->GetHits(0),part.GetPhiDg());	
			}
		}		
		
		// Do basic pi0 analysis on events with 3 NaI cluster and 1 PID hit
		// Fill dE_E for charged particle if pi0 reconstructed
		// Fill associated NaI cluster phi against single PID element		
		if( (fCB->GetNParticle() == 3) && (fPID->GetNhits() == 1) )	
		{
			// Get 3 particles, check for a pi0 from possible pairs
			TA2Particle part1 = fCB->GetParticles(0);
			TA2Particle part2 = fCB->GetParticles(1);
			TA2Particle part3 = fCB->GetParticles(2);
			
			TLorentzVector p4_12 = part1.GetP4() + part2.GetP4();
			TLorentzVector p4_13 = part1.GetP4() + part3.GetP4();
			TLorentzVector p4_23 = part2.GetP4() + part3.GetP4();
			
			// Calculated IM width/20 MeV (acts as a rough IM cut)
			Double_t IM_12 = TMath::Abs(135 - p4_12.M())/20;
			Double_t IM_13 = TMath::Abs(135 - p4_13.M())/20;
			Double_t IM_23 = TMath::Abs(135 - p4_23.M())/20;
			
			if((IM_12 <= 1.0) && (part3.GetVetoEnergy() > 0) && (part3.GetVetoEnergy() < 10))
			{
				Check_CBdE_E_pi0->Fill(part3.GetT(),part3.GetVetoEnergy());
				Check_CBPhiCorr_pi0->Fill(fPID->GetHits(0),part3.GetPhiDg());				
			}
			if((IM_13 <= 1.0) && (part2.GetVetoEnergy() > 0) && (part2.GetVetoEnergy() < 10))
			{
				Check_CBdE_E_pi0->Fill(part2.GetT(),part2.GetVetoEnergy());
				Check_CBPhiCorr_pi0->Fill(fPID->GetHits(0),part2.GetPhiDg());				
			}			
			if((IM_23 <= 1.0) && (part1.GetVetoEnergy() > 0) && (part1.GetVetoEnergy() < 10))
			{
				Check_CBdE_E_pi0->Fill(part1.GetT(),part1.GetVetoEnergy());
				Check_CBPhiCorr_pi0->Fill(fPID->GetHits(0),part1.GetPhiDg());				
			}			
			
		}
	}

	if((fTAPS) && (fVeto))
	{
		// Find the charged particle (if possible) and fill dE_E for this
		// Fill the charged particle phi against the Veto hit
        for(Int_t i=0; i<fTAPS->GetNParticle(); i++)
		{
			TA2Particle part = fTAPS->GetParticles(i);

			if (part.GetVetoEnergy() < 0)  continue;
			if (part.GetVetoEnergy() > 10) continue;			

			for(UInt_t j=0; j<fVeto->GetNhits(); j++)
			{
				Check_TAPSdE_E->Fill(part.GetT(),part.GetVetoEnergy());
				Check_TAPSPhiCorr->Fill(fVeto->GetHits(j),part.GetPhiDg());
			}

			if(fVeto->GetNhits() == 1)
			{	
				Check_TAPSdE_E_1Veto->Fill(part.GetT(),part.GetVetoEnergy());
				Check_TAPSPhiCorr_1Veto->Fill(fVeto->GetHits(0),part.GetPhiDg());	
			}
		}
	}
	
	// Fill CB stability histograms
	if(fNaI)
	{
		for(UInt_t i=0; i<fNaI->GetNhits(); i++)   
				{ Check_CBHits->Fill(eventNumber,fNaI->GetHits(i)); }
			
		if (fNaI->IsRawHits())
		{
			for (UInt_t i=0; i< fNaI->GetNADChits(); i++)
				{ Check_CBADCHits->Fill(eventNumber,fNaI->GetRawEnergyHits()[i]); }	
				
			for (UInt_t i=0; i< fNaI->GetNTDChits(); i++)
				{ Check_CBTDCHits->Fill(eventNumber,fNaI->GetRawTimeHits()[i]);	}	
							
		}
	}
	
	// Fill PID stability histograms
	if(fPID)
	{
		for(UInt_t i=0; i<fPID->GetNhits(); i++)   
				{ Check_PIDHits->Fill(eventNumber,fPID->GetHits(i)); }
			
		if (fPID->IsRawHits())
		{
			for (UInt_t i=0; i< fPID->GetNADChits(); i++)
				{ Check_PIDADCHits->Fill(eventNumber,fPID->GetRawEnergyHits()[i]); }	
				
			for (UInt_t i=0; i< fPID->GetNTDChits(); i++)
				{ Check_PIDTDCHits->Fill(eventNumber,fPID->GetRawTimeHits()[i]);	}	
							
		}
	}	
	
	// Fill TAPS stability histograms
	if(fBaF2PWO)
	{
		for(UInt_t i=0; i<fBaF2PWO->GetNhits(); i++)   
				{ Check_TAPSHits->Fill(eventNumber,fBaF2PWO->GetHits(i)); }
			
		if (fBaF2PWO->IsRawHits())
		{
			if (!(fADC[0] & 1<<15)) // Ignore TAPS Pulser reads
			{			
				for (UInt_t i=0; i< fBaF2PWO->GetNADChits(); i++)
					{ Check_TAPSADCHits->Fill(eventNumber,fBaF2PWO->GetRawEnergyHits()[i]); }	
					
				for (UInt_t i=0; i< fBaF2PWO->GetNTDChits(); i++)
					{ Check_TAPSTDCHits->Fill(eventNumber,fBaF2PWO->GetRawTimeHits()[i]);	}	
			}
		}
	}	
	
	// Fill PID stability histograms
	if(fVeto)
	{
		for(UInt_t i=0; i<fVeto->GetNhits(); i++)   
				{ Check_VetoHits->Fill(eventNumber,fVeto->GetHits(i)); }
			
		if (fVeto->IsRawHits())
		{
			if (!(fADC[0] & 1<<15)) // Ingore TAPS Pulser reads
			{
				for (UInt_t i=0; i< fVeto->GetNADChits(); i++)
					{ Check_VetoADCHits->Fill(eventNumber,fVeto->GetRawEnergyHits()[i]); }	

				for (UInt_t i=0; i< fVeto->GetNTDChits(); i++)
					{ Check_VetoTDCHits->Fill(eventNumber,fVeto->GetRawTimeHits()[i]);	}	
			}			
		}
	}
			
}




void 	TA2GoAT::DataCheckHistogramsAdlarson()
{
    if( fCB->GetNParticle() == 4 )
    {

        Int_t j = Reconstruct4g();
        if( j != 10 )
        {
            TA2Particle part1 = fCB->GetParticles(perm4g[j][0]);
            TA2Particle part2 = fCB->GetParticles(perm4g[j][1]);
            TA2Particle part3 = fCB->GetParticles(perm4g[j][2]);
            TA2Particle part4 = fCB->GetParticles(perm4g[j][3]);

            TLorentzVector p4_12 = part1.GetP4() + part2.GetP4();
            TLorentzVector p4_34 = part3.GetP4() + part4.GetP4();

            for( Int_t j = 0; j < nTagged; j++ )
            {
                TLorentzVector beam(0.,0., photonbeam_E[j], photonbeam_E[j]);
                TLorentzVector target(0.,0.,0.,938.272);

                MM4g->Fill((beam + target - p4_12 - p4_34).M());

            }


            IMetapi0cand->Fill((p4_12+p4_34).M());

            IMgg_vs_CBnr->Fill( p4_12.M() , part1.GetCentralIndex() );
            IMgg_vs_CBnr->Fill( p4_12.M() , part2.GetCentralIndex() );
            IMgg_vs_CBnr->Fill( p4_34.M() , part3.GetCentralIndex() );
            IMgg_vs_CBnr->Fill( p4_34.M() , part4.GetCentralIndex() );

            IMgg_vs_Eg->Fill( p4_12.M(), part1.GetE());
            IMgg_vs_Eg->Fill( p4_12.M(), part2.GetE());
            IMgg_vs_Eg->Fill( p4_34.M(), part3.GetE());
            IMgg_vs_Eg->Fill( p4_34.M(), part4.GetE());

            Eg_vs_CBnr->Fill( part1.GetE() , part1.GetCentralIndex() );
            Eg_vs_CBnr->Fill( part2.GetE() , part2.GetCentralIndex() );
            Eg_vs_CBnr->Fill( part3.GetE() , part3.GetCentralIndex() );
            Eg_vs_CBnr->Fill( part4.GetE() , part4.GetCentralIndex() );

            IMgg_etapi0->Fill(p4_12.M());
            IMgg_etapi0->Fill(p4_34.M());

            if( fTAPS )
                EdEprTAPS->Fill(part1.GetT(),part1.GetVetoEnergy());


        }
    }

/*    if( fTAPS->GetNParticle() == 2  && fCB->GetNParticle() == 3  )
    {
//        Int_t nTaps = 0;

          TA2Particle part1 = fTAPS->GetParticles(1);
          TA2Particle part2 = fTAPS->GetParticles(2);
          if(part1.GetVetoEnergy() > part2.GetVetoEnergy() )
          {
            EdEprTAPS->Fill(part1.GetT(),part1.GetVetoEnergy());
            EdEnTAPS->Fill(part2.GetT(),part2.GetVetoEnergy());
          }
          else
          {
            EdEprTAPS->Fill(part2.GetT(),part2.GetVetoEnergy());
            EdEnTAPS->Fill(part1.GetT(),part1.GetVetoEnergy());
          }



//        TAPSn->Fill(nTaps);
    }
*/

}


Int_t   TA2GoAT::Reconstruct4g()
{
    Double_t        Chi_Sq, Chi_Sq_etapi, Chi_Sq_2pi;
    Double_t        ChiSq_min = 10000;
    Int_t           best_comb = 10;

    for(int i = 0; i < 6; i++)
    {
        // Run through all possible permutations which can form pi0 and eta.
        // Select the best combination via chi2 test.
        // Check that the best combination better fits the etapi0 hypothesis than 2pi0
        // If good candidate best comb i is returned. If failed candidate i = 10

        TA2Particle part1 = fCB->GetParticles(perm4g[i][0]);
        TA2Particle part2 = fCB->GetParticles(perm4g[i][1]);
        TA2Particle part3 = fCB->GetParticles(perm4g[i][2]);
        TA2Particle part4 = fCB->GetParticles(perm4g[i][3]);

        TLorentzVector p4_12 = part1.GetP4() + part2.GetP4();
        TLorentzVector p4_34 = part3.GetP4() + part4.GetP4();


        Chi_Sq = TMath::Abs(p4_12.M() - 135.0)/20 + TMath::Abs(p4_34.M() - 548.0)/40;
        if( Chi_Sq < ChiSq_min )
        {
            best_comb = i;
            ChiSq_min = Chi_Sq;
        }

    }

    if( best_comb != 10)
    {
        TA2Particle part1 = fCB->GetParticles(perm4g[best_comb][0]);
        TA2Particle part2 = fCB->GetParticles(perm4g[best_comb][1]);
        TA2Particle part3 = fCB->GetParticles(perm4g[best_comb][2]);
        TA2Particle part4 = fCB->GetParticles(perm4g[best_comb][3]);

        TLorentzVector p4_12 = part1.GetP4() + part2.GetP4();
        TLorentzVector p4_34 = part3.GetP4() + part4.GetP4();


        Chi_Sq_etapi    = TMath::Abs(p4_12.M() - 135.0)/20 + TMath::Abs(p4_34.M() - 548.0)/40;
        Chi_Sq_2pi      = TMath::Abs(p4_12.M() - 135.0)/20 + TMath::Abs(p4_34.M() - 135.0)/20;

        if( Chi_Sq_2pi < Chi_Sq_etapi )
            best_comb = 10;

    }

    return best_comb;
}

Int_t		TA2GoAT::perm4g[6][4]=
{
    {0,1,2,3},
    {2,3,0,1},

    {0,2,1,3},
    {1,3,0,2},

    {0,3,1,2},
    {1,2,0,3}
};


void    TA2GoAT::Finish()
{
	printf("------------------\n");
	printf("Write Tree to file\n");
	printf("------------------\n");
	
	file->cd();
	
    if(treeTracks)
	{
        treeTracks->Write();	// Write
        delete treeTracks; 	// Close and delete in memory
	}
	if(treeTagger) 
	{
		treeTagger->Write();	// Write	
		delete treeTagger; 	// Close and delete in memory
	}	
	if(treeLinPol) 
	{
		treeLinPol->Write();	// Write	
		delete treeLinPol; 		// Close and delete in memory
	}		
	if(treeTrigger) 
	{
		treeTrigger->Write();	// Write	
		delete treeTrigger; 	// Close and delete in memory
	}		
	if(treeDetectorHits) 
	{
		treeDetectorHits->Write();// Write	
		delete treeDetectorHits;  // Close and delete in memory
	}		
	if(treeScaler)
	{
		treeScaler->Write();	// Write	
		delete treeScaler; 	// Close and delete in memory
    	}

	WriteHistograms();
	
    	if(file) 
		delete file;		// Close and delete in memory

	
	TA2AccessSQL::Finish();
}

void    TA2GoAT::ParseMisc(char* line)
{
	TA2AccessSQL::ParseMisc(line);
}

void 	TA2GoAT::TriggerReconstruction()
{
    if (fNaI) energySum = fNaI->GetTotalEnergy();
	if(gAR->GetProcessType() == EMCProcess) TriggerMC();
	else TriggerHW();

	
}

void 	TA2GoAT::TriggerMC() 
{
	// really rough, just the basic idea 
	// A good example is done in TA2BasePhysics but requires some extra work
	// Also need some flag for new and old 
	// Set some basic discriminator thresh for now
	Double_t DiscTh = 5.0; 
    multiplicity = 0;
		
	if(fNaI) 
	{
		for (Int_t i = 0; i < 45; i++) 
		{ 
			Bool_t flag = kFALSE;
			for (Int_t j = 0; j < 16; j++) 
			{
				if ((fNaI->GetEnergyAll(i*16 + j)) >= DiscTh) flag = kTRUE;
			}
            if (flag == kTRUE) multiplicity++;
		}
	}
		
	if (fBaF2PWO) {
		for (Int_t i = 0; i < 6; i++) { 
			Bool_t flag = kFALSE;
			for (Int_t j = 12; j < 71; j++) {  
			// really add some check of how many crystals are used, skip PbWO4s
				if ((fBaF2PWO->GetEnergyAll(i*71 + j)) >= DiscTh) flag = kTRUE;
			}
            if (flag == kTRUE) multiplicity++;
		}
	}
		
	// True Hardware multiplicities store only M2, M3, M4+
	//  Reduce MC multiplicity to reflect this limitation
    if (multiplicity == 1) multiplicity = 0;
    if (multiplicity > 4) multiplicity = 4;
}
	
void 	TA2GoAT::TriggerHW() 
{
	nTriggerPattern = 0;
    for (Int_t i= 0; i < 16; i++)
	{
		if (fADC[0] & 1<<i) 
		{ 
            triggerPattern[nTriggerPattern] = i;
			nTriggerPattern++;
		}
		if (fADC[1] & 1<<i) 
		{ 
            triggerPattern[nTriggerPattern] = i+16;
			nTriggerPattern++;
		}
	}
	
    multiplicity = 0;
	
    if (fADC[0] & 1<<11) multiplicity+=2;
    if (fADC[0] & 1<<12) multiplicity++;
    if (fADC[0] & 1<<13) multiplicity++;
 	
}

ClassImp(TA2GoAT)
