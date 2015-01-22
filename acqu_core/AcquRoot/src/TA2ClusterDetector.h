//--Author	JRM Annand   30th Sep 2003  Read MC data
//--Rev 	JRM Annand...15th Sep 2003  Generalise methods  
//--Rev 	JRM Annand....9th Mar 2004  New OR variables
//--Rev 	JRM Annand....9th Mar 2005  protected instead of private vars
//--Rev 	JRM Annand...13th Jul 2005  split offs, time OR
//--Rev 	JRM Annand...25th Jul 2005  SetConfig hierarchy
//--Rev 	JRM Annand...20th Oct 2005  Iterative neighbour find (TAPS)
//--Rev 	JRM Annand... 9th Dec 2005  Change ambigous split-off thresh
//--Rev 	JRM Annand... 6th Feb 2006  Bug in SplitSearch
//--Rev 	JRM Annand...21st Apr 2006  Command-key enum to .h
//--Rev 	JRM Annand...22nd Jan 2007  4v0 update
//--Rev 	JRM Annand...12th May 2007  Central-frac and radius
//--Rev 	JRM Annand...18th Jan 2009  TMath::Sort (Root v5.22)
//--Update	JRM Annand   17th Sep 2011  log energy weighting
//--Description
//                *** Acqu++ <-> Root ***
// Online/Offline Analysis of Sub-Atomic Physics Experimental Data 
//
// TA2ClusterDetector
//
// Decoding and calibration methods for EM calorimeters or similar systems
// where a shower or showers of secondary particles fire a cluster or
// clusters neighbouring calorimeter elements, 
// e.g. Crystal Ball NaI(Tl) array
//

#ifndef __TA2ClusterDetector_h__
#define __TA2ClusterDetector_h__

#include "TA2Detector.h"
#include "HitCluster_t.h"       // hit cluster determination

#ifdef WITH_A2DISPLAY
#include "a2display.h"
#endif

// constants for command-line maps
enum { 
  EClustDetMaxCluster = 100, EClustDetNeighbour,  EClustDetAllNeighbour,
  EClustDetSplitOff, EClustDetIterate, EClustEnergyWeight,
  EClustDetEnergy, EClustDetTime, EClustDetCentFrac, EClustDetRadius,
  EClustDetHits, EClustDetMulti,
};

class TA2ClusterDetector : public TA2Detector {
 protected:
  HitCluster_t** fCluster;              // Clusters of hits
  UInt_t* fClustHit;                    // Cluster indices
  Bool_t* fIsSplit;                     // Indices split-off clusters
  UInt_t* fTempHits;                    // Element-Hit store
  UInt_t fNCluster;                     // # of clusters
  Int_t fClustSizeFactor;               // enlarge factor, hit cluster buffers
  UInt_t fNSplit;                       // # low energy clusters
  UInt_t fNSplitMerged;                 // # split-offs merged
  UInt_t fMaxCluster;                   // Max # of clusters
  UInt_t* fNClustHitOR;                 // OR of #hits in individuyal clusters
  Double_t* fTheta;                     // theta of cluster hit
  Double_t* fPhi;                       // phi of cluster hit
  Double_t* fClEnergyOR;                // OR of cluster energies
  Double_t* fClTimeOR;                  // OR of cluster times
  Double_t* fClCentFracOR;              // OR of energy ratios in central elem.
  Double_t* fClRadiusOR;                // OR E-weighted cluster radii
  Double_t fClEthresh;                  // threshold energy for main cluster
  Double_t fEthresh;                    // generic threshold energy for cluster
  Double_t fEthreshSplit;               // threshold energy for split-off
  Double_t fMaxThetaSplitOff;           // max split-off opening angle 
  Double_t fMinPosDiff;                 // min acceptable neighbour pos diff. 
  Double_t fMaxPosDiff;                 // max acceptable neighbour pos diff. 
  Double_t* fSplitAngle;                // opening angles detween clusters
  Double_t fEWgt;                       // energy weighting factor
  Int_t fLEWgt;                         // energy weighting factor switch
  Int_t* fISplit;                       // for sorting cluster opening angles
  Int_t* fIJSplit;                      // for sorting cluster opening angles
  Int_t fMaxSplitPerm;                  // for sorting cluster opening angles
  Bool_t fIsIterate;                    // cluster member find iteration ON/OFF
#ifdef WITH_A2DISPLAY
  Bool_t fDispClusterEnable;
  TH2Crystals* fDispClusterHitsSingle;
  TH2Crystals* fDispClusterHitsEnergy;  
  void DisplayClusters();
#endif
 public:
  TA2ClusterDetector( const char*, TA2System* );// Normal use
  virtual ~TA2ClusterDetector();
  virtual void SetConfig( char*, int );// decode and load setup info
  virtual void ParseDisplay( char* );  // decode histogram setup lines
  virtual void LoadVariable(  );       // name-variable association
  virtual void PostInit( );            // initialise using setup info
  virtual void Decode( );              // hits -> energy procedure
  virtual void DecodeCluster( );       // determine clusters
  virtual void DecodeSaved( );         // decode previously written data
  virtual void Cleanup( );             // end-of-event cleanup
  virtual void SaveDecoded( ) = 0;     // specialist
  virtual void ReadDecoded( ) = 0;     // specialist
  virtual void DeleteClusterArrays();  // flush cluster-specific arrays
  virtual void SplitSearch();          // look for split-off clusters

  HitCluster_t** GetCluster(){ return fCluster; }
  HitCluster_t* GetCluster( UInt_t i ){ return fCluster[i]; }
  UInt_t* GetClustHit(){ return fClustHit; }
  UInt_t* GetClustHit( UInt_t i ){ return fClustHit + i; }
  Bool_t* GetIsSplit(){ return fIsSplit; }
  Bool_t* GetIsSplit( UInt_t i ){ return fIsSplit + i; }
  UInt_t* GetTempHits(){ return fTempHits; }
  UInt_t GetNCluster(){ return fNCluster; }
  Int_t GetClustSizeFactor(){ return fClustSizeFactor; }
  UInt_t GetNSplit(){ return fNSplit; }
  UInt_t GetNSplitMerged(){ return fNSplitMerged; } 
  UInt_t GetMaxCluster(){ return fMaxCluster; }
  UInt_t* GetNClustHitOR(){ return fNClustHitOR; }
  Double_t* GetTheta(){ return fTheta; }
  Double_t* GetPhi(){ return fPhi; }
  Double_t* GetClEnergyOR(){ return fClEnergyOR; }
  Double_t* GetClTimeOR(){ return fClTimeOR; }
  Double_t* GetClCentFracOR(){ return fClCentFracOR; }
  Double_t* GetClRadiusOR(){ return fClRadiusOR; }
  Double_t GetClEthresh(){ return fClEthresh; }
  Double_t GetEthresh(){ return fEthresh; }
  Double_t GetEthreshSplit(){ return fEthreshSplit; }
  Double_t GetMaxThetaSplitOff(){ return fMaxThetaSplitOff; }
  Double_t GetMinPosDiff(){ return fMinPosDiff; } 
  Double_t GetMaxPosDiff(){ return fMaxPosDiff; }
  Double_t* GetSplitAngle(){ return fSplitAngle; }
  Double_t GetEWgt(){ return fEWgt; }
  Int_t GetLEWgt(){ return fLEWgt; }
  Int_t* GetISplit(){ return fISplit; }
  Int_t* GetIJSplit(){ return fIJSplit; }
  Int_t GetMaxSplitPerm(){ return fMaxSplitPerm; }
  Bool_t IsIterate(){ return fIsIterate; }
 
  ClassDef(TA2ClusterDetector,1)
};


#endif
