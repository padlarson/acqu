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

#include "TA2ClusterDetector.h"

// Command-line key words which determine what to read in
static const Map_t kClustDetKeys[] = {
  {"Max-Cluster:",          EClustDetMaxCluster},
  {"Next-Neighbour:",       EClustDetNeighbour},
  {"Moliere-Radius:",       EClustDetMoliereRadius},
  {NULL,          -1}
};

#include <iostream>
#include <set>
#include <map>
#include <algorithm>
#include <sstream>
#include <TA2Analysis.h>

using namespace std;

//---------------------------------------------------------------------------
TA2ClusterDetector::TA2ClusterDetector( const char* name,
          TA2System* apparatus )
  :TA2Detector(name, apparatus)
{
  // Do not allocate any "new" memory here...Root will wipe
  // Set private variables to zero/false/undefined state
  // Add further commands to the default list given by TA2Detector
  AddCmdList( kClustDetKeys );
  fCluster = NULL;
  fClustHit = NULL;
  fNCluster =  fMaxCluster = 0;
  fNClustHitOR = NULL;
  fTheta = fPhi = fClEnergyOR = fClTimeOR = fClCentFracOR = fClRadiusOR = NULL;
  fEthresh = 0.0;

#ifdef WITH_A2DISPLAY
  fDispClusterEnable = kFALSE; // config stuff missing...
  // will be set by child class
  fDispClusterHitsSingle = NULL;
  fDispClusterHitsEnergy = NULL;
#endif
  fIsPos = ETrue;          // override standard detector, must have position
}


//---------------------------------------------------------------------------
void TA2ClusterDetector::Decode( )
{
  // Run basic TA2Detector decode, then cluster decode
  // and then histogram

  DecodeBasic();                    
  DecodeCluster();
#ifdef WITH_A2DISPLAY
  DisplayClusters();
#endif
}

//---------------------------------------------------------------------------
void TA2ClusterDetector::DecodeSaved( )
{
  // Run basic TA2Detector decode, then cluster decode
  // and then histogram

  ReadDecoded();                    
  DecodeCluster();
#ifdef WITH_A2DISPLAY
  DisplayClusters();
#endif
}

static Double_t calc_total_energy(const vector< crystal_t >& cluster) {
  Double_t energy = 0;
  for(size_t i=0;i<cluster.size();i++) {
    energy += cluster[i].Energy;
  }
  return energy;
}

static Double_t opening_angle(const crystal_t& c1, const crystal_t& c2)  {
  // use TMath::ACos since it catches some NaN cases due to double rounding issues
  return TMath::RadToDeg()*TMath::ACos(c1.Position.Unit()*c2.Position.Unit());
}

typedef struct {
  TVector3 Position;
  vector<Double_t> Weights;
} bump_t;

// use local Moliere radius (NaI=4.8, BaF2=3.4, PbWO4=2.2)
// given by config file
static void calc_bump_weights(const vector<crystal_t>& cluster, bump_t& bump) {
  Double_t w_sum = 0;
  for(size_t i=0;i<cluster.size();i++) {
    Double_t r = (bump.Position - cluster[i].Position).Mag();
    Double_t w = cluster[i].Energy*TMath::Exp(-2.0*r/cluster[i].MoliereRadius); // originally 2.5
    bump.Weights[i] = w;
    w_sum += w;
  }
  for(size_t i=0;i<cluster.size();i++) { 
    bump.Weights[i] /= w_sum; 
  }
}

static void update_bump_position(const vector<crystal_t>& cluster, bump_t& bump) {
  Double_t bump_energy = 0;
  for(size_t i=0;i<cluster.size();i++) {
    bump_energy += bump.Weights[i] * cluster[i].Energy;
  }
  TVector3 position(0,0,0);
  Double_t w_sum = 0;
  for(size_t i=0;i<cluster.size();i++) {
    Double_t energy = bump.Weights[i] * cluster[i].Energy;
    Double_t w = 4.0 + TMath::Log(energy / bump_energy);
    w = w<0 ? 0 : w; // ignore too small energy contributions
    position += w * cluster[i].Position;
    w_sum += w;
  }
  position *= 1.0/w_sum;
  bump.Position = position;
}

static bump_t merge_bumps(const vector< list<bump_t>::iterator > bumps) {
  bump_t bump = *(bumps[0]);
  for(size_t i=1;i<bumps.size();i++) {
    bump_t b = *(bumps[i]);
    for(size_t j=0;j<bump.Weights.size();j++) {
      bump.Weights[j] += b.Weights[j];
    }
  }
  // normalize
  for(size_t j=0;j<bump.Weights.size();j++) {
    bump.Weights[j] /= bumps.size();
  }
  return bump;
}

static void split_cluster(const vector<crystal_t>& cluster,
                          vector< vector<crystal_t> >& clusters) {
  
  // make Voting based on relative distance or energy difference
   
  Double_t totalClusterEnergy = 0;
  vector<UInt_t> votes;
  votes.resize(cluster.size(), 0);
  // start searching at the second highest energy (i>0 case in next for loop)
  // since we know that the highest energy always has a vote
  votes[0]++;
  for(size_t i=0;i<cluster.size();i++) {   
    totalClusterEnergy += cluster[i].Energy; // side calculation in this loop, but include i=0
    if(i==0)
      continue;
    
    // i>0 now...
    // for each crystal walk through cluster 
    // according to energy gradient    
    UInt_t currPos = i;
    bool reachedMaxEnergy = false;
    Double_t maxEnergy = 0;
    while(!reachedMaxEnergy) {
      // find neighbours intersection with actually hit clusters
      reachedMaxEnergy = true;
      vector<UInt_t> neighbours = cluster[currPos].NeighbourIndices;
      for(size_t j=0;j<cluster.size();j++) {                      
        for(UInt_t n=0;n<neighbours.size();n++) {
          if(neighbours[n] != cluster[j].Index)
            continue; // cluster element j not neighbour of element currPos, go to next
          Double_t energy = cluster[j].Energy;
          if(maxEnergy < energy) {
            maxEnergy = energy;
            currPos = j;
            reachedMaxEnergy = false;
          }
          break; // neighbour indices are unique, stop iterating
        }
      }
    }
    // currPos is now at max Energy
    votes[currPos]++;
  }
  
  // all crystals vote for highest energy 
  // so this cluster should not be splitted,
  // just add it to the clusters
  if(votes[0] == cluster.size()) {
    clusters.push_back(cluster);
    return;
  }
  
  // find the bumps (crystals voted for)
  // and init the weights
  list<bump_t> bumps;
  for(size_t i=0;i<votes.size();i++) {
    if(votes[i]==0)
      continue;
    // initialize the weights with the position of the crystal
    bump_t bump;
    bump.Position = cluster[i].Position;
    bump.Weights.resize(cluster.size(), 0);
    calc_bump_weights(cluster, bump);
    bumps.push_back(bump);
  }

  // as long as we have overlapping bumps
  Bool_t haveOverlap = kFALSE;
   
  do {  
    // converge the positions of the bumps
    UInt_t iterations = 0;
    list< bump_t > stable_bumps;
    const Double_t positionEpsilon = 0.01;
    while(!bumps.empty()) {
      for(list<bump_t>::iterator b=bumps.begin(); b != bumps.end(); ++b) {
        // calculate new bump position with current weights
        TVector3 oldPos = (*b).Position;
        update_bump_position(cluster, *b);
        Double_t diff = (oldPos - (*b).Position).Mag();
        // check if position is stable
        if(diff>positionEpsilon) {
          // no, then calc new weights with new position
          calc_bump_weights(cluster, *b);
          continue;
        }
        // yes, then save it and erase it from to-be-stabilized bumps
        stable_bumps.push_back(*b);
        b = bumps.erase(b);
      }
      // check max iterations
      iterations++;
      if(iterations<100)
        continue; // go on iterating...
      // ... or we iterated really long without convergence
      // discard the bumps which havn't converged
      bumps.clear();    
    }
    
    // do we have any stable bumps? 
    // Then just the use cluster as is
    if(stable_bumps.empty()) {
      clusters.push_back(cluster);
      return;
    } 
    
    // stable_bumps are now identified, form clusters out of it  
    // check if two bumps share the same crystal of highest energy
    
    // for each crystal, we distribute its energy according to the
    // relative "pull" from each bump
    const Double_t weightEpsilon = 1.0e-5;
    vector<Double_t> weights_sum(cluster.size(), 0);
    for(size_t i=0;i<cluster.size();i++) {
      for(list<bump_t>::iterator b=stable_bumps.begin(); b != stable_bumps.end(); ++b) {
        if((*b).Weights[i] <  weightEpsilon) // ignore really small weights
          continue;
        weights_sum[i] += (*b).Weights[i];
      }
    }
    
    typedef map<UInt_t, vector< list<bump_t>::iterator > > overlaps_t;
    overlaps_t overlaps; // index of highest energy crystal -> corresponding stable bumps
    haveOverlap = kFALSE;
    vector< vector<crystal_t> > potential_clusters;
    for(list<bump_t>::iterator b=stable_bumps.begin(); b != stable_bumps.end(); ++b) {
      vector<crystal_t> bump_cluster;
      for(size_t i=0;i<cluster.size();i++) {
        if((*b).Weights[i] <  weightEpsilon) // ignore really small weights
          continue;
        crystal_t crys = cluster[i];
        crys.Energy *= (*b).Weights[i] / weights_sum[i];
        bump_cluster.push_back(crys);
      }
      // sort also the bump_cluster by energy, get index of highest energy (first one after sorting)
      sort(bump_cluster.begin(), bump_cluster.end());
      UInt_t kmax = bump_cluster[0].Index;
      // check if overlaps contains something at kmax
      if(overlaps.find(kmax) != overlaps.end()) {
        haveOverlap = kTRUE;
      }
      // remember the bump at its highest energy     
      overlaps[kmax].push_back(b);
      // if there's no overlap detected so far, we save the potential cluster
      if(!haveOverlap)
        potential_clusters.push_back(bump_cluster);
    }
    
    // if no overlaps detected,
    // then potential clusters can be added to clusters
    if(!haveOverlap) {
      clusters.insert(clusters.end(), potential_clusters.begin(), potential_clusters.end());
      break; // = continue since haveOverlap is false
    } 
    
    // if overlaps detected, 
    // merge those overlapping bumps and stabilize position again 
    for(overlaps_t::iterator o=overlaps.begin(); o != overlaps.end(); ++o) {
      bumps.push_back(merge_bumps(o->second));
    }
    
  } while(haveOverlap);
  
}


static void build_cluster(list<crystal_t>& crystals, 
                           vector<crystal_t> &cluster) {
  list<crystal_t>::iterator i = crystals.begin();
  
  // start with initial seed list 
  vector<crystal_t> seeds;
  seeds.push_back(*i);
  
  // save i in the current cluster
  cluster.push_back(*i);
  // remove it from the candidates
  crystals.erase(i);
  
  while(seeds.size()>0) {
    // neighbours of all seeds are next seeds
    vector<crystal_t> next_seeds; 
    
    for(vector<crystal_t>::iterator seed=seeds.begin(); seed != seeds.end(); seed++) {
      // find intersection of neighbours and seed 
      for(list<crystal_t>::iterator j = crystals.begin() ; j != crystals.end() ; ) {
        bool foundNeighbour = false;
        for(UInt_t n=0;n<(*seed).NeighbourIndices.size();n++) {
          if((*seed).NeighbourIndices[n] != (*j).Index) 
            continue;    
          next_seeds.push_back(*j);
          cluster.push_back(*j);
          j = crystals.erase(j);
          foundNeighbour = true;
          // neighbours is a list of unique items, we can stop searching
          break;
        }
        // removal moves iterator one forward
        if(!foundNeighbour)
          ++j;
      }
    }
    // set new seeds, if any new found...
    seeds = next_seeds;    
  }
  
  // sort it by energy
  sort(cluster.begin(), cluster.end());
}



void TA2ClusterDetector::DecodeCluster( )
{
  
  // build the hit crystals from the available hits
  // including the neighbour information...
  list<crystal_t> crystals; 
  for(UInt_t i=0;i<fNhits;i++) {
    UInt_t idx = fHits[i];
    crystals.push_back(
          crystal_t(idx,
                    fEnergy[idx], 
                    *(fPosition[idx]),
                    fCluster[idx]->GetNNeighbour(),
                    fCluster[idx]->GetNeighbour(),
                    fCluster[idx]->GetMoliereRadius()
                    )
          );
      }
  // sort hits with highest energy first
  // makes the cluster building faster hopefully
  crystals.sort();
 
  vector< vector<crystal_t> > clusters; // contains the final truly split clusters
  while(crystals.size()>0) {
    vector<crystal_t> cluster; 
    build_cluster(crystals, cluster); // already sorts it by energy
    
    // check if cluster contains bumps
    // if not, cluster is simply added to clusters
    split_cluster(cluster, clusters);
  } 
  
 
  // calculate cluster properties and fill the many arrays of this detector... 
  // assume that each cluster is sorted by fEnergy
  fNCluster = 0;
  for(size_t i=0 ; i < clusters.size() ; i++) {
    vector<crystal_t> cluster = clusters[i];
    
    // kick out really large clusters already here...
    if(cluster.size()>(UInt_t)fCluster[0]->GetMaxHits())
      continue;
    
    fClEnergyOR[fNCluster] = calc_total_energy(cluster);
    // kick out low energetic clusters...
    if(fClEnergyOR[fNCluster]<fEthresh) 
      continue;
    
    // kmax is crystal index with highest energy in cluster
    // cluster is sorted, so it's the first crystal
    UInt_t kmax = cluster[0].Index; 
    fClustHit[fNCluster] = kmax;
    fClTimeOR[fNCluster] = fIsTime ? fTime[kmax] : 0; 
    fNClustHitOR[fNCluster] = cluster.size();
    fClCentFracOR[fNCluster] = cluster[0].Energy / fClEnergyOR[fNCluster];
         
    // go over cluster once more
    TVector3 weightedPosition(0,0,0);
    Double_t weightedSum = 0;
    Double_t weightedRadius = 0;
    Double_t avgRadius = 0;
    Double_t avgOpeningAngle = 0; 
    for(size_t j=0;j<cluster.size();j++) {
   
      // energy weighting, parameter W0=3.5 (should be checked)
      Double_t energy = cluster[j].Energy;
      Double_t wgtE = 4.0 + TMath::Log(energy / fClEnergyOR[fNCluster]); // log to base e
      wgtE = wgtE<0 ? 0 : wgtE; // no negative weights
      weightedSum += wgtE;
      
      // position/radius weighting
      TVector3 pos = cluster[j].Position;
      TVector3 diff = cluster[0].Position - pos;
      weightedPosition += pos * wgtE;
      weightedRadius += diff.Mag() * wgtE;
      
      avgRadius += diff.Mag();
      avgOpeningAngle += opening_angle(cluster[0], cluster[j]);
      }
    
    avgRadius /= cluster.size();    
    avgOpeningAngle /= cluster.size();
    
    weightedPosition *= 1.0/weightedSum;
  
    fPhi[fNCluster] = TMath::RadToDeg() * weightedPosition.Phi();
    fTheta[fNCluster] = TMath::RadToDeg() * weightedPosition.Theta();
    fClRadiusOR[fNCluster] = weightedRadius / weightedSum;
    
    fCluster[kmax]->SetFields(cluster,
                              fClEnergyOR[fNCluster],
                              fClCentFracOR[fNCluster],
                              fClTimeOR[fNCluster],
                              weightedPosition);
    
    fNCluster++;
    // just discard more clusters
    if(fNCluster==fMaxCluster)
      break;
    }
 
  // mark ends in those stupid c-arrays...
  fClustHit[fNCluster] = EBufferEnd;
  fTheta[fNCluster] = EBufferEnd;
  fPhi[fNCluster] = EBufferEnd;
  fNClustHitOR[fNCluster] = EBufferEnd;
  fClEnergyOR[fNCluster] = EBufferEnd;
  fClTimeOR[fNCluster] = EBufferEnd;
  fClCentFracOR[fNCluster] = EBufferEnd;
  fClRadiusOR[fNCluster] = EBufferEnd;
}

//---------------------------------------------------------------------------
void TA2ClusterDetector::Cleanup( )
{
  // end-of-event cleanup
  TA2Detector::Cleanup();
  //  would also clean up any cluster info here
  for(UInt_t i=0; i<fNCluster; i++) fCluster[ fClustHit[i] ]->Cleanup();
}

//---------------------------------------------------------------------------
TA2ClusterDetector::~TA2ClusterDetector()
{
  // Free up all allocated memory
  // ...arrays created at the initialisation stage
  // Start with common arrays of TA2Detector class
  DeleteClusterArrays();
}

//---------------------------------------------------------------------------
void TA2ClusterDetector::DeleteClusterArrays( )
{
  // Free up allocated memory
  DeleteArrays();
  for( UInt_t i=0; i<fMaxCluster; i++ ){
    if( fCluster[i] ) delete fCluster[i];
  }
  delete fCluster;
}

//---------------------------------------------------------------------------
void TA2ClusterDetector::SetConfig( char* line, int key )
{
  // Load Ladder parameters from file or command line
  // Keywords which specify a type of command can be found in
  // the kLaddKeys array at the top of the source .cc file
  // The following are setup...
  //	1. # X-triggers, #elements in ladder array, global TDC offset
  //	2. X-trigger parameters if any, TDC, TDC calib, TDC cut
  //	3. Ladder parameters, TDC, calib, cut window, Eg calib, Scaler
  //	4. Any post initialisation
  //	5. Histograms...should be done after post-initialisation

  // Cluster specific configuration
  switch( key ){
  case EClustDetMaxCluster:
    // Max number of clusters
    if( sscanf( line, "%d%lf", &fMaxCluster, &fEthresh ) < 1 ){
      PrintError(line,"<Parse # clusters>");
      break;
    }
    fCluster = new HitCluster_t*[fNelement];
    fClustHit = new UInt_t[fMaxCluster+1];
    fNClustHitOR = new UInt_t[fMaxCluster+1];
    fTheta = new Double_t[fMaxCluster+1];
    fPhi      = new Double_t[fMaxCluster+1];
    fClEnergyOR  = new Double_t[fMaxCluster+1];
    fClTimeOR  = new Double_t[fMaxCluster+1];
    fClCentFracOR  = new Double_t[fMaxCluster+1];
    fClRadiusOR  = new Double_t[fMaxCluster+1];
    fNCluster = 0;
    break;
  case EClustDetNeighbour:
    // Nearest neighbout input
    if( fNCluster < fNelement )
      fCluster[fNCluster] =
          new HitCluster_t(line, fNCluster);
    fNCluster++;
    break;
  case EClustDetMoliereRadius:
    UInt_t i1, i2;
    Double_t moliere;
    int n;
    n = sscanf(line, "%lf%d%d", &moliere,&i1,&i2);
    if(n==1) {
      // indices not provided, assume global value
      i1 = 0;
      i2 = fNelement-1;
    }
    else if(n != 3) {
      PrintError(line,"<Moliere Radius format wrong>");
    break;
    }
    for(UInt_t i=i1;i<=i2;i++) {
      fCluster[i]->SetMoliereRadius(moliere);
    }
    break;
  default:
    // Command not found...try standard detector
    TA2Detector::SetConfig( line, key );
    break;;
  }
  return;
}

//---------------------------------------------------------------------------
void TA2ClusterDetector::PostInit()
{
  // Some further initialisation after all setup parameters read in
  // Start with alignment offsets
  // Create space for various output arrays
  TA2Detector::PostInit();
}

//-----------------------------------------------------------------------------
void TA2ClusterDetector::LoadVariable( )
{
  // Input name - variable pointer associations for any subsequent
  // cut or histogram setup
  // LoadVariable( "name", pointer-to-variable, type-spec );
  // NB scaler variable pointers need the preceeding &
  //    array variable pointers do not.
  // type-spec ED prefix for a Double_t variable
  //           EI prefix for an Int_t variable
  // type-spec SingleX for a single-valued variable
  //           MultiX  for a multi-valued variable

  //  TA2DataManager::LoadVariable("ClEnergy",  &fClEnergy,   EDSingleX);
  //  TA2DataManager::LoadVariable("ClTimeOR",  f???,         EDSingleX);
  //  TA2DataManager::LoadVariable("ClHits",    fClHits,           EDSingleX);
  //  TA2DataManager::LoadVariable("ClMulti",   fClMulti,         EDMultiX);
  TA2DataManager::LoadVariable("ClNhits",       &fNCluster,      EISingleX);
  TA2DataManager::LoadVariable("ClTotHits",     fClustHit,       EIMultiX);
  TA2DataManager::LoadVariable("ClTheta",       fTheta,          EDMultiX);
  TA2DataManager::LoadVariable("ClPhi",         fPhi,            EDMultiX);
  TA2DataManager::LoadVariable("ClNhitsOR",     fNClustHitOR,    EIMultiX);
  TA2DataManager::LoadVariable("ClEnergyOR",    fClEnergyOR,     EDMultiX);
  TA2DataManager::LoadVariable("ClTimeOR",      fClTimeOR,       EDMultiX);
  TA2DataManager::LoadVariable("ClCentFracOR",  fClCentFracOR,   EDMultiX);
  TA2DataManager::LoadVariable("ClRadiusOR",    fClRadiusOR,     EDMultiX);
  
  TA2Detector::LoadVariable();
}

//-----------------------------------------------------------------------------
void TA2ClusterDetector::ParseDisplay( char* line )
{
  // Input private histogram spec to general purpose parse
  // and setup routine

  // List of local cluster histograms
  const Map_t kClustDetHist[] = {
    {"ClEnergy",         EClustDetEnergy},
    {"ClTime",           EClustDetTime},
    {"ClCentFrac",       EClustDetCentFrac},
    {"ClRadius",         EClustDetRadius},
    {"ClHits",           EClustDetHits},
    {"ClMulti",          EClustDetMulti},
    {NULL,                  -1}
  };

  UInt_t i,j,k,l,chan;
  Double_t low,high;
  Char_t name[EMaxName];
  Char_t histline[EMaxName];

  // Do the rather specialist cluster display setup.
  // This doesn't fit too well with the Name2Variable_t model as
  // the cluster info is stored in separate cluster classes
  // If not specialist cluster display call standard detector display

  if( sscanf(line, "%s%s", histline, name) != 2 ){
    PrintError(line,"<Cluster display - bad format>");
    return;
  }
  i = Map2Key( name, kClustDetHist );
  switch( i ){
  case EClustDetEnergy:
  case EClustDetHits:
  case EClustDetTime:
  case EClustDetMulti:
    // Cluster spectra
    // cluster  chan bins, range low--high, elements j--k,
    if( sscanf(line, "%*s%*s%d%lf%lf%d%d", &chan,&low,&high,&j,&k) != 5 ){
      PrintError(line,"<Cluster display - bad format>");
      return;
    }
    for( l=j; l<=k; l++ ){
      if( l >= fNCluster ){
  PrintError(line,"<Cluster display - element outwith range>");
  return;
      }
      sprintf( histline, "%s%d %d  %lf %lf",name,l,chan,low,high );
      switch( i ){
      case EClustDetEnergy:
  Setup1D( histline, fCluster[l]->GetEnergyPtr() );
  break;
      case EClustDetHits:
  Setup1D( histline, fCluster[l]->GetHits(), EHistMultiX );
  break;
      case EClustDetTime:
  Setup1D( histline, fCluster[l]->GetTimePtr() );
  break;
      case EClustDetMulti:
  Setup1D( histline, fCluster[l]->GetNhitsPtr() );
  break;
      }
    }
    break;
  default:
    // Try the standard display if not specialist cluster stuff
    TA2HistManager::ParseDisplay(line);
    break;
  }

  fIsDisplay = ETrue;
  return;
}

#ifdef WITH_A2DISPLAY
//-----------------------------------------------------------------------------
void TA2ClusterDetector::DisplayClusters() {
  if(!fDispClusterEnable)
    return;

  for(UInt_t i=0;i<fNelement;i++) {
    fDispClusterHitsSingle->SetElement(i,0);    
    fDispClusterHitsEnergy->SetElement(i,0);       
  }
  
  for(UInt_t i=0;i<fNhits;i++) {
    fDispClusterHitsEnergy->SetElement(fHits[i],fEnergy[fHits[i]]);
  }
  
  for(UInt_t i=0;i<fNCluster;i++) {
    HitCluster_t* cl = fCluster[fClustHit[i]];
    UInt_t* hits = cl->GetHits();
    UInt_t nHits = cl->GetNhits();
    for(UInt_t j=0;j<nHits;j++) {
      Double_t val = fDispClusterHitsSingle->GetElement(hits[j]);
      val += 1<<i;
      if(j==0)
        val += 0.1*(i+1);
      fDispClusterHitsSingle->SetElement(hits[j],val);
    }
  }
  std::stringstream ss;
  ss << fDispClusterHitsSingle->GetName();
  ss << " Event=" << gAN->GetNEvent();
  ss << " Clusters=" << fNCluster;
  fDispClusterHitsSingle->SetTitle(ss.str().c_str());
  //usleep(5e5);
}
#endif

ClassImp(TA2ClusterDetector)
