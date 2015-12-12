//--Author	JRM Annand   27th Apr 2003
//--Rev 	JRM Annand...30th Sep 2003  New Detector bases
//--Rev 	JRM Annand...15th Oct 2003  ReadDecoded...MC data
//--Rev 	JRM Annand... 5th Feb 2004  3v8 compatible
//--Rev 	JRM Annand... 7th Jun 2005  ReadDecoded...use fEnergyScale
//--Rev 	JRM Annand...25th Oct 2005  ReadDecoded...energy thresh
//--Rev 	D.Glazier ...24th Aug 2007  ReadDecoded...include detector time
//--Update	JRM Annand ..17th Nov 2007  ReadDecoded total energy fix
//--Description
//                *** Acqu++ <-> Root ***
// Online/Offline Analysis of Sub-Atomic Physics Experimental Data
//
// TA2CalArray
//
// Decoding and calibration methods for the Crystall Ball array of NaI(Tl)
// Configured as forward wall to work in conjunction with the CB
// This can use standard TA2ClusterDetector methods or ones redifined here
//

#include "TA2CalArray.h"
#include <string>
#include <sstream>

#include <TH2CB.h>

enum
{
  ECAEnergyResolution=1900, ECATimeResolution, ECAThetaResolution, ECAPhiResolution,
  ECAOffsetTime, ECASmearThetaMin, ECASmearThetaMax, ECASmearCBThetaBoundary,
  ECAThetaSigma, ECASmearMin, ECASmearMax, ECASmearBoundaryMin, ECASmearBoundaryMax,
  ECASmearEnergyMax
};

static const Map_t kCalArrayKeys[] =
{
  {"Energy-Resolution:",   ECAEnergyResolution},
  {"Time-Resolution:",     ECATimeResolution},
  {"Theta-Resolution:",    ECAThetaResolution},
  {"Phi-Resolution:",      ECAPhiResolution},
  {"Offset-Time:",         ECAOffsetTime},
  {"MC-Smear-ThetaMin:",   ECASmearThetaMin},
  {"MC-Smear-ThetaMax:",   ECASmearThetaMax},
  {"MC-Smear-CBThetaBoundary:", ECASmearCBThetaBoundary},
  {"MC-Smear-ThetaSigma:", ECAThetaSigma},
  {"MC-SmearMin:",         ECASmearMin},
  {"MC-SmearMax:",         ECASmearMax},
  {"MC-SmearBoundaryMin:", ECASmearBoundaryMin},
  {"MC-SmearBoundaryMax:", ECASmearBoundaryMax},
  {"MC-Smear-EnergyMax:",  ECASmearEnergyMax},
  {NULL,            -1}
};

//---------------------------------------------------------------------------
TA2CalArray::TA2CalArray(const char* name, TA2System* apparatus)
            :TA2ClusterDetector(name, apparatus)
{
  // Do not allocate any "new" memory here...Root will wipe
  // Set private variables to zero/false/undefined state
  // Pass kLaddKeys (command-control keywords) and
  // kLaddHist (names of histograms) to progenitor classes

  fUseSigmaEnergy       = 0;
  fUseSigmaTime         = 0;
  fSigmaEnergyFactor    = -1.0;
  fSigmaEnergyPower     = -1.0;
  fSigmaTime            = -1.0;
  fOffsetTime           = 0.0;
  fSigmaTheta           = -1.0;
  fSigmaPhi             = -1.0;
  fEthresh              = 0.0;
  fSmearThetaMin        = 0.0;
  fSmearThetaMax        = 0.0;
  fSmearCBThetaBoundary = 0.0;
  fSmearThetaSigma      = 0.0;
  fSmearMin             = 0.0;
  fSmearMax             = 0.0;
  fSmearBoundaryMin     = 0.0;
  fSmearBoundaryMax     = 0.0;
  fSmearEnergyMax       = 0.0;

  fRandom = new TRandom();

  // defined in base class TA2ClusterDetector
  std::string s_name(GetName());
  std::string s_all = s_name + "_ClustersAll";
  fDispClusterHitsAll = new TH2CB(s_all, s_all);
  std::string s_energy = s_name + "_ClustersEnergy";
  fDispClusterHitsEnergy = new TH2CB(s_energy, s_energy);
  fDispClusterHitsSingle = new TH2Crystals*[MAX_DISP_CLUSTERS];
  for(int i=0;i<MAX_DISP_CLUSTERS;i++) {
    std::stringstream s_single; 
    s_single << s_name << "_ClustersSingle_" << i;
    fDispClusterHitsSingle[i] = new TH2CB(s_single.str(), s_single.str());
  }


  AddCmdList(kCalArrayKeys);                  // for SetConfig()
}

//---------------------------------------------------------------------------

TA2CalArray::~TA2CalArray()
{
  // Free up all allocated memory
  // ...arrays created at the initialisation stage
  // Start with common arrays of TA2Detector class
  DeleteClusterArrays();
}

//---------------------------------------------------------------------------

void TA2CalArray::SetConfig(char* line, int key)
{
  // Load CalArray parameters from file or command line
  // CalArray specific configuration
  switch(key)
  {
  case ECAEnergyResolution:
    // Energy Resolution Read-in Line
    if(sscanf(line, "%lf%lf%d", &fSigmaEnergyFactor, &fSigmaEnergyPower, &fUseSigmaEnergy) < 3)
      PrintError(line,"<TA2CalArray Energy Resolution>");
    break;
  case ECATimeResolution:
    // Time resolution read-in line
    if(sscanf(line, "%lf%d", &fSigmaTime, &fUseSigmaTime) < 2)
      PrintError(line,"<TA2CalArray Time Resolution>");
    break;
  case ECAThetaResolution:
    // Time resolution read-in line
    if(sscanf(line, "%lf", &fSigmaTheta) < 1)
      PrintError(line,"<TA2CalArray Theta Resolution>");
    break;
  case ECAPhiResolution:
    // Time resolution read-in line
    if(sscanf(line, "%lf", &fSigmaPhi) < 1)
      PrintError(line,"<TA2CalArray Phi Resolution>");
    break;
  case ECAOffsetTime:
    // Time offset read-in line
    if(sscanf(line, "%lf", &fOffsetTime) < 1)
      PrintError(line,"<TA2CalArray Time Offset>");
    break;
  case ECASmearThetaMin:
     if(sscanf(line, "%lf", &fSmearThetaMin) < 1)
       PrintError(line,"<TA2CalArray MC-Smear-ThetaMin>");
     break;
  case ECASmearThetaMax:
     if(sscanf(line, "%lf", &fSmearThetaMax) < 1)
       PrintError(line,"<TA2CalArray MC-Smear-ThetaMax>");
     break;
  case ECASmearCBThetaBoundary:
     if(sscanf(line, "%lf", &fSmearCBThetaBoundary) < 1)
       PrintError(line,"<TA2CalArray MC-Smear-CBThetaBoundary>");
     break;
  case ECAThetaSigma:
     if(sscanf(line, "%lf", &fSmearThetaSigma) < 1)
       PrintError(line,"<TA2CalArray MC-Smear-ThetaSigma>");
     break;
  case ECASmearMin:
      if(sscanf(line, "%lf", &fSmearMin) < 1)
        PrintError(line,"<TA2CalArray MC-SmearMin>");
      break;
  case ECASmearMax:
      if(sscanf(line, "%lf", &fSmearMax) < 1)
        PrintError(line,"<TA2CalArray MC-SmearMax>");
      break;
  case ECASmearBoundaryMin:
      if(sscanf(line, "%lf", &fSmearBoundaryMin) < 1)
        PrintError(line,"<TA2CalArray MC-SmearBoundaryMin>");
      break;
  case ECASmearBoundaryMax:
      if(sscanf(line, "%lf", &fSmearBoundaryMax) < 1)
        PrintError(line,"<TA2CalArray MC-SmearBoundaryMax>");
      break;
  case ECASmearEnergyMax:
      if(sscanf(line, "%lf", &fSmearEnergyMax) < 1)
        PrintError(line,"<TA2CalArray MC-Smear-EnergyMax>");
      break;

  default:
    // Command not found...possible pass to next config
    TA2ClusterDetector::SetConfig(line, key);
    break;
  }
  return;
}

//-----------------------------------------------------------------------------

void TA2CalArray::ParseDisplay(char* line)
{
  // Input private histogram spec to general purpose parse
  // and setup routine

  //  const Name2Variable_t hist[] = {
  // Name          ->variable  single/mult    Fill-condition
  //    {"Nphoton_Minv", fM_Nphoton, EHistSingleX},
  //    {NULL,          0,         0}
  //  };
  // Do not remove the final NULL line

  TA2ClusterDetector::ParseDisplay(line);
  return;
}

//---------------------------------------------------------------------------

void TA2CalArray::PostInit()
{
  // Some further initialisation after all setup parameters read in
  // Start with alignment offsets
  // Create space for various output arrays

  fEnergyAll = new Double_t[fNelem+1];

  TA2ClusterDetector::PostInit();
}

//---------------------------------------------------------------------------

void TA2CalArray::SaveDecoded()
{
  // Save decoded info to Root Tree file
}

//---------------------------------------------------------------------------

Double_t TA2CalArray::SmearClusterEnergy(double energy)
{
    return (energy += fRandom->Gaus(0.0, GetSigmaEnergyGeV(energy)));
}

Double_t TA2CalArray::SmearClusterEnergy(Double_t energy, std::vector<crystal_t> cluster)
{
    crystal_t center = cluster.front();
    double theta = center.Position.Theta()*TMath::RadToDeg();
//    const double theta_min = 20., theta_max = 160., theta_boundary = 28., theta_sigma = 2.5;
//    const double smear_min = 0.01, smear_max = 0.06; // max and min smearing value to be applied.
//    const double smear_boundary_min = 0.001, smear_boundary_max = 0.1, E_max = 1.2;
    double sigma, c; // smearing to be applied, decay constant

    // convert to GeV
    energy /= 1000.;

    // smear theta angle
    theta += fRandom->Gaus(0.0, fSmearThetaSigma);

    // "decay" constant to mimic the experimental resolution
    c = ( TMath::Log(fSmearMax/fSmearMin) )/( fSmearThetaMax-fSmearThetaMin );
    // calculate smearing value, mutliply by the minimum smearing value
    sigma = fSmearMin*TMath::Exp(c*(fSmearThetaMax-theta));

    // increase smearing closer to the CB boundary region
    if( (theta < fSmearCBThetaBoundary) && (energy < fSmearEnergyMax)){
        // linearly decreasing in E
        sigma += fSmearBoundaryMax - energy*(fSmearBoundaryMax - fSmearBoundaryMin)/fSmearEnergyMax;
    }

    sigma *= TMath::Power(energy, fSigmaEnergyPower);
    energy += fRandom->Gaus(0.0, sigma);
    energy *= 1000.;

    return energy;
}


//---------------------------------------------------------------------------

ClassImp(TA2CalArray)
