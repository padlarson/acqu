/*************************************************************************
 * Author: Dominik Werthmueller, 2008-2011
 *************************************************************************/

//////////////////////////////////////////////////////////////////////////
//                                                                      //
// MyHitClusterTAPS_t                                                   //
//                                                                      //
// Cluster reconstruction class for TAPS.                               //
// Based on HitCluster_t by JRM Annand.                                 //
//                                                                      //
//////////////////////////////////////////////////////////////////////////


#ifndef __MyHitClusterTAPS_t_h__
#define __MyHitClusterTAPS_t_h__

#include "TA2ClusterDetector.h"
#include "HitCluster_t.h"


class MyHitClusterTAPS_t : public HitCluster_t
{

public:
    MyHitClusterTAPS_t(Char_t* line, UInt_t index);
    virtual ~MyHitClusterTAPS_t() { }
    
    void ClusterDetermine(TA2ClusterDetector* cl);
};

#endif

