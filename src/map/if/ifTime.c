/**CFile****************************************************************

  FileName    [ifTime.c]

  SystemName  [ABC: Logic synthesis and verification system.]

  PackageName [FPGA mapping based on priority cuts.]

  Synopsis    [Computation of delay paramters depending on the library.]

  Author      [Alan Mishchenko]
  
  Affiliation [UC Berkeley]

  Date        [Ver. 1.0. Started - November 21, 2006.]

  Revision    [$Id: ifTime.c,v 1.00 2006/11/21 00:00:00 alanmi Exp $]

***********************************************************************/

#include "if.h"
#include "base/io/ioAbc.h"
#include "base/abc/abc.h"
#include "map/scl/sclLib.h"
#include <math.h>

ABC_NAMESPACE_IMPL_START

////////////////////////////////////////////////////////////////////////
///                        DECLARATIONS                              ///
////////////////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////////////////
///                     FUNCTION DEFINITIONS                         ///
////////////////////////////////////////////////////////////////////////

/**Function*************************************************************

  Synopsis    [Get coordinates by node name string using pNtkOrig + pNtkCoords.]

  Description [Wire-aware mapping: looks up the node with given blif name
  in the original (pre-strash) ABC network (pNtkOrig), then finds its
  coordinates from the coordinate network (pNtkCoords).

  Returns 1 on success, 0 if not found.]

  SideEffects []

  SeeAlso     []

***********************************************************************/
// Helper: get x coordinate of If_Obj by leaf ID
static float If_CutDelayGetX( If_Man_t * p, int LeafId, void * pCoords )
{
    if ( pCoords == NULL || p == NULL || p->vNodeNameMap == NULL )
        return 0.0f;
    void * pBlifNodePtr = p->vNodeNameMap[LeafId];
    if ( pBlifNodePtr == NULL )
        return 0.0f;
    Abc_Obj_t * pBlifNode = (Abc_Obj_t *)pBlifNodePtr;
    float x, y;
    if ( Io_ReadCoordsGetCoordByName( pCoords, Abc_ObjName(pBlifNode), &x, &y ) )
        return x;
    return 0.0f;
}

// Helper: get y coordinate of If_Obj by leaf ID
static float If_CutDelayGetY( If_Man_t * p, int LeafId, void * pCoords )
{
    if ( pCoords == NULL || p == NULL || p->vNodeNameMap == NULL )
        return 0.0f;
    void * pBlifNodePtr = p->vNodeNameMap[LeafId];
    if ( pBlifNodePtr == NULL )
        return 0.0f;
    Abc_Obj_t * pBlifNode = (Abc_Obj_t *)pBlifNodePtr;
    float x, y;
    if ( Io_ReadCoordsGetCoordByName( pCoords, Abc_ObjName(pBlifNode), &x, &y ) )
        return y;
    return 0.0f;
}

/**Function*************************************************************

  Synopsis    [Sorts the pins in the decreasing order of delays.]

  Description []
               
  SideEffects []

  SeeAlso     []

***********************************************************************/
void If_CutSortInputPins( If_Man_t * p, If_Cut_t * pCut, int * pPinPerm, float * pPinDelays )
{
    If_Obj_t * pLeaf;
    int i, j, best_i, temp;
    // start the trivial permutation and collect pin delays
    If_CutForEachLeaf( p, pCut, pLeaf, i )
    {
        pPinPerm[i] = i;
        pPinDelays[i] = If_ObjCutBest(pLeaf)->Delay;
    }
    // selection sort the pins in the decreasible order of delays
    // this order will match the increasing order of LUT input pins
    for ( i = 0; i < (int)pCut->nLeaves-1; i++ )
    {
        best_i = i;
        for ( j = i+1; j < (int)pCut->nLeaves; j++ )
            if ( pPinDelays[pPinPerm[j]] > pPinDelays[pPinPerm[best_i]] )
                best_i = j;
        if ( best_i == i )
            continue;
        temp = pPinPerm[i]; 
        pPinPerm[i] = pPinPerm[best_i]; 
        pPinPerm[best_i] = temp;
    }
/*
    // verify
    assert( pPinPerm[0] < (int)pCut->nLeaves );
    for ( i = 1; i < (int)pCut->nLeaves; i++ )
    {
        assert( pPinPerm[i] < (int)pCut->nLeaves );
        assert( pPinDelays[pPinPerm[i-1]] >= pPinDelays[pPinPerm[i]] );
    }
*/
}


/**Function*************************************************************

  Synopsis    [Computes delay.]

  Description []
               
  SideEffects []

  SeeAlso     []

***********************************************************************/
float If_CutDelay( If_Man_t * p, If_Obj_t * pObj, If_Cut_t * pCut )
{
    static int pPinPerm[IF_MAX_LUTSIZE];
    static float pPinDelays[IF_MAX_LUTSIZE];
    static int print_count = 0;
    static int print_limit = 20;
    char * pPerm = If_CutPerm( pCut );
    If_Obj_t * pLeaf;
    float Delay, DelayCur;
    float * pLutDelays;
    int i, Shift, Pin2PinDelay;//, iLeaf;
    Delay = -IF_FLOAT_LARGE;
    float cell_delay = -IF_FLOAT_LARGE;
    float wire_delay = 0.0f;
    float wirelength_um = 0.0f;
    float r_per_um = 0.0f, c_per_um = 0.0f;

    // ============================================================
    // Wire-aware mapping: cut-to-cut wirelength and Elmore delay
    //
    // Model: Each cut's output sits at the centroid of its fanin cells.
    //        Wire between cut_i and cut_j = Manhattan distance between
    //        their fanin centroids.
    //        wirelength = max over fanins: (fanin.LValue + dist(fanin_pos, centroid))
    //        wire_delay = Elmore: R * C * L^2 / 2
    //
    // vNodeNameMap[If_ObjId] -> blif node pointer (stable, built before DFS)
    // pObj->LValue  -> accumulated wirelength of this node's cut output (um)
    // pLeaf->LValue -> fanin's accumulated wirelength
    // ============================================================
    if ( p->pPars->WireDelay > 0.0 && p->pPars->pNtkCoords && p->vNodeNameMap )
    {
        void * pCoords = p->pPars->pNtkCoords;

        // Read wire RC
        Io_ReadCoordsGetWireRCFromCoords( pCoords, &r_per_um, &c_per_um );

        // ----- Phase 1: compute cut centroid from fanin cell positions -----
        // (cut internal wire is handled by cell_delay / tdelay)
        float sum_x = 0.0f, sum_y = 0.0f;
        int n_fanins = 0;

        If_CutForEachLeaf( p, pCut, pLeaf, i )
        {
            // Prefer precomputed If-node -> output-net name (scheme A).
            // Fallback to old logic for compatibility.
            void * pBlifNodePtr = p->vNodeNameMap ? p->vNodeNameMap[pLeaf->Id] : NULL;
            const char * pName = NULL;
            if ( p->vNodeOutputNetName )
                pName = (const char *)p->vNodeOutputNetName[pLeaf->Id];
            if ( (pName == NULL || pName[0] == '\0') && pBlifNodePtr )
            {
                Abc_Obj_t * pBlifNode = (Abc_Obj_t *)pBlifNodePtr;
                int origNodeId = Abc_ObjId( pBlifNode );
                if ( p->vNodeDrivingPoName && origNodeId >= 0 )
                    pName = (const char *)p->vNodeDrivingPoName[origNodeId];
            }
            // Fallback: use node name (may have backslash escaping)
            if ( (pName == NULL || pName[0] == '\0') && pBlifNodePtr )
            {
                Abc_Obj_t * pBlifNode = (Abc_Obj_t *)pBlifNodePtr;
                pName = Abc_ObjName( pBlifNode );
            }
            if ( pName == NULL || pName[0] == '\0' )
                continue;
            // Handle backslash escaping: "ctrl.state.out\[2\]" -> "ctrl.state.out[2]"
            if ( pName && strchr( pName, '\\' ) != NULL )
            {
                static char unesc[1000];
                int jj = 0;
                for ( int kk = 0; pName[kk] != '\0' && jj < 999; kk++ )
                    if ( pName[kk] != '\\' )
                        unesc[jj++] = pName[kk];
                unesc[jj] = '\0';
                pName = unesc;
            }
            float x, y;
            if ( Io_ReadCoordsGetCoordByName( pCoords, pName, &x, &y ) )
            {
                sum_x += x;
                sum_y += y;
                n_fanins++;
            }
        }

        if ( n_fanins == 0 )
        {
            // No coordinates: wirelength stays 0
            wirelength_um = 0.0f;
        }
        else
        {
            float centroid_x = sum_x / n_fanins;
            float centroid_y = sum_y / n_fanins;

            // ----- Phase 2: wirelength = max(fanin.LValue + dist(fanin, centroid)) -----
            // PhyLS formula: wirelength = max over fanins of
            //   (fanin_accumulated_wirelength + Manhattan distance)
            float worst_wl = 0.0f;

            If_CutForEachLeaf( p, pCut, pLeaf, i )
            {
                // Prefer precomputed If-node -> output-net name (scheme A).
                // Fallback to old logic for compatibility.
                void * pBlifNodePtr = p->vNodeNameMap ? p->vNodeNameMap[pLeaf->Id] : NULL;
                const char * pName = NULL;
                if ( p->vNodeOutputNetName )
                    pName = (const char *)p->vNodeOutputNetName[pLeaf->Id];
                if ( (pName == NULL || pName[0] == '\0') && pBlifNodePtr )
                {
                    Abc_Obj_t * pBlifNode = (Abc_Obj_t *)pBlifNodePtr;
                    int origNodeId = Abc_ObjId( pBlifNode );
                    if ( p->vNodeDrivingPoName && origNodeId >= 0 )
                        pName = (const char *)p->vNodeDrivingPoName[origNodeId];
                }
                if ( (pName == NULL || pName[0] == '\0') && pBlifNodePtr )
                {
                    Abc_Obj_t * pBlifNode = (Abc_Obj_t *)pBlifNodePtr;
                    pName = Abc_ObjName( pBlifNode );
                }
                if ( pName == NULL || pName[0] == '\0' )
                    continue;
                if ( pName && strchr( pName, '\\' ) != NULL )
                {
                    static char unesc2[1000];
                    int jj = 0;
                    for ( int kk = 0; pName[kk] != '\0' && jj < 999; kk++ )
                        if ( pName[kk] != '\\' )
                            unesc2[jj++] = pName[kk];
                    unesc2[jj] = '\0';
                    pName = unesc2;
                }
                float x, y;
                if ( Io_ReadCoordsGetCoordByName( pCoords, pName, &x, &y ) )
                {
                    float dist = fabsf( x - centroid_x ) + fabsf( y - centroid_y );
                    float wl = pLeaf->LValue + dist;
                    if ( wl > worst_wl )
                        worst_wl = wl;
                }
            }
            wirelength_um = worst_wl;

            // ----- Phase 3: propagate centroid position to pObj->LValue -----
            // pObj->LValue = wirelength of this cut's output = worst_wl
            pObj->LValue = wirelength_um;

            // ----- Phase 4: compute Elmore wire delay -----
            // wire_delay_ps = R(ohm/um) * C(fF/um) * L(um)^2 / 2 * 1e-3
            if ( r_per_um > 0.0f && c_per_um > 0.0f )
            {
                wire_delay = r_per_um * c_per_um * wirelength_um * wirelength_um * 0.5f * 1e-3f;
            }
            else
            {
                // Fallback: use fixed WireDelay coefficient (unit delay * length)
                wire_delay = wirelength_um * p->pPars->WireDelay;
            }
        }
    }

    // ============================================================
    // Cell delay computation (accumulate into cell_delay)
    // ============================================================
    if ( pCut->fAndCut )
    {
        If_CutForEachLeaf( p, pCut, pLeaf, i )
        {
            DelayCur = If_ObjCutBest(pLeaf)->Delay + p->pPars->nAndDelay;
            cell_delay = IF_MAX( cell_delay, DelayCur );
        }
    }
    else if ( p->pPars->pLutLib )
    {
        assert( !p->pPars->fLiftLeaves );
        pLutDelays = p->pPars->pLutLib->pLutDelays[pCut->nLeaves];
        if ( p->pPars->pLutLib->fVarPinDelays )
        {
            If_CutSortInputPins( p, pCut, pPinPerm, pPinDelays );
            for ( i = 0; i < (int)pCut->nLeaves; i++ )
            {
                DelayCur = pPinDelays[pPinPerm[i]] + pLutDelays[i];
                cell_delay = IF_MAX( cell_delay, DelayCur );
            }
        }
        else
        {
            If_CutForEachLeaf( p, pCut, pLeaf, i )
            {
                DelayCur = If_ObjCutBest(pLeaf)->Delay + pLutDelays[0];
                cell_delay = IF_MAX( cell_delay, DelayCur );
            }
        }
    }
    else if ( p->pPars->pSclLib && p->pPars->WireDelay > 0.0 )
    {
        // Wire-aware mapping with SC_Lib (tdelay-based cell delay)
        SC_Lib * pSclLib = (SC_Lib *)p->pPars->pSclLib;
        SC_Cell * pCellBest = NULL;
        SC_Cell * pCell = NULL;
        int i_temp;
        SC_LibForEachCell( pSclLib, pCell, i_temp )
        {
            if ( pCell->n_inputs == (int)pCut->nLeaves )
            {
                pCellBest = pCell->pRepr ? pCell->pRepr : pCell;
                break;
            }
        }
        if ( pCellBest == NULL )
        {
            If_CutForEachLeaf( p, pCut, pLeaf, i )
            {
                DelayCur = If_ObjCutBest(pLeaf)->Delay + 1.0;
                cell_delay = IF_MAX( cell_delay, DelayCur );
            }
        }
        else
        {
            float Slew = Abc_SclComputeAverageSlew( pSclLib );
            float Load = SC_CellPinCapAve( pCellBest );
            If_CutForEachLeaf( p, pCut, pLeaf, i )
            {
                float tdelay_pin = Scl_LibPinArrivalEstimate( pCellBest, i, Slew, Load );
                DelayCur = If_ObjCutBest(pLeaf)->Delay + tdelay_pin;
                cell_delay = IF_MAX( cell_delay, DelayCur );
            }
        }
    }
    else if ( pCut->fUser )
    {
        assert( !p->pPars->fLiftLeaves );
        If_CutForEachLeaf( p, pCut, pLeaf, i )
        {
            Pin2PinDelay = pPerm ? (pPerm[i] == IF_BIG_CHAR ? -IF_BIG_CHAR : pPerm[i]) : 1;
            DelayCur = If_ObjCutBest(pLeaf)->Delay + (float)Pin2PinDelay;
            cell_delay = IF_MAX( cell_delay, DelayCur );
        }
    }
    else
    {
        if ( p->pPars->fLiftLeaves )
        {
            If_CutForEachLeafSeq( p, pCut, pLeaf, Shift, i )
            {
                DelayCur = If_ObjCutBest(pLeaf)->Delay - Shift * p->Period;
                cell_delay = IF_MAX( cell_delay, DelayCur + 1.0 );
            }
        }
        else
        {
            If_CutForEachLeaf( p, pCut, pLeaf, i )
            {
                DelayCur = If_ObjCutBest(pLeaf)->Delay + 1.0;
                cell_delay = IF_MAX( cell_delay, DelayCur );
            }
        }
    }

    // Total delay = cell_delay + wire_delay
    Delay = cell_delay + wire_delay;

    // ============================================================
    // Debug output (first N cuts)
    // ============================================================
    if ( print_count < print_limit )
    {
        print_count++;
        printf( "[DEBUG If_CutDelay #%d] obj=%d nLeaves=%d  "
                "wirelength=%.2fum  R=%.4f ohm/um  C=%.4f fF/um  "
                "cell_delay=%.4f  wire_delay=%.4fps  total=%.4f\n",
                print_count, pObj->Id, pCut->nLeaves,
                wirelength_um, r_per_um, c_per_um,
                cell_delay, wire_delay, Delay );
    }

    return Delay;
}

/**Function*************************************************************

  Synopsis    []

  Description []
               
  SideEffects []

  SeeAlso     []

***********************************************************************/
void If_CutPropagateRequired( If_Man_t * p, If_Obj_t * pObj, If_Cut_t * pCut, float ObjRequired )
{
    static int pPinPerm[IF_MAX_LUTSIZE];
    static float pPinDelays[IF_MAX_LUTSIZE];
    If_Obj_t * pLeaf;
    float * pLutDelays;
    float Required;
    int i, Pin2PinDelay;//, iLeaf;
    assert( !p->pPars->fLiftLeaves );
    // compute the pins
    if ( pCut->fAndCut )
    {
        If_CutForEachLeaf( p, pCut, pLeaf, i )
            pLeaf->Required = IF_MIN( pLeaf->Required, ObjRequired - p->pPars->nAndDelay );
    }
    else if ( p->pPars->pLutLib )
    {
        pLutDelays = p->pPars->pLutLib->pLutDelays[pCut->nLeaves];
        if ( p->pPars->pLutLib->fVarPinDelays )
        {
            // compute the delay using sorted pins
            If_CutSortInputPins( p, pCut, pPinPerm, pPinDelays );
            for ( i = 0; i < (int)pCut->nLeaves; i++ )
            {
                Required = ObjRequired - pLutDelays[i];
                pLeaf = If_ManObj( p, pCut->pLeaves[pPinPerm[i]] );
                pLeaf->Required = IF_MIN( pLeaf->Required, Required );
            }
        }
        else
        {
            Required = ObjRequired;
            If_CutForEachLeaf( p, pCut, pLeaf, i )
                pLeaf->Required = IF_MIN( pLeaf->Required, Required - pLutDelays[0] );
        }
    }
    // Wire-aware mapping with SC_Lib (tdelay-based backward required time)
    // Mirror the forward cell delay computation: subtract tdelay[pin] from ObjRequired
    else if ( p->pPars->pSclLib && p->pPars->WireDelay > 0.0 )
    {
        SC_Lib * pSclLib = (SC_Lib *)p->pPars->pSclLib;
        SC_Cell * pCellBest = NULL;
        SC_Cell * pCell = NULL;
        int i_temp;
        SC_LibForEachCellClass( pSclLib, pCell, i_temp )
        {
            if ( pCell->n_inputs == (int)pCut->nLeaves )
            {
                pCellBest = pCell->pRepr ? pCell->pRepr : pCell;
                break;
            }
        }
        if ( pCellBest == NULL )
        {
            // Fallback: use unit delay
            Required = ObjRequired;
            If_CutForEachLeaf( p, pCut, pLeaf, i )
                pLeaf->Required = IF_MIN( pLeaf->Required, Required - (float)1.0 );
        }
        else
        {
            float Slew = Abc_SclComputeAverageSlew( pSclLib );
            float Load = SC_CellPinCapAve( pCellBest );
            If_CutForEachLeaf( p, pCut, pLeaf, i )
            {
                float tdelay_pin = Scl_LibPinArrivalEstimate( pCellBest, i, Slew, Load );
                pLeaf->Required = IF_MIN( pLeaf->Required, ObjRequired - tdelay_pin );
            }
        }
    }
    else if ( p->pPars->fUserLutDec || p->pPars->fUserLut2D )
    {
        Required = ObjRequired;
        If_CutForEachLeaf( p, pCut, pLeaf, i )
            pLeaf->Required = IF_MIN( pLeaf->Required, Required - If_LutDecPinRequired( p, pCut, i, ObjRequired ) );
    }
    else
    {
        if ( pCut->fUser )
        {
            char Perm[IF_MAX_FUNC_LUTSIZE], * pPerm = Perm;
            if ( p->pPars->fDelayOpt )
            {
                int Delay = If_CutSopBalancePinDelays( p, pCut, pPerm );
                assert( -Delay > IF_INFINITY/2 || Delay > IF_INFINITY/2 || Delay == (int)pCut->Delay );
            }
            else if ( p->pPars->fDelayOptLut )
            {
                int Delay = If_CutLutBalancePinDelays( p, pCut, pPerm );
                assert( -Delay > IF_INFINITY/2 || Delay > IF_INFINITY/2 || Delay == (int)pCut->Delay );
            }
            else if ( p->pPars->fDsdBalance )
            {
                int Delay = If_CutDsdBalancePinDelays( p, pCut, pPerm );
                assert( -Delay > IF_INFINITY/2 || Delay > IF_INFINITY/2 || Delay == (int)pCut->Delay );
            }
            else
                pPerm = If_CutPerm(pCut);
            If_CutForEachLeaf( p, pCut, pLeaf, i )
            {
                Pin2PinDelay = pPerm ? (pPerm[i] == IF_BIG_CHAR ? -IF_BIG_CHAR : pPerm[i]) : 1;
                Required = ObjRequired - (float)Pin2PinDelay;
                pLeaf->Required = IF_MIN( pLeaf->Required, Required );
            }
        }
        else
        {
            Required = ObjRequired;
            If_CutForEachLeaf( p, pCut, pLeaf, i )
                pLeaf->Required = IF_MIN( pLeaf->Required, Required - (float)1.0 );
        }
    }
}

/**Function*************************************************************

  Synopsis    [Returns the max delay of the POs.]

  Description []
               
  SideEffects []

  SeeAlso     []

***********************************************************************/
float If_ManDelayMax( If_Man_t * p, int fSeq )
{
    If_Obj_t * pObj;
    float DelayBest;
    int i;
    if ( p->pPars->fLatchPaths && (p->pPars->nLatchesCi == 0 || p->pPars->nLatchesCo == 0) )
    {
        Abc_Print( 0, "Delay optimization of latch path is not performed because there is no latches.\n" );
        p->pPars->fLatchPaths = 0;
    }
    DelayBest = -IF_FLOAT_LARGE;
    if ( fSeq )
    {
        assert( p->pPars->nLatchesCi > 0 );
        If_ManForEachPo( p, pObj, i )
            if ( DelayBest < If_ObjArrTime(If_ObjFanin0(pObj)) )
                 DelayBest = If_ObjArrTime(If_ObjFanin0(pObj));
    }
    else if ( p->pPars->fLatchPaths )
    {
        If_ManForEachLatchInput( p, pObj, i )
            if ( DelayBest < If_ObjArrTime(If_ObjFanin0(pObj)) )
                 DelayBest = If_ObjArrTime(If_ObjFanin0(pObj));
    }
    else 
    {
        If_ManForEachCo( p, pObj, i )
            if ( DelayBest < If_ObjArrTime(If_ObjFanin0(pObj)) )
                 DelayBest = If_ObjArrTime(If_ObjFanin0(pObj));
    }
    return DelayBest;
}

/**Function*************************************************************

  Synopsis    [Computes the required times of all nodes.]

  Description []
               
  SideEffects []

  SeeAlso     []

***********************************************************************/
void If_ManComputeRequired( If_Man_t * p )
{
    If_Obj_t * pObj;
    int i, Counter;
    float reqTime;

    // compute area, clean required times, collect nodes used in the mapping
//    p->AreaGlo = If_ManScanMapping( p );
    If_ManMarkMapping( p );
    if ( p->pManTim == NULL )
    {
        // get the global required times
        p->RequiredGlo = If_ManDelayMax( p, 0 );

        // consider the case when the required times are given
        if ( p->pPars->pTimesReq && !p->pPars->fAreaOnly )
        {
            // make sure that the required time hold
            Counter = 0;
            If_ManForEachCo( p, pObj, i )
            {
                if ( If_ObjArrTime(If_ObjFanin0(pObj)) > p->pPars->pTimesReq[i] + p->fEpsilon )
                {
                    If_ObjFanin0(pObj)->Required = If_ObjArrTime(If_ObjFanin0(pObj));
                    Counter++;
    //                Abc_Print( 0, "Required times are violated for output %d (arr = %d; req = %d).\n", 
    //                    i, (int)If_ObjArrTime(If_ObjFanin0(pObj)), (int)p->pPars->pTimesReq[i] );
                }
                else
                    If_ObjFanin0(pObj)->Required = p->pPars->pTimesReq[i];
            }
            if ( Counter && !p->fReqTimeWarn )
            {
                Abc_Print( 0, "Required times are exceeded at %d output%s. The earliest arrival times are used.\n", Counter, Counter > 1 ? "s":"" );
                p->fReqTimeWarn = 1;
            }
        }
        else
        {
            // find new delay target
            if ( p->pPars->nRelaxRatio && p->pPars->DelayTargetNew == 0 )
                p->pPars->DelayTargetNew = p->RequiredGlo * (100.0 + p->pPars->nRelaxRatio) / 100.0; 

            // update the required times according to the target
            if ( p->pPars->DelayTarget != -1 )
            {
                if ( p->RequiredGlo > p->pPars->DelayTarget + p->fEpsilon )
                {
                    if ( p->fNextRound == 0 )
                    {
                        p->fNextRound = 1;
                        Abc_Print( 0, "Cannot meet the target required times (%4.2f). Mapping continues anyway.\n", p->pPars->DelayTarget );
                    }
                }
                else if ( p->RequiredGlo < p->pPars->DelayTarget - p->fEpsilon )
                {
                    if ( p->fNextRound == 0 )
                    {
                        p->fNextRound = 1;
//                        Abc_Print( 0, "Relaxing the required times from (%4.2f) to the target (%4.2f).\n", p->RequiredGlo, p->pPars->DelayTarget );
                    }
                    p->RequiredGlo = p->pPars->DelayTarget;
                }
            }
            else if ( p->pPars->DelayTargetNew > 0 ) // relax the required times 
                p->RequiredGlo = p->pPars->DelayTargetNew;
            // do not propagate required times if area minimization is requested
            if ( p->pPars->fAreaOnly ) 
                return;
            // set the required times for the POs
            if ( p->pPars->fDoAverage )
            {
                if ( p->pPars->nRelaxRatio )
                {
                    If_ManForEachCo( p, pObj, i )
                        If_ObjFanin0(pObj)->Required = If_ObjArrTime(If_ObjFanin0(pObj)) * (100.0 + p->pPars->nRelaxRatio) / 100.0;
                }
                else
                {
                    If_ManForEachCo( p, pObj, i )
                        If_ObjFanin0(pObj)->Required = If_ObjArrTime(If_ObjFanin0(pObj));
                }
            }
            else if ( p->pPars->fLatchPaths )
            {
                If_ManForEachLatchInput( p, pObj, i )
                    If_ObjFanin0(pObj)->Required = p->RequiredGlo;
            }
            else 
            {
                If_ManForEachCo( p, pObj, i )
                    If_ObjFanin0(pObj)->Required = p->RequiredGlo;
            }
        }
        // go through the nodes in the reverse topological order
    //    Vec_PtrForEachEntry( If_Obj_t *, p->vMapped, pObj, i )
    //        If_CutPropagateRequired( p, pObj, If_ObjCutBest(pObj), pObj->Required );
        If_ManForEachObjReverse( p, pObj, i )
        {
            if ( pObj->nRefs == 0 )
                continue;
            If_CutPropagateRequired( p, pObj, If_ObjCutBest(pObj), pObj->Required );
        }
    }
    else
    {
        // get the global required times
        p->RequiredGlo = If_ManDelayMax( p, 0 );

        // find new delay target
        if ( p->pPars->nRelaxRatio && p->pPars->DelayTargetNew == 0 )
            p->pPars->DelayTargetNew = p->RequiredGlo * (100.0 + p->pPars->nRelaxRatio) / 100.0; 

        // update the required times according to the target
        if ( p->pPars->DelayTarget != -1 )
        {
            if ( p->RequiredGlo > p->pPars->DelayTarget + p->fEpsilon )
            {
                if ( p->fNextRound == 0 )
                {
                    p->fNextRound = 1;
                    Abc_Print( 0, "Cannot meet the target required times (%4.2f). Mapping continues anyway.\n", p->pPars->DelayTarget );
                }
            }
            else if ( p->RequiredGlo < p->pPars->DelayTarget - p->fEpsilon )
            {
                if ( p->fNextRound == 0 )
                {
                    p->fNextRound = 1;
//                    Abc_Print( 0, "Relaxing the required times from (%4.2f) to the target (%4.2f).\n", p->RequiredGlo, p->pPars->DelayTarget );
                }
                p->RequiredGlo = p->pPars->DelayTarget;
            }
        }
        else if ( p->pPars->DelayTargetNew > 0 ) // relax the required times 
            p->RequiredGlo = p->pPars->DelayTargetNew;

        // do not propagate required times if area minimization is requested
        if ( p->pPars->fAreaOnly ) 
            return;
        // set the required times for the POs
        Tim_ManIncrementTravId( p->pManTim );
        if ( p->vCoAttrs )
        {
            assert( If_ManCoNum(p) == Vec_IntSize(p->vCoAttrs) );
            If_ManForEachCo( p, pObj, i )
            { 
                if ( Vec_IntEntry(p->vCoAttrs, i) == -1 )       // -1=internal
                    continue;
                if ( Vec_IntEntry(p->vCoAttrs, i) == 0 )        //  0=optimize
                    Tim_ManSetCoRequired( p->pManTim, i, p->RequiredGlo );
                else if ( Vec_IntEntry(p->vCoAttrs, i) == 1 )   //  1=keep
                    Tim_ManSetCoRequired( p->pManTim, i, If_ObjArrTime(If_ObjFanin0(pObj)) );
                else if ( Vec_IntEntry(p->vCoAttrs, i) == 2 )   //  2=relax
                    Tim_ManSetCoRequired( p->pManTim, i, IF_FLOAT_LARGE );
                else assert( 0 );
            }
        }
        else if ( p->pPars->fDoAverage )
        {
            if ( p->pPars->nRelaxRatio )
            {
                If_ManForEachCo( p, pObj, i )
                    Tim_ManSetCoRequired( p->pManTim, i, If_ObjArrTime(If_ObjFanin0(pObj)) * (100.0 + p->pPars->nRelaxRatio) / 100.0 );
            }
            else
            {
                If_ManForEachCo( p, pObj, i )
                    Tim_ManSetCoRequired( p->pManTim, i, If_ObjArrTime(If_ObjFanin0(pObj)) );
            }
        }
        else if ( p->pPars->fLatchPaths )
        {
            If_ManForEachPo( p, pObj, i )
                Tim_ManSetCoRequired( p->pManTim, i, IF_FLOAT_LARGE );
            If_ManForEachLatchInput( p, pObj, i )
                Tim_ManSetCoRequired( p->pManTim, i, p->RequiredGlo );
        }
        else  
        {
            Tim_ManInitPoRequiredAll( p->pManTim, p->RequiredGlo );
//            If_ManForEachCo( p, pObj, i )
//                Tim_ManSetCoRequired( p->pManTim, pObj->IdPio, p->RequiredGlo );
        }
        // go through the nodes in the reverse topological order
        If_ManForEachObjReverse( p, pObj, i )
        {
            if ( If_ObjIsAnd(pObj) )
            {
                if ( pObj->nRefs == 0 )
                    continue;
                If_CutPropagateRequired( p, pObj, If_ObjCutBest(pObj), pObj->Required );
            }
            else if ( If_ObjIsCi(pObj) )
            {
                reqTime = pObj->Required;
                Tim_ManSetCiRequired( p->pManTim, pObj->IdPio, reqTime );
            }
            else if ( If_ObjIsCo(pObj) )
            {
                reqTime = Tim_ManGetCoRequired( p->pManTim, pObj->IdPio );
                If_ObjFanin0(pObj)->Required = IF_MIN( reqTime, If_ObjFanin0(pObj)->Required );
            }
            else if ( If_ObjIsConst1(pObj) )
            {
            }
            else // add the node to the mapper
                assert( 0 );
        }
    }
}

////////////////////////////////////////////////////////////////////////
///                       END OF FILE                                ///
////////////////////////////////////////////////////////////////////////


ABC_NAMESPACE_IMPL_END
