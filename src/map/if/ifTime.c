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
    float wirelength = 0.0f;
    float r_per_um = 0.0f, c_per_um = 0.0f;

    // ============================================================
    // Wire-aware mapping: compute wire delay using centroid distance
    // Based on phyLS approach: centroid of cut leaves = gate output position,
    // worst leaf-to-centroid Manhattan distance = wire length.
    // If_Obj_t::pCopy points to the corresponding Abc_Obj_t in pNtkOrig.
    // ============================================================
    if ( p->pPars->WireDelay > 0.0 && p->pPars->pNtkCoords )
    {
        Abc_Ntk_t * pNtkCoords = (Abc_Ntk_t *)p->pPars->pNtkCoords;
        float pin_x_sum = 0.0f, pin_y_sum = 0.0f;
        int nPins = 0;

        // Phase 1: collect unique pin coordinates (deduplicate by instance name)
        // Output pin
        {
            Abc_Obj_t * pAbcObj = (Abc_Obj_t *)pObj->pCopy;
            if ( pAbcObj != NULL )
            {
                char * pName = Abc_ObjName( pAbcObj );
                float x, y;
                if ( Io_ReadCoordsGetCoord( pNtkCoords, pName, &x, &y ) )
                {
                    pin_x_sum += x;
                    pin_y_sum += y;
                    nPins++;
                }
            }
        }

        // Input pins (cut leaves) - collect unique coordinates only
        If_CutForEachLeaf( p, pCut, pLeaf, i )
        {
            Abc_Obj_t * pAbcObj = (Abc_Obj_t *)pLeaf->pCopy;
            if ( pAbcObj == NULL )
                continue;
            char * pName = Abc_ObjName( pAbcObj );

            // Dedup: skip if same instance as output pin
            {
                Abc_Obj_t * pOutObj = (Abc_Obj_t *)pObj->pCopy;
                if ( pOutObj != NULL )
                {
                    char * pOutName = Abc_ObjName( pOutObj );
                    if ( pOutName && pName && strcmp( pOutName, pName ) == 0 )
                        continue;
                }
            }

            // Dedup against previously processed leaves
            {
                int j, is_dup = 0;
                If_Obj_t * pOther;
                If_CutForEachLeaf( p, pCut, pOther, j )
                {
                    if ( j >= i )
                        break;
                    Abc_Obj_t * pPrevObj = (Abc_Obj_t *)pOther->pCopy;
                    if ( pPrevObj == NULL )
                        continue;
                    char * pPrevName = Abc_ObjName( pPrevObj );
                    if ( pPrevName && pName && strcmp( pPrevName, pName ) == 0 )
                    {
                        is_dup = 1;
                        break;
                    }
                }
                if ( is_dup )
                    continue;
            }

            float x, y;
            if ( Io_ReadCoordsGetCoord( pNtkCoords, pName, &x, &y ) )
            {
                pin_x_sum += x;
                pin_y_sum += y;
                nPins++;
            }
        }

        if ( nPins == 0 )
        {
            // No coordinates found; fall back to zero wire delay
            wirelength = 0.0f;
        }
        else
        {
            // Phase 2: worst (max) Manhattan distance from centroid to any unique pin
            float centroid_x = pin_x_sum / nPins;
            float centroid_y = pin_y_sum / nPins;
            float worst_wirelength = 0.0f;

            // Check output pin
            {
                Abc_Obj_t * pAbcObj = (Abc_Obj_t *)pObj->pCopy;
                if ( pAbcObj != NULL )
                {
                    char * pName = Abc_ObjName( pAbcObj );
                    float x, y;
                    if ( Io_ReadCoordsGetCoord( pNtkCoords, pName, &x, &y ) )
                    {
                        float wl = fabsf( centroid_x - x ) + fabsf( centroid_y - y );
                        if ( wl > worst_wirelength )
                            worst_wirelength = wl;
                    }
                }
            }

            // Check each unique leaf pin (already deduplicated in Phase 1)
            If_CutForEachLeaf( p, pCut, pLeaf, i )
            {
                Abc_Obj_t * pAbcObj = (Abc_Obj_t *)pLeaf->pCopy;
                if ( pAbcObj == NULL )
                    continue;
                char * pName = Abc_ObjName( pAbcObj );
                float x, y;
                if ( Io_ReadCoordsGetCoord( pNtkCoords, pName, &x, &y ) )
                {
                    float wl = fabsf( centroid_x - x ) + fabsf( centroid_y - y );
                    if ( wl > worst_wirelength )
                        worst_wirelength = wl;
                }
            }
            wirelength = worst_wirelength;
        }

        // Phase 3: compute wire delay from wirelength and RC
        Io_ReadCoordsGetWireRCFromCoords( pNtkCoords, &r_per_um, &c_per_um );
        if ( r_per_um > 0.0f && c_per_um > 0.0f )
        {
            // wire_delay_ps = R(ohm/um) * C(fF/um) * L(um) * 1e-3
            wire_delay = r_per_um * c_per_um * wirelength * 1e-3f;
        }
        else
        {
            // Fallback: use fixed WireDelay coefficient
            wire_delay = wirelength * p->pPars->WireDelay;
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
                "wirelength=%.2f  R=%.4f ohm/um  C=%.4f fF/um  "
                "cell_delay=%.4f  wire_delay=%.4f  total=%.4f\n",
                print_count, pObj->Id, pCut->nLeaves,
                wirelength, r_per_um, c_per_um,
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
