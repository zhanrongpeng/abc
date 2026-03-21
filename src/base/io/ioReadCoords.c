/**CFile****************************************************************

  FileName    [ioReadCoords.c]

  SystemName  [ABC: Logic synthesis and verification system.]

  PackageName [Command processing package.]

  Synopsis    [Procedures to read coordinate files for ABC-aware remapping.]

  Author      [OpenROAD]
  
  Affiliation [UC Berkeley]

  Date        [Ver. 1.0. Started - March 2026]

  Revision    [$Id: ioReadCoords.c,v 1.00 2026/03/07 openroad Exp $]

***********************************************************************/

#include "ioAbc.h"
#include "base/main/main.h"
#include "base/main/mainInt.h"
#include "aig/aig/aig.h"

ABC_NAMESPACE_IMPL_START

////////////////////////////////////////////////////////////////////////
///                        DECLARATIONS                              ///
////////////////////////////////////////////////////////////////////////

typedef struct Io_ReadCoords_t_        Io_ReadCoords_t;
struct Io_ReadCoords_t_
{
    // file info
    char *               pFileName;
    FILE *               pFile;
    int                  LineCur;
    // coordinate data
    Vec_Ptr_t *          vNames;
    Vec_Flt_t *          vCoordsX;
    Vec_Flt_t *          vCoordsY;
    // wire RC data (set by read_coords command)
    float                wire_r_per_um;
    float                wire_c_per_um;
    // error handling
    char                 sError[1000];
    int                  fError;
};

static Io_ReadCoords_t * Io_ReadCoordsFile( char * pFileName );
static void Io_ReadCoordsFree( Io_ReadCoords_t * p );
static void Io_ReadCoordsPrintErrorMessage( Io_ReadCoords_t * p );
static int Io_ReadCoordsNetwork( Io_ReadCoords_t * p );

////////////////////////////////////////////////////////////////////////
///                     FUNCTION DEFINITIONS                         ///
////////////////////////////////////////////////////////////////////////

/**Function*************************************************************

  Synopsis    [Reads the coordinate file.]

  Description [Format:
    # Instance coordinates for ABC-aware remapping
    # Format: instance_name x y
    N
    inst_name1 x1 y1
    inst_name2 x2 y2
    ...
  ]
               
  SideEffects []

  SeeAlso     []

***********************************************************************/
Abc_Ntk_t * Io_ReadCoords( char * pFileName )
{
    Io_ReadCoords_t * p;
    Abc_Ntk_t * pNtk = NULL;

    // start the file
    p = Io_ReadCoordsFile( pFileName );
    if ( p == NULL )
        return NULL;

    // read coordinates
    if ( !Io_ReadCoordsNetwork( p ) )
    {
        Io_ReadCoordsPrintErrorMessage( p );
        Io_ReadCoordsFree( p );
        return NULL;
    }

    // print statistics
    printf( "Read %d instance coordinates from file \"%s\".\n",
            Vec_PtrSize(p->vNames), pFileName );
    printf( "  [DEBUG] Wire RC from file: R=%.6e ohm/um, C=%.6e fF/um\n",
            p->wire_r_per_um, p->wire_c_per_um );

    // Set wire RC in ABC frame (if specified in file)
    if ( p->wire_r_per_um > 0.0 || p->wire_c_per_um > 0.0 )
    {
        Abc_FrameSetWireRC( p->wire_r_per_um, p->wire_c_per_um );
        printf( "Wire RC set from file: R=%.6f ohm/um, C=%.6f fF/um\n",
                p->wire_r_per_um, p->wire_c_per_um );
    }
    else
    {
        printf( "  [DEBUG] Wire RC is zero, will use fallback WireDelay coefficient.\n" );
    }

    // create a dummy network to hold the coordinates
    pNtk = Abc_NtkAlloc( ABC_NTK_NETLIST, ABC_FUNC_SOP, 1 );
    pNtk->pSpec = Extra_UtilStrsav( pFileName );
    
    // store coordinates in the network's pData field
    // This will be accessed later when converting to AIG
    pNtk->pData = p;

    // don't free p here - it's stored in pNtk->pData
    return pNtk;
}

/**Function*************************************************************

  Synopsis    []

  Description []
               
  SideEffects []

  SeeAlso     []

***********************************************************************/
static Io_ReadCoords_t * Io_ReadCoordsFile( char * pFileName )
{
    Io_ReadCoords_t * p;

    // allocate the reader
    p = ABC_CALLOC( Io_ReadCoords_t, 1 );
    p->pFileName = pFileName;
    p->LineCur = 1;
    p->vNames = Vec_PtrAlloc( 100 );
    p->vCoordsX = Vec_FltAlloc( 100 );
    p->vCoordsY = Vec_FltAlloc( 100 );

    // open the file
    p->pFile = fopen( pFileName, "r" );
    if ( p->pFile == NULL )
    {
        sprintf( p->sError, "Cannot open file \"%s\".", pFileName );
        p->fError = 1;
        return p;
    }

    return p;
}

/**Function*************************************************************

  Synopsis    []

  Description []
               
  SideEffects []

  SeeAlso     []

***********************************************************************/
static void Io_ReadCoordsFree( Io_ReadCoords_t * p )
{
    if ( p->vNames ) Vec_PtrFree( p->vNames );
    if ( p->vCoordsX ) Vec_FltFree( p->vCoordsX );
    if ( p->vCoordsY ) Vec_FltFree( p->vCoordsY );
    if ( p->pFile ) fclose( p->pFile );
    ABC_FREE( p );
}

/**Function*************************************************************

  Synopsis    []

  Description []
               
  SideEffects []

  SeeAlso     []

***********************************************************************/
static void Io_ReadCoordsPrintErrorMessage( Io_ReadCoords_t * p )
{
    if ( p->sError[0] == '\0' )
        return;
    fprintf( stdout, "Error in line %d of file \"%s\": %s\n", 
             p->LineCur, p->pFileName, p->sError );
}

/**Function*************************************************************

  Synopsis    []

  Description []
               
  SideEffects []

  SeeAlso     []

***********************************************************************/
static int Io_ReadCoordsNetwork( Io_ReadCoords_t * p )
{
    char sBuffer[1000];
    int nInstances = 0;
    int i;

    // Initialize wire RC defaults
    p->wire_r_per_um = 0.0f;
    p->wire_c_per_um = 0.0f;

    // skip comments
    while ( fgets(sBuffer, 1000, p->pFile) != NULL )
    {
        p->LineCur++;
        // skip empty lines and comments (except wire_rc line)
        if ( sBuffer[0] == '\n' )
            continue;
        if ( sBuffer[0] == '#' )
        {
            // Check for wire_rc directive: "wire_rc R C"
            float r_val = 0.0f, c_val = 0.0f;
            if ( sscanf( sBuffer + 1, "wire_rc %f %f", &r_val, &c_val ) == 2 )
            {
                p->wire_r_per_um = r_val;
                p->wire_c_per_um = c_val;
                printf( "Wire RC: R=%.6f ohm/um, C=%.6f fF/um\n", r_val, c_val );
            }
            continue;
        }
        break;
    }

    // read the number of instances
    if ( sscanf( sBuffer, "%d", &nInstances ) != 1 )
    {
        sprintf( p->sError, "Cannot read number of instances." );
        return 0;
    }

    // read each instance coordinate
    for ( i = 0; i < nInstances; i++ )
    {
        char sName[500];
        float x, y;

        if ( fgets(sBuffer, 1000, p->pFile) == NULL )
        {
            sprintf( p->sError, "Unexpected end of file." );
            return 0;
        }
        p->LineCur++;

        if ( sscanf( sBuffer, "%s %f %f", sName, &x, &y ) != 3 )
        {
            sprintf( p->sError, "Cannot parse coordinate line." );
            return 0;
        }

        Vec_PtrPush( p->vNames, Extra_UtilStrsav(sName) );
        Vec_FltPush( p->vCoordsX, x );
        Vec_FltPush( p->vCoordsY, y );
    }

    return 1;
}

/**Function*************************************************************

  Synopsis    [Get coordinate by instance name.]

  Description []
               
  SideEffects []

  SeeAlso     []

***********************************************************************/
int Io_ReadCoordsGetCoord( Abc_Ntk_t * pNtk, char * pName, float * px, float * py )
{
    Io_ReadCoords_t * p;
    int i;

    if ( pNtk == NULL || pNtk->pData == NULL )
        return 0;

    p = (Io_ReadCoords_t *)pNtk->pData;
    if ( p->vNames == NULL )
        return 0;

    // find the instance
    Vec_PtrForEachEntry( char *, p->vNames, pName, i )
    {
        if ( !strcmp( (char *)Vec_PtrEntry(p->vNames, i), pName ) )
        {
            *px = Vec_FltEntry( p->vCoordsX, i );
            *py = Vec_FltEntry( p->vCoordsY, i );
            return 1;
        }
    }

    return 0;
}

/**Function*************************************************************

  Synopsis    [Get coordinate by index.]

  Description []
               
  SideEffects []

  SeeAlso     []

***********************************************************************/
void Io_ReadCoordsGetCoordByIndex( Abc_Ntk_t * pNtk, int i, float * px, float * py )
{
    Io_ReadCoords_t * p;

    if ( pNtk == NULL || pNtk->pData == NULL )
    {
        *px = *py = 0.0f;
        return;
    }

    p = (Io_ReadCoords_t *)pNtk->pData;
    if ( i < 0 || i >= Vec_FltSize(p->vCoordsX) )
    {
        *px = *py = 0.0f;
        return;
    }

    *px = Vec_FltEntry( p->vCoordsX, i );
    *py = Vec_FltEntry( p->vCoordsY, i );
}

/**Function*************************************************************

  Synopsis    [Get number of instances.]

  Description []
               
  SideEffects []

  SeeAlso     []

***********************************************************************/
int Io_ReadCoordsGetCount( Abc_Ntk_t * pNtk )
{
    Io_ReadCoords_t * p;

    if ( pNtk == NULL || pNtk->pData == NULL )
        return 0;

    p = (Io_ReadCoords_t *)pNtk->pData;
    return Vec_PtrSize( p->vNames );
}

/**Function*************************************************************

  Synopsis    [Free coordinate data.]

  Description []
               
  SideEffects []

  SeeAlso     []

***********************************************************************/
void Io_ReadCoordsFreeData( Abc_Ntk_t * pNtk )
{
    if ( pNtk && pNtk->pData )
    {
        Io_ReadCoordsFree( (Io_ReadCoords_t *)pNtk->pData );
        pNtk->pData = NULL;
    }
}

/**Function*************************************************************

  Synopsis    [Get wire RC from Io_ReadCoords_t structure.]

  Description [Extracts R and C per micron from the Io_ReadCoords_t stored in pData.]

  SideEffects []

  SeeAlso     []

***********************************************************************/
void Io_ReadCoordsGetWireRCFromCoords( Abc_Ntk_t * pNtkCoords, float * pR, float * pC )
{
    Io_ReadCoords_t * p;

    if ( pNtkCoords == NULL || pNtkCoords->pData == NULL )
    {
        *pR = 0.0f;
        *pC = 0.0f;
        return;
    }

    p = (Io_ReadCoords_t *)pNtkCoords->pData;
    *pR = p->wire_r_per_um;
    *pC = p->wire_c_per_um;
}

/**Function*************************************************************

  Synopsis    [Get coordinate data from ABC frame.]

  Description []

  SideEffects []

  SeeAlso     []

***********************************************************************/
Abc_Ntk_t * Io_ReadCoordsGetFromFrame( Abc_Frame_t * pAbc )
{
    return pAbc ? pAbc->pNtkCoords : NULL;
}

/**Function*************************************************************

  Synopsis    [Get wire RC from ABC frame.]

  Description [Returns R and C per micron from the frame's wire_rc settings.]

  SideEffects []

  SeeAlso     []

***********************************************************************/
void Io_ReadCoordsGetWireRC( Abc_Frame_t * pAbc, float * pR, float * pC )
{
    if ( pAbc )
    {
        *pR = pAbc->wire_r_per_um;
        *pC = pAbc->wire_c_per_um;
    }
    else
    {
        *pR = 0.0f;
        *pC = 0.0f;
    }
}

ABC_NAMESPACE_IMPL_END
