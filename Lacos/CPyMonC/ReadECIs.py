"""This module is all about to read the ECIs and the cluster file obtained from PyCLEX module."""
import copy
import numpy as np
from Lacos.CPyMonC.Utility import Cartesian2Fractions, Fractions2Cartesian
from Lacos.CPyMonC.miscfunctions import *

NORM_TOLERANCE = 1.e-3
COORD_TOLERANCE = 1.e-5
VERY_SMALL_NUMBER = 1.0e-14


class Cluster:
    """
    The cluster class consists of Parent cluster and all symmetry equivalent Child clusters. The symmetry
    equivalent clusters are stored rather than generating them through operation on-demand. The ECIs and the
    correlation index is also included in the class.
    """
    Ksites = None
    Parent = None
    Child = None
    ClusterID = None
    ECIs = []
    CIndexList = []

    def __init__(self, inCluster, id_num, eci, correl_index):
        self.Ksites = len( inCluster )
        self.ClusterID = str( self.Ksites ) + '-' + str( id_num ).zfill( 5 )
        if self.Ksites > 0:
            self.Parent = SMALL_NUMBER_2_ZERO( inCluster )

        self.ECIs = copy.copy( eci )
        self.CIndexList = copy.copy( correl_index )

    def Equivalent_Clusters(self, base_lattice, rot_symmetry):
        self.Child = {j: [] for j in range( len( base_lattice[ 'Coordinates'] ) )}
        self.Parent = [Fractions2Cartesian( base_lattice['Lattice_Vectors'], _ ) for _ in self.Parent]
        if self.Ksites == 1:
            for (j, _) in enumerate( base_lattice['Coordinates'] ): self.Child[j] = [self.Parent]
        else:
            if rot_symmetry is None: exit_script(
                "Error :: No symmetry operators found for cluster generation. Check."
            )
            self.__Clusters_From_Parent( base_lattice, rot_symmetry )

        return self

    def __Clusters_From_Parent(self, base_lattice, rot_symmetry):
        clust = [r for r in self.Parent]
        rotclust_List = [[list( np.array( symm ).dot( np.array( _ ) ) ) for _ in clust] for symm in
                         rot_symmetry]

        for (j, R) in enumerate( base_lattice['Real_Coordinates'] ):
            dispclust = [(np.array( _ ) + np.array( R[:3] )).tolist() for _ in rotclust_List]

            for _, cluster in zip( dispclust, rotclust_List ):
                if Valid_Cluster( _, base_lattice ) and Cluster_Not_In_List( cluster, self.Child[j] ):
                    fcluster = Cartesian2Fractions( base_lattice['Lattice_Vectors'], cluster )
                    self.Child[j].append( SMALL_NUMBER_2_ZERO( fcluster ) )

        if sum( len( self.Child[_] ) for _ in self.Child ) == 0:
            exit_script(
                "Error :: The cluster does not exist in the system. If it is read from a file, re-check cluster "
                "co-ordinates (if not, problem with lattice symmetry)."
            )
        self.Parent = Cartesian2Fractions( base_lattice['Lattice_Vectors'], self.Parent )

        return self


def Load_Clusters_and_ECIs(cluster_filename):
    """
    This method will Read Clusters from file and equivalent clusters are also read. this method is having some
    private methods. :param cluster_filename: :return:
    """
    def Splitted_Line(nline, sr):
        return [r.lstrip().rstrip() for r in nline.split( sr ) if r != '']

    def Find_Cluster_Header(Lines, sline):
        if sline == len( Lines ): return -1, -1
        nstart = sline
        while True:
            if all( keywords in Lines[nstart] for keywords in ["Cluster Type", "k-site"] ): break
            nstart += 1
            if nstart == len( Lines ): return -1, -1
        nend = nstart + 1
        while True:
            if all( keywords in Lines[nend] for keywords in ["Cluster Type", "k-site"] ): break
            nend += 1
            if nend == len( Lines ): return nstart, nend
        return nstart, nend

    def Find_ClusterNumber(Lines, Start, End, n):
        pstart = -1
        for j, nl in enumerate( Lines[Start:End] ):
            if all( _ in nl for _ in ["Cluster", "Number", str( n )] ): pstart = Start + j
            if all( _ in nl for _ in ["Cluster", "Number", str( n + 1 )] ): return pstart, Start + j

        if pstart < 0:
            exit_script( "Error in reading ECI file :: Unable to find %d cluster" % n )
        else:
            return pstart, End

    def Find_Number_of_kClusters(Lines, Start, End, k):
        for line in Lines[Start:End]:
            if all( _ in line.lower() for _ in ['clusters', 'number', str( k ) + '-site'] ):
                return int( re.split( r'[;,\s]\s*', line.lstrip() )[-1] )
        exit_script(
          "Error in reading ECI file :: Unable to find number of clusters for %d-site" % k
        )

    def Find_ECIs_and_Correlation_Index(subLines, ksite):
        for (j, nl) in enumerate( subLines ):
            if all( r in nl for r in ["ECIs", "Correlation", "Index"] ): break
        if j == len( subLines ): exit_script(
            "Error in reading Clusters and ECIs: Unable to locate the start of ECIs. Exit.\n" )
        for (jj, nl) in enumerate( subLines ):
            if all( r in nl for r in ["Parent", "Cluster"] ): break

        if ksite == 0:
            return [float( subLines[j + 1].lstrip().rstrip() )]
        else:
            eci = [float( nl.split( ':' )[0] ) for nl in subLines[j + 1:jj]]
            correl_index = [[int( r ) for r in Splitted_Line( nl.split( ':' )[1], ' ' )] for nl in subLines[j + 1:jj]]
            return eci, correl_index

    def Find_Parent_Cluster(subLines, ksite):
        for (j, nl) in enumerate( subLines ):
            if all( r in nl for r in ["Parent", "Cluster"] ): break
        coord = [[float( r ) for r in Splitted_Line( nl, ' ' )] for nl in subLines[j + 1:j + ksite + 1]]
        return coord

    try:
        fp = open( cluster_filename, "r" )
    except IOError as e:
        exit_script( "Error: Cluster file <%s> does not exist. Exit" % cluster_filename )
    All_Lines = [(r.rstrip( '\n' )).lstrip() for r in fp]
    All_Lines = [re.sub( r'[\x00-\x1f\x7f-\x9f]', '', r ) for r in All_Lines]
    TotalClusters = [int( r ) for r in All_Lines[0].split( ' ' ) if r.isdigit()][0]
    print( "\nFile %s contains %d clusters" % (cluster_filename, TotalClusters) )
    elmlist = [elm.lstrip() for elm in (All_Lines[1].split( ':' )[1]).split( ',' )]
    ElementRep = {_: None for _ in elmlist}
    if all( _ in All_Lines[2] for _ in ['Point', 'function', 'elements'] ):
        for k in range( len( elmlist ) ):
            tmp = [_ for _ in re.split( r'[;,:\s]\s*', All_Lines[3 + k] ) if _ != '']
            ElementRep[tmp[0]] = [float( _ ) for _ in tmp[1:]]
    ClusterList, Cluster_Id = [], '00001X'
    start_line = 0
    while True:
        Start, End = Find_Cluster_Header( All_Lines, start_line )
        if Start == -1: break
        ksite = int( All_Lines[Start].split( ":" )[1] )
        if ksite == 0:
            eci = Find_ECIs_and_Correlation_Index( All_Lines[Start:End], ksite )
            ClusterList.append( Cluster( [], Cluster_Id, eci, [] ) )
            print( "    Collected 1 0-site clusters" )
            Cluster_Id = str( int( Cluster_Id[:-1] ) + 1 ).zfill( 5 ) + 'X'
        else:
            Nkcluster = Find_Number_of_kClusters( All_Lines, Start, End, ksite )
            for j in range( 1, Nkcluster + 1 ):
                pStart, pEnd = Find_ClusterNumber( All_Lines, Start, End, j )
                eci, correl_index = Find_ECIs_and_Correlation_Index( All_Lines[pStart:pEnd], ksite )
                pcluster = Find_Parent_Cluster( All_Lines[pStart:pEnd], ksite )
                ClusterList.append( Cluster( pcluster, Cluster_Id, eci, correl_index ) )

                Cluster_Id = str( int( Cluster_Id[:-1] ) + 1 ).zfill( 5 ) + 'X'
            print( "    Collected %d %d-site clusters" % (Nkcluster, ksite) )

        start_line = End

    print( "\n    Read cluster file %s and collected %d clusters." % (cluster_filename, len( ClusterList )) )
    print(
        "\n The ECI file was prepared for following elements: %s\n If it is not what you are intending, Quit now."
        " " % ",".join( _ for _ in list( ElementRep.keys() ) )
    )

    return ClusterList, ElementRep


def SMALL_NUMBER_2_ZERO(a):
    """
    This method will convert small number to zero.
    :param a:
    :return: A list type.
    """
    _ = np.array( a )
    _[np.abs( _ ) < VERY_SMALL_NUMBER] = 0.0
    return _.tolist()


def SameCoord(Ra, Rb):
    if all( abs( pa - pb ) < NORM_TOLERANCE for (pa, pb) in zip( Ra, Rb ) ): return True
    return False


def SameCluster_Unordered(dRl1, dRl2):
    """
    This method will check similar clusters order.
    :param dRl1:
    :param dRl2:
    :return: Boolean type.
    """
    tmp = copy.copy( dRl2 )
    for R1 in dRl1:
        for (j, R2) in enumerate( tmp ):
            if SameCoord( R1, R2 ):
                tmp.pop( j )
                break
    if tmp == []:
        return True
    else:
        False


def SameCluster_Ordered(dRl1, dRl2):
    for (R1, R2) in zip( dRl1, dRl2 ):
        if not SameCoord( R1, R2 ): return False
    return True


def Cluster_Not_In_List(newclust, inList):
    """
    This method will chek the clusters are in the list or not.
    :param newclust:
    :param inList:
    :return: Boolean type.
    """
    if any( SameCluster_Ordered( newclust, _ ) for _ in inList ): return False
    return True


def Valid_Cluster(cls, base_lattice):
    """
    This method will find the valid clusters.
    :param cls:
    :param base_lattice:
    :return: Boolean value.
    """
    for R in cls:
        r0 = Cartesian2Fractions( base_lattice['Lattice_Vectors'], R )
        _ = np.round( np.array( r0[:3] ) % np.array( [1.0, 1.0, 1.0] ), 6 )
        _[abs( _ - np.array(
            [1.0, 1.0, 1.0] ) ) < COORD_TOLERANCE] = 0.0
        if not Coord_In_List( _, base_lattice['Coordinates'] ): return False
    return True


def Coord_In_List(R, inList):
    if any( all( abs( pa - pb ) < NORM_TOLERANCE for (pa, pb) in zip( R[:3], _[:3] ) ) for _ in inList ): return True
    return False
