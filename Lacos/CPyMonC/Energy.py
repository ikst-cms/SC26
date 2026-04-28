import numpy as np
from math import *
from operator import mul
from Lacos.CPyMonC.miscfunctions import *
from Lacos.CPyMonC.Symmetry import Load_Symmetry_Operators
from Lacos.CPyMonC.ReadECIs import Load_Clusters_and_ECIs
from functools import reduce

COORD_TOLERANCE = 1.e-3


def ConfigurationEnergy(config, ClusterList, ECIList, ThetaRep):
    """
    This method used to calculate the energy of a configurations.
    :param config:
    :param ClusterList:
    :param ECIList:
    :param ThetaRep:
    :return: Multi Dimension array type.
    """
    def __Element_at_Site(rin, pList, ncell, COORD_KEYS):
        R = np.round( rin % ncell, 3 )
        R[abs( R - ncell ) < COORD_TOLERANCE] = 0.0
        try:
            return pList[int( np.round( np.sum( COORD_KEYS * R ) ) )][3]
        except:
            for p in pList:
                if all( abs( pa - pb ) < COORD_TOLERANCE for (pa, pb) in zip( R, pList[p][:3] ) ): return pList[p][
                    3]
            exit_script( "Error :: Unable to locate coordinate %s in list. Exit." % ' '.join(
                str( r ) for r in rin ) )

    def __Get_Cluster_Occupations(cls, R, Rlist, ncell, COORD_KEYS):
        return [R[3]] + [__Element_at_Site( np.array( r[:3] ) + np.array( R[:3] ), Rlist, ncell, COORD_KEYS ) for r in
                         cls[1:]]
    correl = np.array( [1.0] )
    for cls in ClusterList:
        Nindex = len( cls.CIndexList )
        TupleIndex = [tuple( np.array( _ ) ) for _ in cls.CIndexList]
        AvgCorrel = [0.0] * Nindex
        for key, R in config['Coordinates'].items():
            for cluster in cls.Child[R[4]]:
                cls_occupation = __Get_Cluster_Occupations( cluster, R, config['Coordinates'],
                                                            np.array( config['nCell'] ), np.array( config['Keys'] ) )
                pp = np.array(
                  reduce( np.multiply.outer, [ThetaRep[cls_occupation[_]] for _ in range( cls.Ksites )] )
                )
                for j in range( Nindex ): AvgCorrel[j] += pp[TupleIndex[j]]
        ncls = sum( len( cls.Child[_] ) for _ in cls.Child )
        ncls *= reduce( mul, config['nCell'] ) * cls.Ksites
        for j in range( Nindex ): AvgCorrel[
            j] /= ncls
        correl = np.append( correl, AvgCorrel )
    return np.dot( correl, ECIList )


def ConfigurationEnergy_Fast(config, ClusterList, ECIList, ThetaRep, SiteBasedClusterList):
    """
    This method will calculate energy of configurations
    :param config:
    :param ClusterList:
    :param ECIList:
    :param ThetaRep: 
    :param SiteBasedClusterList:
    :return: ndarray type.
    """
    correl = np.array( [1.0] )
    for cls in ClusterList:
        Nindex = len( cls.CIndexList )
        TupleIndex = [tuple( np.array( _ ) ) for _ in cls.CIndexList]
        AvgCorrel = [0.0] * Nindex
        for key, R in config['Coordinates'].items():
            for j in range( len( cls.Child[R[4]] ) ):
                cls_occupation = [R[3]] + [config['Coordinates'][_][3] for _ in
                                           SiteBasedClusterList[key][cls.ClusterID][j]]
                pp = np.array( reduce( np.multiply.outer, [ThetaRep[cls_occupation[_]] for _ in range( cls.Ksites )] ) )
                for jj in range( Nindex ): AvgCorrel[jj] += pp[TupleIndex[jj]]
        ncls = sum( len( cls.Child[_] ) for _ in cls.Child )
        ncls *= reduce( mul, config['nCell'] ) * cls.Ksites
        for j in range( Nindex ): AvgCorrel[
            j] /= ncls
        correl = np.append( correl, AvgCorrel )
    return np.dot( correl, ECIList )


def Convert_Clusters_to_IndexList(config, ClusterList):

    def __SiteKey(rin, pList, ncell, COORD_KEYS):
        R = np.round( rin % ncell, 3 )
        R[abs( R - ncell ) < COORD_TOLERANCE] = 0.0
        if int( np.sum( COORD_KEYS * R ) ) in list( pList.keys() ):
            return int( np.round( np.sum( COORD_KEYS * R ) ) )
        else:
            for (
            key, p) in pList.items():
                if all( abs( pa - pb ) < COORD_TOLERANCE for (pa, pb) in zip( R, p[:3] ) ): return key
            exit_script(
                " Error :: Unable to locate coordinate %s in list. Exit." % ' '.join( str( r ) for r in rin )
            )

    def __Get_Cluster_Sites(cls, R, Rlist, ncell, COORD_KEYS):
        return [__SiteKey( np.array( r[:3] ) + np.array( R[:3] ), Rlist, ncell, COORD_KEYS ) for r in cls[1:]]

    def __Clusters_at_R(cls, R, Rlist, ncell, COORD_KEYS):
        tmp = []
        for cluster in cls.Child[R[4]]: tmp.append( __Get_Cluster_Sites( cluster, R, Rlist, ncell, COORD_KEYS ) )
        return tmp

    def __Remove_Duplicate_Clusters(sClusterList):
        tmp = list( sClusterList.keys() )
        for (j, key) in enumerate( tmp ):
            for clsID in sClusterList[key]:
                for q in sClusterList[key][clsID]:
                    a = set( [key] + q )
                    for key1 in tmp[j + 1:]:
                        if clsID in list( sClusterList[key1].keys() ):
                            for (indx, p) in enumerate( sClusterList[key1][clsID] ):
                                if a == set( [key1] + p ): sClusterList[key1][clsID].pop( indx )

        return sClusterList
    sClusterList = {}
    for key, R in config['Coordinates'].items():
        sClusterList[key] = {
            cls.ClusterID: __Clusters_at_R(
                cls, R, config['Coordinates'], config['nCell'], config['Keys'] ) for cls in ClusterList}
    return sClusterList


def Initalize_Correlation_Calculator(Base_Lattice):
    """
    This method used for connecting to drivers and main function.
    :param Base_Lattice:
    :return: ClusterList, ECIList, ThetaRep.
    """
    Translation_Vector, Rotation_Matrix = Load_Symmetry_Operators( Base_Lattice )
    ClusterList, ThetaRep = Load_Clusters_and_ECIs( read_string_input(
        "\nGive the file containing clusters and ECIs (should have the format of output from PyLCE code): " )
    )
    print(
        "\nCreating equivalent clusters in base lattice....might take time (depends on number of clusters)"
    )
    for cluster in ClusterList:
        if cluster.Ksites == 0:
            continue
        else:
            cluster.Equivalent_Clusters( Base_Lattice, Rotation_Matrix )
    ECIList = np.array( [eci for cls in ClusterList for eci in cls.ECIs] )
    for _ in ClusterList:
        if _.Ksites == 0: ClusterList.pop( ClusterList.index( _ ) )

    return ClusterList, ECIList, ThetaRep
