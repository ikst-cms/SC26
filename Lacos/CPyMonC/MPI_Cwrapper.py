import ctypes
import operator
import sys
from functools import reduce
from Lacos.CPyMonC.Utility import *
import psutil
from Lacos import CPyMonC
import os
import glob

try:
    MPI_C_MCFunc = ctypes.cdll.LoadLibrary(glob.glob(os.path.join(os.path.dirname(CPyMonC.__file__),"lib/MPI_C_MCinterface.*.so"))[0])
except ImportError:
    exit_script( "Error in Cwrapper.py. Unable to load C_MCinterface.so. Exit." )


class C_Coordinate( ctypes.Structure ):
    _fields_ = [('x', ctypes.c_float), ('y', ctypes.c_float), ('z', ctypes.c_float)]


class C_LabeledCoordinate( ctypes.Structure ):
    _fields_ = [('R', C_Coordinate), ('Label', ctypes.c_char_p)]


class C_ElementRep( ctypes.Structure ):
    _fields_ = [('Name', ctypes.c_char_p), ('Theta', ctypes.POINTER( ctypes.c_double ))]


class C_AtomLabeledCoordinate( ctypes.Structure ):
    _fields_ = [('R', C_Coordinate), ('Atom', ctypes.c_char_p), ('iAtom', ctypes.c_int), ('iBaseLat', ctypes.c_int)]


class C_CharList( ctypes.Structure ):
    _fields_ = [('p', ctypes.c_char_p)]


class C_Cluster( ctypes.Structure ):
    _fields_ = [('R', ctypes.POINTER( C_Coordinate ))]


class C_ChildGroup( ctypes.Structure ):
    _fields_ = [('Child', ctypes.POINTER( C_Cluster )), ('nChild', ctypes.c_int), ('BL_Index', ctypes.c_int)]


class C_Indexed_ECI( ctypes.Structure ):
    _fields_ = [('ECI', ctypes.c_double), ('Index', ctypes.POINTER( ctypes.c_int ))]


class C_ClusterFamily( ctypes.Structure ):
    _fields_ = [('Ksites', ctypes.c_int), ('ClusterID', ctypes.c_char_p), ('Parent', C_Cluster),
                ('ChildGroup_List', ctypes.POINTER( C_ChildGroup )), ('nChildGroups', ctypes.c_int),
                ('ECI_List', ctypes.POINTER( C_Indexed_ECI )), ('nECIs', ctypes.c_int)]


def MPI_C_MonteCarlo(LogFile_Name, OutputFile_Name, Base_Lattice, ElementRep, SuperLattice, Frozen_Elements,
                 Temperature_Sequence, Equilibrium_Steps, Equilibrium_Tolerance, nSTEP_MOVAVG,
                 MinEquil_Steps, Calculation_Steps, Block_Size, ClusterList, ECIList, numList, listCap):
    """
    This is the method used for the main function calling the C-code for Monte Carlo simulation.First pass all
    parameters and inputs to the C-code to create global arrays. The Call_C_MonteCarlo is where the actual
    Monte-Carlo simulation is performed.
    :param LogFile_Name:
    :param OutputFile_Name:
    :param Base_Lattice:
    :param ElementRep:
    :param SuperLattice:
    :param Frozen_Elements:
    :param Temperature_Sequence:
    :param Equilibrium_Steps:
    :param Equilibrium_Tolerance:
    :param nSTEP_MOVAVG:
    :param MinEquil_Steps:
    :param Calculation_Steps:
    :param Block_Size:
    :param ClusterList:
    :param ECIList:
    :return: None.
    """
    Enter_C_MonteCarlo_and_Send_Elements( ElementRep, LogFile_Name )
    Send_Base_Lattice_to_C( Base_Lattice )
    Send_SuperLattice_to_C( SuperLattice, Frozen_Elements )
    Send_MC_Parameters_to_C( Temperature_Sequence, Equilibrium_Steps, Equilibrium_Tolerance, nSTEP_MOVAVG,
                             Calculation_Steps, Block_Size, MinEquil_Steps )
    Send_ClusterInfo_to_C( ClusterList, ECIList[0] )
    Send_Exchange_Dist_Output_File_to_C( OutputFile_Name, numList, listCap )
    return


def Enter_C_MonteCarlo_and_Send_Elements(ElementRep, LogFile_Name):
    """
    This method will enter to C code of monte carlo and send the elements.
    :param ElementRep:
    :param LogFile_Name:
    :return: None.
    """
    MPI_C_MCFunc.First_Entry_and_Receive_Elements.argtypes = (
      ctypes.POINTER( C_ElementRep ), ctypes.c_int, ctypes.c_char_p
    )
    MPI_C_MCFunc.First_Entry_and_Receive_Elements.restype = ctypes.c_void_p
    p2c_elementrep = (C_ElementRep * len( ElementRep ))()
    for (j, elm) in enumerate( ElementRep ):
        p2c_elementrep[j].Name = ctypes.c_char_p( elm.encode( 'utf8' ) )
        p2c_elementrep[j].Theta = (ctypes.c_double * len( ElementRep[elm] ))( *ElementRep[elm] )
    p2c_logfilename = ctypes.c_char_p( LogFile_Name.encode( 'utf8' ) )

    MPI_C_MCFunc.First_Entry_and_Receive_Elements( p2c_elementrep, len( ElementRep ), p2c_logfilename )

    return


def Send_Base_Lattice_to_C(Base_Lattice):
    """
    This method will send the base lattice to c code.
    :param Base_Lattice:
    :return: None.
    """
    MPI_C_MCFunc.Receive_BaseLattice.argtypes = (
    ctypes.POINTER( C_LabeledCoordinate ), ctypes.c_int, ctypes.c_char_p, ctypes.c_int, ctypes.POINTER( C_Coordinate )
    )
    MPI_C_MCFunc.Receive_BaseLattice.restype = ctypes.c_void_p
    p2c_coordarray = (C_LabeledCoordinate * len( Base_Lattice['Coordinates'] ))()
    for j, r in enumerate( Base_Lattice['Coordinates'] ):
        p2c_coordarray[j].R.x, p2c_coordarray[j].R.y, p2c_coordarray[j].R.z = r[0], r[1], r[2]
        p2c_coordarray[j].Label = ctypes.c_char_p( r[3].encode( 'utf8' ) )  # r[3]

    labels = list( set( [r[3] for r in Base_Lattice['Coordinates']] ) )
    labels = ''.join( _ for _ in labels )
    p2c_labels = ctypes.c_char_p( labels.encode( 'utf8' ) )

    p2c_latvec = (C_Coordinate * 3)()
    for (j, lv) in enumerate( Base_Lattice['Lattice_Vectors'] ):
        p2c_latvec[j].x, p2c_latvec[j].y, p2c_latvec[j].z = lv[0], lv[1], lv[2]

    MPI_C_MCFunc.Receive_BaseLattice( p2c_coordarray, len( Base_Lattice['Coordinates'] ), p2c_labels, len( labels ),
                                  p2c_latvec )

    return


def Send_SuperLattice_to_C(SuperLattice, Frozen_Elements):
    """
    This method used to send the super lattice to c code.
    :param SuperLattice:
    :param Frozen_Elements:
    :return: None.
    """
    MPI_C_MCFunc.Receive_SuperLattice.argtypes = (
    ctypes.POINTER( C_AtomLabeledCoordinate ), ctypes.c_int, ctypes.POINTER( C_Coordinate ),
    ctypes.POINTER( ctypes.c_int ), ctypes.POINTER( C_CharList ), ctypes.c_int)
    MPI_C_MCFunc.Receive_SuperLattice.restype = ctypes.c_void_p

    p2c_coordarray = (C_AtomLabeledCoordinate * len( SuperLattice['Coordinates'] ))()
    for j, r in enumerate( SuperLattice['Coordinates'] ):
        p2c_coordarray[j].R.x, p2c_coordarray[j].R.y, p2c_coordarray[j].R.z = r[0], r[1], r[2]
        p2c_coordarray[j].Atom = ctypes.c_char_p( r[3].encode( 'utf8' ) )
        p2c_coordarray[j].iAtom = -1
        p2c_coordarray[j].iBaseLat = r[4]

    p2c_latvec = (C_Coordinate * 3)()
    for (j, lv) in enumerate( SuperLattice['Lattice_Vectors'] ):
        p2c_latvec[j].x, p2c_latvec[j].y, p2c_latvec[j].z = lv[0], lv[1], lv[2]

    p2c_ncell = (ctypes.c_int * len( SuperLattice['nCell'] ))( *SuperLattice['nCell'] )
    p2c_frozen_elements = (C_CharList * len( Frozen_Elements ))()
    for j, elm in enumerate( Frozen_Elements ):
        p2c_frozen_elements[j].p = ctypes.c_char_p( elm.encode( 'utf8' ) )

    MPI_C_MCFunc.Receive_SuperLattice( p2c_coordarray, len( SuperLattice['Coordinates'] ), p2c_latvec, p2c_ncell,
                                   p2c_frozen_elements, len( Frozen_Elements ) )

    return


def Send_MC_Parameters_to_C(Temperature_Sequence, Equilibrium_Steps, Equilibrium_Tolerance, nSTEP_MOVAVG,
                            Calculation_Steps, Block_Size, MinEquil_Steps):
    """
    This method will send monte carlo parameters to c- code.
    :param Temperature_Sequence:
    :param Equilibrium_Steps:
    :param Equilibrium_Tolerance:
    :param nSTEP_MOVAVG:
    :param Calculation_Steps:
    :param Block_Size:
    :param MinEquil_Steps:
    :return: None.
    """
    MPI_C_MCFunc.Receive_MC_Parameters.argtypes = (
    ctypes.POINTER( ctypes.c_float ), ctypes.c_int, ctypes.c_int, ctypes.c_float, ctypes.c_int, ctypes.c_int,
    ctypes.c_int, ctypes.c_int
    )
    MPI_C_MCFunc.Receive_MC_Parameters.restype = ctypes.c_void_p
    p2c_tlist = (ctypes.c_float * len( Temperature_Sequence ))( *Temperature_Sequence )
    if Equilibrium_Steps is None:
        p2c_eqsteps = -1
    else:
        p2c_eqsteps = Equilibrium_Steps

    if Equilibrium_Tolerance is None:
        p2c_eqtol = -1.0
    else:
        p2c_eqtol = Equilibrium_Tolerance

    p2c_nh = nSTEP_MOVAVG
    p2c_calcsteps = Calculation_Steps
    p2c_blocksize = Block_Size
    p2c_minequilsteps = MinEquil_Steps

    MPI_C_MCFunc.Receive_MC_Parameters( p2c_tlist, len( Temperature_Sequence ), p2c_eqsteps, p2c_eqtol, p2c_nh,
                                    p2c_calcsteps, p2c_blocksize, p2c_minequilsteps )

    return


def Send_ClusterInfo_to_C(ClusterList, ECI0):
    """
    This method will send cluster information to C code.
    :param ClusterList:
    :param ECI0:
    :return: None.
    """
    MPI_C_MCFunc.Receive_ClusterInfo.argtypes = (ctypes.POINTER( C_ClusterFamily ), ctypes.c_int, ctypes.c_double)
    MPI_C_MCFunc.Receive_ClusterInfo.restype = ctypes.c_void_p
    p2c_cflist = (C_ClusterFamily * len( ClusterList ))()
    for j, cls in enumerate( ClusterList ):
        p2c_cflist[j].Ksites = cls.Ksites
        p2c_cflist[j].ClusterID = ctypes.c_char_p( cls.ClusterID.encode( 'utf8' ) )
        p2c_cflist[j].nECIs = len( cls.ECIs )
        p2c_cflist[j].ECI_List = (C_Indexed_ECI * len( cls.ECIs ))()
        for k in range( len( cls.ECIs ) ):
            p2c_cflist[j].ECI_List[k].ECI = cls.ECIs[k]
            p2c_cflist[j].ECI_List[k].Index = (ctypes.c_int * len( cls.CIndexList[k] ))( *cls.CIndexList[k] )

        p2c_cflist[j].Parent.R = (C_Coordinate * cls.Ksites)()
        for k, R in enumerate( cls.Parent ):
            p2c_cflist[j].Parent.R[k].x, p2c_cflist[j].Parent.R[k].y, p2c_cflist[j].Parent.R[k].z = R[0], R[1], R[2]

        p2c_cflist[j].nChildGroups = len( cls.Child )
        p2c_cflist[j].ChildGroup_List = (
                    C_ChildGroup * len( cls.Child ))()
        for k, jsite in enumerate( cls.Child ):
            p2c_cflist[j].ChildGroup_List[k].Child = (C_Cluster * len(
                cls.Child[jsite] ))()
            p2c_cflist[j].ChildGroup_List[k].BL_Index = jsite
            p2c_cflist[j].ChildGroup_List[k].nChild = len( cls.Child[jsite] )
            for kk, child_cls in enumerate( cls.Child[jsite] ):
                p2c_cflist[j].ChildGroup_List[k].Child[kk].R = (C_Coordinate * cls.Ksites)()
                for t, R in enumerate( child_cls ):
                    p2c_cflist[j].ChildGroup_List[k].Child[kk].R[t].x = R[0]
                    p2c_cflist[j].ChildGroup_List[k].Child[kk].R[t].y = R[1]
                    p2c_cflist[j].ChildGroup_List[k].Child[kk].R[t].z = R[2]

    MPI_C_MCFunc.Receive_ClusterInfo( p2c_cflist, len( ClusterList ), ECI0 )
    return


def Send_Exchange_Dist_Output_File_to_C( OutputFile_Name, numList, listCap):
    """
    This method will send neighbour list info to c code.
    :return: None.
    """
    MPI_C_MCFunc.Get_Exchange_Dist_and_Output_File.argtypes = (ctypes.c_float, ctypes.c_char_p, ctypes.c_int, ctypes.c_int)
    MPI_C_MCFunc.Get_Exchange_Dist_and_Output_File.restype = ctypes.c_void_p
    ext_minDist = read_float_input(
        "\n\nGive distance within which to exchange atoms, calculated using fractional coordinates (for default "
        "Nearest-neighbor, just \"Enter\"):", "Allow Blank"
    )
    if ext_minDist is None: ext_minDist = -1.0
    p2c_outfilename = ctypes.c_char_p( OutputFile_Name.encode( 'utf8' ) )
    MPI_C_MCFunc.Get_Exchange_Dist_and_Output_File( ext_minDist, p2c_outfilename, numList, listCap )
    return
