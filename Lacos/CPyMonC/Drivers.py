import sys
import os
from os import path
import numpy as np
import string
import time
import datetime
import copy
from random import sample, shuffle
import bisect
from shutil import copy2
import shutil
from Lacos.CPyMonC.miscfunctions import *
from Lacos.CPyMonC.Energy import ConfigurationEnergy, ConfigurationEnergy_Fast, Convert_Clusters_to_IndexList
from Lacos.CPyMonC.Utility import *
from Lacos.CPyMonC.POSCAR_handler import *
from Lacos.CPyMonC.Composition import *
from Lacos.CPyMonC.Plots import plot_histogram

BASE_DIR__gl = None
EnergyDiff = 1.0e-8


def DirectoryCalc(Base_Lattice, ClusterList, ECIList, ElementRep):
    """
    This method used to calculate the energy of POSCARs in a directory.
    :param Base_Lattice:
    :param ClusterList:
    :param ECIList:
    :param ElementRep:
    :return: None.
    """
    global BASE_DIR__gl
    write_header( 'Calculating energy of structures in directory' )
    base_dir = read_string_input( "Give the base directory path (Default is ./): ", "Allow Blank" )
    BASE_DIR__gl = base_dir
    if base_dir == '': base_dir = './'
    base_dir = get_full_path( base_dir )
    if not os.path.isdir( base_dir ): exit_script( "Directory does not exist. Exit." )

    filebase = read_string_input( "Give base file name to search (Default is POSCAR): ", "Allow Blank" )
    if filebase == '': filebase = 'POSCAR'

    print( "\n Energy of the structures in the database %s:\n" % base_dir )
    MASTER_ID = 0
    counter = 0
    for root, dirs, files in os.walk( base_dir ):
        for f in files:
            if filebase in f and f[0] != '.':
                poscar = Read_POSCAR( os.path.abspath( root ) + '/' + f )
                inConfig = POSCAR_2_Config( poscar, Base_Lattice )
                inConfig['Coordinates'], inConfig['Keys'] = Create_Coordinate_Dictionary( inConfig['Coordinates'] )
                print( os.path.abspath( root ) + '/' + f,
                       ConfigurationEnergy( inConfig, ClusterList, ECIList, ElementRep ) )
                did = os.path.abspath( root )
                MASTER_ID = did[-5:]
                energy = ConfigurationEnergy( inConfig, ClusterList, ECIList, ElementRep )
                idstr = str( datetime.date.today() ) + '-' + str( MASTER_ID ).zfill( 5 )
                lvec = inConfig['Lattice_Vectors']
                atomNumber = inConfig['AtomNumbers']
                Natoms = sum( atomNumber )
                elmlist = inConfig['ElementList']
                cood = poscar['Coordinates']
                cood1 = [item[0:3] for item in cood]

                fp = open( os.path.abspath( root ) + '/' + 'DBinfo.' + idstr, 'a+' )
                fp.write( 'DB Id = ' + idstr + '\n' )
                fp.write( "Run Dir = " + os.path.abspath( root ) + '/' + '\n' )
                fp.write( 'Element List = ' + ' '.join( r for r in elmlist ) + '\n' )
                fp.write( 'Number of atoms = ' + ' '.join( str( r ) for r in atomNumber ) + '\n' )
                fp.write(
                    'Total Energy & energy/atom = ' + ''.join( "{0:14.11f}".format( float( energy * Natoms ) ) ) + str(
                        '  ' ) + ''.join( "{0:14.11f}".format( float( energy ) ) ) + '\n' )
                fp.write( 'System Volume = ' + ''.join( "{0:14.11f}".format( float( energy * 1000 ) ) ) + '\n' )
                fp.write( 'Fermi Energy  = ' + ''.join( "{0:14.11f}".format( float( energy * 0 ) ) ) + '\n' )
                fp.write( 'Band Gap  = ' + ''.join( "{0:14.11f}".format( float( energy * 0 ) ) ) + '\n' )
                fp.write( '\nNumber of Configuration iteration = 1 \n' )
                fp.write( '\n#--------Iteration = 1 ----------\n' )
                fp.write(
                    'Total Energy & energy/atom = ' + ''.join( "{0:14.11f}".format( float( energy * Natoms ) ) ) + str(
                        '  ' ) + ''.join( "{0:14.11f}".format( float( energy ) ) ) + '\n' )
                fp.write( 'System Volume = ' + ''.join( "{0:14.11f}".format( float( energy * 1000 ) ) ) + '\n' )
                fp.write( 'Lattice Vector' + '\n' )
                for kk in range( 3 ):
                    fp.write( '   '.join( "{0:14.11f}".format( float( r ) ) for r in lvec[kk] ) + '\n' )
                fp.write( 'Atom Positions: \n' )
                for kk in range( Natoms ):
                    fp.write( '   '.join( "{0:10.7f}".format( float( r ) ) for r in cood1[kk] ) + '\n' )
                fp.close

    return


def ConfigurationSearch(Base_Lattice, ClusterList, ECIList, ElementRep):
    """
    Calculate energy of random configurations for given system size and composition and
    generate P(E)
    :param Base_Lattice:
    :param ClusterList:
    :param ECIList:
    :param ElementRep:
    :return: None
    """
    def __Shuffle_Atom_Labels(AtomLabels):
        for s in AtomLabels:
            pp = [_[1] for _ in AtomLabels[s]]
            shuffle( pp )
            for (_, elm) in zip( AtomLabels[s], pp ): _[1] = elm
        return AtomLabels

    write_header(
      'Energy Distribution (random walk in configuration space)'
    )
    Ncell = read_int_list_input( "\nSize of the supercell: ", 3 )
    if any( r == 0 for r in Ncell ): exit_script(
      "No dimension can be less than 1. Exit"
    )
    SuperLattice = Create_Supercell_Fractional( Base_Lattice, Ncell, 'rel' )
    SuperLattice['nCell'] = Ncell
    Composition = Read_Composition( Base_Lattice, list( ElementRep.keys() ) )
    SuperLattice = Fix_Composition_in_Supercell( SuperLattice, Base_Lattice,
                                                 Composition )

    MaxConfig = read_int_input(
      "\nNumber of configurations to span (Default is 1e4): ", "Allow Blank"
    )
    if not MaxConfig: MaxConfig = 10000
    nLowEnergyCollect = read_int_input(
        "\nNumber of low-energy configurations (from the lowest) to collect (Default is 10): ", "Allow Blank"
    )
    if not nLowEnergyCollect: nLowEnergyCollect = 10
    SuperLattice['Coordinates'], SuperLattice['Keys'] = Create_Coordinate_Dictionary( SuperLattice['Coordinates'] )
    AtomLabels = {s: [[_, SuperLattice['Coordinates'][_][3]] for _ in SuperLattice['Coordinates'] if
                      SuperLattice['Coordinates'][_][3] in list( Composition[s].keys() )] for s in Composition}

    EnergyList, minEnergy, LowEnergyConfig_List = [], 0.0, []

    SiteBasedClusterList = Convert_Clusters_to_IndexList( SuperLattice, ClusterList )

    for counter in range( MaxConfig ):
        AtomLabels = __Shuffle_Atom_Labels( AtomLabels )
        for s in AtomLabels:
            for _ in AtomLabels[s]: SuperLattice['Coordinates'][_[0]][3] = _[1]
        EnergyList.append(
            ConfigurationEnergy_Fast( SuperLattice, ClusterList, ECIList, ElementRep, SiteBasedClusterList )
        )

        if all( abs( EnergyList[-1] - e[0] ) > EnergyDiff for e in LowEnergyConfig_List ):
            bisect.insort( LowEnergyConfig_List, [EnergyList[-1], copy.deepcopy( SuperLattice )] )
            if len( LowEnergyConfig_List ) > 10: LowEnergyConfig_List.pop( 10 )

        if EnergyList[-1] < minEnergy: minEnergy = copy.deepcopy( EnergyList[-1] )
    try:
        plot_histogram( EnergyList, 'CE-Calculated' )
    except:
        print( "Unable to plot distribution/historgram.\n" )
    try:
        with open( "EnergyList.log", "w" ) as ftmp:
            ftmp.write(
              "List of all energies calculated during configuration search (eV/atom)\n"
            )
            for _ in EnergyList: ftmp.write( "%f " % _ + '\n' )
            ftmp.write( "---------------------------------------------------------------------\n" )
            ftmp.write(
              "Lowest %d energy/energies during the search (eV/atom)\n" % len( LowEnergyConfig_List )
            )
            for _ in LowEnergyConfig_List: ftmp.write( "%f " % _[0] + '\n' )
            ftmp.close()
        print(
            "\n  Energies of all configurations generated, alongwith lowest %d energy/energies, are listed in "
            "EnergyList.log file.\n" % len(LowEnergyConfig_List )
        )
    except IOError as e:
        print( e )

    if len( LowEnergyConfig_List ) > 0:
        print(
            "%d Configurations written to POSCAR.x files, with x increasing with energy (POSCAR.00 for the lowest "
            "energy).\n" % len(LowEnergyConfig_List )
        )
        for (j, _) in enumerate( LowEnergyConfig_List ):
            Write_POSCAR( Config_2_POSCAR( _[1], Base_Lattice ), "POSCAR with " + str( _[0] ),
                          "./POSCAR." + str( j ).zfill( 2 ) )
    return
