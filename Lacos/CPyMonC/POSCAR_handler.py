import copy
import sys

from Lacos.CPyMonC.Utility import (Cartesian2Fractions,
                     RelFractional_2_Cartesian,
                     Fractions2Cartesian,
                     AtomPosition_Rel2_Base_Lattice)
from Lacos.CPyMonC.Utility import Get_Supercell_Size, Shift_Coordinates_Positive_Quad
from Lacos.CPyMonC.miscfunctions import *


def POSCAR_2_Config(poscar, base_lattice):
    """
    This method used to convert poscar to configuration.
    :param poscar:
    :param base_lattice:
    :return: config: A dictionary type.
    """
    config = {'ElementList': copy.copy( poscar['ElementList'] ), 'AtomNumbers': copy.copy( poscar['AtomNumbers'] ),
              'Lattice_Vectors': copy.copy( poscar['Lattice_Vectors'] ),
              'nCell': Get_Supercell_Size( poscar['Lattice_Vectors'], base_lattice['Lattice_Vectors'] ),
              'Coordinates': AtomPosition_Rel2_Base_Lattice( poscar['Coordinates'], base_lattice )}
    config['Coordinates'] = Shift_Coordinates_Positive_Quad( config['Coordinates'],
                                                             config['nCell'] )

    return config


def Config_2_POSCAR(config, base_lattice):
    """
    This method will  convert config to poscar.
    :param config:
    :param base_lattice:
    :return: poscar: A dictionary type.
    """
    poscar = {}
    if isinstance( config['Coordinates'], dict ):
        pp = [config['Coordinates'][_][3] for _ in config['Coordinates']]
    else:
        pp = [_[3] for _ in config['Coordinates']]
    elmlist = list( set( pp ) )
    atomnumber = [pp.count( _ ) for _ in elmlist]
    poscar['ElementList'], poscar['AtomNumbers'] = list(
      zip( *sorted( r for r in zip( elmlist, atomnumber ) ) )
    )
    poscar['Lattice_Vectors'] = copy.copy( config['Lattice_Vectors'] )
    coords = {_: [] for _ in poscar['ElementList']}

    if isinstance( config['Coordinates'], list ):
        for _ in config['Coordinates']: coords[_[3]].append( _ )
    else:
        for _ in config['Coordinates']: coords[config['Coordinates'][_][3]].append( config['Coordinates'][_] )

    poscar['Coordinates'] = Cartesian2Fractions( config['Lattice_Vectors'], RelFractional_2_Cartesian(
        [_ for elm in poscar['ElementList'] for _ in coords[elm]], base_lattice ) )
    return poscar


def Read_POSCAR(infile):
    """
    This method will Read POSCAR as a dictionary with keys: Lattice_Vectors, ElementList, AtomNumbers and Coordinates
    from input file. :param infile: :return: poscar: A dictionary type.
    """
    inposcar = []
    try:
        ftmp = open( infile, 'r' )
        inposcar = [line.rstrip( '\n' ) for line in ftmp]
    except IOError as e:
        exit_script( e )
    poscar = {'Title': copy.copy( inposcar[0] )}
    pp = [[float( _ ) for _ in inposcar[2].split()], [float( _ ) for _ in inposcar[3].split()],
          [float( _ ) for _ in inposcar[4].split()]]
    poscar['Lattice_Vectors'] = pp
    poscar['ElementList'] = [_ for _ in inposcar[5].split()]
    if any( not ValidElementSymbol( _ ) for _ in poscar['ElementList'] ):
        exit_script( "Error in POSCAR :: unknow Element symbol. Exit" )
    poscar['AtomNumbers'] = [int( _ ) for _ in inposcar[6].split()]

    if 's' in inposcar[8] or 'S' in inposcar[8]:
        start_index = 9
    else:
        start_index = 8
    nelm, j, coord = len( poscar['AtomNumbers'] ), 0, []
    for i in range( nelm ):
        for k in range( poscar['AtomNumbers'][i] ):
            coord.append( [float( _ ) for _ in inposcar[start_index + j].split()] + [poscar['ElementList'][i]] )
            j += 1
    if 'direct' in inposcar[7].lower():
        poscar['Coordinates'] = Fractions2Cartesian( poscar['Lattice_Vectors'], coord )
    else:
        poscar['Coordinates'] = coord

    return poscar


def Write_POSCAR(poscar, title, filename):
    """
    This method is used to write poscar to file.
    :param poscar:
    :param title:
    :param filename:
    :return: None.
    """
    sposcar = [title + '\n', '1.000000' + '\n']
    if (sys.version_info[0] == 2 and sys.version_info[1] == 6) or (
            sys.version_info[0] == 3 and sys.version_info[1] == 0):
        for r in poscar['Lattice_Vectors']:
            sposcar.append( ' '.join( "{0:16.13f}".format( float( _ ) ) for _ in r ) + '\n' )
    else:
        for r in poscar['Lattice_Vectors']:
            sposcar.append( ' '.join( "{:16.13f}".format( float( _ ) ) for _ in r ) + '\n' )

    sposcar.append( ' '.join( r for r in poscar['ElementList'] ) + '\n' )
    sposcar.append( ' '.join( str( _ ) for _ in poscar['AtomNumbers'] ) + '\n' )
    sposcar.append( 'Direct\n' )
    if (sys.version_info[0] == 2 and sys.version_info[1] == 6) or (
            sys.version_info[0] == 3 and sys.version_info[1] == 0):
        for coord in poscar['Coordinates']: sposcar.append(
            ' '.join( "{0:12.9f}".format( float( _ ) ) for _ in coord[:3] ) + '\n' )
    else:
        for coord in poscar['Coordinates']: sposcar.append(
            ' '.join( "{:12.9f}".format( float( _ ) ) for _ in coord[:3] ) + '\n' )

    try:
        with open( filename, 'w' ) as ftmp:
            for r in sposcar: ftmp.write( r )
            ftmp.close()
    except IOError as e:
        print( e )

    return
