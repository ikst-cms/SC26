""" This module is to Finds symmetry of the user defined lattice using Cogue and identifies all SG transformation matrices. """

import sys
import os
import os.path
import string
from math import *
import numpy as np
from Lacos.CPyMonC.miscfunctions import *
sys.path.append (sys.path[0 ] +'/.. ')
import cogue.interface.vasp_io as vio
import cogue.crystal.symmetry as symm


def Load_Symmetry_Operators (lattice):
    """
    This method used to load the symmetry operators and convert the base lattice to POSCAR format which will be
    required for cogue.
    :param lattice:
    :return: translations, rotations: Both are list type.
    """
    write_header ("Picking Symmetry operators" )
    poscar = convert_lattice_to_poscar( lattice )
    cell = vio.parse_poscar( poscar )
    fullinfo = symm.get_symmetry_dataset( cell )
    translations, rotations = [], []
    for r in fullinfo['translations']: translations.append( list( r ) )
    for r in fullinfo['rotations']: rotations.append( r.tolist() )
    lattice['ElementList'] = list( set( fullinfo['wyckoffs'] ) )
    lattice['AtomNumbers'] = [fullinfo['wyckoffs'].count( r ) for r in lattice['ElementList']]
    tmp = []
    for _ in lattice['Coordinates']:
        if len( _ ) > 3: tmp.append( _[3] )
    if tmp == []:
        for k, row in enumerate( lattice['Coordinates'] ): lattice['Coordinates'][k] += fullinfo['wyckoffs'][k]
    return translations, rotations


def convert_lattice_to_poscar(lattice):
    """
    Thsi method will convert lattice to POSCAR.
    :param lattice:
    :return: poscar: A list type.
    """
    poscar = ['Tempfile', '1.0\n']
    for j in range( 3 ):
        poscar.append( ' '.join( str( r ) for r in lattice['Lattice_Vectors'][j] ) + '\n' )
    poscar.append( 'Lr' + '\n' )
    poscar.append( str( len(
        lattice['Coordinates'] ) ) + '\n' )
    poscar.append( 'Direct' )
    for row in lattice['Coordinates']:
        poscar.append( ' '.join( str( r ) for r in row ) + '\n' )

    return poscar
      
