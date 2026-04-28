"""This module is about the method  related to fix the compositions in a system"""

from random import shuffle
from Lacos.CPyMonC.Utility import Label_Coordinates
from Lacos.CPyMonC.miscfunctions import *


def Read_Composition(Base_Lattice, ElementList):
    """
  This method reads the composition from the command line.
    :param Base_Lattice:
    :param ElementList:
    :return: A dictionary type.
    """
    write_header( 'Input Composition' )
    Site_Labels = list( set( [_[3] for _ in Base_Lattice['Coordinates']] ) )
    if len( Site_Labels ) == 1:
        print(
            "Following site-labels found (if not in template file, then this is Wyckoff labels): %s" % Site_Labels[0]
        )
    else:
        print( "Following site-labels identified (if not in template, then these are Wyckoff labels): %s" % ", ".join(
            _ for _ in Site_Labels )
               )

    print(
        "\nGive composition (site fraction) for each label (in template file) or Wyckoff label (each element to be "
        "separated by , or ;) \n(e.g., For element a: Al:0.3, Ni:0.7)"
    )
    tmpComp = {}
    for s in Site_Labels:
        while True:
            pp = read_string_input( "  For site-label %s: " % s )
            pp = [_ for _ in re.split( r'[;,\s]\s*', pp.lstrip() ) if _ != '' and _ != ':']
            if len( pp ) >= 1: break
        tmp = [_ for x in pp for _ in x.split( ':' ) if _ is not '']
        if any( not re.match( "[0-9.]+", tmp[_] ) for _ in range( 1, len( tmp ), 2 ) ): exit_script(
            "Error in giving composition. Check." )
        pp = [tmp[_] + ':' + tmp[_ + 1] for _ in range( 0, len( tmp ), 2 )]
        if any( _.split( ':' )[0] not in ElementList for _ in pp ):
            exit_script(
                "  Error :: Cannot handle elements that are not in ECI file obtained from Cluster Expansion."
            )
        if len( pp ) == 1:
            print(
                "   WARNING :: Element %s will be frozen during calculation.\n" % pp[0].split( ':' )[0]
            )
            if abs( float( pp[0].split( ':' )[1] ) - 1.0 ) > 1e-3:
                exit_script(
                    "   Error :: Composition for each site should add to 1.0 (within 0.001 tolerance). Exit"
                )
            else:
                tmpComp[s] = {pp[0].split( ':' )[0]: float( pp[0].split( ':' )[1] )}
        else:
            tmpComp[s] = {_.split( ':' )[0]: float( _.split( ':' )[1] ) for _ in pp}
            if abs( sum( tmpComp[s][_] for _ in tmpComp[s] ) - 1.0 ) > 1e-3:
                exit_script(
                    "   Error :: Composition for each site should add to 1.0 (within 0.001 tolerance). Exit"
                )
    return tmpComp


def Fix_Composition_in_Supercell(inLattice, base_lattice, Composition):
    """
  The method fixes the composition in the input lattice.
    :param inLattice:
    :param base_lattice:
    :param Composition:
    :return: inLattice
    """
    tmpComp = Fraction_2_AtomNumber( Composition, inLattice )
    inLattice['Coordinates'] = Label_Coordinates( Assign_Atom_to_Site( inLattice['Coordinates'], tmpComp ),
                                                  base_lattice )
    inLattice['ElementList'] = list( set( _[3] for _ in inLattice['Coordinates'] ) )
    return inLattice


def Fraction_2_AtomNumber(Composition, inLattice):
    """
  This method onverts the atomic fraction in composition to atom number.
    :param Composition:
    :param inLattice:
    :return: A dictionary type.
    """
    import numpy as np
    p = [_[3] for _ in inLattice['Coordinates']]
    Site_Count = {k: p.count( k ) for k in Composition}
    NatomComp = {}
    for s in Composition:
        NatomComp[s] = {}
        for elm in Composition[s]:
            NatomComp[s][elm] = int( round( Composition[s][elm] * Site_Count[s], 2 ) )
        diff = Site_Count[s] - sum( NatomComp[s][_] for _ in NatomComp[s] )
        if abs( diff ) != 0:
            for elm in NatomComp[s]:
                NatomComp[s][elm] += 1 * np.sign( diff )
                diff -= np.sign( diff )
                if diff == 0: break
    print( '\nDistributing following number of atoms on the lattice:' )
    for _ in NatomComp:
        print(
            '  at site %s: ' % _, ", ".join( elm + ':' + str( NatomComp[_][elm] ) for elm in NatomComp[_] )
        )
    return NatomComp


def Assign_Atom_to_Site(Coord_List, tmpComp):
    """
  This method used to assign the atoms to lattice sites
    :param Coord_List:
    :param tmpComp:
    :return: A list type.
    """
    for s in tmpComp:
        AtomLabels = []
        for elm in tmpComp[s]:
            AtomLabels += [elm] * tmpComp[s][elm]
        for _ in range( 10 ): shuffle( AtomLabels )
        k = 0
        for _ in Coord_List:
            if _[3] == s:
                _[3] = AtomLabels[k]
                k += 1
    return Coord_List


def List_Frozen_Elements(Composition):
    """
   This method Takes the list of elements that needs to be frozen/defrozen during Monte Carlo simulations
   based on the condition if certain sites are active or not.
    :param Composition:
    :return: A list type.
    """
    tmp = []
    for s in Composition:
        if len( Composition[s] ) == 1: tmp.append( list( Composition[s].keys() )[0] )

    if len( tmp ) > 0:
        if len( tmp ) == 1:
            print(
              "\n %s is frozen on the lattice during MonteCarlo (since at. fraction requested is 1.0)" % tmp[0]
            )
        elif len( tmp ) > 1:
            print(
                "\n Elements %s are frozen on the lattice during MonteCarlo (since at. fraction requested for these "
                "is 1.0)" % ", ".join( _ for _ in tmp )
            )
        defreeze_list = read_string_input(
            " If you want to freeze the frozen elements, list them now, else \"Enter\": ", "Allow Blank"
        )
        pp = [_ for _ in re.split( r'[;,\s]\s*', defreeze_list.lstrip() ) if _ != '']
        for _ in pp:
            if _ in tmp: del tmp[tmp.index( _ )]

    freeze_list = read_string_input(
        "\nInput any element to be frozen during MonteCarlo simulations? (\"Enter\" for none): ", "Allow Blank"
    )
    if freeze_list == "":
        return tmp
    else:
        pp = [_ for _ in re.split( r'[;,\s]\s*', freeze_list.lstrip() ) if _ != '']
        for _ in pp:
            if _ not in tmp: tmp.append( _ )
    return tmp
