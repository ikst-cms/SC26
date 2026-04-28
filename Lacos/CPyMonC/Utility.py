import copy
import numpy as np
from functools import reduce
from operator import mul
from Lacos.CPyMonC.miscfunctions import *

COORD_TOLERANCE = 1.e-3


def Create_Supercell_Fractional(inLattice, N, Flag=None):
    """
    This method will create super cell for fractional.
    :param inLattice:
    :param N:
    :param Flag:
    :return: A dictionary type.
    """
    newcoords = copy.copy( inLattice['Coordinates'] )

    for k in range( 1, N[0] ): newcoords += [[r[0] + k] + r[1:] for r in inLattice['Coordinates']]
    newcoords1 = copy.copy( newcoords )
    for k in range( 1, N[1] ): newcoords += [[r[0]] + [r[1] + k] + r[2:] for r in newcoords1]
    newcoords1 = copy.copy( newcoords )
    for k in range( 1, N[2] ): newcoords += [[r[0], r[1]] + [r[2] + k] + r[3:] for r in
                                             newcoords1]
    if Flag == None:
        AA = [[float( N[0] ), 0, 0], [0, float( N[1] ), 0], [0, 0, float( N[2] )]]
        newcoords = Cartesian2Fractions( AA, newcoords )
    elif 'rel' in Flag.lower():
        pass
    else:
        exit_script( "Error in Create Supercell Fractional." )
    lvec = [list( n * np.array( r ) ) for (n, r) in zip( N, inLattice['Lattice_Vectors'] )]

    if all( _ in list( inLattice.keys() ) for _ in ['AtomNumbers', 'ElementList'] ):
        ntot = reduce( mul, N, 1 )
        newatomnumber_list = [ntot * r for r in inLattice['AtomNumbers']]
        return {'Lattice_Vectors': lvec, 'ElementList': inLattice['ElementList'], 'AtomNumbers': newatomnumber_list,
                'Coordinates': copy.deepcopy( newcoords )}
    else:
        return {'Lattice_Vectors': lvec, 'Coordinates': copy.deepcopy( newcoords )}


def Create_Supercell_Cartesian(inLattice, N):
    """
    This method will create supper cell for cartesian.
    :param inLattice:
    :param N:
    :return: A dictionary type.
    """
    newcoords = [r for r in inLattice['Coordinates']]
    for k in range( 1, N[0] ):
        x = np.array( [k * r for r in inLattice['Lattice_Vectors'][0]] )
        newcoords += [list( np.array( r[:3] ) + x ) + r[3:] for r in inLattice['Coordinates']]

    coords = copy.copy( newcoords )
    for k in range( 1, N[1] ):
        y = np.array( [k * r for r in inLattice['Lattice_Vectors'][1]] )
        newcoords += [list( np.array( r[:3] ) + y ) + r[3:] for r in coords]
    coords = copy.copy( newcoords )
    for k in range( 1, N[2] ):
        z = np.array( [k * r for r in inLattice['Lattice_Vectors'][2]] )
        newcoords += [list( np.array( r[:3] ) + z ) + r[3:] for r in coords]
    lvec = [list( n * np.array( r ) ) for (n, r) in zip( N, inLattice['Lattice_Vectors'] )]

    if all( _ in list( inLattice.keys() ) for _ in ['AtomNumbers', 'ElementList'] ):
        ntot = reduce( mul, N, 1 )
        newatomnumber_list = [ntot * r for r in inLattice['AtomNumbers']]
        return {'Lattice_Vectors': lvec, 'ElementList': inLattice['ElementList'], 'AtomNumbers': newatomnumber_list,
                'Coordinates': newcoords}
    else:
        return {'Lattice_Vectors': lvec, 'Coordinates': newcoords}


def Cartesian2Fractions(latvec, coords):
    """
    This method will transfer the cartesian to fractions coordinates system.
    :param latvec:
    :param coords:
    :return: A list type.
    """
    A = np.transpose( np.array( latvec ) )
    Ainv = np.linalg.inv( A )
    if isinstance( coords[0], list ):
        xcoords = []
        for r in coords: xcoords.append( [round( p, 7 ) for p in Ainv.dot( r[:3] )] + r[3:] )
    else:
        xcoords = list( Ainv.dot( coords[:3] ) ) + coords[3:]
    return xcoords


def Fractions2Cartesian(latvec, coords):
    """
    This method will convert fractions to cartesian coordinate system.
    :param latvec:
    :param coords:
    :return: A dictionary type.
    """
    A = np.transpose( np.array( latvec ) )
    if isinstance( coords[0], list ):
        xcoords = []
        for r in coords: xcoords.append( [round( p, 7 ) for p in A.dot( r[:3] )] + r[3:] )
    else:
        xcoords = list( A.dot( coords[:3] ) ) + coords[3:]
    return xcoords


def Get_Supercell_Size(LatVec, BaseVec):
    """
    This method will return the super cell size for lattice vector.
    :param LatVec:
    :param BaseVec:
    :return:
    """
    p = [np.linalg.norm( np.array( r1 ) ) / np.linalg.norm( np.array( r2 ) ) for (r1, r2) in zip( LatVec, BaseVec )]
    q = [int( round( r ) ) for r in p]
    if any( abs( r1 - r2 ) > 0.1 for (r1, r2) in zip( p, q ) ):
        exit_script(
            "Check the Lattice Vector in input lattice --> it is not a integer multiple of base lattice. Exit."
        )
    if any( r == 0 for r in q ):
        exit_script(
            "Input lattice has a smaller lattice size than base lattice in template file. Exit."
        )
    return q


def AtomPosition_Rel2_Base_Lattice(istruct, base_lattice):
    """
    This method will Converts atomic position in "istruct" relative to reference "base_lattice".
    :param istruct:
    :param base_lattice:
    :return: xcoords: A list type.
    """
    A = np.transpose( np.array( base_lattice['Lattice_Vectors'] ) )
    Ainv = np.linalg.inv( A )
    xcoords = []
    for r in istruct: xcoords.append(
        [round( p, 6 ) for p in Ainv.dot( r[:3] )] + [r[3]] )
    xcoords = Label_Coordinates( xcoords, base_lattice )
    return xcoords


def Create_Coordinate_Dictionary(coords):
    """
    This method will create coordinate dictionary.
    :param coords:
    :return: rlist, COORD_KEY.
    """
    COORD_KEY = [1000, 1000000, 1000000000]
    rlist = {}
    for j, r in enumerate( coords ):
        r0 = [round( _, 3 ) for _ in r[:3]]
        indx = int( sum( np.array( COORD_KEY ) * np.array( r0[:3] ) ) )
        if indx in list( rlist.keys() ): print( "Error in Creating Hash Keys. Collision between ", r0, ' and ',
                                                rlist[indx][:3] )
        rlist[indx] = r
    return rlist, COORD_KEY


def RelFractional_2_Cartesian(istruct, base_lattice):
    """
    This method converts real fractional to cartesian and uses base lattice vectors rather actual lattice vector
    since the coordinates are relative fractional.
    :param istruct:
    :param base_lattice:
    :return: xcoords: A list type.
    """
    A = np.transpose( np.array( base_lattice['Lattice_Vectors'] ) )
    xcoords = []
    for r in istruct:
        xcoords.append( [round( p, 6 ) for p in A.dot( r[:3] )] + [r[3]] )
    return xcoords


def Label_Coordinates(istruct, base_lattice):
    """
    This method used to label the coordinates. This adds an additional label to each coordinate which is index of
    its image in base_lattice
    :param istruct:
    :param base_lattice:
    :return: A list type.
    """

    def __Site_Label(rin, base_lattice_coords):
        r1 = np.round( np.array( rin ) % np.array( [1, 1, 1] ), 6 )
        r1[abs( r1 - np.array( [1, 1, 1] ) ) < COORD_TOLERANCE] = 0.0
        for j, r2 in enumerate( base_lattice_coords ):
            if np.linalg.norm( np.array( r1 ) - np.array(
                    r2[:3] ) ) < COORD_TOLERANCE:
                return [j]
        exit_script( "Unable to identify the coordinate --> Check Symmetry of template structure. Exit." )

    tmp = []
    for r in istruct:
        tmp.append( r + __Site_Label( r[:3], base_lattice[
            'Coordinates'] ) )
    return tmp


def Shift_Coordinates_Positive_Quad(coords, ncell):
    """
    This works for fractional or relative-fractional coordinates (between 0 and 1).
    :param coords:
    :param ncell:
    :return: A list type.
    """
    tmp = []
    for r in coords: tmp.append( list( np.array( r[:3] ) % np.array( ncell ) ) + r[3:] )
    return tmp


def Write_CheckPoint_File(SuperLattice, Frozen_Elements, T):
    """
    This method will write a check point in the file.
    :param SuperLattice:
    :param Frozen_Elements:
    :param T:
    :return: None.
    """
    try:
        fp = open( "CheckPoint.mc", "w" )
    except IOError as e:
        exit_script( "Error: Unable to open CheckPoint File to write. Exit" )
    fp.write( "Temperature = %f" % T + "\n" )
    fp.write( "Lattice_Vectors = \n" )
    for i in range( 3 ):
        fp.write( "  %14.10f %14.10f %14.10f" % (
            SuperLattice['Lattice_Vectors'][i][0], SuperLattice['Lattice_Vectors'][i][1],
            SuperLattice['Lattice_Vectors'][i][2]) + "\n" )

    fp.write(
        "nCell = %d %d %d" % (SuperLattice['nCell'][0], SuperLattice['nCell'][1], SuperLattice['nCell'][2]) + '\n' )
    fp.write( "Element_List = " + " ".join( _ for _ in SuperLattice['ElementList'] ) + "\n" )
    fp.write( "Frozen_Elements = " + " ".join( _ for _ in Frozen_Elements ) + "\n" )
    fp.write( "Atom_Numbers = %d" % SuperLattice['AtomNumbers'][0] + "\n" )
    fp.write( "Coordinates = \n" )
    tmp = list( SuperLattice['Coordinates'].values() )
    for r in tmp:
        fp.write( "  %14.10f %14.10f %14.10f  " % (r[0], r[1], r[2]) + r[3] + "  " + str( r[4] ) + "\n" )

    return


def Load_Template_Structure():
    """
    This method will read the template structure file.
    :return: A dictionary type.
    """
    infile = read_string_input(
        "\nGive the structure template file (MUST be the unit cell and can contain labeled coordinates\n- "
        "alphabetically -  to specify sites): "
    )
    try:
        fp = open( infile, "r" )
    except IOError as e:
        exit_script( "Error: template file <%s> does not exist. Exit" % infile )
    lines = [r for r in fp]
    k, lat = 0, []
    while k < len( lines ):
        if all( substr in lines[k].lower() for substr in ['lattice', 'vectors'] ):
            for kk in range( 1, 4 ): lat.append( [float( r ) for r in (lines[k + kk].rstrip( '\n' )).split()] )
            break
        k += 1
    if not lat: exit_script( "Error in template file %s. Unable to find Lattice vectors." % infile )
    k, coords = 0, []
    while k < len( lines ):
        if 'coordinates' in lines[k].lower():
            for l in lines[k + 1:]:
                if l[0] == '#': break
                if l == "": continue
                p = [r for r in (l.rstrip( '\n' )).split()]
                coords.append( [float( r ) for r in p[:3]] + p[3:] )
            break
        k += 1
    if not coords: exit_script( "Error in template file %s. Unable to find coordinates." % infile )

    tmp = []
    for _ in coords:
        if len( _ ) > 3: tmp.append( _[3] )
    if tmp != []:
        tmp = list( set( tmp ) )
        _, k = list( string.ascii_lowercase ), 0
        while True:
            if _[k] in tmp:
                k += 1
            else:
                tmp.append( _[k] )
                break
        coords = [_[3].append( tmp[-1] ) if len( _ ) == 3 else _ for _ in coords]

    Shift_Coordinates_Positive_Quad( coords, [1, 1, 1] )
    rcoords = Fractions2Cartesian( lat, coords )

    return {'Lattice_Vectors': lat, 'Coordinates': coords, 'Real_Coordinates': rcoords}
