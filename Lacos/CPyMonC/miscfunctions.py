""" This module is abut some misc methods used by CPyMonC module."""
import os
import os.path
import re
import string

from terminaltables import AsciiTable

List_of_Characters = list( string.ascii_lowercase )
List_of_Characters += list( string.ascii_uppercase )
List_of_Characters += list( string.digits )
List_of_Characters += ['-', '_']

PeriodicTable = ('H', 'He', 'Li', 'Be', 'B', 'C', 'N', 'O', 'F', 'Ne', 'Na', 'Mg', 'Al', 'Si', 'P', 'S', 'Cl',
                 'Ar', 'K', 'Ca', 'Sc', 'Ti', 'V', 'Cr', 'Mn', 'Fe', 'Co', 'Ni', 'Cu', 'Zn', 'Ga', 'Ge', 'As',
                 'Se', 'Br', 'Kr', 'Rb', 'Sr', 'Y', 'Zr', 'Nb', 'Mo', 'Tc', 'Ru', 'Rh', 'Pd', 'Ag', 'Cd', 'In',
                 'Sn', 'Sb', 'Te', 'I', 'Xe', 'Cs', 'Ba', 'Hf', 'Ta', 'W', 'Re', 'Os', 'Ir', 'Pt', 'Au', 'Tl',
                 'Pb', 'Bi', 'Po', 'At', 'Rn', 'La', 'Ce', 'Pr', 'Nd', 'Sm', 'Eu', 'Gd', 'Tb', 'Dy', 'Ho', 'Er',
                 'Tm', 'Yb', 'Lu', 'Ac', 'Th', 'Pa', 'U', 'Zz')


def get_keyed_string(instr):
    try:
        return input( instr )
    except KeyboardInterrupt:
        exit_script( "Quit" )


def create_subdirectory(base_dir, subdir_name, overule_flag):
    """
    This method used to create sub directory.
    :param base_dir:
    :param subdir_name:
    :param overule_flag:
    :return: Subdir
    """
    if base_dir[-1] == '/':
        subdir = base_dir + subdir_name
    else:
        subdir = base_dir + '/' + subdir_name

    if os.path.exists( subdir ):
        if overule_flag == 1:
            return subdir, -1
        else:
            tmpstr = read_yesno_input(
                "\n Found %s dir in the given path %s.\n  Use the same dir (yes|no): " % ( subdir_name, base_dir)
            )
            if tmpstr == 'n': exit_script( "Give a new name to the dir....till then ciao. Bye." )
            return subdir, -1
    else:
        try:
            os.makedirs( subdir )
        except:
            exit_script(
              "Unable to create base dir %s...trouble in paradise. Exiting Script." % subdir
            )
        return subdir, 1


def list_of_sublists(rin, m):
    """
    This method finds list of sub-lists. This also groups object of rin into list of len(m) sublist
    with k-th sublist containing m[k] objets sequentially.
    :param rin:
    :param m:
    :return: P : A list type.
    """

    if len( rin ) != sum( m ):
        print( "Error: length of sublist does not add to length of lists. Exit." )
        exit_script( ' ' )
    p = []
    s = 0
    for i in range( len( m ) ):
        e = s + m[i]
        p.append( rin[s:e] )
        s = e
    return p


def any_special_character(instr):
    """
    This method used to find any special character in the input string.
    :param instr:
    :return: Boolean type.
    """
    if any( r not in List_of_Characters for r in instr ):
        return True
    else:
        return False


def get_full_path(tmpname):
    if '~' in tmpname:
        return os.path.expanduser( tmpname )
    elif './' in tmpname:
        return os.path.realpath( tmpname )
    else:
        return os.path.abspath( tmpname )


def read_yesno_input(instr, Blank=None):
    """
    This method used to read yes no input string. it also allows black space as in put.
    :param instr:
    :param Blank:
    :return: Boolean type
    """
    n = 0
    while True:
        Input_Str = get_keyed_string( instr )
        Input_Str = Input_Str.rstrip()
        if Blank != None and 'blank' in Blank.lower() and Input_Str == '':
            return Input_Str
        if Input_Str.isdigit():
            if int( Input_Str ) == 0:
                return 'n'
            elif int( Input_Str ) == 1:
                return 'y'
            else:
                n += 1
                if n <= 2:
                    print( "   Input either yes|no or 1|0." )
                else:
                    exit_script( "Error: Don't know what you want. Exit" )
        else:
            if 'y' in Input_Str.lower() and 'n' in Input_Str.lower():
                print( "   Input either yes|no or 1|0." )
                n += 1
                continue
            if 'y' in Input_Str.lower():
                return 'y'
            elif 'n' in Input_Str.lower():
                return 'n'
            else:
                n += 1
                if n <= 3:
                    print( "   Input either yes|no or 1|0." )
                else:
                    exit_script( "Error: Don't know what you want. Exit" )


def read_string_input(instr, Blank=None):
    """
    This method will read a string input by the user.
    :param instr:
    :param Blank:
    :return: String type..
    """
    n = 0
    while True:
        Input_Str = get_keyed_string( instr )
        r = Input_Str.rstrip()
        if r != '':
            return r
        else:
            if Blank != None and 'blank' in Blank.lower(): return r
            n += 1
            if n <= 2:
                print( "    Error: Blank space not allowed." )
            else:
                exit_script( "Error: Could not read any string. Exit." )


def read_float_input(instr, Blank=None):
    """
    This method will read float values input by the user.
    :param instr:
    :param Blank:
    :return: A float value.
    """
    n = 0
    while True:
        Input_Str = get_keyed_string( instr )
        if Blank != None and 'blank' in Blank.lower() and Input_Str == '': return None
        try:
            r = float( Input_Str.rstrip() )
            return r
        except:
            n += 1
            if n <= 2:
                print( "    Error: Input digits only." )
            else:
                exit_script( "Error: Unable to read digits. Exit." )


def read_int_input(instr, Blank=None):
    """
    This method will read an integer value given by he user.
    :param instr:
    :param Blank:
    :return: An integer value.
    """
    n = 0
    while True:
        Input_Str = get_keyed_string( instr )
        if Blank != None and 'blank' in Blank.lower() and Input_Str == '': return None
        try:
            r = int( Input_Str.rstrip() )
            return r
        except:
            n += 1
            if n <= 2:
                print( "    Error: Input integer only." )
            else:
                exit_script( "Error: Unable to read integer. Exit." )


def read_string_list_input(instr, m):
    """
    This method will read a string list from the user.
    :param instr:
    :param m:
    :return: A list type.
    """
    n, flag, ps = 0, 0, []
    while True:
        if flag == 0:
            Input_Str = get_keyed_string( instr )
        else:
            Input_Str = get_keyed_string( "  Give %d more inputs: " % (m - len( ps )) )
        Input_Str = [r for r in re.split( r'[;,\s]\s*', Input_Str ) if r]

        flag = 1
        for r in Input_Str:
            try:
                ps.append( r )
            except:
                flag = 0
        if flag == 1:
            if len( ps ) > m:
                print(
                  "  Too many inputs. Picking first %d: %s" % (m, ' '.join( str( r ) for r in ps[:m] ))
                )
                if read_yesno_input( "  Acceptable? (Y|N): " ) == 'y':
                    return ps[:m]
                else:
                    flag, n, ps = 0, n + 1, []
                    print( "  Try once more" )
                    if n == 3: exit_script( "\n  Error: Too many attempts. Exit." )
            elif len( ps ) == m:
                return ps
            else:
                print( "  Need %d inputs. But got %d only. " % (m, len( ps )) )
        else:
            n += 1
            if n <= 2:
                print( "    Error: string only." )
            else:
                exit_script( "\n  Error: non-valid inputs or too many attempts. Exit." )
            flag = 0


def read_int_list_input(instr, m):
    """
    This method will read a list of integers input by the user.
    :param instr:
    :param m:
    :return: A list type.
    """
    n, flag, ps = 0, 0, []
    while True:
        if flag == 0:
            Input_Str = get_keyed_string( instr )
        else:
            Input_Str = get_keyed_string( "  Input %d more integers: " % (m - len( ps )) )
        Input_Str = [r for r in re.split( r'[;,\s]\s*', Input_Str ) if r]
        flag = 1
        for r in Input_Str:
            try:
                ps.append( int( r ) )
            except:
                flag = 0
        if flag == 1:
            if len( ps ) > m:
                print(
                  "  Too many inputs. Picking first %d numbers: %s" % (m, ' '.join( str( r ) for r in ps[:m] ))
                )
                if read_yesno_input( "  Acceptable? (Y|N): " ) == 'y':
                    return ps[:m]
                else:
                    flag, n, ps = 0, n + 1, []
                    print( "  Try once more" )
                    if n == 3: exit_script( "\n  Error: Too many attempts. Exit." )
            elif len( ps ) == m:
                return ps
            else:
                print( "  Need %d integers, got %d. " % (m, len( ps )) )
        else:
            n += 1
            if n <= 2:
                print( "    Error: Input integers only." )
            else:
                exit_script( "\n  Error: no integer inputs or too many attempts. Exit." )
            flag = 0


def read_float_list_input(instr, m):
    """
    This method will read a list of float values.
    :param instr:
    :param m:
    :return: A list type.
    """
    n, flag, ps = 0, 0, []
    while True:
        if flag == 0:
            Input_Str = get_keyed_string( instr )
        else:
            Input_Str = get_keyed_string( "  Input %d more numbers: " % (m - len( ps )) )
        Input_Str = [r for r in re.split( r'[;,\s]\s*', Input_Str ) if r]
        flag = 1
        for r in Input_Str:
            try:
                ps.append( float( r ) )
            except:
                flag = 0
        if flag == 1:
            if len( ps ) > m:
                print( "  Too many inputs. Picking first %d numbers: %s" % (m, ' '.join( str( r ) for r in ps[:m] )) )
                if read_yesno_input( "  Acceptable? (Y|N): " ) == 'y':
                    return ps[:m]
                else:
                    flag, n, ps = 0, n + 1, []
                    print( "  Try once more" )
                    if n == 3: exit_script( "\n  Error: Too many attempts. Exit." )
            elif len( ps ) == m:
                return ps
            else:
                print( "  Need %d numbers, got %d. " % (m, len( ps )) )
        else:
            n += 1
            if n <= 2:
                print( "    Error: Input numbers only." )
            else:
                exit_script( "\n  Error: non-numbers in input field or too many attempts. Exit." )
            flag = 0


def float_sequence(start, end, move):
    """
    This method will form a sequence of float values.
    :param start:
    :param end:
    :param move:
    :return: A list type.
    """
    tmp = [start]
    if end < start:
        m = abs( move )
        while tmp[-1] > end:
            tmp.append( tmp[-1] - m )
        if (tmp[-1] - end) < -1: del tmp[-1]
    else:
        m = abs( move )
        while tmp[-1] < end:
            tmp.append( tmp[-1] + m )
        if (tmp[-1] - end) > 1: del tmp[-1]
    return tmp


def exit_script(instr):
    print( '\n\n>>>>>>  ' + instr + '  <<<<<<\n' )
    raise SystemExit


def write_header(instr, ftmp=None):
    """
    This method used to write a customised header in the log or in terminal.
    :param instr:
    :param ftmp:
    :return: None.
    """
    if ftmp == None:
        print( "\n" )
        print( '-' * 90 )
        print( instr.center( 90 ) )
        print( '-' * 90 + '\n' )
    else:
        ftmp.write( '-' * 90 + '\n' )
        ftmp.write( instr.center( 90 ) + '\n' )
        ftmp.write( '-' * 90 + '\n' )
    return


def print_warning(instr):
    print( '\n>>>>> WARNING :: ' + instr + '  <<<<< ' )
    return


def print_tabular_format(p):
    """
    This method will display the Ascii values in tabular format.
    :param p:
    :return: None.
    """
    table = AsciiTable( p )
    print( table.table )
    return


def print_message(instr):
    print( instr + '\n' )
    return


def write_LogAndDisplay(ftmp, instr, Warn=False):
    """
    This method will write customised log as well as display during execution in the terminal.
    :param ftmp:
    :param instr:
    :param Warn:
    :return: None.
    """
    if Warn:
        print_warning( instr )
    else:
        print( instr )
    ftmp.write( instr + '\n' )
    return


def ValidElementSymbol(x):
    if x in PeriodicTable:
        return True
    else:
        False
