"""
/******************************************************************************
  This source file is part of the Avogadro project.
  This source code is released under the New BSD License, (the "License").
******************************************************************************/
"""

import argparse
import json
import sys
from orient import GeometryXYZ,OperationList,StaticReflect
import numpy as np
# Some globals:
debug = True
name = {1: 'H', 2: 'He', 3: 'Li', 4: 'Be', 5: 'B', 6: 'C', 7: 'N', 8: 'O', 9: 'F', 10: 'Ne',
        11: 'Na', 12: 'Mg', 13: 'Al', 14: 'Si', 15: 'P', 16: 'S', 17: 'Cl', 18: 'Ar', 19: 'K', 20: 'Ca',
        21: 'Sc', 22: 'Ti', 23: 'V', 24: 'Cr', 25: 'Mn', 26: 'Fe', 27: 'Ni', 28: 'Co', 29: 'Cu', 30: 'Zn',
        31: 'Ga', 32: 'Ge', 33: 'As', 34: 'Se', 35: 'Br', 36: 'Kr', 37: 'Rb', 38: 'Sr', 39: 'Y', 40: 'Zr',
        41: 'Nb', 42: 'Mo', 43: 'Tc', 44: 'Ru', 45: 'Rh', 46: 'Pd', 47: 'Ag', 48: 'Cd', 49: 'In', 50: 'Sn',
        51: 'Sb', 52: 'Te', 53: 'I', 54: 'Xe', 55: 'Cs', 56: 'Ba', 57: 'La', 58: 'Ce', 59: 'Pr', 60: 'Nd',
        61: 'Pm', 62: 'Sm', 63: 'Eu', 64: 'Gd', 65: 'Tb', 66: 'Dy', 67: 'Ho', 68: 'Er', 69: 'Tm', 70: 'Yb',
        71: 'Lu', 72: 'Hf', 73: 'Ta', 74: 'W', 75: 'Re', 76: 'Os', 77: 'Ir', 78: 'Pt', 79: 'Au', 80: 'Hg',
        81: 'Tl', 82: 'Pb', 83: 'Bi', 84: 'Po', 85: 'At', 86: 'Rn', 87: 'Fr', 88: 'Ra', 89: 'Ac', 90: 'Th',
        91: 'Pa', 92: 'U', 93: 'Np', 94: 'Pu', 95: 'Am', 96: 'Cm', 97: 'Bk', 98: 'Cf', 99: 'Es', 100: 'Fm',
        101: 'Md', 102: 'No', 103: 'Lr', 104: 'Rf', 105: 'Db', 106: 'Sg', 107: 'Bh', 108: 'Hs', 109: 'Mt',
        110: 'Ds', 111: 'Rg', 112: 'Cn', 113: 'Nh', 114: 'Fl', 115: 'Mc', 116: 'Lv', 117: 'Ts', 118: 'Og'}


def getOptions():
    userOptions = {}
        
    userOptions['axis'] = {}
    userOptions['axis']['label'] = 'Reflection across'
    userOptions['axis']['type'] = 'stringList'
    userOptions['axis']['values'] = ['yz plane (x axis)','xz plane (y axis)','xy plane (z axis)']
    userOptions['axis']['default'] = 'yz plane (x axis)'
    opts = {'userOptions': userOptions}
    return opts


def Reflect(opts,mol):
    ax=opts['axis']
    if(ax=='yz plane (x axis)'):
        ax='x'
    elif(ax=='xz plane (y axis)'):
        ax='y'
    elif(ax=='xy plane (z axis)'):
        ax='z'
    atomic_numbers = mol['atoms']['elements']['number']
    selected_atoms  = []
    atom_names=[]
    
    for i in range(len(atomic_numbers)):
        # if mol['atoms']['selected'][i]:
        selected_atoms.append(i)
        atom_names.append(name.get(atomic_numbers[i]))


    # print(selected_atoms)
    # print(selected_coordinates)
    # print(atom_names)
    operation_list = OperationList()
    coordinates_array = np.array(mol['atoms']['coords']['3d']).reshape(-1, 3)
    # print(coordinates_array)
    geometry = GeometryXYZ(
    names=atom_names,
    coordinates=coordinates_array,
)

    # translate_to_origin = CentroidTranslate(selected_atoms, fac=-1.0)
    # operation_list.append(translate_to_origin)
    normal = np.zeros(3)
    normal["xyz".index(ax)] = 1.0

    reflect = StaticReflect(normal)
    # reflect = PlaneReflect(selected_atoms)
    operation_list.append(reflect)
    
    # shift = ShiftedOperation(translate_to_origin, reflect)
    # operation_list.append(shift)

    # print(selected_geometry.coordinates)
    for operation in operation_list:
        operation(geometry.coordinates)
    # Print the modified molecule
    # print(geometry.coordinates)
    mol['atoms']['coords']['3d']=[]
    for i in geometry.coordinates:
        mol['atoms']['coords']['3d'].extend(i)
    return mol


def runCommand():
    # Read options from stdin
    stdinStr = sys.stdin.read()

    # Parse the JSON strings
    opts = json.loads(stdinStr)

    # Prepare the result
    result = {}
    result['cjson'] = Reflect(opts,opts['cjson'])
    return result


if __name__ == "__main__":
    parser = argparse.ArgumentParser('Reflect across plane defined by chosen axis as normal')
    parser.add_argument('--debug', action='store_true')
    parser.add_argument('--print-options', action='store_true')
    parser.add_argument('--run-command', action='store_true')
    parser.add_argument('--display-name', action='store_true')
    parser.add_argument('--menu-path', action='store_true')
    parser.add_argument('--lang', nargs='?', default='en')
    args = vars(parser.parse_args())

    debug = args['debug']

    if args['display_name']:
        print("Reflect across plane defined by chosen axis as normal")
    if args['menu_path']:
        print("&Build")
    if args['print_options']:
        print(json.dumps(getOptions()))
    elif args['run_command']:
        print(json.dumps(runCommand()))
