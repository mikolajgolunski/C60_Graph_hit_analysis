import csv
import time
import copy
import gzip
import numpy as np

from tqdm import tqdm

from useful_functions import *


def readLammpstrj(data_file):
    # read data from lammpstrj file

    print('Reading from .lammpstrj file')

    csv.register_dialect('lammpstrj', delimiter=' ', strict=True)

    lammpstrj = []

    with gzip.open(data_file, 'rt', newline='') as infile:
        reader = csv.reader(infile, dialect='lammpstrj')

        # initialization code
        frame_nr = 0
        lammpstrj.append({'frame_nr': frame_nr, 'atoms': {}})
        atoms = lammpstrj[-1]['atoms']
        part = None
        header = None

        pbar = tqdm(leave=False, unit='line', desc='lammpstrj file')

        for line_previous, line, line_next in previousAndNext(reader):
            pbar.update()
            if part == 'NUMBER OF ATOMS':
                lammpstrj[-1]['atoms_nr'] = int(line[0])
                part = None
            elif part == 'ATOMS':
                atom = {}
                if 'id' in header:
                    atom['id'] = int(line[header.index('id')])
                if 'element' in header:
                    atom['element'] = line[header.index('element')]
                if 'q' in header:
                    atom['q'] = float(line[header.index('q')])
                if 'mass' in header:
                    atom['mass'] = float(line[header.index('mass')])
                if 'c_sput_ke' in header:
                    atom['ke_lammps'] = float(line[header.index('c_sput_ke')])
                if 'c_sput_pe' in header:
                    atom['pe_lammps'] = float(line[header.index('c_sput_pe')])
                if 'x' in header and 'y' in header and 'z' in header:
                    atom['position'] = np.matrix(
                        [
                            line[header.index('x')],
                            line[header.index('y')],
                            line[header.index('z')]
                        ], dtype=np.dtype('f8')
                    )
                if 'vx' in header and 'vy' in header and 'vz' in header:
                    atom['velocity'] = np.matrix(
                        [
                            line[header.index('vx')],
                            line[header.index('vy')],
                            line[header.index('vz')]
                        ], dtype=np.dtype('f8')
                    )

                atoms[atom['id']] = atom

                if line_next is None or line_next[0] == 'ITEM:':
                    frame_nr += 1
                    lammpstrj.append({'frame': frame_nr, 'atoms': {}})
                    atoms = lammpstrj[-1]['atoms']
                    part = None
            else:
                if line[0] == 'ITEM:':
                    part = ' '.join([token for token in line[1:] if token.isupper()])
                    header = [token for token in line[1:] if token.islower()]
        if len(lammpstrj[-1]['atoms']) == 0:
            del (lammpstrj[-1])

        pbar.close()
        time.sleep(0.1)
    return lammpstrj


def readRdx(data_file):
    print('Reading from .rdx file')

    csv.register_dialect('rdx', delimiter=' ', strict=True)

    with gzip.open(data_file, 'rt', newline='') as infile:
        reader = csv.reader(infile, dialect='rdx')

        # initialization code
        part = None
        search_for_content_q = False
        masses = {}
        atoms = {}
        velocities = {}

        pbar = tqdm(leave=False, unit='line', desc='rdx file')

        for line in reader:
            pbar.update()

            if len(line) == 0:
                if search_for_content_q:
                    continue
                else:
                    part = None
            else:
                if part is not None:
                    search_for_content_q = False
                    if part == 'MASSES':
                        masses[int(line[0])] = float(line[1])
                    elif part == 'ATOMS':
                        atom = {
                            'id': int(line[0]),
                            'type': int(line[1]),
                            'q': float(line[2]),
                            'position': np.matrix(line[3:6], dtype=np.dtype('f8'))
                        }
                        atoms[atom['id']] = atom
                    elif part == 'VELOCITIES':
                        velocities[int(line[0])] = np.matrix(line[1:], dtype=np.dtype('f8'))
                else:
                    if line[0] == 'Masses':
                        part = 'MASSES'
                        search_for_content_q = True
                    elif line[0] == 'Atoms':
                        part = 'ATOMS'
                        search_for_content_q = True
                    elif line[0] == 'Velocities':
                        part = 'VELOCITIES'
                        search_for_content_q = True

        pbar.close()
        time.sleep(0.1)

    rdx = atoms
    if len(masses) != 0:
        if 'mass' not in next(iter(rdx.values())).keys():
            for atom in rdx.values():
                atom['mass'] = masses[atom['type']]
    if len(velocities) != 0:
        if 'velocity' not in next(iter(rdx.values())).keys():
            for atom in rdx.values():
                atom['velocity'] = velocities[atom['id']]

    return rdx


def readBonds(data_file):
    print('Reading from bonds file')

    csv.register_dialect('bonds', delimiter=' ', strict=True, skipinitialspace=True)

    with gzip.open(data_file, 'rt', newline='') as infile:
        reader = csv.reader(infile, dialect='bonds')

        # initialization code
        bonds = {}

        pbar = tqdm(leave=False, unit='line', desc='bonds file')

        for line in reader:
            pbar.update()
            if line[0] == '#':
                continue
            else:
                bond = {
                    'from_id': int(line[0]),
                    'from_type': int(line[1]),
                    'to_nr': int(line[2]),
                    'to_ids': [int(to_id) for to_id in line[3:3+int(line[2])]]
                }
                bonds[bond['from_id']] = bond

        pbar.close()
        time.sleep(0.1)

    return bonds


def createMolecules(bonds):
    print('Creating molecules\' list from bonds file')

    molecules = []

    pbar = tqdm(total=len(bonds), leave=False, unit='atom', desc='creating molecules')

    for bond in bonds.values():
        pbar.update()

        ids = {bond['from_id']}.union(bond['to_ids'])
        for molecule in molecules:
            if len(ids.intersection(molecule)) != 0:
                molecule.update(ids)
                break
        else:
            molecules.append(ids)

    pbar.close()
    time.sleep(0.1)

    return molecules


def unionAtoms(atoms_first, atoms_second):
    print('Unifying atoms\' descriptions')

    atoms = copy.deepcopy(atoms_first)

    keys_second = set(next(iter(atoms_second.values())).keys())
    keys_first = set(next(iter(atoms_first.values())).keys())
    diff = keys_second.difference(keys_first)

    for atom in atoms.values():
        for key in diff:
            atom[key] = atoms_second[atom['id']][key]

    return atoms


def calculateAtomEnergies(atoms_in):
    print('Calculating atoms\' energies')

    atoms = copy.deepcopy(atoms_in)

    for atom in atoms.values():
        atom['ke'] = (atom['mass'] * np.linalg.norm(atom['velocity']))/2

    return atoms


def calculateCOMProperties(molecules, atoms):
    print('Calculating COM properties')

    molecules = copy.deepcopy(molecules)
    molecules_out = []

    pbar = tqdm(total=len(molecules), leave=False, unit='molecule', desc='calculating COM properties')

    for i, molecule in enumerate(molecules):
        pbar.update()

        mass = 0
        com_position = np.matrix([0., 0., 0.])
        com_velocity = np.matrix([0., 0., 0.])
        for mol_id in molecule:
            mass += atoms[mol_id]['mass']
            com_position = com_position + atoms[mol_id]['position'] * atoms[mol_id]['mass']
            com_velocity = com_velocity + atoms[mol_id]['velocity'] * atoms[mol_id]['mass']
        com_position = com_position / mass
        com_velocity = com_velocity / mass

        molecule_out = {
            'id': i+1,
            'atoms': molecule,
            'mass': mass,
            'com_position': com_position,
            'com_velocity': com_velocity,
            'com_v': np.linalg.norm(com_velocity),
            'com_ke': (mass * np.linalg.norm(com_velocity)**2) / 2
        }
        molecules_out.append(molecule_out)

    pbar.close()
    time.sleep(0.1)

    return molecules_out


def calculateYields(molecules, atoms, ejected_limit = 10, reflected_limit = -10):
    print('Calculating yields')

    yields = {'molecules': {'ejected': [], 'reflected': []}, 'atoms': {'ejected': [], 'reflected': []}}

    for molecule in molecules:
        if molecule['com_position'].item(2) > 0 and molecule['com_velocity'].item(2) > 0:
            for atom_id in molecule['atoms']:
                if atoms[atom_id]['position'].item(2) <= ejected_limit:
                    break
            else:
                yields['molecules']['ejected'].append(molecule)
        elif molecule['com_position'].item(2) < 0 and molecule['com_velocity'].item(2) < 0:
            for atom_id in molecule['atoms']:
                if atoms[atom_id]['position'].item(2) >= reflected_limit:
                    break
            else:
                yields['molecules']['reflected'].append(molecule)

    ejected = {}
    for mol in yields['molecules']['ejected']:
        for atom in mol['atoms']:
            if atoms[atom]['type'] not in ejected.keys():
                ejected[atoms[atom]['type']] = 1
            else:
                ejected[atoms[atom]['type']] += 1
    reflected = {}
    for mol in yields['molecules']['reflected']:
        for atom in mol['atoms']:
            if atoms[atom]['type'] not in reflected.keys():
                reflected[atoms[atom]['type']] = 1
            else:
                reflected[atoms[atom]['type']] += 1

    yields['atoms']['ejected'] = ejected
    yields['atoms']['reflected'] = reflected

    return yields


def calculateMassSpectrum(yields):
    print('Calculating mass spectrum')

    mass_spectrum = {'ejected': [], 'reflected': []}
    for mol in yields['molecules']['ejected']:
        mass_spectrum['ejected'].append(mol['mass'])
    for mol in yields['molecules']['reflected']:
        mass_spectrum['reflected'].append(mol['mass'])
    return mass_spectrum
