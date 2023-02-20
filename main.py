from functions import *

from pprint import pprint
from datetime import datetime

from pathlib import Path

path = Path('./data')

files = []
for file in path.glob('**/*'):
    if file.is_file() and file.suffix == '.gz':

        fileparts = file.parts
        filetype = file.suffixes[-2]
        filetype = filetype.lower()
        coords = [filepart.lower() for filepart in fileparts[-1].split('_')]
        coords = [coords[-2], coords[-1].split(filetype)[0]]
        layers = fileparts[-2]
        energy = fileparts[-3]

        if fileparts[1] == 'impact_angle':
            continue

        for member in files:
            if member['id'] == [energy, layers, coords]:
                member['paths'][filetype] = file
                break
        else:
            files.append({'id': [energy, layers, coords], 'paths': {filetype: file}})

files_tmp = []
for fileset in files:
    if len(fileset['paths']) == 3:
        files_tmp.append(fileset)
files = files_tmp

fileslen = len(files)
for i, fileset in enumerate(files):
    print('---')
    print(str(datetime.now()))
    print('File ' + str(i+1) + ' out of ' + str(fileslen))
    pprint(fileset['paths'])
    lammpstrj = readLammpstrj(fileset['paths']['.lammpstrj'])
    rdx = readRdx(fileset['paths']['.rdx'])
    bonds = readBonds(fileset['paths']['.dat'])
    molecules = createMolecules(bonds)

    if len(lammpstrj) == 0:
        continue

    try:
        atoms = unionAtoms(rdx, lammpstrj[-1]['atoms'])
        atoms = calculateAtomEnergies(atoms)

        molecules = calculateCOMProperties(molecules, atoms)

        yields = calculateYields(molecules, atoms)
        mass_spectrum = calculateMassSpectrum(yields)
    except:
        continue

    fileset['data'] = {'yields': yields, 'mass_spectrum': mass_spectrum}

    with open('results.csv', 'a', newline='') as csvfile:
        writer = csv.writer(csvfile, dialect='excel-tab')
        writer.writerow([
            '#',
            'Energy:', fileset['id'][0],
            'Layers:', fileset['id'][1],
            'Coordinates:', fileset['id'][2][0], fileset['id'][2][1]
        ])
        writer.writerow([''])
        writer.writerow(['Mass spectrum'])
        writer.writerow(['Ejected:'] + fileset['data']['mass_spectrum']['ejected'])
        writer.writerow(['Reflected:'] + fileset['data']['mass_spectrum']['reflected'])
        writer.writerow([''])
        writer.writerow(['Yields'])
        writer.writerow(['Ejected'])
        writer.writerow(['Type', 'Number'])
        writer.writerows(fileset['data']['yields']['atoms']['ejected'].items())
        writer.writerow(['Reflected'])
        writer.writerow(['Type', 'Number'])
        writer.writerows(fileset['data']['yields']['atoms']['reflected'].items())
        writer.writerow([''])
        writer.writerow(['#'*10])

exit()

lammpstrj_file = 'data/10keV/8L/10keVC60_8L_OK_00deg_10ps_X_0_0.lammpstrj.gz'
rdx_file = 'data/10keV/8L/final_0_0.rdx.gz'
bonds_file = 'data/10keV/8L/bonds_0_0.dat.gz'

lammpstrj = readLammpstrj(lammpstrj_file)
rdx = readRdx(rdx_file)
bonds = readBonds(bonds_file)
molecules = createMolecules(bonds)

atoms = unionAtoms(rdx, lammpstrj[-1]['atoms'])
atoms = calculateAtomEnergies(atoms)

molecules = calculateCOMProperties(molecules, atoms)

yields = calculateYields(molecules, atoms)
mass_spectrum = calculateMassSpectrum(yields)


