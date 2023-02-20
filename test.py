from pathlib import Path
import gzip
import shutil
from tqdm import tqdm

energies_folders = ['05keV', '10keV', '20keV', '40keV']

path = Path('./data')

print(len([f for f in path.glob('**/*') if f.is_file()]))

pbar = tqdm(total=1856)

for file in path.glob('**/*'):
    pbar.update()
    if file.is_file() and file.suffix != '.gz':
        print(file)

        with open(file, 'rb') as f_in:
            with gzip.open(file.parent / (file.name+'.gz'), 'wb') as f_out:
                shutil.copyfileobj(f_in, f_out)
        file.unlink()
pbar.close()
