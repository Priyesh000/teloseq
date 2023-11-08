# import pandas as pd 
from pathlib import Path
from snakemake import shell
from urllib.parse import urlparse
from os.path import exists 
import re
from loguru import logger
snakemake = snakemake # noqa

# logger.add(snakemake.log[0])

# https://stackoverflow.com/questions/68626097/pythonic-way-to-identify-a-local-file-or-a-url
def is_local(url):
    url = urlparse(url)
    if url.scheme in ('file', ''):
        return Path(url.path).exists()
    return False

def get_extension(fn):
    fn = Path(fn)
    r = re.compile(r'(\.tar\.gz|\.tar\.bz2|\.tar|\.gz)$').search(str(fn.name))
    if r:
        exten = r.group(1)
        logger.debug(f'Found extension: {exten} for {fn.name }')
        return exten
    return
    # logger.exception(f'Don\'t know how to handle this files {fn.name}')
    # raise ValueError(f'Don\'t know how to handle this files {fn.name}')

def is_dir(fn):
    # ext = get_extension(fn)
    # tmp = Path(str(fn).replace(ext,''))
    return Path(fn).suffix in ['.txt', '.fasta', '.fa', '.fastq', '.fna', '.fq']
    
    

def decompress(input_file, output_file):
    decompressor = {
                '.tar':      f'tar -xvf "{input_file}" -C {output_file}',
                '.tar.gz':   f'tar -xvf "{input_file}"  -C {output_file}',
                '.tar.bz2':  f'bzip2 -dc "{input_file}" | tar -xv -C {output_file}',
                '.gz':       f'gzip -dc "{input_file}" > {output_file}',

                # 'is_dir': 'ln -s "{input_file}" {output_file}'
        }
    extension = get_extension(input_file)
    cmd  = decompressor.get(extension)
    logger.info(f'Decompression cmd:\n{cmd}')
    shell(cmd)
    

input_file = snakemake.params.inputfile
output_file = snakemake.output[0]



if is_local(input_file):
    logger.info(f'Local file: {input} ')
    shell(f'cp {input_file} {output_file}')
else:
    output = Path(output_file)
    output_dir  = output.parent
    logger.info(f'Input URL: {input_file}')
    if get_extension(input_file):
        tmp_file  = output_dir / Path(input_file).name
        logger.info(f'Downloading from: {input_file}')
        if not is_dir(output):
            output = output.parent ## Output dir is created during decompression  
            output.mkdir(parents=True, exist_ok=True)
        shell(f'wget -O {tmp_file} {input_file}')
        decompress(tmp_file, output)
    else:
        logger.info(f'Downloading from: {input_file}')
        shell(f'wget -O {output} {input_file}')

    ## TODO: If in the future there is a requirement for decompressing multiple file this can be achieved by looking at the extensions 
    
    # if Path(input_file).suffix in ['.gz']: 
        # shell(f'wget -O {output}.gz {input_file}')
    #     shell(f'gzip -d {output}.gz')
    # else:
    #     shell(f'wget -O {output} {input_file}')
        