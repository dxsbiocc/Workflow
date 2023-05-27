import os

# -------------------- Global Parameters -------------------- #
# database path
if config['data']['db'] == 'local':
    DATABASE = os.path.join(PATH, 'data')
else:
    DATABASE = config['data']['db']
# paired
if config['control']['paired']:
    RUN = ['R1', 'R2']
else:
    RUN = ['R1']

# reference
REFERENCE = config['data']['ref']
INDEX = config['data']['index']
GTF = config['data']['gtf']
TRIMMING = config['control']['trimming']
MAPPING = config['control']['mapping']
DEDUP = config['control']['dedup']
# ------------------------- End ------------------------ #