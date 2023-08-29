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
# --------------------- Reference Data ---------------------- #
GENOME = config['data']['genome']
REFERENCE = config['data']['ref']
INDEX = config['data']['index']
GTF = config['data']['gtf']
TRIMMING = config['control']['trimming']
MAPPING = config['control']['mapping']
DEDUP = config['control']['dedup']
# --------------------- Local Database ---------------------- #
REF_GC = os.path.join(DATABASE, 'GC', f'{GENOME}.gc')
# ---------------------- Samples Info ----------------------- #
# read the sample file using pandas lib (sample names+ fastq names) and 
# create index using the sample name
DATA = pd.read_csv(config['data']['sample_file'], sep='\t', index_col=0)
SAMPLES = DATA.index.to_list()
if config['data']['sample_info']:
    SAMPLE_MAP = json.load(open(config['data']['sample_info']))
    PAIRS = list(SAMPLE_MAP.keys())
else:
    SAMPLE_MAP = None
    PAIRS = SAMPLES
# ---------------------- Workflow Envs ---------------------- #
PIPELINE = config['pipeline']

if not config.get('workdir'):
    config['workdir'] = os.getcwd()
# workdir
workdir: config['workdir']
# ------------------- Wildcard constraints ------------------ #
wildcard_constraints:
    sample = "|".join(SAMPLES),
    pair = "|".join(PAIRS)

# ---------------------- Include rules ---------------------- #
include: "utils.smk"