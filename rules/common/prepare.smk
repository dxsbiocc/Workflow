import os


# -------------------- Global Parameters -------------------- #
# output directory
OUTDIR = os.path.abspath(config['outdir'])
# database path (local or remote)
if config['data']['db'] == 'local':
    DATABASE = os.path.join(PATH, 'data')
else:
    DATABASE = config['data']['db']
# paired end or single end sequencing
PAIRED = config['control']['paired']
if PAIRED:
    RUN = ['R1', 'R2']
else:
    RUN = ['R1']
# --------------------- Reference Data ---------------------- #
# genome name
GENOME = config['data']['genome']
# reference genome
REFERENCE = config['data']['ref']
# genome index
INDEX = config['data'].get('index', None)
if not INDEX:
    include: "index.smk"
# annotation file
GTF = config['data']['gtf']
# --------------------- Control Flags ----------------------- #
# trimming tool
TRIMMING = config['control']['trimming']
# mapping tool
MAPPING = config['control']['mapping']
# deduplication tool
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
# pipeline name
PIPELINE = config['pipeline']
# workdir
if not config.get('workdir'):
    config['workdir'] = os.getcwd()
# workdir
workdir: config['workdir']
# set output directory if not specified
if not OUTDIR:
    OUTDIR = os.path.abspath(config['workdir'])
# ------------------- Wildcard constraints ------------------ #
wildcard_constraints:
    sample = "|".join(SAMPLES),
    pair = "|".join(PAIRS)

# ---------------------- Include rules ---------------------- #
include: "utils.smk"