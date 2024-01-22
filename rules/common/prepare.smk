import os


# pipeline name
PIPELINE = config['pipeline']
# overwrite configurations
if PIPELINE in ['ATAC-seq', 'Cut&Tag', 'ChIP-seq']:     # ATAC-seq, Cut&Tag, ChIP-seq
    configfile: os.path.join(PATH, "config/ATAC-seq.yaml")
    PIPENAME = 'ATAC'
elif PIPELINE in ['WGS', 'WES']:  # DNA-seq, targeted sequencing
    configfile: os.path.join(PATH, "config/DNA-seq.yaml")
    PIPENAME = 'DNA'
elif PIPELINE == 'RNA-seq':       # RNA-seq
    configfile: os.path.join(PATH, "config/RNA-seq.yaml")
    PIPENAME = 'RNA'
elif PIPELINE == 'WGBS':
    configfile: os.path.join(PATH, "config/WGBS.yaml")
    PIPENAME = 'WGBS'
elif PIPELINE == 'HIC':
    configfile: os.path.join(PATH, "config/HIC.yaml")
    PIPENAME = 'HIC'
elif PIPELINE == 'Assembly':
    configfile: os.path.join(PATH, "config/Assembly.yaml")
    PIPENAME = 'Assembly'
else:
    raise ValueError(f"Unknown pipeline: {PIPELINE}")
# -------------------- Global Parameters -------------------- #
# output directory
OUTDIR = config.get('outdir')
# set output directory if not specified
if not OUTDIR:
    OUTDIR = os.getcwd()
else:
    OUTDIR = os.path.abspath(OUTDIR)
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
# --------------------- Control Flags ----------------------- #
# trimming tool
TRIMMING = config['control']['trimming']
# mapping tool
MAPPING = config['control']['mapping']
# deduplication tool
DEDUP = config['control']['dedup']
# --------------------- Reference Data ---------------------- #
# genome name
GENOME = config['data']['genome']
# reference genome
REFERENCE = config['data']['ref']
# genome index
INDEX = config['data'].get('index', None)
# annotation file
GTF = config['data']['gtf']
# RNA sequence
RNA = config['data'].get('rna', None)
# ----------------------- Util rules ------------------------ #
include: "utils.smk"
if (not INDEX) or (not os.path.exists(INDEX)) or (os.listdir(os.path.dirname(INDEX)) == []):
    include: "index.smk"
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
# workdir
if not config.get('workdir'):
    config['workdir'] = os.getcwd()
# workdir
workdir: config['workdir']
# ------------------- Wildcard constraints ------------------ #
wildcard_constraints:
    sample = "|".join(SAMPLES),
    pair = "|".join(PAIRS)