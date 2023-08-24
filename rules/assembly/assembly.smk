# configfile
configfile: os.path.join(PATH, "config/ATAC-seq.yaml")

# set config
TRIMMING == config["control"]["trimming"]

# include rules
include: os.path.join(PATH, "rules/common/trimmed.smk")
include: os.path.join(PATH, "rules/assembly/polishing.smk")
include: os.path.join(PATH, "rules/assembly/annotation.smk")
include: os.path.join(PATH, "rules/assembly/quality_control.smk")