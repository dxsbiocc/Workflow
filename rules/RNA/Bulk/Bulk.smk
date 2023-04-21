include: os.path.join(PATH, "rules/common/utils.smk")


# ------------------------- common rules ---------------------- #
# trimming
get_trimmed(config['control']['trimmed_tool'])
# mapping
get_mapping(config['control']['mapping_tool'])
# ------------------------- special rules --------------------- #