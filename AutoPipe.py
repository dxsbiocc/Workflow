import os
import yaml
import argparse


PWD = os.path.dirname(os.path.abspath(__file__))

def parser():
    pass


if __name__ == '__main__':
    config_path = os.path.join(PWD, 'config/config.yaml')
    with open(config_path, 'r') as f:
        config = yaml.load(f, Loader=yaml.FullLoader)