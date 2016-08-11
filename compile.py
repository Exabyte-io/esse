#!/usr/bin/env python
"""
Creates final json files by resolving `include` statements for all files within a folder
"""
import os
import json
import shutil
import argparse
import json_include

parser = argparse.ArgumentParser(description='compile json example resolving inclusions')
group = parser.add_argument_group('Required arguments')
# group.add_argument('-u', '--username', action="store", required=True, help="job owner username")
parser.add_argument('-m', '--minify', default=False, action="store_true", help="Minify json.")
args = parser.parse_args()

THISFILE_DIR = os.path.dirname(os.path.realpath(__file__))
COMPILED_DIR = os.path.join(THISFILE_DIR, 'compiled')
EXAMPLES_DIR = os.path.join(THISFILE_DIR, 'example')


def create_dir_tree(source, destination):
    """
    Creates directory tree by copying source to destination
    """
    if os.path.isdir(destination):
        shutil.rmtree(destination)
    shutil.copytree(source, destination)


def compile_json_files(dir):
    """
    Go through all json files in a directory and resolve json inclusions
    """
    for root, dirs, files in os.walk(dir):
        for file in files:
            if file.endswith(".json"):
                print root, file
                json_data = json_include.build_json_include(root, file)
                with open(os.path.join(root, file), 'w') as f:
                    if args.minify:
                        json_data = json.dumps(json.loads(json_data), separators=(',', ':'))
                    f.write(json_data)


if __name__ == '__main__':

    create_dir_tree(EXAMPLES_DIR, COMPILED_DIR)
    compile_json_files(COMPILED_DIR)
