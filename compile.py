#!/usr/bin/env python
"""
Creates final json files by resolving `include` statements for all files within a folder
"""
import os
import json
import shutil
import argparse
import json_include
from itertools import islice

parser = argparse.ArgumentParser(description='compile json example resolving inclusions')
group = parser.add_argument_group('Required arguments')
# group.add_argument('-u', '--username', action="store", required=True, help="job owner username")
parser.add_argument('-m', '--minify', default=False, action="store_true", help="Minify json.")
parser.add_argument('-f', '--file', default=False, action="store", dest='file_mode', help="Single json file")
parser.add_argument('-d', '--dir', default=True, action="store", dest='dir_mode', help="Top of directory tree containing json files")
parser.add_argument('-e', '--element', default="Si", action="store", dest='element', help="Element to be used in the calculation")

args = parser.parse_args()

THISFILE_DIR = os.path.dirname(os.path.realpath(__file__))
COMPILED_DIR = os.path.join(THISFILE_DIR, 'compiled')
EXAMPLES_DIR = os.path.join(THISFILE_DIR, 'example')
SCHEMA_DIR = os.path.join(THISFILE_DIR, 'schema')

def create_dir_tree(source, destination):
    """
    Creates directory tree by copying source to destination
    """
    if os.path.isdir(destination):
        shutil.rmtree(destination)
    shutil.copytree(source, destination)


def compile_json_dir(dir):
    """
    Go through all json files in a directory and resolve json inclusions
    """
    for root, dirs, files in os.walk(dir):
        for file in files:
            if file.endswith(".json"):
                if file.endswith("_template.json"):
                    with open(file, 'r') as f2:
                        newText = f2.read().replace('USER_SPECIFIED_ELEMENT', args.element)
                    with open(file, 'w') as f3:
                        f3.write(newText)
                    compile_json_file_template(os.path.join(root.replace("compiled\/examples","example"),file))
                    file = file.replace("_template","")
                print root, file
                json_data = json_include.build_json_include(root, file)
                with open(os.path.join(root, file), 'w') as f:
                    if args.minify:
                        json_data = json.dumps(json.loads(json_data), separators=(',', ':'))
                    f.write(json_data)

def compile_json_file_template(file):
    """
    Go through the json file and resolve json inclusions
    """
    if file.endswith(".json"):
        absFile = os.path.abspath(file)
        with open(absFile, 'r') as f:
            for line in f:
                if '\"ELEMENT\"' in line:
                    stub = ''.join(islice(f, 1)).strip('\n')
                    if 'value' in stub:
                        value=((((stub.strip()).strip('"')).strip('value\"\:')).strip()).strip('\"')
        with open(file, 'r') as f2:
            newText = f2.read().replace('{{ELEMENT}}', value)
        with open(absFile.replace("_template",""), "w") as f3:
            f3.write(newText)

if __name__ == '__main__':

    if (args.file_mode):
        compile_json_file_template(args.file_mode)
        exit()
    if (args.dir_mode):
        create_dir_tree(SCHEMA_DIR, os.path.join(COMPILED_DIR, 'schema'))
        create_dir_tree(EXAMPLES_DIR, os.path.join(COMPILED_DIR, 'examples'))
        compile_json_dir(COMPILED_DIR)
        exit()
