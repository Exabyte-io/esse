import os

from esse import ESSE, SCHEMA_DR, EXAMPLE_DR

if __name__ == "__main__":
    es = ESSE()
    for root, dirs, files in os.walk(EXAMPLE_DR):
        for file_ in files:
            example_path = os.path.join(root, file_)
            schema_path = example_path.replace(EXAMPLE_DR, SCHEMA_DR)
            try:
                schema = es._get_json(schema_path)
                example = es._get_json(example_path)
                es.validate(example, schema)
            except:
                print("{} is invalid!".format(example_path))
                raise
