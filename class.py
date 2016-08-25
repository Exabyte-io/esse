import os
import json
# import python_jsonschema_objects as pjs
import warlock

THISFILE_DIR = os.path.dirname(os.path.realpath(__file__))

def read_file(file):
    filename = os.path.join(THISFILE_DIR, file)
    with open(filename) as f:
        text = f.read()
    return text

Unit = warlock.model_factory(json.loads(read_file('compiled/schema/job/model/method/workflow/unit.json')))

assignment1 = Unit(**{
    "type": "assignment",
    "name": "initialize_previous_energy",
    "flowchartId": "initialize_previous_energy",
    "next": "initialize_PARAMETER",
    "head": True,
    "assignment": {
        "name": "PARAMETER",
        "value": "0"
    }
})

print json.dumps(assignment1)

assignment2 = Unit(**{
    "type": "assignment",
    "name": "initialize_previous_energy",
    "head": False,
    "flowchartId": "initialize_previous_energy",
    "next": "initialize_tolerance",
    "assignment": {
        "name": "TOLERANCE",
        "value": "0"
    }
})
#
Workflow = warlock.model_factory(json.loads(read_file('compiled/schema/job/model/method/workflow.json')))

sample_workflow = Workflow(**{
    "name": "sample_workflow",
    "units": []
})

print json.dumps(sample_workflow)

