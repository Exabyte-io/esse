import json
from itertools import chain
from typing import Tuple, Union, List


class LabeledNode:
    def __init__(self, value, label):
        self.label = label
        self.value = value
        self.nodes = []

    def add_node(self, lab_node: "LabeledNode") -> "LabeledNode":
        if isinstance(lab_node, LabeledNode) and self.has_correct_label(lab_node):
            self.nodes.append(lab_node)
            return lab_node

    def has_correct_label(self, lab_node: "LabeledNode") -> bool:
        if len(self.nodes) == 0:
            return True
        return all([lab_node.label == node.label for node in self.nodes])

    @property
    def has_children(self) -> bool:
        return len(self.nodes) > 0

    def items(self) -> Tuple[str, Union[str, List[str]]]:
        for node in self.nodes:
            yield node.label, node.value

    @property
    def options(self) -> List[str]:
        return [node.value for node in self.nodes]

    @property
    def child_label(self):
        if not self.has_children:
            return None
        return self.nodes[0].label


def build_dependency(parent: str, child: str, option_map: dict) -> dict:
    cases = [build_case(parent, child, option={p: c}) for p, c in option_map.items()]
    return {
        "dependencies": {
            parent: {
                "oneOf": cases
            }
        }
    }


def build_case(parent: str, child: str, option: dict, extra_fields: dict = None) -> dict:
    if extra_fields is None:
        extra_fields = {}
    return {
        "properties": {
            parent: {
                "enum": list(option.keys()),
                **(extra_fields.get(parent, {}))
            },
            child: {
                "enum": list(chain.from_iterable(option.values())),
                **(extra_fields.get(child, {}))
            }
        }
    }


def recursive_dependency(nodelist: List[LabeledNode]) -> dict:
    if len(nodelist) == 0 or not any([node.has_children for node in nodelist]):
        return {}
    parent_label = nodelist[0].label
    return {
        "dependencies": {
            parent_label: {
                # only parent_values, child_values, and dependencies will change here
                "oneOf": [build_case_v2(
                    parent_name=node.label,
                    parent_value=node.value,
                    child_name=node.child_label,
                    child_values=node.options,
                    dependencies=recursive_dependency(node.nodes)) for node in nodelist]
            }
        }
    }


def build_case_v2(parent_name: str,
                  parent_value: str,
                  child_name: str,
                  child_values: list,
                  dependencies: dict = None,
                  extra_fields: dict = None
                  ) -> dict:
    if extra_fields is None:
        extra_fields = {}
    if dependencies is None:
        dependencies = {}
    has_dependency = len(dependencies) > 0
    return {
        "properties": {
            parent_name: {
                "enum": [parent_value],
                **(extra_fields.get(parent_name, {}))
            },
            child_name: {
                "enum": child_values,
                **(extra_fields.get(child_name, {}))
            }
        },
        **(dependencies if has_dependency else {})
    }


def non_recursive() -> None:
    o_map = {
        "1": ["ABC", "BCD"],
        "2": ["CDE", "DEF"],
    }
    deps = build_dependency("a", "b", o_map)
    print(deps)


def recursive() -> None:
    root = LabeledNode("root", None)
    dft = root.add_node(LabeledNode("DFT", "type"))
    lda = dft.add_node(LabeledNode("LDA", "subtype"))
    gga = dft.add_node(LabeledNode("GGA", "subtype"))
    svwn = lda.add_node(LabeledNode("SVWN", "functional"))
    pz = lda.add_node(LabeledNode("PZ", "functional"))
    pbe = gga.add_node(LabeledNode("PBE", "functional"))
    pw91 = gga.add_node(LabeledNode("PW91", "functional"))

    deps = recursive_dependency(root.nodes)
    print(json.dumps(deps, indent=2))


if __name__ == "__main__":
    recursive()
