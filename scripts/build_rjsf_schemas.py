import json
from itertools import chain
from typing import Tuple, Union, List


class LabeledNode:
    def __init__(self, value: str, label: str) -> None:
        self.label = label
        self.value = value
        self.nodes = []
        self.default = False

    def __repr__(self):
        return f"LabeledNode(label={self.label}, value={self.value}, default={self.default})"

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
    def child_label(self) -> Union[str, None]:
        if not self.has_children:
            return None
        return self.nodes[0].label

    def set_default(self, is_default) -> "LabeledNode":
        self.default = is_default
        return self

    @property
    def default_option(self) -> str:
        default_node = LabeledNode.get_default_node(self.nodes)
        if default_node:
            return default_node.value

    @staticmethod
    def get_default_node(nodelist: List["LabeledNode"]) -> Union["LabeledNode", None]:
        filtered = list(filter(lambda n: n.default, nodelist))
        if len(filtered) == 1:
            return filtered[0]

    @staticmethod
    def get_default_option(nodelist: List["LabeledNode"]) -> str:
        default_node = LabeledNode.get_default_node(nodelist)
        if default_node:
            return default_node.value


def recursive_dependency(nodelist: List[LabeledNode]) -> dict:
    if len(nodelist) == 0 or not any([node.has_children for node in nodelist]):
        return {}
    parent_label = nodelist[0].label
    return {
        "dependencies": {
            parent_label: {
                # only parent_values, child_values, and dependencies will change here
                "oneOf": [build_case(
                    parent_name=node.label,
                    parent_value=node.value,
                    child_name=node.child_label,
                    child_values=node.options,
                    dependencies=recursive_dependency(node.nodes)) for node in nodelist]
            }
        }
    }


def build_case(parent_name: str,
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


def build_properties(tree: LabeledNode) -> dict:
    if tree.value != "root":
        return {}
    level1_label = tree.child_label
    level1 = tree.nodes
    level2_label = tree.nodes[0].child_label
    level2 = list(chain.from_iterable([c.nodes for c in level1]))
    level3_label = tree.nodes[0].nodes[0].child_label
    level3 = list(chain.from_iterable([c.nodes for c in level2]))
    string_type = {"type": "string"}
    return {
        "type": "object",
        "properties": {
            level1_label: {
                **string_type,
                "enum": tree.options,
                "default": LabeledNode.get_default_option(level1)
            },
            level2_label: {
                **string_type,
                "default": LabeledNode.get_default_option(level2)
            },
            level3_label: {
                **string_type,
                "default": LabeledNode.get_default_option(level3)
            },
        }
    }


def build_rjsf_schema(tree: LabeledNode) -> dict:
    if tree.value != "root" or not tree.has_children:
        return {}
    return {
        **(build_properties(tree)),
        **(recursive_dependency(tree.nodes))
    }


def recursive() -> None:
    root = LabeledNode("root", "root")
    dft = root.add_node(LabeledNode("density functional theory", "type").set_default(True))
    lda = dft.add_node(LabeledNode("local density approximation", "subtype"))
    gga = dft.add_node(LabeledNode("generalized gradient approximation", "subtype").set_default(True))
    other = dft.add_node(LabeledNode("other", "subtype"))

    lda.add_node(LabeledNode("pz", "functional"))
    lda.add_node(LabeledNode("pw", "functional"))
    lda.add_node(LabeledNode("vwn", "functional"))
    lda.add_node(LabeledNode("other", "functional"))

    gga.add_node(LabeledNode("pbe", "functional").set_default(True))
    gga.add_node(LabeledNode("pbesol", "functional"))
    gga.add_node(LabeledNode("pw91", "functional"))
    gga.add_node(LabeledNode("other", "functional"))

    other.add_node(LabeledNode("other", "functional"))

    schema = build_rjsf_schema(tree=root)
    print(json.dumps(schema, indent=2))


if __name__ == "__main__":
    recursive()
