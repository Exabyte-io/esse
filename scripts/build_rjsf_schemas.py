import json
from itertools import chain
from typing import Tuple, Union, List
from pathlib import Path
from anytree import Node, RenderTree


class ModelNode(Node):

    def __init__(self, name, default: bool = False, label: str = None, value: str = None, **kwargs):
        super().__init__(name, **kwargs)
        self.default = default
        self.label = label
        self.value = value

    def items(self) -> Tuple[str, Union[str, List[str]]]:
        for node in self.children:
            yield node.label, node.value

    def has_correct_label(self, lab_node: "ModelNode") -> bool:
        if not self.has_children:
            return True
        return all([lab_node.label == node.label for node in self.children])

    def add_node(self, new_node: "ModelNode") -> "ModelNode":
        if self.has_correct_label(new_node):
            self.children = [*self.children, new_node]
            return new_node

    @property
    def has_children(self) -> bool:
        return len(self.children) > 0

    @property
    def child_label(self) -> Union[str, None]:
        if not self.has_children:
            return None
        return self.children[0].label

    @property
    def options(self) -> List[str]:
        return [node.value for node in self.children]

    def set_default(self, is_default) -> "ModelNode":
        self.default = is_default
        return self

    @staticmethod
    def get_default_node(nodelist: List["ModelNode"]) -> Union["ModelNode", None]:
        filtered = list(filter(lambda n: n.default, nodelist))
        if len(filtered) == 1:
            return filtered[0]

    @staticmethod
    def get_default_option(nodelist: List["ModelNode"]) -> str:
        default_node = ModelNode.get_default_node(nodelist)
        if default_node:
            return default_node.value


def recursive_dependency(nodelist: List[ModelNode]) -> dict:
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
                    dependencies=recursive_dependency(node.children)) for node in nodelist]
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


def build_properties(tree: ModelNode) -> dict:
    if not tree.is_root:
        return {}
    level1_label = tree.child_label
    level1_value = tree.children
    level2_label = tree.children[0].child_label
    level2_value = list(chain.from_iterable([c.children for c in level1_value]))
    level3_label = tree.children[0].children[0].child_label
    level3_value = list(chain.from_iterable([c.children for c in level2_value]))
    string_type = {"type": "string"}
    return {
        "type": "object",
        "properties": {
            level1_label: {
                **string_type,
                "enum": tree.options,
                "default": ModelNode.get_default_option(level1_value)
            },
            level2_label: {
                **string_type,
                "default": ModelNode.get_default_option(level2_value)
            },
            level3_label: {
                **string_type,
                "default": ModelNode.get_default_option(level3_value)
            },
        }
    }


def build_rjsf_schema(tree: ModelNode) -> dict:
    if not tree.is_root or not tree.has_children:
        return {}
    return {
        **(build_properties(tree)),
        **(recursive_dependency(tree.children))
    }


def recursive() -> None:
    root = ModelNode("root", label="root")
    dft = ModelNode("dft", parent=root, label="type", value="density functional theory").set_default(True)
    lda = ModelNode("lda", parent=dft, label="subtype", value="local density approximation")
    gga = ModelNode("gga", parent=dft, label="subtype", value="generalized gradient approximation").set_default(True)
    other = ModelNode("other", parent=dft, label="subtype", value="other")

    ModelNode("pz", parent=lda, label="functional", value="pz")
    ModelNode("pw", parent=lda, label="functional", value="pw")
    ModelNode("vwn", parent=lda, label="functional", value="vwn")
    ModelNode("other", parent=lda, label="functional", value="other")

    ModelNode("pbe", parent=gga, label="functional", value="pbe").set_default(True)
    ModelNode("pbesol", parent=gga, label="functional", value="pbesol")
    ModelNode("pw91", parent=gga, label="functional", value="pw91")
    ModelNode("other", parent=gga, label="functional", value="other")

    ModelNode("other", parent=other, label="functional", value="other")

    print(RenderTree(root))

    schema = build_rjsf_schema(tree=root)
    print(json.dumps(schema, indent=2))


def extract_slug_from_entry(schema, key) -> Union[str, None]:
    if "enum" in schema[key]:
        return schema[key]["enum"][0]["slug"]


def resolve_schema_path(schema_dir: Path, ref: str) -> Path:
    return (schema_dir / ref).resolve()


def follow_all_of(store, schema_path: Path, max_depth=3) -> dict:
    if max_depth == 0:
        return store
    current_schema_dir = schema_path.parent.resolve()
    with open(schema_path) as f:
        schema = json.load(f)

    if "allOf" not in schema:
        return store
    n_refs = len(schema["allOf"])
    if n_refs >= 1:
        ref_paths = [resolve_schema_path(current_schema_dir, schema["allOf"][i]["$ref"]) for i in range(n_refs)]
        store.update({schema_path.name: ref_paths})
        for ref_path in ref_paths:
            follow_all_of(store=store, schema_path=ref_path, max_depth=max_depth-1)
    return store



def test():
    schema_file = Path("../schema/model_units/pb/qm/dft/ksdft/lda.json")
    print(follow_all_of({}, schema_file))

if __name__ == "__main__":
    recursive()
