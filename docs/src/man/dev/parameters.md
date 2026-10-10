# Parameters

The input deck is read into typed structs. Every section and every model declares
its YAML keys as an `@params` struct (`src/Support/Parameters/Spec`): type,
required or default, `min` / `max`, allowed values, a quantity (documentation only;
PeriLab has no fixed units) and a description.

`PeriLab.InputDeck.read_input` validates a deck against these declarations and
reports every problem at once (unknown keys with a suggestion, wrong types, values
out of range). With `Strict Validation: false` (or `--no_strict`) unknown keys are
warnings instead of errors.

A model module declares its own keys and registers them under its model name, e.g.

    @params struct MyMaterialParams
        "Stress at which the material starts to yield"
        yield_stress::Dependent = req("Yield Stress"; min = 0, quantity = :stress)
    end
    __init__() = register_material("My Material", MyMaterialParams)

A docstring above a field is its description (the same as
`description = "..."`; giving both is an error).

See the templates under `src/Models/*/…_template` for every category.

The [input reference](../../generated/input_sections.md) lists every section and model with its
keys; it is generated from the declarations (`PeriLab.generate_parameter_docs`). In a
Julia session, `PeriLab.describe("Correspondence Plastic")` prints the keys of one
model and `PeriLab.describe("Correspondence Plastic"; template = true)` a YAML block
to start from. `PeriLab.to_json_schema(PeriLab.InputDeck.PeriLabInput)` returns a
JSON Schema of the whole deck for editors.

!!! note "Good start"
    Please check some of the full scale tests. There are several yaml files with parameter definitions.
