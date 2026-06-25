using SciMLTesting, ProcessSimulator, Test
using JET

run_qa(
    ProcessSimulator;
    explicit_imports = true,
    jet_kwargs = (; target_defined_modules = true),
    ei_kwargs = (;
        # scalarize is owned by SymbolicUtils and re-exported (non-public) via
        # ModelingToolkit, which is where ProcessSimulator explicitly imports it.
        all_explicit_imports_via_owners = (; ignore = (:scalarize,)),
        all_explicit_imports_are_public = (; ignore = (:scalarize,)),
    ),
)
