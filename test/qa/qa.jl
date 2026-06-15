using ProcessSimulator, Aqua, JET
using SciMLTesting

run_qa(ProcessSimulator; Aqua = Aqua, JET = JET, jet = true,
    jet_kwargs = (; target_defined_modules = true))
