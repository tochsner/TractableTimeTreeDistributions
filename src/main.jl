using TractableTimeTreeDistributions
using Distributions

include("interactive_workflow.jl")

function julia_main()::Cint
    if length(ARGS) == 1
        run_interactive_workflow(ARGS[1])
    end

    return 0
end

julia_main()
