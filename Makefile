format:
	julia --project=. -e 'using JuliaFormatter; format(".")'

tests:
	julia --project=. -e 'using Pkg; Pkg.test()'

build:
	julia --project=. scripts/build.jl
	zip -r bin.zip bin

PHONY: format tests build