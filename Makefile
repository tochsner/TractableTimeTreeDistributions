format:
	julia --project=. -e 'using JuliaFormatter; format(".")'

build:
	julia --project=. scripts/build.jl
	zip -r bin.zip bin