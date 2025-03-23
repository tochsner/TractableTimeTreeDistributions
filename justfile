format:
	julia --project=. -e 'using JuliaFormatter; format(".")'

test:
	julia --project=. -e 'using Pkg; Pkg.test()'

build:
	rm -rf bin
	julia --project=. scripts/build.jl
	zip -r bin.zip bin
