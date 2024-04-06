JL = julia

fmt:
	$(JL) --project=. -e 'using JuliaFormatter; format("src", BlueStyle())'

fmtex:
	$(JL) --project=. -e 'using JuliaFormatter; format("examples", BlueStyle())'