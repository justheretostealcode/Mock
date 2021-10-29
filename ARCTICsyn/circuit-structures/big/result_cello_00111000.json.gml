Creator "JGraphT GML Exporter"
Version 1
graph
[
	label ""
	directed 1
	node
	[
		id "a"
		label "A"
	]
	node
	[
		id "b"
		label "B"
	]
	node
	[
		id "c"
		label "C"
	]
	node
	[
		id "YFP"
		label "OR2"
	]
	node
	[
		id "A1_AmtR"
		label "NOT,E1_BetI"
	]
	node
	[
		id "B3_BM3R1"
		label "NOT,H1_HlyIIR"
	]
	node
	[
		id "H1_HlyIIR"
		label "NOR2,A1_AmtR"
	]
	node
	[
		id "E1_BetI"
		label "NOR2,B3_BM3R1"
	]
	node
	[
		id "R1_PsrA"
		label "NOR2,S2_SrpR"
	]
	node
	[
		id "P3_PhlF"
		label "NOR2,P1_PhlF"
	]
	node
	[
		id "O"
		label "X"
	]
	edge
	[
		id 1
		source "E1_BetI"
		target "R1_PsrA"
	]
	edge
	[
		id 2
		source "b"
		target "R1_PsrA"
	]
	edge
	[
		id 3
		source "c"
		target "B3_BM3R1"
	]
	edge
	[
		id 4
		source "R1_PsrA"
		target "YFP"
	]
	edge
	[
		id 5
		source "P3_PhlF"
		target "YFP"
	]
	edge
	[
		id 6
		source "A1_AmtR"
		target "H1_HlyIIR"
	]
	edge
	[
		id 7
		source "B3_BM3R1"
		target "H1_HlyIIR"
	]
	edge
	[
		id 8
		source "b"
		target "A1_AmtR"
	]
	edge
	[
		id 9
		source "H1_HlyIIR"
		target "E1_BetI"
	]
	edge
	[
		id 10
		source "a"
		target "E1_BetI"
	]
	edge
	[
		id 11
		source "E1_BetI"
		target "P3_PhlF"
	]
	edge
	[
		id 12
		source "a"
		target "P3_PhlF"
	]
	edge
	[
		id 13
		source "YFP"
		target "O"
	]
]
