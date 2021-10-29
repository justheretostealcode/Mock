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
		id "S4_SrpR"
		label "NOT,A1_AmtR"
	]
	node
	[
		id "H1_HlyIIR"
		label "NOT,E1_BetI"
	]
	node
	[
		id "R1_PsrA"
		label "NOR2,B3_BM3R1"
	]
	node
	[
		id "P2_PhlF"
		label "NOR2,P3_PhlF"
	]
	node
	[
		id "E1_BetI"
		label "NOR2,S2_SrpR"
	]
	node
	[
		id "O"
		label "X"
	]
	edge
	[
		id 1
		source "P2_PhlF"
		target "YFP"
	]
	edge
	[
		id 2
		source "E1_BetI"
		target "YFP"
	]
	edge
	[
		id 3
		source "S4_SrpR"
		target "P2_PhlF"
	]
	edge
	[
		id 4
		source "H1_HlyIIR"
		target "P2_PhlF"
	]
	edge
	[
		id 5
		source "a"
		target "R1_PsrA"
	]
	edge
	[
		id 6
		source "c"
		target "R1_PsrA"
	]
	edge
	[
		id 7
		source "a"
		target "S4_SrpR"
	]
	edge
	[
		id 8
		source "c"
		target "H1_HlyIIR"
	]
	edge
	[
		id 9
		source "b"
		target "E1_BetI"
	]
	edge
	[
		id 10
		source "R1_PsrA"
		target "E1_BetI"
	]
	edge
	[
		id 11
		source "YFP"
		target "O"
	]
]
