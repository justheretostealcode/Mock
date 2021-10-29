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
		id "E1_BetI"
		label "NOT,N1_LmrA"
	]
	node
	[
		id "B3_BM3R1"
		label "NOT,H1_HlyIIR"
	]
	node
	[
		id "S4_SrpR"
		label "NOR2,A1_AmtR"
	]
	node
	[
		id "A1_AmtR"
		label "NOT,S2_SrpR"
	]
	node
	[
		id "H1_HlyIIR"
		label "NOR2,E1_BetI"
	]
	node
	[
		id "P3_PhlF"
		label "NOR2,P3_PhlF"
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
		target "H1_HlyIIR"
	]
	edge
	[
		id 2
		source "B3_BM3R1"
		target "H1_HlyIIR"
	]
	edge
	[
		id 3
		source "A1_AmtR"
		target "YFP"
	]
	edge
	[
		id 4
		source "P3_PhlF"
		target "YFP"
	]
	edge
	[
		id 5
		source "b"
		target "E1_BetI"
	]
	edge
	[
		id 6
		source "b"
		target "S4_SrpR"
	]
	edge
	[
		id 7
		source "c"
		target "S4_SrpR"
	]
	edge
	[
		id 8
		source "H1_HlyIIR"
		target "P3_PhlF"
	]
	edge
	[
		id 9
		source "S4_SrpR"
		target "P3_PhlF"
	]
	edge
	[
		id 10
		source "c"
		target "B3_BM3R1"
	]
	edge
	[
		id 11
		source "a"
		target "A1_AmtR"
	]
	edge
	[
		id 12
		source "YFP"
		target "O"
	]
]
