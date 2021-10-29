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
		id "H1_HlyIIR"
		label "NOR2,H1_HlyIIR"
	]
	node
	[
		id "S4_SrpR"
		label "NOR2,A1_AmtR"
	]
	node
	[
		id "A1_AmtR"
		label "NOT,B3_BM3R1"
	]
	node
	[
		id "E1_BetI"
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
		target "YFP"
	]
	edge
	[
		id 2
		source "P3_PhlF"
		target "YFP"
	]
	edge
	[
		id 3
		source "b"
		target "H1_HlyIIR"
	]
	edge
	[
		id 4
		source "c"
		target "H1_HlyIIR"
	]
	edge
	[
		id 5
		source "S4_SrpR"
		target "A1_AmtR"
	]
	edge
	[
		id 6
		source "H1_HlyIIR"
		target "S4_SrpR"
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
		source "A1_AmtR"
		target "E1_BetI"
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
		source "A1_AmtR"
		target "P3_PhlF"
	]
	edge
	[
		id 11
		source "c"
		target "P3_PhlF"
	]
	edge
	[
		id 12
		source "YFP"
		target "O"
	]
]
