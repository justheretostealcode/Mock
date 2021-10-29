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
		label "NOT,H1_HlyIIR"
	]
	node
	[
		id "S4_SrpR"
		label "NOR2,A1_AmtR"
	]
	node
	[
		id "E1_BetI"
		label "NOR2,S2_SrpR"
	]
	node
	[
		id "N1_LmrA"
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
		source "H1_HlyIIR"
		target "N1_LmrA"
	]
	edge
	[
		id 2
		source "c"
		target "N1_LmrA"
	]
	edge
	[
		id 3
		source "N1_LmrA"
		target "P3_PhlF"
	]
	edge
	[
		id 4
		source "S4_SrpR"
		target "P3_PhlF"
	]
	edge
	[
		id 5
		source "b"
		target "H1_HlyIIR"
	]
	edge
	[
		id 6
		source "P3_PhlF"
		target "YFP"
	]
	edge
	[
		id 7
		source "E1_BetI"
		target "YFP"
	]
	edge
	[
		id 8
		source "a"
		target "E1_BetI"
	]
	edge
	[
		id 9
		source "c"
		target "E1_BetI"
	]
	edge
	[
		id 10
		source "b"
		target "S4_SrpR"
	]
	edge
	[
		id 11
		source "a"
		target "S4_SrpR"
	]
	edge
	[
		id 12
		source "YFP"
		target "O"
	]
]
