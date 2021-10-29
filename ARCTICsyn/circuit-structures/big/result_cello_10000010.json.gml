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
		id "A1_AmtR"
		label "NOT,N1_LmrA"
	]
	node
	[
		id "B3_BM3R1"
		label "NOT,E1_BetI"
	]
	node
	[
		id "F1_AmeR"
		label "NOR2,A1_AmtR"
	]
	node
	[
		id "S2_SrpR"
		label "NOT,S2_SrpR"
	]
	node
	[
		id "H1_HlyIIR"
		label "NOR2,Q2_QacR"
	]
	node
	[
		id "E1_BetI"
		label "NOR2,H1_HlyIIR"
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
		source "A1_AmtR"
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
		source "a"
		target "A1_AmtR"
	]
	edge
	[
		id 4
		source "b"
		target "B3_BM3R1"
	]
	edge
	[
		id 5
		source "H1_HlyIIR"
		target "E1_BetI"
	]
	edge
	[
		id 6
		source "F1_AmeR"
		target "E1_BetI"
	]
	edge
	[
		id 7
		source "a"
		target "F1_AmeR"
	]
	edge
	[
		id 8
		source "b"
		target "F1_AmeR"
	]
	edge
	[
		id 9
		source "c"
		target "S2_SrpR"
	]
	edge
	[
		id 10
		source "S2_SrpR"
		target "P3_PhlF"
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
		source "P3_PhlF"
		target "O"
	]
]
