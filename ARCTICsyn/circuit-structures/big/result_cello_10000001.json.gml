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
		id "F1_AmeR"
		label "NOT,A1_AmtR"
	]
	node
	[
		id "A1_AmtR"
		label "NOT,B3_BM3R1"
	]
	node
	[
		id "B3_BM3R1"
		label "NOR2,F1_AmeR"
	]
	node
	[
		id "E1_BetI"
		label "NOR2,H1_HlyIIR"
	]
	node
	[
		id "S4_SrpR"
		label "NOR2,S3_SrpR"
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
		source "b"
		target "F1_AmeR"
	]
	edge
	[
		id 2
		source "E1_BetI"
		target "H1_HlyIIR"
	]
	edge
	[
		id 3
		source "B3_BM3R1"
		target "H1_HlyIIR"
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
		source "c"
		target "B3_BM3R1"
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
		source "A1_AmtR"
		target "E1_BetI"
	]
	edge
	[
		id 8
		source "a"
		target "A1_AmtR"
	]
	edge
	[
		id 9
		source "A1_AmtR"
		target "S4_SrpR"
	]
	edge
	[
		id 10
		source "c"
		target "S4_SrpR"
	]
	edge
	[
		id 11
		source "H1_HlyIIR"
		target "P3_PhlF"
	]
	edge
	[
		id 12
		source "S4_SrpR"
		target "P3_PhlF"
	]
	edge
	[
		id 13
		source "P3_PhlF"
		target "O"
	]
]
