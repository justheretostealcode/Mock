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
		id "H1_HlyIIR"
		label "NOT,A1_AmtR"
	]
	node
	[
		id "B3_BM3R1"
		label "NOT,F1_AmeR"
	]
	node
	[
		id "A1_AmtR"
		label "NOT,S3_SrpR"
	]
	node
	[
		id "S4_SrpR"
		label "NOR2,E1_BetI"
	]
	node
	[
		id "Q1_QacR"
		label "NOT,H1_HlyIIR"
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
		source "c"
		target "A1_AmtR"
	]
	edge
	[
		id 2
		source "a"
		target "H1_HlyIIR"
	]
	edge
	[
		id 3
		source "S4_SrpR"
		target "Q1_QacR"
	]
	edge
	[
		id 4
		source "A1_AmtR"
		target "P3_PhlF"
	]
	edge
	[
		id 5
		source "Q1_QacR"
		target "P3_PhlF"
	]
	edge
	[
		id 6
		source "b"
		target "B3_BM3R1"
	]
	edge
	[
		id 7
		source "H1_HlyIIR"
		target "S4_SrpR"
	]
	edge
	[
		id 8
		source "B3_BM3R1"
		target "S4_SrpR"
	]
	edge
	[
		id 9
		source "P3_PhlF"
		target "O"
	]
]
