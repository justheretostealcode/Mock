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
		id "B1_BM3R1"
		label "NOT,A1_AmtR"
	]
	node
	[
		id "A1_AmtR"
		label "NOT,S3_SrpR"
	]
	node
	[
		id "S1_SrpR"
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
		source "B1_BM3R1"
		target "S1_SrpR"
	]
	edge
	[
		id 2
		source "b"
		target "S1_SrpR"
	]
	edge
	[
		id 3
		source "a"
		target "B1_BM3R1"
	]
	edge
	[
		id 4
		source "S1_SrpR"
		target "Q1_QacR"
	]
	edge
	[
		id 5
		source "c"
		target "A1_AmtR"
	]
	edge
	[
		id 6
		source "A1_AmtR"
		target "P3_PhlF"
	]
	edge
	[
		id 7
		source "Q1_QacR"
		target "P3_PhlF"
	]
	edge
	[
		id 8
		source "P3_PhlF"
		target "O"
	]
]
