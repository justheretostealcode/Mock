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
		id "F1_AmeR"
		label "NOT,B1_BM3R1"
	]
	node
	[
		id "P1_PhlF"
		label "NOT,A1_AmtR"
	]
	node
	[
		id "S3_SrpR"
		label "NOR2,S2_SrpR"
	]
	node
	[
		id "E1_BetI"
		label "NOR2,Q1_QacR"
	]
	node
	[
		id "H1_HlyIIR"
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
		source "b"
		target "F1_AmeR"
	]
	edge
	[
		id 2
		source "F1_AmeR"
		target "E1_BetI"
	]
	edge
	[
		id 3
		source "P1_PhlF"
		target "E1_BetI"
	]
	edge
	[
		id 4
		source "b"
		target "S3_SrpR"
	]
	edge
	[
		id 5
		source "a"
		target "S3_SrpR"
	]
	edge
	[
		id 6
		source "S3_SrpR"
		target "YFP"
	]
	edge
	[
		id 7
		source "H1_HlyIIR"
		target "YFP"
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
		source "E1_BetI"
		target "H1_HlyIIR"
	]
	edge
	[
		id 10
		source "a"
		target "P1_PhlF"
	]
	edge
	[
		id 11
		source "YFP"
		target "O"
	]
]
