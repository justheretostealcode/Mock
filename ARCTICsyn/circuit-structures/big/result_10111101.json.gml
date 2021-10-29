Creator "JGraphT GML Exporter"
Version 1
graph
[
	label ""
	directed 1
	node
	[
		id "O"
		label "X"
	]
	node
	[
		id "OR2_0"
		label "OR2"
	]
	node
	[
		id "NOR2_1"
		label "NOR2,S3_SrpR"
	]
	node
	[
		id "c"
		label "C"
	]
	node
	[
		id "NOR2_2"
		label "NOR2,B3_BM3R1"
	]
	node
	[
		id "NOT_3"
		label "NOT,E1_BetI"
	]
	node
	[
		id "b"
		label "B"
	]
	node
	[
		id "NOR2_4"
		label "NOR2,Q2_QacR"
	]
	node
	[
		id "NOR2_5"
		label "NOR2,N1_LmrA"
	]
	node
	[
		id "a"
		label "A"
	]
	node
	[
		id "NOR2_6"
		label "NOR2,P3_PhlF"
	]
	node
	[
		id "NOT_7"
		label "NOT,H1_HlyIIR"
	]
	edge
	[
		id 1
		source "c"
		target "NOR2_1"
	]
	edge
	[
		id 2
		source "b"
		target "NOT_3"
	]
	edge
	[
		id 3
		source "NOT_3"
		target "NOR2_2"
	]
	edge
	[
		id 4
		source "a"
		target "NOR2_5"
	]
	edge
	[
		id 5
		source "b"
		target "NOR2_5"
	]
	edge
	[
		id 6
		source "NOR2_5"
		target "NOR2_4"
	]
	edge
	[
		id 7
		source "a"
		target "NOR2_4"
	]
	edge
	[
		id 8
		source "NOR2_4"
		target "NOR2_2"
	]
	edge
	[
		id 9
		source "NOR2_2"
		target "NOR2_1"
	]
	edge
	[
		id 10
		source "NOR2_1"
		target "OR2_0"
	]
	edge
	[
		id 11
		source "c"
		target "NOT_7"
	]
	edge
	[
		id 12
		source "NOT_7"
		target "NOR2_6"
	]
	edge
	[
		id 13
		source "NOR2_6"
		target "OR2_0"
	]
	edge
	[
		id 14
		source "OR2_0"
		target "O"
	]
	edge
	[
		id 15
		source "NOR2_5"
		target "NOR2_6"
	]
]
