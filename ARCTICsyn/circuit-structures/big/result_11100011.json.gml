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
		id "NOR2_0"
		label "NOR2,P3_PhlF"
	]
	node
	[
		id "NOR2_1"
		label "NOR2,S4_SrpR"
	]
	node
	[
		id "b"
		label "B"
	]
	node
	[
		id "NOR2_3"
		label "NOR2,H1_HlyIIR"
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
		id "NOT_6"
		label "NOT,E1_BetI"
	]
	node
	[
		id "a"
		label "A"
	]
	node
	[
		id "c"
		label "C"
	]
	edge
	[
		id 1
		source "a"
		target "NOT_6"
	]
	edge
	[
		id 2
		source "NOT_6"
		target "NOR2_5"
	]
	edge
	[
		id 3
		source "c"
		target "NOR2_5"
	]
	edge
	[
		id 4
		source "NOR2_5"
		target "NOR2_4"
	]
	edge
	[
		id 5
		source "b"
		target "NOR2_4"
	]
	edge
	[
		id 6
		source "NOR2_4"
		target "NOR2_3"
	]
	edge
	[
		id 7
		source "b"
		target "NOR2_3"
	]
	edge
	[
		id 8
		source "a"
		target "NOR2_1"
	]
	edge
	[
		id 9
		source "NOR2_1"
		target "NOR2_0"
	]
	edge
	[
		id 10
		source "NOR2_3"
		target "NOR2_0"
	]
	edge
	[
		id 11
		source "NOR2_0"
		target "O"
	]
	edge
	[
		id 12
		source "NOR2_4"
		target "NOR2_1"
	]
]
