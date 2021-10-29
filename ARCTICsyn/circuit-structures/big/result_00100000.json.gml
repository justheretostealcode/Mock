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
		label "NOR2,S2_SrpR"
	]
	node
	[
		id "NOT_1"
		label "NOT,P1_PhlF"
	]
	node
	[
		id "NOR2_2"
		label "NOR2,B3_BM3R1"
	]
	node
	[
		id "NOT_3"
		label "NOT,H1_HlyIIR"
	]
	node
	[
		id "c"
		label "C"
	]
	node
	[
		id "NOT_4"
		label "NOT,A1_AmtR"
	]
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
	edge
	[
		id 1
		source "c"
		target "NOT_3"
	]
	edge
	[
		id 2
		source "NOT_3"
		target "NOR2_2"
	]
	edge
	[
		id 3
		source "a"
		target "NOT_4"
	]
	edge
	[
		id 4
		source "NOT_4"
		target "NOR2_2"
	]
	edge
	[
		id 5
		source "NOR2_2"
		target "NOT_1"
	]
	edge
	[
		id 6
		source "NOT_1"
		target "NOR2_0"
	]
	edge
	[
		id 7
		source "b"
		target "NOR2_0"
	]
	edge
	[
		id 8
		source "NOR2_0"
		target "O"
	]
]
