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
		id "NOT_0"
		label "NOT,S2_SrpR"
	]
	node
	[
		id "NOR2_1"
		label "NOR2,P1_PhlF"
	]
	node
	[
		id "NOT_2"
		label "NOT,B3_BM3R1"
	]
	node
	[
		id "NOR2_3"
		label "NOR2,E1_BetI"
	]
	node
	[
		id "NOR2_4"
		label "NOR2,A1_AmtR"
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
	node
	[
		id "c"
		label "C"
	]
	edge
	[
		id 1
		source "a"
		target "NOR2_4"
	]
	edge
	[
		id 2
		source "b"
		target "NOR2_4"
	]
	edge
	[
		id 3
		source "NOR2_4"
		target "NOR2_3"
	]
	edge
	[
		id 4
		source "c"
		target "NOR2_3"
	]
	edge
	[
		id 5
		source "NOR2_3"
		target "NOT_2"
	]
	edge
	[
		id 6
		source "NOT_2"
		target "NOR2_1"
	]
	edge
	[
		id 7
		source "b"
		target "NOR2_1"
	]
	edge
	[
		id 8
		source "NOR2_1"
		target "NOT_0"
	]
	edge
	[
		id 9
		source "NOT_0"
		target "O"
	]
]
