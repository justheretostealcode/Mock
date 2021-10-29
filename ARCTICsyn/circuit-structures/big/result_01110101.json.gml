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
		label "NOR2,P3_PhlF"
	]
	node
	[
		id "NOT_2"
		label "NOT,A1_AmtR"
	]
	node
	[
		id "c"
		label "C"
	]
	node
	[
		id "NOR2_3"
		label "NOR2,E1_BetI"
	]
	node
	[
		id "NOT_4"
		label "NOT,H1_HlyIIR"
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
		target "NOT_2"
	]
	edge
	[
		id 2
		source "NOT_2"
		target "NOR2_1"
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
		target "NOR2_3"
	]
	edge
	[
		id 5
		source "b"
		target "NOR2_3"
	]
	edge
	[
		id 6
		source "NOR2_3"
		target "NOR2_1"
	]
	edge
	[
		id 7
		source "NOR2_1"
		target "NOT_0"
	]
	edge
	[
		id 8
		source "NOT_0"
		target "O"
	]
]
