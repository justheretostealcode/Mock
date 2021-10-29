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
		id "NOR2_2"
		label "NOR2,H1_HlyIIR"
	]
	node
	[
		id "a"
		label "A"
	]
	node
	[
		id "NOR2_5"
		label "NOR2,E1_BetI"
	]
	node
	[
		id "NOR2_6"
		label "NOR2,A1_AmtR"
	]
	node
	[
		id "NOR2_7"
		label "NOR2,B3_BM3R1"
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
		target "NOR2_6"
	]
	edge
	[
		id 2
		source "b"
		target "NOR2_7"
	]
	edge
	[
		id 3
		source "c"
		target "NOR2_7"
	]
	edge
	[
		id 4
		source "NOR2_7"
		target "NOR2_6"
	]
	edge
	[
		id 5
		source "NOR2_6"
		target "NOR2_5"
	]
	edge
	[
		id 6
		source "a"
		target "NOR2_5"
	]
	edge
	[
		id 7
		source "NOR2_7"
		target "NOR2_2"
	]
	edge
	[
		id 8
		source "NOR2_2"
		target "NOR2_1"
	]
	edge
	[
		id 9
		source "NOR2_5"
		target "NOR2_1"
	]
	edge
	[
		id 10
		source "NOR2_1"
		target "NOT_0"
	]
	edge
	[
		id 11
		source "NOT_0"
		target "O"
	]
	edge
	[
		id 12
		source "NOR2_6"
		target "NOR2_2"
	]
]
