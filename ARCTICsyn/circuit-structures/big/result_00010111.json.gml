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
		id "NOR2_2"
		label "NOR2,P3_PhlF"
	]
	node
	[
		id "NOR2_3"
		label "NOR2,A1_AmtR"
	]
	node
	[
		id "a"
		label "A"
	]
	node
	[
		id "NOR2_4"
		label "NOR2,Q2_QacR"
	]
	node
	[
		id "NOT_5"
		label "NOT,H1_HlyIIR"
	]
	node
	[
		id "c"
		label "C"
	]
	node
	[
		id "NOT_6"
		label "NOT,N1_LmrA"
	]
	node
	[
		id "b"
		label "B"
	]
	node
	[
		id "NOR2_7"
		label "NOR2,E1_BetI"
	]
	edge
	[
		id 1
		source "a"
		target "NOR2_3"
	]
	edge
	[
		id 2
		source "c"
		target "NOT_5"
	]
	edge
	[
		id 3
		source "NOT_5"
		target "NOR2_4"
	]
	edge
	[
		id 4
		source "b"
		target "NOT_6"
	]
	edge
	[
		id 5
		source "NOT_6"
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
		source "NOR2_3"
		target "NOR2_2"
	]
	edge
	[
		id 8
		source "b"
		target "NOR2_7"
	]
	edge
	[
		id 9
		source "c"
		target "NOR2_7"
	]
	edge
	[
		id 10
		source "NOR2_7"
		target "NOR2_2"
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
		source "NOR2_2"
		target "NOT_0"
	]
]
