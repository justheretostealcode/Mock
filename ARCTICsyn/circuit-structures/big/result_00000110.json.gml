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
		id "a"
		label "A"
	]
	node
	[
		id "NOR2_1"
		label "NOR2,P3_PhlF"
	]
	node
	[
		id "NOR2_2"
		label "NOR2,R1_PsrA"
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
		id "NOT_5"
		label "NOT,A1_AmtR"
	]
	node
	[
		id "b"
		label "B"
	]
	node
	[
		id "NOR2_6"
		label "NOR2,Q2_QacR"
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
		target "NOR2_0"
	]
	edge
	[
		id 2
		source "c"
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
		source "b"
		target "NOT_5"
	]
	edge
	[
		id 5
		source "NOR2_7"
		target "NOR2_6"
	]
	edge
	[
		id 6
		source "NOR2_2"
		target "NOR2_1"
	]
	edge
	[
		id 7
		source "NOR2_6"
		target "NOR2_1"
	]
	edge
	[
		id 8
		source "NOR2_1"
		target "NOR2_0"
	]
	edge
	[
		id 9
		source "NOR2_0"
		target "O"
	]
	edge
	[
		id 10
		source "NOT_3"
		target "NOR2_7"
	]
	edge
	[
		id 11
		source "NOT_5"
		target "NOR2_7"
	]
	edge
	[
		id 12
		source "NOT_5"
		target "NOR2_6"
	]
	edge
	[
		id 13
		source "NOR2_7"
		target "NOR2_2"
	]
]
