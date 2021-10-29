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
		id "NOT_1"
		label "NOT,S2_SrpR"
	]
	node
	[
		id "a"
		label "A"
	]
	node
	[
		id "NOT_2"
		label "NOT,P1_PhlF"
	]
	node
	[
		id "NOR2_3"
		label "NOR2,B3_BM3R1"
	]
	node
	[
		id "NOR2_4"
		label "NOR2,Q2_QacR"
	]
	node
	[
		id "NOT_5"
		label "NOT,R1_PsrA"
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
		id "NOR2_6"
		label "NOR2,A1_AmtR"
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
		target "NOT_1"
	]
	edge
	[
		id 2
		source "NOT_1"
		target "OR2_0"
	]
	edge
	[
		id 3
		source "b"
		target "NOT_5"
	]
	edge
	[
		id 4
		source "NOT_5"
		target "NOR2_4"
	]
	edge
	[
		id 5
		source "c"
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
		target "NOR2_7"
	]
	edge
	[
		id 8
		source "c"
		target "NOR2_7"
	]
	edge
	[
		id 9
		source "NOR2_7"
		target "NOR2_6"
	]
	edge
	[
		id 10
		source "b"
		target "NOR2_6"
	]
	edge
	[
		id 11
		source "NOR2_6"
		target "NOR2_3"
	]
	edge
	[
		id 12
		source "NOR2_3"
		target "NOT_2"
	]
	edge
	[
		id 13
		source "NOT_2"
		target "OR2_0"
	]
	edge
	[
		id 14
		source "OR2_0"
		target "O"
	]
]
