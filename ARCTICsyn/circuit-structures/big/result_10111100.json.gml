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
		label "NOR2,H1_HlyIIR"
	]
	node
	[
		id "NOT_2"
		label "NOT,A1_AmtR"
	]
	node
	[
		id "NOR2_3"
		label "NOR2,E1_BetI"
	]
	node
	[
		id "NOT_4"
		label "NOT,F1_AmeR"
	]
	node
	[
		id "b"
		label "B"
	]
	node
	[
		id "NOT_5"
		label "NOT,R1_PsrA"
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
	node
	[
		id "NOR2_6"
		label "NOR2,S2_SrpR"
	]
	edge
	[
		id 1
		source "b"
		target "NOT_4"
	]
	edge
	[
		id 2
		source "NOT_4"
		target "NOR2_3"
	]
	edge
	[
		id 3
		source "a"
		target "NOT_5"
	]
	edge
	[
		id 4
		source "NOT_5"
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
		source "c"
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
		source "a"
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
		target "NOR2_0"
	]
	edge
	[
		id 12
		source "NOR2_0"
		target "O"
	]
]
