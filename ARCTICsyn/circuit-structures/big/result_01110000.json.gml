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
		id "NOR2_2"
		label "NOR2,E1_BetI"
	]
	node
	[
		id "NOT_3"
		label "NOT,R1_PsrA"
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
		id "b"
		label "B"
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
		target "NOR2_0"
	]
	edge
	[
		id 3
		source "c"
		target "NOT_3"
	]
	edge
	[
		id 4
		source "NOT_3"
		target "NOR2_2"
	]
	edge
	[
		id 5
		source "b"
		target "NOT_4"
	]
	edge
	[
		id 6
		source "NOT_4"
		target "NOR2_2"
	]
	edge
	[
		id 7
		source "NOR2_2"
		target "NOR2_0"
	]
	edge
	[
		id 8
		source "NOR2_0"
		target "O"
	]
]
