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
		label "NOT,H1_HlyIIR"
	]
	node
	[
		id "NOR2_2"
		label "NOR2,E1_BetI"
	]
	node
	[
		id "NOT_3"
		label "NOT,F1_AmeR"
	]
	node
	[
		id "b"
		label "B"
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
		id "NOT_5"
		label "NOT,S3_SrpR"
	]
	node
	[
		id "c"
		label "C"
	]
	edge
	[
		id 1
		source "b"
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
		source "c"
		target "NOT_5"
	]
	edge
	[
		id 8
		source "NOT_5"
		target "NOR2_0"
	]
	edge
	[
		id 9
		source "NOR2_0"
		target "O"
	]
]
