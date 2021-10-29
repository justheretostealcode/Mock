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
		label "NOR2,S3_SrpR"
	]
	node
	[
		id "NOT_2"
		label "NOT,B3_BM3R1"
	]
	node
	[
		id "c"
		label "C"
	]
	node
	[
		id "b"
		label "B"
	]
	node
	[
		id "NOR2_3"
		label "NOR2,E1_BetI"
	]
	node
	[
		id "NOR2_4"
		label "NOR2,F1_AmeR"
	]
	node
	[
		id "a"
		label "A"
	]
	node
	[
		id "NOR2_5"
		label "NOR2,H1_HlyIIR"
	]
	node
	[
		id "NOT_6"
		label "NOT,A1_AmtR"
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
		source "b"
		target "NOR2_1"
	]
	edge
	[
		id 4
		source "NOR2_1"
		target "NOR2_0"
	]
	edge
	[
		id 5
		source "a"
		target "NOR2_4"
	]
	edge
	[
		id 6
		source "b"
		target "NOR2_4"
	]
	edge
	[
		id 7
		source "NOR2_4"
		target "NOR2_3"
	]
	edge
	[
		id 8
		source "a"
		target "NOT_6"
	]
	edge
	[
		id 9
		source "NOT_6"
		target "NOR2_5"
	]
	edge
	[
		id 10
		source "NOT_2"
		target "NOR2_5"
	]
	edge
	[
		id 11
		source "NOR2_5"
		target "NOR2_3"
	]
	edge
	[
		id 12
		source "NOR2_3"
		target "NOR2_0"
	]
	edge
	[
		id 13
		source "NOR2_0"
		target "O"
	]
]
