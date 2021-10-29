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
		id "NOR2_2"
		label "NOR2,B3_BM3R1"
	]
	node
	[
		id "a"
		label "A"
	]
	node
	[
		id "NOR2_3"
		label "NOR2,A1_AmtR"
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
		id "NOR2_4"
		label "NOR2,S3_SrpR"
	]
	edge
	[
		id 1
		source "a"
		target "NOR2_2"
	]
	edge
	[
		id 2
		source "b"
		target "NOR2_3"
	]
	edge
	[
		id 3
		source "c"
		target "NOR2_3"
	]
	edge
	[
		id 4
		source "NOR2_3"
		target "NOR2_2"
	]
	edge
	[
		id 5
		source "NOR2_2"
		target "NOR2_1"
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
		target "NOR2_0"
	]
	edge
	[
		id 8
		source "a"
		target "NOR2_4"
	]
	edge
	[
		id 9
		source "NOR2_2"
		target "NOR2_4"
	]
	edge
	[
		id 10
		source "NOR2_4"
		target "NOR2_0"
	]
	edge
	[
		id 11
		source "NOR2_0"
		target "O"
	]
]
