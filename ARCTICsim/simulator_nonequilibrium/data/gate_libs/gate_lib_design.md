# Structure of the new GateLib

Each entry in the gate lib has the following standard entries:

```
identifier
collection
name
```

Entries specifying the parameters of a model additionally include

```
model_info
```

## Promoters

Additional entries of promoters

```
cognate_transcription_factors
sequence_ids
```

Exemplary promoter entry

```
{
        "identifier": "promoter_P1_PhlF",
        "name": "P1_PhlF",
        "collection": "promoters",
        "cognate_transcription_factors": [PhlF, PhlF_a, PhlF_b],
        "technology_mapping": {
            "inflection_point": {
                "x_m": 401.50904740377604,
                "y_m": 7299.533329199669,
                "d_y_m": -16.363749316592973
            }
        },
        "model_info": {
                    "N_PROMOTER_STATES": 6,
                    "N_PARAMETERS": 14,
                    "INPUT_SCALING_FACTOR": 0.0001,
                    "PROMOTER_ACTIVITY": [0, 0, 0, 1, 1, 1],
                    "INFINITESIMAL_GENERATOR_FUNCTION": {"matrices": {"1": [[]],
                    "c": [[]]
                    }}
        },
        "sequence_ids": []
}       
```

## Coding Sequence

```
sequence_ids
```

Exemplary coding sequence entry entry

```
{
        "identifier": "cds_PhlF",
        "name": "PhlF",
        "collection": "coding_sequence",
        "sequence_ids": ["sequence_PhlF_1", "sequence_PhlF_2"]
}
```

## Protein

This class is especially relevant as it includes the transcription factors for te promoters.

```
{
    "identifier": "protein_PhlF",
    "name": "PhlF",
    "collection": "proteins",
    "sequence_ids": ["sequence_PhlF_1", "sequence_PhlF_2"],
    "model_info": {"rna": {
                           "transcription_rate": 0.066006600660066,
                           "degradation_rate": 0.000583,
                           "e": 16,
                           "e_const": 0,
                           "length": 606.0,
                        },
                    "protein": {
                            "translation_rate": 0.0297029702970297,
                            "degradation_rate": 0.000283,
                            "e": 42,
                            "e_const": 52,                       
                            "length": 202.0,
                       }
                   }
}
```

## TranscriptionfactorInputs

This class represents the inputs of the genetic circuit as it maps the Boolean value to a transcription factor Level.

```
{
    "identifier": "tf_input_1",
    "name": "input_1",    
    "collection": "tf_inputs",    
    "sequence_ids": [],
    "model_info": {
                    "LUT": {
                            0: 100,
                            1: 1
                            }
                   }
}

```

## Sequence

```
{
    "identifier": "sequence_PhlF_a",
    "name": "PhlF_a",
    "collection": "sequences",
    "sequence": "ATGGCCCGTACGCCGAGTCGTAGCTCCATAGGATCTCTACGTTCTCCACACACACACAAGGCGATACTTACCAGCACCATTGAGATCCTGAAAGAATGCGGTTACAGCGGTCTTTCAATTGAGTCAGTCGCAAGAAGAGCAGGTGCCTCAAAACCAACAATTTACAGATGGTGGACGAATAAGGCTGCTCTTATTGCTGAGGTCTACGAAAACGAATCTGAGCAAGTAAGAAAATTCCCCGATTTGGGTTCCTTCAAAGCGGACCTTGACTTTTTGCTGAGAAACCTATGGAAGGTTTGGAGGGAAACTATATGCGGAGAGGCGTTCAGATGCGTCATAGCTGAAGCCCAGCTTGATCCAGCTACTCTTACACAGCTTAAAGATCAGTTCATGGAGAGGAGGAGAGAAATGCCCAAGAAATTGGTTGAGAACGCCATCTCCAACGGTGAACTTCCTAAGGATACCAACAGAGAGCTTCTTTTGGACATGATATTTGGATTCTGCTGGTACAGGCTACTGACGGAGCAATTAACTGTAGAACAAGACATCGAGGAGTTTACGTTTCTGCTAATCAATGGGGTTTGTCCTGGTACACAAAGAGGGAGCCCCAAAAAAAAAAGGAAGGTGTAA",
}
```

## UTR

Currently only a dummy class linking to a sequence

```
{
    "identifier": "utr_a",
    "name": "utr_a",
    "collection": "utrs",
    "sequence_id": "sequence_utr_a",
}
```

## Terminator

Currently only a dummy class linking to a sequence

```
{
    "identifier": "utr_a",
    "name": "utr_a",
    "collection": "terminators",
    "sequence_id": "sequence_utr_a",
}
```

## Device

This class represents the genetic gates available
A standard genetic gate includes the transcription factor internally signalling the information to the corresponding
promoter.
The cds_id is also used for TranscriptionfactorInputs.

Additional Fields:

Exemplary genetic gate

```
{
    "identifier": "device_a",
    "name": "device_a",
    "group": "a",
    "collection": "devices",
    "primitive_identifier": [
            "NOT",
            "NOR2"
        ],
    "utr_id": "utr_a",
    "cds_id": "protein_PhlF",
    "terminator_id": nan,
    "promoter_id": "promoter_P1_PhlF"
}
```

Exemplary Input

```
{
    "identifier": "input_1",
    "name": "input_1",
    "group": "input_1",
    "collection": "devices",
    "primitive_identifier": [
            "INPUT",
        ],
    "promoter_id": "promoter_xyl",
    "utr_id": nan,
    "cds_id": "cds_xylR",
    "terminator_id": nan
}
```

Output OR und Output Buffer

```
{
    "identifier": "output_YFP",
    "name": "output_YFP",
    "group": "YFP",
    "collection": "devices",
    "primitive_identifier": [
            "OUTPUT_OR2",
            "OUTPUT_BUFFER",
        ],
    "promoter_id": nan,
    "utr_id": nan,
    "cds_id": "protein_YFP",
    "terminator_id": nan
}
```

