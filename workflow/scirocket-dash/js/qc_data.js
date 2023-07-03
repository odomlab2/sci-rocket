var data = {
    "version": "0.1",
    "sequencing_name": "test_run",
    "n_pairs": 0,
    "n_pairs_success": 0,
    "n_pairs_failure": 0,
    "n_corrected_p5": 0,
    "n_corrected_p7": 0,
    "n_corrected_ligation": 0,
    "n_corrected_rt": 0,
    "p5_index_counts": [
        {
            "row": "A",
            "col": "2",
            "frequency": 0
        },
        {
            "row": "B",
            "col": "2",
            "frequency": 0
        }
    ],
    "ligation_barcode_counts": [
        {
            "barcode": "LIG1",
            "frequency": 1
        },
        {
            "barcode": "LIG2",
            "frequency": 2
        },
        {
            "barcode": "LIG3",
            "frequency": 3
        }
    ],
    "rt_barcode_counts": {
        "P01": [
            {
                "row": "C",
                "col": "4",
                "frequency": 1
            },
            {
                "row": "F",
                "col": "1",
                "frequency": 2
            }
        ],
    },
    "uncorrectables_sankey": [
        {
            "source": "(True, True, True, True)",
            "value": 0
        },
        {
            "source": "(True, True, True, False)",
            "value": 1
        },
        {
            "source": "(True, True, False, True)",
            "value": 2
        }
    ],
    "n_uncorrectable_p5": 1,
    "n_uncorrectable_p7": 2,
    "n_uncorrectable_ligation": 3,
    "n_uncorrectable_rt": 4,
    "n_cells": 0,
    "n_cells_umi_100": 0,
    "n_cells_umi_1000": 0,
    "top_uncorrectables": {
        "p5": [
            {
                "barcode": "GGGGGGGGGG",
                "frequency": 1
            }
        ],
        "p7": [
            {
                "barcode": "GGGGGGGGGG",
                "frequency": 1
            }
        ],
        "ligation": [
            {
                "barcode": "AGATCCTGTC",
                "frequency": 1
            }
        ],
        "rt": [
            {
                "barcode": "ACATCTCCGA",
                "frequency": 1
            }
        ]
    },
    "samples_qc": {
        "sample1": {
            "n_pairs_success": 1,
            "n_cells": 2,
            "n_UMIs": 3,
            "n_cells_umi_100": 4,
            "n_cells_umi_1000": 5,
            "dup_rate": 0.9
        },
        "sample2": {
            "n_pairs_success": 1,
            "n_cells": 2,
            "n_UMIs": 3,
            "n_cells_umi_100": 4,
            "n_cells_umi_1000": 5,
            "dup_rate": 0.5
        }
    }
}