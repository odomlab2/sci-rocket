from collections import Counter
import pandas as pd


def sanity_samples(log, samples, barcodes, config):
    """
    Checks the provided sample sheet for sanity and inconsistensies with the supplied config.

    Args:
        log (logging.Logger): Logger object.
        samples (pandas.DataFrame): Imported sample-sheet.
        barcodes (pandas.DataFrame): Imported barcodes.
        config (dict): Imported config.

    Returns:
        bool: True if the sample sheet is valid, False otherwise.
    """

    samples = pd.read_csv("~/jvanriet/git/snakemake-sciseq/workflow/examples/example_samplesheet.tsv", sep="\t") 

    log.info("Checking sanity: Sample sheet.")

    # Check if the sample sheet is empty.
    if samples.empty:
        log.error("Sanity check (Sample sheet) - Sample sheet is empty!")
        return False

    # Check if the sample sheet contains the required columns.
    required_columns = set(["path_bcl", "sequencing_name", "experiment_name", "well_start", "well_end", "barcode_rt", "sample_name", "species"])
    if not required_columns.issubset(samples.columns):
        log.error("Sanity check (Sample sheet) - Missing required column(s): {}".format(", ".join(required_columns.difference(samples.columns))))
        return False

    # Check if the sample sheet contains species which are not defined in the config.
    if not set(samples["species"].unique()).issubset(config["species"].keys()):
        log.error("Sanity check (Sample sheet) - Contains species which are not defined in the config: {}".format(", ".join(set(samples["species"].unique()).difference(config["species"].keys()))))
        return False

    # Check if different species are assigned to the same sample name.
    if samples.groupby("sample_name")["species"].apply(lambda x: len(x.unique()) > 1).any():
        # Check which samples have different species assigned.
        samples_with_different_species = samples.groupby("sample_name")["species"].apply(lambda x: ", ".join(x.unique())).reset_index()

        # Print the samples with different species assigned.
        for _, row in samples_with_different_species.iterrows():
            if len(row["species"].split(", ")) > 1:
                log.error("Sanity check (Sample sheet) - Different species assigned to {}: {}".format(row["sample_name"], row["species"]))

        return False

    # Retrieve the RT barcodes.
    barcodes_rt = barcodes.query("type == 'rt'")["barcode"].unique()

    # Check if the sample sheet contains RT barcodes which are not defined in the config.
    if not set(samples["barcode_rt"].unique()).issubset(barcodes_rt):
        log.error("Sanity check (Sample sheet) - Contains RT barcodes which are not defined in the config: {}".format(", ".join(set(samples["barcode_rt"].unique()).difference(barcodes_rt))))
        return False

    # region Sanity of 96-well coordinates --------------------------------------------------------------------------------
    # Check if well_start and well_end fall within 96-well plate (between A:01 and H:12)
    if not samples["well_start"].str.contains("^[A-H][0-1][0-2]$").all():
        log.error(
            "Sanity check (Sample sheet) - Well start is not a valid 96-well plate well coordinate (between A01 and H12): {}".format(
                ", ".join(samples[~samples["well_start"].str.contains("^[A-H][0-1][0-2]$")]["well_start"].unique())
            )
        )
        return False

    if not samples["well_end"].str.contains("^[A-H][0-1][0-2]$").all():
        log.error(
            "Sanity check (Sample sheet) - Well end is not a valid 96-well plate well coordinate (between A:01 and H:12): {}".format(
                ", ".join(samples[~samples["well_end"].str.contains("^[A-H][0-1][0-2]$")]["well_end"].unique())
            )
        )
        return False

    # Check if well_start and well_end are in the correct order.
    if not (samples["well_start"] <= samples["well_end"]).all():
        log.error("Sanity check (Sample sheet) - Well start is not smaller than well end: {}".format(", ".join(samples[samples["well_start"] > samples["well_end"]]["sample_name"].unique())))
        return False

    # Generate a list of all 96-wells combinations (A01 to H12).
    wells = ["{}{}".format(chr(65 + i), "{:02d}".format(j)) for i in range(8) for j in range(1, 13)]

    # Generate a list of samples + barcode_rt per well.
    wells = pd.DataFrame(wells, columns=["well"])

    # Generate a list of barcodes_rt per well, including duplicates.
    wells["barcode_rt"] = wells["well"].apply(lambda x: samples[(samples["well_start"] <= x) & (samples["well_end"] >= x)]["barcode_rt"].values.tolist())

    # Identify duplicate RT barcodes per well.
    wells["barcode_rt_duplicates"] = wells["barcode_rt"].apply(lambda x: [item for item, count in Counter(x).items() if count > 1])

    # Check if the sample sheet contains duplicate RT barcodes per well.
    if wells["barcode_rt_duplicates"].apply(lambda x: len(x) > 0).any():
        # Report the duplicate RT barcodes per well and add the sample_name(s) of the duplicate RT barcodes.
        for _, row in wells[wells["barcode_rt_duplicates"].apply(lambda x: len(x) > 0)].iterrows():
            log.error(
                "Sanity check (Sample sheet) - Contains duplicate RT barcodes in the same 96-well for different samples: {} ({} - {})".format(
                    row["well"], ",".join(row["barcode_rt_duplicates"]), ", ".join(samples[samples["barcode_rt"].isin(row["barcode_rt_duplicates"])]["sample_name"].unique())
                )
            )

        return False

    # endregion ----------------------------------------------------------------------------------------------------------------

    # Otherwise, the sample sheet is valid.
    return True


def sanity_barcodes(log, barcodes):
    """
    Checks the provided barcodes file for sanity.

    Args:
        log (logging.Logger): Logger object.
        barcodes (pandas.DataFrame): Imported barcodes.

    Returns:
        bool: True if the barcodes file is valid, False otherwise.
    """

    log.info("Checking sanity: Barcode files.")

    # Check if the barcodes file is empty.
    if barcodes.empty:
        log.error("Sanity check (Barcodes) - Barcodes file is empty!")
        return False

    # Check if the barcodes file contains the required columns.
    required_columns = set(["type", "barcode", "sequence"])
    if not set(required_columns).issubset(barcodes.columns):
        log.error("Sanity check (Barcodes) - Barcodes file is missing required column(s): {}".format(", ".join(required_columns.difference(barcodes.columns))))
        return False

    # Check if all barcode types are defined.
    required_columns = set(["rt", "p5", "p7", "ligation"])
    if not set(required_columns).issubset(barcodes["type"].unique()):
        log.error("Sanity check (Barcodes) - Barcodes file is missing barcode type(s): {}".format(", ".join(required_columns.difference(barcodes["type"].unique()))))
        return False

    # Check if the barcodes file contains duplicate barcodes per type.
    if barcodes.groupby("type")["barcode"].apply(lambda x: x.duplicated().any()).any():
        # Print the duplicate barcodes per type with barcode (type)
        duplicated_barcodes = barcodes.groupby("type")["barcode"].apply(lambda x: x[x.duplicated()]).unique()
        for barcode in duplicated_barcodes:
            log.error("Sanity check (Barcodes) - Barcode {} is duplicated in type {}".format(barcode, barcodes.query("barcode == @barcode")["type"].unique()[0]))

        return False

    # Check if barcode sequences contain only valid nucleotides (ATCG).
    if not barcodes["sequence"].str.contains("^[ATCG]+$").all():
        # Print the invalid barcodes
        invalid_barcodes = barcodes[~barcodes["sequence"].str.contains("^[ATCG]+$")]
        log.error("Sanity check (Barcodes) - Barcodes file contains invalid nucleotides (only ATCG allowed):\n{}".format(invalid_barcodes))

        return False

    # If all checks pass, the barcodes file is valid.
    return True
