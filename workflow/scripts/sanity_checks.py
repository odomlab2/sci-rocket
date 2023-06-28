import pandas as pd


def retrieve_p5_barcodes(log, p5, barcodes_p5):
    """
    Retrieves all p5 barcodes which are defined by the p5 start/end indexes.

    Parameters:
        log (logging.Logger): Logger object.
        p5 (list): List of all unique p5 start/end coordinates (e.g., A09:H09,B09:H09).
        barcodes_p5 (set): Set of all p5 names (e.g., A09, A10, A11, ..., H09).

    Returns:
        list or booleaan: List of p5 indexes (e.g., A09, A10, A11, ..., H09). False if the sanity check failed.
    """

    # Split the start and end indexes in case of multiple strips (, separated; e.g., A09:H09,B09:H09).
    p5 = [i.split(",") for i in p5]

    # Split the start/end indexes or p5 indexes (: separated).
    p5 = [[j.split(":") for j in i] for i in p5]

    # Retrieve the unique p5 indexes.
    p5 = list(set([tuple(sorted(j)) for i in p5 for j in i]))
    p5_indexes = []

    # Per p5 index, check if the start index is smaller than the end index.
    for i in p5:
        if i[0] > i[1]:
            log.error("Sanity check (Sample sheet) - P5 start index is larger than the end index: {} > {}".format(i[0], i[1]))
            return False
        else:
            # Generate the p5 indexes in which the first letter is variable and the second number is constant.
            letter_start = i[0][0]
            letter_end = i[1][0]
            number_start = int(i[0][1:])
            number_end = int(i[1][1:])

            if number_start != number_end:
                log.error("Sanity check (Sample sheet) - P5 start and end index are not on the same row (should be {}: {} and {}".format(number_start, i[0], i[1]))
                return False

            # Add a zero to the number if it is a single digit.
            if number_start < 10:
                number_start = "0{}".format(number_start)

            # Generate the p5 indexes.
            p5_indexes += ["{}{}".format(letter, number_start) for letter in list(map(chr, range(ord(letter_start), ord(letter_end) + 1)))]

    # Check if the p5 indexes are defined in the barcodes file.
    if not set(p5_indexes).issubset(barcodes_p5):
        log.error("Sanity check (Sample sheet) - Contains p5 indexes which are not defined in the barcodes file: {}".format(", ".join(set(p5_indexes).difference(barcodes_p5))))
        return False

    # Return the p5 indexes if all checks passed.
    return p5_indexes


def retrieve_p7_barcodes(log, p7, barcodes_p7):
    """
    Retrieves all p7 barcodes which are defined by the p7 start/end indexes.

    Parameters:
        log (logging.Logger): Logger object.
        p5 (list): List of all unique p7 start/end coordinates (e.g., A01:A12,B01:B12).
        barcodes_p5 (set): Set of all p5 names (e.g., A01, B02 etc)

    Returns:
        list or booleaan: List of p7 indexes (e.g., A09, A10, A11, ..., H09). False if the sanity check failed.
    """

    # Split the start and end indexes in case of multiple strips (, separated; e.g., A01:A12,B01:B12).
    p7 = [i.split(",") for i in p7]

    # Split the start/end indexes or p5 indexes (: separated).
    p7 = [[j.split(":") for j in i] for i in p7]

    # Retrieve the unique p5 indexes.
    p7 = list(set([tuple(sorted(j)) for i in p7 for j in i]))
    p7_indexes = []

    # Per p7 index, check if the start index is smaller than the end index.
    for i in p7:
        if i[0] > i[1]:
            log.error("Sanity check (Sample sheet) - P7 start index is larger than the end index: {} > {}".format(i[0], i[1]))
            return False
        else:
            # Generate the p7 indexes in which the first letter is constant and the second number is variable.
            letter_start = i[0][0]
            letter_end = i[1][0]
            number_start = int(i[0][1:])
            number_end = int(i[1][1:])

            if letter_start != letter_end:
                log.error("Sanity check (Sample sheet) - P7 start and end index are not on the same row (should be {}: {} and {}".format(letter_start, i[0], i[1]))
                return False

            # Generate the p7 indexes.
            p7_indexes += ["{}{}".format(letter_start, "0{}".format(number) if number < 10 else number) for number in range(number_start, number_end + 1)]

    # Check if the p5 indexes are defined in the barcodes file.
    if not set(p7_indexes).issubset(barcodes_p7):
        log.error("Sanity check (Sample sheet) - Contains p7 indexes which are not defined in the barcodes file: {}".format(", ".join(set(p7_indexes).difference(barcodes_p7))))
        return False

    # Return the p7 indexes if all checks passed.
    return p7_indexes


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
    log.info("Checking sanity: Sample sheet.")

    # Check if the sample sheet is empty.
    if samples.empty:
        log.error("Sanity check (Sample sheet) - Sample sheet is empty!")
        return False

    # Check if the sample sheet contains the required columns.
    required_columns = set(["path_bcl", "sequencing_name", "experiment_name", "p5", "p7", "barcode_rt", "sample_name", "species"])
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
        log.error("Sanity check (Sample sheet) - Contains RT barcodes which are not defined barcodes file: {}".format(", ".join(set(samples["barcode_rt"].unique()).difference(barcodes_rt))))
        return False

    # region Sanity of p5/p7 indexes ---------------------------------------------------------------

    barcodes_p5 = barcodes.query("type == 'p5'")["barcode"].unique()
    indexes_p5 = retrieve_p5_barcodes(log, samples["p5"].unique(), barcodes_p5)

    if not indexes_p5:
        return False

    barcodes_p7 = barcodes.query("type == 'p7'")["barcode"].unique()
    indexes_p7 = retrieve_p7_barcodes(log, samples["p7"].unique(), barcodes_p7)

    if not indexes_p7:
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
