#!/usr/bin/env python3
import argparse 
import typing as ty

from _utils import connect_to_database, setup_logger


class_list = [
    "Open-chain polyketides",
    "Macrolide lactones",
    "Polyene macrolides",
    "Erythromycins",
    "Epothilones",
    "Tylosins",
    "Bafilomycins",
    "Antimycins",
    "Lactam bearing macrolide lactones",
    "Cyclic peptides",
    "Linear peptides",
    "Lipopeptides",
    "Microcystins",
    "Cyanopeptolins",
]


def cli() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--mibig", type=str, required=True, help="path to mibig jsons folder")
    parser.add_argument("--output", type=str, required=True, help="path to output csv file")
    parser.add_argument("--log-level", required=False, default="INFO", help="log level for logging", choices=["DEBUG", "INFO", "WARNING", "ERROR", "CRITICAL"])
    return parser.parse_args()


def retrieve_compounds(cur, logger, return_inchikey=False) -> ty.List[ty.Tuple[str, str]]:
    """Retrieve NPAtlas compounds from the database"""
    if return_inchikey:
        return_type = "c.inchikey"
    else:
        return_type = "c.smiles"


    cur.execute(
        f"""
        SELECT DISTINCT c.npa_id, {return_type}
        FROM compounds AS c 
        JOIN compounds_biosynthetic_classes AS cbc 
        ON c.npa_id = cbc.compound_id 
        WHERE cbc.biosynthetic_class_level = 'class' 
        AND cbc.biosynthetic_class_classifier = 'npclassifier' 
        AND cbc.biosynthetic_class_name IN %s
        ORDER BY c.npa_id
        """,
        (tuple(class_list),),
    )
    npatlas_compounds = cur.fetchall()
    logger.info(f"Selected {len(npatlas_compounds)} NPAtlas compounds")

    return npatlas_compounds


def main() -> None:
    args = cli()
    log_level = args.log_level

    logger = setup_logger(log_level)
    logger.info("Selecting NPAtlas compounds")

    conn, cur = connect_to_database()
    npatlas_compounds = retrieve_compounds(cur, logger)

    # write to file
    with open(args.output, "w") as f:
        f.write("npa_id,smiles\n")
        for npa_id, smiles in npatlas_compounds:
            f.write(f"{npa_id},{smiles}\n")

    conn.close()
    cur.close()


if __name__ == "__main__":
    main()
