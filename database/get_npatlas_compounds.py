#!/usr/bin/env python3
import argparse 
import typing as ty
from populate_database import connect_to_database


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
    parser.add_argument("--output", type=str, required=True, help="path to output csv file")
    return parser.parse_args()


def retrieve_compounds(cur, return_inchikey=False) -> ty.List[ty.Tuple[str, str]]:
    """Retrieve NPAtlas compounds from the database"""
    if return_inchikey: return_type = "c.inchikey"
    else: return_type = "c.smiles"
    cur.execute(
        f"""
        SELECT DISTINCT c.compound_id, {return_type}
        FROM compounds c 
        JOIN compounds_biosynthetic_classes cbc 
        ON c.compound_id = cbc.compound_id 
        JOIN biosynthetic_classes bc 
        ON cbc.biosynthetic_class_id = bc.id 
        WHERE bc.level = 'class' 
        AND bc.classifier = 'npclassifier' 
        AND bc.name IN %s;
        """,
        (tuple(class_list),),
    )
    npatlas_compounds = cur.fetchall()
    return npatlas_compounds


def main() -> None:
    args = cli()
    cur, conn = connect_to_database()
    npatlas_compounds = retrieve_compounds(cur)

    # write to file
    with open(args.output, "w") as f:
        f.write("npa_id,smiles\n")
        for npa_id, smiles in npatlas_compounds:
            f.write(f"{npa_id},{smiles}\n")

    conn.close()
    cur.close()


if __name__ == "__main__":
    main()
