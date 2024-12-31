import os
import json
import joblib
from parasect.api import run_parasect # pip install git+https://github.com/BTheDragonMaster/parasect.git@webapp

from step3_add_retromol_results import RETROMOL_VERSION


TEMP_DIR = os.path.join(os.path.dirname(__file__), "../temp")
MODEL_DIR = os.path.join(os.path.dirname(__file__), "../models")
PARASECT_MODEL = joblib.load(os.path.join(MODEL_DIR, "model.parasect"))
PARASECT_MODEL.set_params(n_jobs=1)


def add_asdb(cur, path_asdb):
    """Add ASDB records to database."""
    # get all json files in asdb directory
    asdb_files = [f for f in os.listdir(path_asdb) if f.endswith('.json')]
    for file_path in asdb_files:
        with open(os.path.join(path_asdb, file_path), "r") as file:
            data = json.load(file)
        
        # only keep records where primary sequence has length 7 or more (includes start and end, so 5 units in total)
        min_target_length = 7
        record_length = len(data["primary_sequence"])
        if record_length < min_target_length:
            continue  # skip record

        # if 'amino acid' in data, predict substrate specificities and exchange 'amino acid' with predicted specificities
        if 'amino acid' in data["primary_sequence"] and len(data["a_domain_fasta"]) > 0:
            try:
                primary_sequence = data["primary_sequence"]
                to_predict = data["a_domain_fasta"] 
                results = run_parasect(selected_input=to_predict, selected_input_type="fasta", path_temp_dir=TEMP_DIR, model=PARASECT_MODEL)
                predicted_specificities = []
                for j, result in enumerate(results):
                    result_json = result.to_json()
                    predictions = result_json["predictions"]
                    predictions = sorted(predictions, key=lambda x: x["probability"], reverse=True)
                    top_prediction = predictions[0]
                    top_prediction_name = top_prediction["substrate_name"]
                    predicted_specificities.append(top_prediction_name)

                # noticed that if predicted is +1 in length, usually extra alanine at -2 position
                if len(predicted_specificities) == primary_sequence.count('amino acid') + 1 and predicted_specificities[-2] == 'alanine':
                    predicted_specificities.pop(-2)

                # if predicted length is not same size as number of amino acid, check if between start and end
                # only amino acid, in that case just replace the whole sequence with predicted
                if len(predicted_specificities) != primary_sequence.count('amino acid') and set(primary_sequence[1:-1]) == {'amino acid'}:
                    primary_sequence = ["start"] + predicted_specificities + ["end"]
                else:
                    assert len(predicted_specificities) == primary_sequence.count('amino acid'), f"{data['record_id']} -- predicted_specificities: {predicted_specificities}, primary_sequence: {primary_sequence}"

                new_primary_sequence = []
                for motif in primary_sequence:
                    if motif == 'amino acid': new_primary_sequence.append(predicted_specificities.pop(0))
                    else: new_primary_sequence.append(motif)
                primary_sequence = new_primary_sequence
            except Exception as e:
                print(f"Error predicting specificities for {data['record_id']}: {e}")
                continue
        else:
            primary_sequence = data["primary_sequence"]







        # TODO: parse out data from parsed file so that can be uploaded to database
        exit("TMP")








        # add record to as protocluster
        cur.execute(
            """
            INSERT INTO protoclusters (input_file, start_pos, end_pos)
            VALUES (%s, %s, %s)
            ON CONFLICT (input_file, start_pos, end_pos) DO NOTHING
            DO UPDATE SET input_file = protoclusters.input_file
            RETURNING id
            """,
            (input_file, start_pos, end_pos)
        )
        protocluster_id = cur.fetchone()[0]

        # add producing organism
        cur.execute(
            """
            INSERT INTO organisms (type, genus, species)
            VALUES (%s, %s, %s)
            ON CONFLICT (type, genus, species)
            DO UPDATE SET type = organisms.type
            RETURNING id
            """,
            (type_, genus, species)
        )
        organism_id = cur.fetchone()[0]
        
        # link protocluster to organism
        cur.execute(
            """
            INSERT INTO protoclusters_organisms (protocluster_id, organism_id)
            VALUES (%s, %s)
            ON CONFLICT (protocluster_id, organism_id) DO NOTHING
            """,
            (protocluster_id, organism_id)
        )

        # add primary sequence and link to protocluster
        cur.execute(
            """
            INSERT INTO primary_sequences (retromol_version)
            VALUES (%s)
            RETURNING id;
            """,
            (RETROMOL_VERSION,)
        )
        primary_sequence_id = cur.fetchone()[0]

        # link primary sequence to protocluster it originates from
        cur.execute(
            """
            INSERT INTO protoclusters_primary_sequences (protocluster_id, primary_sequence_id)
            VALUES (%s, %s);
            """,
            (protocluster_id, primary_sequence_id)
        )

        # add primary sequence motifs
        for position, motif in enumerate(primary_sequence):
            cur.execute(
                """
                INSERT INTO primary_sequences_motifs (primary_sequence_id, position, motif_id)
                VALUES (%s, %s, %s);
                """,
                (primary_sequence_id, position, motif)
            )
