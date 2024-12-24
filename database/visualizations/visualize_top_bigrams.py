import argparse
import os
import psycopg2
import matplotlib.pyplot as plt

parser = argparse.ArgumentParser()
parser.add_argument("--out", required=True, help="output file path svg")
args = parser.parse_args()

# Arial font
plt.rcParams["font.family"] = "Arial"

def connect_to_database():
    """Connect to the database and return the connection and cursor."""
    conn = psycopg2.connect(
        dbname=os.getenv("DB_NAME", "retromol"),
        user=os.getenv("DB_USER", "davidmeijer"),
        password=os.getenv("DB_PASSWORD", "postgres"),
        host=os.getenv("DB_HOST", "localhost"),
        port=os.getenv("DB_PORT", "5432")
    )
    return conn.cursor(), conn

cur, conn = connect_to_database()

cur.execute(
    """
    SELECT CONCAT(m1.motif_name, '-', m2.motif_name) AS bigram, COUNT(*) AS occurrence 
    FROM compounds_primary_sequences cps 
    JOIN primary_sequences_motifs psm1 ON cps.primary_sequence_id = psm1.primary_sequence_id 
    JOIN primary_sequences_motifs psm2 ON cps.primary_sequence_id = psm2.primary_sequence_id AND psm2.position = psm1.position + 1 
    JOIN motifs m1 ON psm1.motif_id = m1.motif_name JOIN motifs m2 ON psm2.motif_id = m2.motif_name 
    GROUP BY bigram 
    ORDER BY occurrence DESC;
    """
)
results = cur.fetchall()
top_bigrams = results[:10][::-1]
print(top_bigrams)
cur.close()
conn.close()

# plot the top 20 bigrams in a bar chart
bigrams, occurrences = zip(*top_bigrams)
plt.bar(bigrams, occurrences, color="#56b4e9", zorder=3, edgecolor="black")
plt.grid(axis="y", linestyle="--", zorder=0)
plt.xticks(rotation=90)
plt.xticks(fontweight="bold")
plt.yticks(range(0, max(occurrences)+1, 2000))
plt.ylabel("occurrences")
plt.xlabel("bigram")
plt.tight_layout()
plt.savefig(args.out)