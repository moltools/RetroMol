-- *********************************************************************************************************************
-- *********************************************************************************************************************
-- Compounds, organisms, bioactivities, biosynthetic classes, primary sequences, and motifs.
-- *********************************************************************************************************************
-- *********************************************************************************************************************

-- =====================================
-- Create table for compounds.
-- =====================================
CREATE TABLE compounds (
    id SERIAL PRIMARY KEY,
    name VARCHAR(255) NOT NULL,
    smiles VARCHAR(1000)
);

-- =====================================
-- Create table for organisms.
-- =====================================
CREATE TABLE organisms (
    id SERIAL PRIMARY KEY,
    genus VARCHAR(255) NOT NULL,
    species VARCHAR(255) NOT NULL,
    ncbi_id VARCHAR(50)
);

-- =====================================
-- Many-to-many: compounds_organisms.
-- =====================================
CREATE TABLE compounds_organisms (
    compound_id INT NOT NULL,
    organism_id INT NOT NULL,
    PRIMARY KEY (compound_id, organism_id),
    FOREIGN KEY (compound_id) REFERENCES compounds(id) ON DELETE CASCADE,
    FOREIGN KEY (organism_id) REFERENCES organisms(id) ON DELETE CASCADE
);

-- =====================================
-- Create table for bioactivities.
-- =====================================
CREATE TABLE bioactivities (
    id SERIAL PRIMARY KEY,
    bioactivity_type VARCHAR(255) NOT NULL
);

-- =====================================
-- Many-to-many: compounds_bioactivities.
-- =====================================
CREATE TABLE compounds_bioactivities (
    compound_id INT NOT NULL,
    bioactivity_id INT NOT NULL,
    PRIMARY KEY (compound_id, bioactivity_id),
    FOREIGN KEY (compound_id) REFERENCES compounds(id) ON DELETE CASCADE,
    FOREIGN KEY (bioactivity_id) REFERENCES bioactivities(id) ON DELETE CASCADE
);

-- =====================================
-- Create table for biosynthetic_classes.
-- =====================================
CREATE TABLE biosynthetic_classes (
    id SERIAL PRIMARY KEY,
    name VARCHAR(255) NOT NULL
);

-- =====================================
-- Many-to-many: compounds_biosynthetic_classes.
-- =====================================
CREATE TABLE compounds_biosynthetic_classes (
    compound_id INT NOT NULL,
    biosynthetic_class_id INT NOT NULL,
    PRIMARY KEY (compound_id, biosynthetic_class_id),
    FOREIGN KEY (compound_id) REFERENCES compounds(id) ON DELETE CASCADE,
    FOREIGN KEY (biosynthetic_class_id) REFERENCES biosynthetic_classes(id) ON DELETE CASCADE
);

-- =====================================
-- Create table for primary_sequences.
-- =====================================
CREATE TABLE primary_sequences (
    id SERIAL PRIMARY KEY,
    name VARCHAR(255) NOT NULL
);

-- =====================================
-- Many-to-many: compounds_primary_sequences.
-- =====================================
CREATE TABLE compounds_primary_sequences (
    compound_id INT NOT NULL,
    primary_sequence_id INT NOT NULL,
    PRIMARY KEY (compound_id, primary_sequence_id),
    FOREIGN KEY (compound_id) REFERENCES compounds(id) ON DELETE CASCADE,
    FOREIGN KEY (primary_sequence_id) REFERENCES primary_sequences(id) ON DELETE CASCADE
);

-- =====================================
-- Create table for motifs.
-- =====================================
CREATE TABLE motifs (
    id SERIAL PRIMARY KEY,
    motif_name VARCHAR(255) NOT NULL
);

-- =====================================
-- Primary_sequences_motifs.
-- =====================================
CREATE TABLE primary_sequences_motifs (
    primary_sequence_id INT NOT NULL,
    position INT NOT NULL,
    motif_id INT NOT NULL,
    PRIMARY KEY (primary_sequence_id, position),
    FOREIGN KEY (primary_sequence_id) REFERENCES primary_sequences(id) ON DELETE CASCADE,
    FOREIGN KEY (motif_id) REFERENCES motifs(id),
    CHECK (position >= 0)
);

-- =====================================
-- Trigger function to enforce continuity of positions.
-- This function checks, on INSERT, that the position is either:
-- * 0 if no motifs exist for that primary_sequence, or
-- * max_position + 1 otherwise.
-- =====================================
CREATE OR REPLACE FUNCTION check_continuity() 
RETURNS TRIGGER AS $$
DECLARE
    max_pos INT;
BEGIN
    SELECT MAX(position) INTO max_pos FROM primary_sequences_motifs
    WHERE primary_sequence_id = NEW.primary_sequence_id;

    IF max_pos IS NULL THEN
        -- no existing motifs, position must start at 0
        IF NEW.position <> 0 THEN
            RAISE EXCEPTION 'first motif in a sequence must have position 0, attempted: %', NEW.position;
        END IF;
    ELSE
        -- sequence exists, next position should be max_pos + 1
        IF NEW.position <> (max_pos + 1) THEN
            RAISE EXCEPTION 'positions must be continuous with no gaps. Next position should be % but got %', max_pos + 1, NEW.position;
        END IF;
    END IF;
    RETURN NEW;
END;
$$ LANGUAGE plpgsql;

-- Create the trigger on primary_sequences_motifs for continuity check
CREATE TRIGGER trg_continuity
BEFORE INSERT ON primary_sequences_motifs
FOR EACH ROW EXECUTE FUNCTION check_continuity();

-- *********************************************************************************************************************
-- *********************************************************************************************************************
-- Linking antiSMASH database records to compounds through primary_sequences.
-- *********************************************************************************************************************
-- *********************************************************************************************************************

-- =====================================
-- Create the antismash_database_records table.
-- =====================================
CREATE TABLE antismash_database_records (
    id SERIAL PRIMARY KEY,
    filename VARCHAR(255) NOT NULL,
    mibig_accession VARCHAR(50)
);

-- =====================================
-- Many-to-many: antismash_database_records and organisms.
-- =====================================
CREATE TABLE antismash_database_records_organisms (
    antismash_record_id INT NOT NULL,
    organism_id INT NOT NULL,
    PRIMARY KEY (antismash_record_id, organism_id),
    FOREIGN KEY (antismash_record_id) REFERENCES antismash_database_records(id) ON DELETE CASCADE,
    FOREIGN KEY (organism_id) REFERENCES organisms(id) ON DELETE CASCADE
);

-- =====================================
-- Many-to-many: antismash_database_records and primary_sequences.
-- =====================================
CREATE TABLE antismash_database_records_primary_sequences (
    antismash_record_id INT NOT NULL,
    primary_sequence_id INT NOT NULL,
    PRIMARY KEY (antismash_record_id, primary_sequence_id),
    FOREIGN KEY (antismash_record_id) REFERENCES antismash_database_records(id) ON DELETE CASCADE,
    FOREIGN KEY (primary_sequence_id) REFERENCES primary_sequences(id) ON DELETE CASCADE
);