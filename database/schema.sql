-- *********************************************************************************************************************
-- *********************************************************************************************************************
-- Compounds, organisms, bioactivities, biosynthetic classes, primary sequences, and motifs.
-- *********************************************************************************************************************
-- *********************************************************************************************************************

-- =====================================
-- Create table for compounds.
-- =====================================
CREATE TABLE compounds (
    npa_id VARCHAR(9) PRIMARY KEY,
    name VARCHAR(1000) NOT NULL,
    smiles VARCHAR(1000) NOT NULL,
    inchikey VARCHAR(27) NOT NULL,
    inchikey_prefix VARCHAR(14) NOT NULL
);

-- =====================================
-- Create table for organisms.
-- =====================================
CREATE TABLE organisms (
    type VARCHAR(255) NOT NULL
        CHECK (type IN ('bacterium', 'fungus')),
    genus VARCHAR(255) NOT NULL,
    species VARCHAR(255) NOT NULL,
    PRIMARY KEY (genus, species)
);

-- =====================================
-- Many-to-many: compounds_organisms.
-- =====================================
CREATE TABLE compounds_organisms (
    compound_id VARCHAR(9) NOT NULL,
    organism_genus VARCHAR(255) NOT NULL,
    organism_species VARCHAR(255) NOT NULL,
    PRIMARY KEY (compound_id, organism_genus, organism_species),
    FOREIGN KEY (compound_id) REFERENCES compounds(npa_id) ON DELETE CASCADE,
    FOREIGN KEY (organism_genus, organism_species) 
        REFERENCES organisms(genus, species) ON DELETE CASCADE
);

-- =====================================
-- Create table for bioactivities.
-- =====================================
CREATE TABLE bioactivities (
    bioactivity_type VARCHAR(255) PRIMARY KEY
);

-- Insert predefined list of bioactivities
INSERT INTO bioactivities (bioactivity_type) VALUES
('antibacterial'),
('antifungal'),
('antigramnegative'),
('antigrampositive'),
('antioxidant'),
('anticoronaviral'),
('antiretroviral'),
('antiviral'),
('cytotoxic'),
('inhibitor'),
('not_antibacterial'),
('not_antifungal'),
('not_antiviral'),
('siderophore'),
('signalling'),
('surfactant');

CREATE OR REPLACE FUNCTION prevent_modifications_bioactivities()
RETURNS TRIGGER AS $$
BEGIN
    RAISE EXCEPTION 'modifications to the bioactivities table are not allowed';
END;
$$ LANGUAGE plpgsql;

CREATE TRIGGER no_modifications
BEFORE INSERT OR UPDATE OR DELETE ON bioactivities
FOR EACH ROW
EXECUTE FUNCTION prevent_modifications_bioactivities();

-- =====================================
-- Many-to-many: compounds_bioactivities.
-- =====================================
CREATE TABLE compounds_bioactivities (
    compound_id VARCHAR(9) NOT NULL,
    bioactivity_type VARCHAR(255) NOT NULL,
    PRIMARY KEY (compound_id, bioactivity_type),
    FOREIGN KEY (compound_id) REFERENCES compounds(npa_id) ON DELETE CASCADE,
    FOREIGN KEY (bioactivity_type) REFERENCES bioactivities(bioactivity_type) ON DELETE CASCADE
);

-- =====================================
-- Create table for biosynthetic_classes.
-- =====================================
CREATE TABLE biosynthetic_classes (
    name VARCHAR(255) NOT NULL,
    level VARCHAR(255) NOT NULL,
    classifier VARCHAR(255) NOT NULL
        CHECK (classifier IN ('npclassifier', 'classyfire', 'manual')),
    PRIMARY KEY (name, level, classifier)
);

-- =====================================
-- Many-to-many: compounds_biosynthetic_classes.
-- =====================================
CREATE TABLE compounds_biosynthetic_classes (
    compound_id VARCHAR(9) NOT NULL,
    biosynthetic_class_name VARCHAR(255) NOT NULL,
    biosynthetic_class_level VARCHAR(255) NOT NULL,
    biosynthetic_class_classifier VARCHAR(255) NOT NULL,
    PRIMARY KEY (compound_id, biosynthetic_class_name, biosynthetic_class_level, biosynthetic_class_classifier),
    FOREIGN KEY (compound_id) REFERENCES compounds(npa_id) ON DELETE CASCADE,
    FOREIGN KEY (biosynthetic_class_name, biosynthetic_class_level, biosynthetic_class_classifier) 
        REFERENCES biosynthetic_classes(name, level, classifier) ON DELETE CASCADE
);

-- =====================================
-- Create table for primary_sequences.
-- =====================================
CREATE TABLE primary_sequences (
    id SERIAL PRIMARY KEY,
    retromol_version VARCHAR(50) NOT NULL
);

-- =====================================
-- Many-to-many: compounds_primary_sequences.
-- =====================================
CREATE TABLE compounds_primary_sequences (
    compound_id VARCHAR(9) NOT NULL,
    primary_sequence_id INT NOT NULL,
    PRIMARY KEY (compound_id, primary_sequence_id),
    FOREIGN KEY (compound_id) REFERENCES compounds(npa_id) ON DELETE CASCADE,
    FOREIGN KEY (primary_sequence_id) REFERENCES primary_sequences(id) ON DELETE CASCADE
);

-- =====================================
-- Create table for modifications.
-- =====================================
CREATE TABLE modifications (
    id SERIAL PRIMARY KEY,
    epoxidation BOOLEAN NOT NULL,
    glycosylation BOOLEAN NOT NULL,
    halogenation BOOLEAN NOT NULL,
    macrocyclization BOOLEAN NOT NULL,
    methylation BOOLEAN NOT NULL
);

-- =====================================
-- Many-to-many: primary_sequences_modifications.
-- =====================================
CREATE TABLE primary_sequences_modifications (
    primary_sequence_id INT NOT NULL,
    modification_id INT NOT NULL,
    PRIMARY KEY (primary_sequence_id, modification_id),
    FOREIGN KEY (primary_sequence_id) REFERENCES primary_sequences(id) ON DELETE CASCADE,
    FOREIGN KEY (modification_id) REFERENCES modifications(id) ON DELETE CASCADE
);

-- =====================================
-- Create table for motifs.
-- =====================================
CREATE TABLE motifs (
    motif_name VARCHAR(255) PRIMARY KEY
);

-- Insert predefined list of motifs
INSERT INTO motifs (motif_name) VALUES
('start'),
('end'),
('other'),
('amino acid'),
('polyketide'),
('2-aminoadipic acid'),
('2-aminoisobutyric acid'),
('2,4-diaminobutyric acid'),
('3,5-dihydroxyphenylglycine'),
('4-hydroxyphenylglycine'),
('alanine'),
('arginine'),
('asparagine'),
('aspartic acid'),
('beta-alanine'),
('beta-hydroxytyrosine'),
('cysteine'),
('D-alanine'),
('glutamic acid'),
('glutamine'),
('glycine'),
('histidine'),
('isoleucine'),
('leucine'),
('lysine'),
('N5-formyl-N5-hydroxyornithine'),
('N5-hydroxyornithine'),
('ornithine'),
('phenylalanine'),
('pipecolic acid'),
('proline'),
('serine'),
('threonine'),
('tryptophan'),
('tyrosine'),
('valine'),
('2,3-dihydroxybenzoic acid'),
('anthranillic acid'),
('salicylic acid'),
('A'),
('A1'),
('A2'),
('A3'),
('A4'),
('A5'),
('A6'),
('A7'),
('A8'),
('A9'),
('A10'),
('A11'),
('B'),
('B1'),
('B2'),
('B3'),
('B4'),
('B5'),
('B6'),
('B7'),
('B8'),
('B9'),
('B10'),
('B11'),
('C'),
('C1'),
('C2'),
('C4'),
('D'),
('D1'),
('D2'),
('D3'),
('D4'),
('D5'),
('D6'),
('D7'),
('D8'),
('D9'),
('D10'),
('D11'),
('D12');

CREATE OR REPLACE FUNCTION prevent_modifications_motifs()
RETURNS TRIGGER AS $$
BEGIN
    RAISE EXCEPTION 'modifications to the motifs table are not allowed';
END;
$$ LANGUAGE plpgsql;

CREATE TRIGGER no_modifications
BEFORE INSERT OR UPDATE OR DELETE ON motifs
FOR EACH ROW
EXECUTE FUNCTION prevent_modifications_motifs();

-- =====================================
-- Primary_sequences_motifs.
-- =====================================
CREATE TABLE primary_sequences_motifs (
    primary_sequence_id INT NOT NULL,
    position INT NOT NULL,
    motif_id VARCHAR(255) NOT NULL,
    PRIMARY KEY (primary_sequence_id, position),
    FOREIGN KEY (primary_sequence_id) REFERENCES primary_sequences(id) ON DELETE CASCADE,
    FOREIGN KEY (motif_id) REFERENCES motifs(motif_name),
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
    organism_genus VARCHAR(255) NOT NULL,
    organism_species VARCHAR(255) NOT NULL,
    PRIMARY KEY (antismash_record_id, organism_genus, organism_species),
    FOREIGN KEY (antismash_record_id) REFERENCES antismash_database_records(id) ON DELETE CASCADE,
    FOREIGN KEY (organism_genus, organism_species) 
        REFERENCES organisms(genus, species) ON DELETE CASCADE
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

-- *********************************************************************************************************************
-- *********************************************************************************************************************
-- Views.
-- *********************************************************************************************************************
-- *********************************************************************************************************************

-- =====================================
-- Create a view to display compound details.
-- =====================================
CREATE OR REPLACE VIEW compound_details_view AS
SELECT
    c.npa_id AS compound_id,
    c.name AS compound_name,
    o.type AS organism_type,
    o.genus AS organism_genus,
    o.species AS organism_species,
    bc.name AS biosynthetic_class_name,
    bc.level AS biosynthetic_class_level,
    bc.classifier AS biosynthetic_class_classifier
FROM
    compounds AS c
LEFT JOIN
    compounds_organisms AS co ON c.npa_id = co.compound_id
LEFT JOIN
    organisms AS o ON co.organism_genus = o.genus AND co.organism_species = o.species
LEFT JOIN
    compounds_biosynthetic_classes AS cbc ON c.npa_id = cbc.compound_id
LEFT JOIN
    biosynthetic_classes AS bc ON cbc.biosynthetic_class_name = bc.name 
                                AND cbc.biosynthetic_class_level = bc.level 
                                AND cbc.biosynthetic_class_classifier = bc.classifier;

-- =====================================
-- Create a view to display compound bioactivities.
-- =====================================
CREATE OR REPLACE VIEW compound_bioactivities_view AS
SELECT
    c.npa_id AS compound_id,
    c.name AS compound_name,
    b.bioactivity_type AS bioactivity
FROM
    compounds AS c
LEFT JOIN
    compounds_bioactivities AS cb ON c.npa_id = cb.compound_id
LEFT JOIN
    bioactivities AS b ON cb.bioactivity_type = b.bioactivity_type;