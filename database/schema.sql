-- *********************************************************************************************************************
-- *********************************************************************************************************************
-- Compounds, organisms, bioactivities, biosynthetic classes, primary sequences, and motifs.
-- *********************************************************************************************************************
-- *********************************************************************************************************************

-- =====================================
-- Create table for compounds.
-- =====================================
CREATE TABLE compounds (
    compound_id      VARCHAR(20) NOT NULL,
    source           VARCHAR(50) NOT NULL 
        CHECK (source IN ('npatlas', 'manual')),
    name             VARCHAR(1000) NOT NULL,
    smiles           VARCHAR(1000) NOT NULL,
    inchikey         VARCHAR(27)   NOT NULL,
    inchikey_prefix  VARCHAR(14)   NOT NULL,
    -- If source = 'npatlas', enforce: compound_id is 9 chars AND starts with 'npa'.
    CHECK (
       (source <> 'npatlas')
       OR (LENGTH(compound_id) = 9 AND compound_id LIKE 'NPA%')
    ),
    PRIMARY KEY (compound_id)
);

-- =====================================
-- Create table for organisms.
-- =====================================
CREATE TABLE organisms (
    id SERIAL PRIMARY KEY,
    type VARCHAR(255) NOT NULL
        CHECK (type IN ('bacterium', 'fungus', 'plant', 'unknown')),
    genus VARCHAR(255) NOT NULL,
    species VARCHAR(255) NOT NULL,
    UNIQUE (type, genus, species)
);

-- =====================================
-- Many-to-many: compounds_organisms.
-- =====================================
CREATE TABLE compounds_organisms (
    compound_id VARCHAR(20) NOT NULL,
    organism_id INT NOT NULL,
    PRIMARY KEY (compound_id, organism_id),
    FOREIGN KEY (compound_id) REFERENCES compounds(compound_id) ON DELETE CASCADE,
    FOREIGN KEY (organism_id) REFERENCES organisms(id) ON DELETE CASCADE
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
    compound_id VARCHAR(20) NOT NULL,
    bioactivity_type VARCHAR(255) NOT NULL,
    PRIMARY KEY (compound_id, bioactivity_type),
    FOREIGN KEY (compound_id) REFERENCES compounds(compound_id) ON DELETE CASCADE,
    FOREIGN KEY (bioactivity_type) REFERENCES bioactivities(bioactivity_type) ON DELETE CASCADE
);

-- =====================================
-- Create table for biosynthetic_classes.
-- =====================================
CREATE TABLE biosynthetic_classes (
    id SERIAL PRIMARY KEY,
    name VARCHAR(255) NOT NULL,
    level VARCHAR(255) NOT NULL,
    classifier VARCHAR(255) NOT NULL 
        CHECK (classifier IN ('npclassifier', 'classyfire', 'manual')),
    UNIQUE (name, level, classifier)
);

-- =====================================
-- Many-to-many: compounds_biosynthetic_classes.
-- =====================================
CREATE TABLE compounds_biosynthetic_classes (
    compound_id VARCHAR(20) NOT NULL,
    biosynthetic_class_id INT NOT NULL,
    PRIMARY KEY (compound_id, biosynthetic_class_id),
    FOREIGN KEY (compound_id) REFERENCES compounds(compound_id) ON DELETE CASCADE,
    FOREIGN KEY (biosynthetic_class_id) REFERENCES biosynthetic_classes(id) ON DELETE CASCADE
);

-- =====================================
-- Create table for primary_sequences.
-- =====================================
CREATE TABLE primary_sequences (
    id SERIAL PRIMARY KEY,
    retromol_version VARCHAR(50) NOT NULL,
    created_at TIMESTAMP NOT NULL DEFAULT NOW()
);

-- =====================================
-- Many-to-many: compounds_primary_sequences.
-- =====================================
CREATE TABLE compounds_primary_sequences (
    compound_id VARCHAR(20) NOT NULL,
    primary_sequence_id INT NOT NULL,
    PRIMARY KEY (compound_id, primary_sequence_id),
    FOREIGN KEY (compound_id) REFERENCES compounds(compound_id) ON DELETE CASCADE,
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
    motif_name VARCHAR(255) PRIMARY KEY,
    motif_type VARCHAR(10) NOT NULL
        CHECK (motif_type IN ('nrp', 'pk', 'other', 'start', 'end'))
);

-- Insert predefined list of motifs
INSERT INTO motifs (motif_name, motif_type) VALUES
('start', 'start'),
('end', 'end'),
('other', 'other'),
('amino acid', 'nrp'),
('polyketide', 'pk'),
('2-aminoadipic acid', 'nrp'),
('2-aminoisobutyric acid', 'nrp'),
('2,4-diaminobutyric acid', 'nrp'),
('3,5-dihydroxyphenylglycine', 'nrp'),
('4-hydroxyphenylglycine', 'nrp'),
('alanine', 'nrp'),
('arginine', 'nrp'),
('asparagine', 'nrp'),
('aspartic acid', 'nrp'),
('beta-alanine', 'nrp'),
('R-beta-hydroxytyrosine', 'nrp'),
('cysteine', 'nrp'),
('D-alanine', 'nrp'),
('glutamic acid', 'nrp'),
('glutamine', 'nrp'),
('glycine', 'nrp'),
('histidine', 'nrp'),
('isoleucine', 'nrp'),
('leucine', 'nrp'),
('lysine', 'nrp'),
('N5-formyl-N5-hydroxyornithine', 'nrp'),
('N5-hydroxyornithine', 'nrp'),
('ornithine', 'nrp'),
('phenylalanine', 'nrp'),
('pipecolic acid', 'nrp'),
('proline', 'nrp'),
('serine', 'nrp'),
('threonine', 'nrp'),
('tryptophan', 'nrp'),
('tyrosine', 'nrp'),
('valine', 'nrp'),
('2,3-dihydroxybenzoic acid', 'nrp'),
('anthranilic acid', 'nrp'),
('salicylic acid', 'nrp'),
('A', 'pk'),
('A1', 'pk'),
('A2', 'pk'),
('A3', 'pk'),
('A4', 'pk'),
('A5', 'pk'),
('A6', 'pk'),
('A7', 'pk'),
('A8', 'pk'),
('A9', 'pk'),
('A10', 'pk'),
('A11', 'pk'),
('B', 'pk'),
('B1', 'pk'),
('B2', 'pk'),
('B3', 'pk'),
('B4', 'pk'),
('B5', 'pk'),
('B6', 'pk'),
('B7', 'pk'),
('B8', 'pk'),
('B9', 'pk'),
('B10', 'pk'),
('B11', 'pk'),
('C', 'pk'),
('C1', 'pk'),
('C2', 'pk'),
('C4', 'pk'),
('D', 'pk'),
('D1', 'pk'),
('D2', 'pk'),
('D3', 'pk'),
('D4', 'pk'),
('D5', 'pk'),
('D6', 'pk'),
('D7', 'pk'),
('D8', 'pk'),
('D9', 'pk'),
('D10', 'pk'),
('D11', 'pk'),
('D12', 'pk');

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
-- Linking antiSMASH database and MIBiG protocluster records to compounds through primary_sequences.
-- *********************************************************************************************************************
-- *********************************************************************************************************************

-- =====================================
-- Create the protoclusters table.
-- =====================================
CREATE TABLE protoclusters (
    id SERIAL PRIMARY KEY,
    input_file VARCHAR(255) NOT NULL,
    mibig_accession VARCHAR(50), -- NULL if not associated to a MIBiG entry
    start_pos INT NOT NULL,
    end_pos INT NOT NULL,
    UNIQUE (input_file, start_pos, end_pos)
);

-- =====================================
-- Many-to-many: protoclusters_organisms, as MIBiG entries and antiSMASH records can be associated with multiple organisms.
-- =====================================
CREATE TABLE protoclusters_organisms (
    protocluster_id INT NOT NULL,
    organism_id INT NOT NULL,
    PRIMARY KEY (protocluster_id, organism_id),
    FOREIGN KEY (protocluster_id) REFERENCES protoclusters(id) ON DELETE CASCADE,
    FOREIGN KEY (organism_id) REFERENCES organisms(id) ON DELETE CASCADE
);

-- =====================================
-- Many-to-many: protoclusters_compounds, for when a MIBiG entry is directly linked to a NPAtlas ID.
-- =====================================
CREATE TABLE protoclusters_compounds (
    protocluster_id INT NOT NULL,
    compound_id VARCHAR(20) NOT NULL,
    PRIMARY KEY (protocluster_id, compound_id),
    FOREIGN KEY (protocluster_id) REFERENCES protoclusters(id) ON DELETE CASCADE,
    FOREIGN KEY (compound_id) REFERENCES compounds(compound_id) ON DELETE CASCADE
);

-- =====================================
-- Many-to-many: protoclusters_primary_sequences.
-- =====================================
CREATE TABLE protoclusters_primary_sequences (
    protocluster_id INT NOT NULL,
    primary_sequence_id INT NOT NULL,
    PRIMARY KEY (protocluster_id, primary_sequence_id),
    FOREIGN KEY (protocluster_id) REFERENCES protoclusters(id) ON DELETE CASCADE,
    FOREIGN KEY (primary_sequence_id) REFERENCES primary_sequences(id) ON DELETE CASCADE
);

-- *********************************************************************************************************************
-- *********************************************************************************************************************
-- Views.
-- *********************************************************************************************************************
-- *********************************************************************************************************************

-- =====================================
-- View for retrieving all primary_sequences associated with a compound_id.
-- Example usage: 
--      SELECT compound_id, compound_name, primary_sequence
--      FROM compound_primary_sequences
--      WHERE compound_id IN ('NPA009987', 'NPA015585');
-- =====================================
CREATE VIEW compound_primary_sequences
AS SELECT 
    c.compound_id, 
    c.name AS compound_name, 
    STRING_AGG(m.motif_name, '|' ORDER BY psm.position) AS primary_sequence 
FROM compounds c 
JOIN compounds_primary_sequences cps ON c.compound_id = cps.compound_id 
JOIN primary_sequences_motifs psm ON cps.primary_sequence_id = psm.primary_sequence_id 
JOIN motifs m ON psm.motif_id = m.motif_name GROUP BY c.compound_id, c.name;
