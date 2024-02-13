-- Create the ligands_smiles table

CREATE TABLE ligands_smiles (
    "Name" VARCHAR(255),
    "SMILES" VARCHAR(255)
);

COPY ligands_smiles FROM PROGRAM 'zcat /ligands_smiles.csv.gz' DELIMITER ',' CSV HEADER;
