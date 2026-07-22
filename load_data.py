from pathlib import Path
import sqlite3

import pandas as pd


ROOT = Path(__file__).resolve().parent
CSV_PATH = ROOT / "data" / "cell-count.csv"
DB_PATH = ROOT / "clinical_trial_patients.db"

METADATA_COLUMNS = {
    "project",
    "subject",
    "condition",
    "age",
    "sex",
    "treatment",
    "response",
    "sample",
    "sample_type",
    "time_from_treatment_start",
}

DROP_TABLES_SQL = """
DROP TABLE IF EXISTS cell_counts;
DROP TABLE IF EXISTS cell_populations;
DROP TABLE IF EXISTS samples;
DROP TABLE IF EXISTS subjects;
DROP TABLE IF EXISTS projects;
"""

CREATE_TABLES_SQL = """
CREATE TABLE projects (
    project_id TEXT PRIMARY KEY
);

CREATE TABLE subjects (
    subject_id TEXT PRIMARY KEY,
    project_id TEXT NOT NULL,
    condition TEXT NOT NULL,
    age INTEGER NOT NULL CHECK (age >= 0),
    sex TEXT NOT NULL CHECK (sex IN ('M', 'F')),
    treatment TEXT NOT NULL,
    response TEXT CHECK (response IN ('yes', 'no') OR response IS NULL),
    FOREIGN KEY (project_id) REFERENCES projects(project_id)
);

CREATE TABLE samples (
    sample_id TEXT PRIMARY KEY,
    subject_id TEXT NOT NULL,
    sample_type TEXT NOT NULL,
    time_from_treatment_start INTEGER NOT NULL,
    FOREIGN KEY (subject_id) REFERENCES subjects(subject_id)
);

CREATE TABLE cell_populations (
    population TEXT PRIMARY KEY
);

CREATE TABLE cell_counts (
    sample_id TEXT NOT NULL,
    population TEXT NOT NULL,
    count INTEGER NOT NULL CHECK (count >= 0),
    PRIMARY KEY (sample_id, population),
    FOREIGN KEY (sample_id) REFERENCES samples(sample_id),
    FOREIGN KEY (population)
        REFERENCES cell_populations(population)
);

CREATE INDEX idx_subjects_project ON subjects(project_id);
CREATE INDEX idx_subjects_analysis
    ON subjects(condition, treatment, response, sex);
CREATE INDEX idx_samples_subject ON samples(subject_id);
CREATE INDEX idx_samples_analysis
    ON samples(sample_type, time_from_treatment_start);
CREATE INDEX idx_cell_counts_population ON cell_counts(population);
"""

INSERT_PROJECT_SQL = """
INSERT INTO projects (project_id)
VALUES (?);
"""

INSERT_SUBJECT_SQL = """
INSERT INTO subjects (
    subject_id,
    project_id,
    condition,
    age,
    sex,
    treatment,
    response
)
VALUES (?, ?, ?, ?, ?, ?, ?);
"""

INSERT_SAMPLE_SQL = """
INSERT INTO samples (
    sample_id,
    subject_id,
    sample_type,
    time_from_treatment_start
)
VALUES (?, ?, ?, ?);
"""

INSERT_CELL_POPULATION_SQL = """
INSERT INTO cell_populations (population)
VALUES (?);
"""

INSERT_CELL_COUNT_SQL = """
INSERT INTO cell_counts (sample_id, population, count)
VALUES (?, ?, ?);
"""


def optional_text(value: object) -> str | None:
    """Convert a CSV value to text, preserving missing values as SQL NULL."""
    return None if pd.isna(value) else str(value)


def load_data() -> None:
    """Recreate the SQLite schema and load all rows from the source CSV."""
    data = pd.read_csv(CSV_PATH)
    population_names = [
        column for column in data.columns if column not in METADATA_COLUMNS
    ]

    if not population_names:
        raise ValueError("No cell-population columns were found in the CSV.")

    invalid_count_columns = [
        column
        for column in population_names
        if not pd.api.types.is_numeric_dtype(data[column])
    ]
    if invalid_count_columns:
        raise ValueError(
            "Cell-count columns must be numeric: "
            + ", ".join(invalid_count_columns)
        )

    projects = [(project,) for project in data["project"].drop_duplicates()]
    cell_populations = [(population,) for population in population_names]

    subject_columns = [
        "subject",
        "project",
        "condition",
        "age",
        "sex",
        "treatment",
        "response",
    ]
    subjects = [
        (
            row.subject,
            row.project,
            row.condition,
            int(row.age),
            row.sex,
            row.treatment,
            optional_text(row.response),
        )
        for row in data[subject_columns].drop_duplicates("subject").itertuples(
            index=False
        )
    ]

    samples = [
        (
            row.sample,
            row.subject,
            row.sample_type,
            int(row.time_from_treatment_start),
        )
        for row in data[
            ["sample", "subject", "sample_type", "time_from_treatment_start"]
        ].itertuples(index=False)
    ]

    cell_counts = [
        (row.sample, population, int(getattr(row, population)))
        for row in data.itertuples(index=False)
        for population in population_names
    ]

    with sqlite3.connect(DB_PATH) as connection:
        connection.execute("PRAGMA foreign_keys = ON;")
        connection.executescript(DROP_TABLES_SQL)
        connection.executescript(CREATE_TABLES_SQL)
        connection.executemany(INSERT_PROJECT_SQL, projects)
        connection.executemany(INSERT_SUBJECT_SQL, subjects)
        connection.executemany(INSERT_SAMPLE_SQL, samples)
        connection.executemany(
            INSERT_CELL_POPULATION_SQL, cell_populations
        )
        connection.executemany(INSERT_CELL_COUNT_SQL, cell_counts)

    print(
        f"Loaded {len(projects)} projects, {len(subjects)} subjects, "
        f"{len(samples)} samples, {len(cell_populations)} cell populations, "
        f"and {len(cell_counts)} cell counts into {DB_PATH}"
    )


if __name__ == "__main__":
    load_data()
