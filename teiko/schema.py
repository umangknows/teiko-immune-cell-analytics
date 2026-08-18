SCHEMA_SQL = """
DROP VIEW IF EXISTS sample_cell_frequencies;
DROP VIEW IF EXISTS sample_metadata;
DROP TABLE IF EXISTS cell_counts;
DROP TABLE IF EXISTS cell_populations;
DROP TABLE IF EXISTS samples;
DROP TABLE IF EXISTS subjects;
DROP TABLE IF EXISTS projects;

CREATE TABLE projects (
    project_id TEXT PRIMARY KEY
);

CREATE TABLE subjects (
    subject_id TEXT PRIMARY KEY,
    project_id TEXT NOT NULL REFERENCES projects(project_id),
    condition TEXT NOT NULL,
    age INTEGER NOT NULL,
    sex TEXT NOT NULL,
    treatment TEXT NOT NULL,
    response TEXT
);

CREATE TABLE samples (
    sample_id TEXT PRIMARY KEY,
    subject_id TEXT NOT NULL REFERENCES subjects(subject_id),
    sample_type TEXT NOT NULL,
    time_from_treatment_start INTEGER NOT NULL,
    UNIQUE(subject_id, sample_type, time_from_treatment_start)
);

CREATE TABLE cell_populations (
    population_id INTEGER PRIMARY KEY,
    population_name TEXT NOT NULL UNIQUE
);

CREATE TABLE cell_counts (
    sample_id TEXT NOT NULL REFERENCES samples(sample_id),
    population_id INTEGER NOT NULL REFERENCES cell_populations(population_id),
    count INTEGER NOT NULL CHECK(count >= 0),
    PRIMARY KEY (sample_id, population_id)
);

CREATE INDEX idx_subjects_project ON subjects(project_id);
CREATE INDEX idx_subjects_condition_treatment_response
    ON subjects(condition, treatment, response);
CREATE INDEX idx_samples_subject_type_time
    ON samples(subject_id, sample_type, time_from_treatment_start);
CREATE INDEX idx_counts_population ON cell_counts(population_id);

CREATE VIEW sample_metadata AS
SELECT
    p.project_id AS project,
    s.subject_id AS subject,
    s.condition,
    s.age,
    s.sex,
    s.treatment,
    s.response,
    sm.sample_id AS sample,
    sm.sample_type,
    sm.time_from_treatment_start
FROM samples sm
JOIN subjects s ON s.subject_id = sm.subject_id
JOIN projects p ON p.project_id = s.project_id;

CREATE VIEW sample_cell_frequencies AS
WITH totals AS (
    SELECT sample_id, SUM(count) AS total_count
    FROM cell_counts
    GROUP BY sample_id
)
SELECT
    sm.project,
    sm.subject,
    sm.condition,
    sm.age,
    sm.sex,
    sm.treatment,
    sm.response,
    sm.sample,
    sm.sample_type,
    sm.time_from_treatment_start,
    t.total_count,
    cp.population_name AS population,
    cc.count,
    CASE
        WHEN t.total_count = 0 THEN NULL
        ELSE 100.0 * cc.count / t.total_count
    END AS percentage
FROM cell_counts cc
JOIN totals t ON t.sample_id = cc.sample_id
JOIN cell_populations cp ON cp.population_id = cc.population_id
JOIN sample_metadata sm ON sm.sample = cc.sample_id;
"""

