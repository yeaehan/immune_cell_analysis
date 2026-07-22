# Clinical Trial Immune Cell Analysis

This project analyzes immune-cell populations from a clinical trial. It loads
the supplied CSV into a normalized SQLite database, calculates cell-type
frequencies, compares miraclib responders with non-responders, performs the
requested baseline subset analysis, and presents the results in an interactive
Streamlit dashboard.

## Dashboard

**Deployed dashboard link:** https://immunecell.streamlit.app/

The dashboard can also be run locally or in GitHub Codespaces using
`make dashboard`. In Codespaces, open the forwarded URL for port `8501` from
the **Ports** tab after starting the dashboard.

## Reproduce the analysis

Python 3.11 or newer is required. From the repository root, run:

```bash
make setup
make pipeline
```

To start the interactive dashboard, run:

```bash
make dashboard
```

The Makefile targets are:

- `make setup`: installs the project and all dependencies from
  `requirements.txt`.
- `make pipeline`: recreates the SQLite database from the input CSV and
  generates every analysis table and plot.
- `make dashboard`: starts the Streamlit dashboard on port `8501`.

The two pipeline stages can also be run directly:

```bash
python load_data.py
python analysis.py
```

`load_data.py` requires no command-line arguments. It reads
`data/cell-count.csv` and creates `clinical_trial_patients.db` in the
repository root. Running the pipeline again safely rebuilds the database and
outputs, which prevents stale or duplicated results.

## Database schema

The source CSV is transformed from a wide file into five related tables:

| Table | Purpose | Key relationships |
|---|---|---|
| `projects` | One row per clinical project | `project_id` is the primary key |
| `subjects` | Condition, age, sex, treatment, and response for each subject | References `projects` |
| `samples` | Sample type and treatment-start time for each collected sample | References `subjects` |
| `cell_populations` | Valid immune-cell population names | `population` is the primary key |
| `cell_counts` | One count per sample and cell population | References `samples` and `cell_populations` |

The schema is normalized so subject metadata is stored once rather than
repeated for every sample and every cell population. Primary keys prevent
duplicate records, foreign keys prevent orphaned records, and constraints
reject invalid values such as negative cell counts. Indexes cover fields used
frequently by the analysis, including condition, treatment, response, sample
type, time point, and population.

Cell populations are stored as rows instead of fixed database columns. A new
population can therefore be added to `cell_populations` and measured through
new `cell_counts` rows without altering the table structure. The loader also
discovers numeric cell-population columns from the input CSV automatically.

### Scaling

This structure supports hundreds of projects and thousands of samples by
adding rows to the existing tables. It also supports different analytics by
joining subject, sample, and measurement data without duplicating metadata.

## Code structure

```text
.
├── data/cell-count.csv              # Supplied input data
├── load_data.py                     # Creates and loads SQLite
├── analysis.py                      # Parts 2–4 and generated outputs
├── dashboard.py                     # Interactive Streamlit application
├── clinical_trial_patients.db       # Generated SQLite database
├── results/                         # Generated tables and plot
├── .streamlit/config.toml           # Dashboard server configuration
├── requirements.txt                 # Python dependencies
└── Makefile                         # Automated grading commands
```

The code separates data loading, analysis, and presentation so each stage has
one responsibility. `load_data.py` owns the database schema and ingestion.
`analysis.py` queries the database, generates reproducible results, and saves
them under `results/`. `dashboard.py` reads the database and generated outputs
to provide interactive filters and result displays.

## Analysis methods

For each sample, relative frequency is calculated as:

```text
percentage = 100 × population count / total sample count
```

The response comparison includes only PBMC samples from melanoma subjects
treated with miraclib. The boxplot displays relative frequencies for all
eligible samples. For statistical testing, frequencies are averaged per
subject across treatment time points, and responders are compared with
non-responders using a two-sided Mann–Whitney U test. P-values are adjusted
across the five cell populations using the Benjamini–Hochberg procedure.

The baseline subset contains melanoma PBMC samples collected at time `0` from
miraclib-treated subjects. The pipeline reports sample counts by project and
unique-subject counts by response and sex. Age is included as descriptive context because immune-cell distributions can vary with age.

The requested B-cell average uses male melanoma responders at time `0` across all treatments and sample types.

## Generated outputs

Running `make pipeline` creates:

- `results/frequency_summary.csv`
- `results/statistical_results.csv`
- `results/response_boxplot.png`
- `results/baseline_samples.csv`
- `results/baseline_project_counts.csv`
- `results/baseline_response_counts.csv`
- `results/baseline_sex_counts.csv`
- `results/average_b_cell_count.csv`
