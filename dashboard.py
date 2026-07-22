from pathlib import Path
import sqlite3

import pandas as pd
import streamlit as st


ROOT = Path(__file__).resolve().parent
DB_PATH = ROOT / "clinical_trial_patients.db"
RESULTS_DIR = ROOT / "results"


@st.cache_data
def load_frequency_data(database_path: Path) -> pd.DataFrame:
    """Return one row per sample and cell population."""
    query = """
    SELECT
        s.sample_id AS sample,
        sub.project_id AS project,
        sub.subject_id AS subject,
        sub.condition,
        sub.age,
        sub.sex,
        sub.treatment,
        sub.response,
        s.sample_type,
        s.time_from_treatment_start,
        cc.population,
        cc.count,
        SUM(cc.count) OVER (PARTITION BY s.sample_id) AS total_count,
        100.0 * cc.count
            / SUM(cc.count) OVER (PARTITION BY s.sample_id) AS percentage
    FROM samples AS s
    JOIN subjects AS sub
        ON sub.subject_id = s.subject_id
    JOIN cell_counts AS cc
        ON cc.sample_id = s.sample_id
    ORDER BY s.sample_id, cc.population;
    """

    with sqlite3.connect(database_path) as connection:
        return pd.read_sql_query(query, connection)


st.set_page_config(
    page_title="Clinical Trial Cell Type Analysis",
    page_icon="🧬",
    layout="wide",
)

st.title("Clinical Trial Immune Cell Analysis")

if not DB_PATH.exists():
    st.error(
        "The database does not exist. Run `python load_data.py` before "
        "starting the dashboard."
    )
    st.stop()

frequency_data = load_frequency_data(DB_PATH)

cell_types = ", ".join(
    population.replace("_", " ")
    for population in sorted(frequency_data["population"].unique())
)
treatments = ", ".join(sorted(frequency_data["treatment"].unique()))
conditions = ", ".join(sorted(frequency_data["condition"].unique()))
sex_names = {"F": "female", "M": "male"}
sexes = ", ".join(
    sex_names.get(sex, sex) for sex in sorted(frequency_data["sex"].unique())
)
response_names = {"no": "non-responder", "yes": "responder"}
responses = ", ".join(
    response_names.get(response, response)
    for response in sorted(frequency_data["response"].dropna().unique())
)
minimum_age = int(frequency_data["age"].min())
maximum_age = int(frequency_data["age"].max())
minimum_day = int(frequency_data["time_from_treatment_start"].min())
maximum_day = int(frequency_data["time_from_treatment_start"].max())

st.markdown(
    f"""
This dataset contains immune-cell measurements for **{cell_types}**. Subjects
belong to the **{conditions}** condition groups and received **{treatments}**.
The dataset includes **{sexes}** subjects, with treatment response recorded as
**{responses}** where applicable. Ages range from **{minimum_age} to
{maximum_age} years**, and samples were collected from treatment day
**{minimum_day} through {maximum_day}**.
"""
)

sample_count, subject_count= st.columns(2)
sample_count.metric("Samples", frequency_data["sample"].nunique())
subject_count.metric("Subjects", frequency_data["subject"].nunique())


sample_metadata = frequency_data[
    [
        "sample",
        "project",
        "subject",
        "condition",
        "age",
        "sex",
        "treatment",
        "response",
        "sample_type",
        "time_from_treatment_start",
    ]
].drop_duplicates("sample")

st.divider()
st.subheader("Immune cell type frequencies by sample")
selection_method = st.radio(
    "Find a sample by",
    ["Sample ID", "Metadata"],
    horizontal=True,
)
selected_sample = None

if selection_method == "Sample ID":
    selected_sample = st.selectbox(
        "Sample ID",
        sorted(sample_metadata["sample"].unique()),
    )
else:
    metadata_filters = sample_metadata.copy()
    first_row = st.columns(4)
    selected_project = first_row[0].selectbox(
        "Project",
        ["All"] + sorted(sample_metadata["project"].unique()),
        key="sample_project",
    )
    selected_condition = first_row[1].selectbox(
        "Condition",
        ["All"] + sorted(sample_metadata["condition"].unique()),
        key="sample_condition",
    )
    selected_treatment = first_row[2].selectbox(
        "Treatment",
        ["All"] + sorted(sample_metadata["treatment"].unique()),
        key="sample_treatment",
    )
    selected_type = first_row[3].selectbox(
        "Sample type",
        ["All"] + sorted(sample_metadata["sample_type"].unique()),
        key="sample_type",
    )

    second_row = st.columns(3)
    selected_day = second_row[0].selectbox(
        "Treatment day",
        ["All"]
        + sorted(sample_metadata["time_from_treatment_start"].unique()),
        key="sample_day",
    )
    selected_sex = second_row[1].selectbox(
        "Sex",
        ["All"] + sorted(sample_metadata["sex"].unique()),
        key="sample_sex",
    )
    selected_response = second_row[2].selectbox(
        "Response",
        ["All", "no", "yes", "not recorded"],
        key="sample_response",
    )

    filters = {
        "project": selected_project,
        "condition": selected_condition,
        "treatment": selected_treatment,
        "sample_type": selected_type,
        "time_from_treatment_start": selected_day,
        "sex": selected_sex,
    }
    for column, value in filters.items():
        if value != "All":
            metadata_filters = metadata_filters[
                metadata_filters[column] == value
            ]
    if selected_response == "not recorded":
        metadata_filters = metadata_filters[
            metadata_filters["response"].isna()
        ]
    elif selected_response != "All":
        metadata_filters = metadata_filters[
            metadata_filters["response"] == selected_response
        ]

    st.caption(f"{len(metadata_filters):,} samples match these criteria.")
    if metadata_filters.empty:
        st.warning("No samples match the selected metadata.")

    selected_sample = st.selectbox(
        "Matching sample",
        [None] + sorted(metadata_filters["sample"].unique()),
        format_func=lambda sample: (
            "Select a sample (optional)" if sample is None else sample
        ),
    )

    matching_frequency_data = frequency_data[
        frequency_data["sample"].isin(metadata_filters["sample"])
    ]
    st.markdown("**All matching sample frequencies**")
    st.dataframe(
        matching_frequency_data[
            [
                "sample",
                "subject",
                "project",
                "population",
                "count",
                "total_count",
                "percentage",
            ]
        ],
        width="stretch",
        hide_index=True,
    )

if selected_sample is not None:
    selected_sample_data = frequency_data[
        frequency_data["sample"] == selected_sample
    ]
    st.markdown(f"**Cell proportions for {selected_sample}**")
    st.bar_chart(
        selected_sample_data.set_index("population")["percentage"],
        y_label="Relative frequency (%)",
    )

    st.dataframe(
        selected_sample_data[
            ["sample", "total_count", "population", "count", "percentage"]
        ],
        width="stretch",
        hide_index=True,
    )

st.divider()
st.header("Immune cell type frequencies by Miraclib response")
st.caption(
    "This analysis includes only PBMC samples from melanoma patients who "
    "received miraclib. Relative cell-population frequencies are compared "
    "between patients classified as responders and non-responders. For the "
    "statistical tests, each subject's frequencies are averaged across "
    "treatment time points."
)

statistics_path = RESULTS_DIR / "statistical_results.csv"
boxplot_path = RESULTS_DIR / "response_boxplot.png"

if statistics_path.exists() and boxplot_path.exists():
    statistical_results = pd.read_csv(statistics_path)
    st.image(
        str(boxplot_path),
        caption=(
            "Immune cell relative frequencies in PBMC from melanoma patients "
            " receiving Miraclib, comparing treatment response"
        ),
        width="stretch",
    )
    st.dataframe(
        statistical_results,
        width="stretch",
        hide_index=True,
    )

    raw_significant = statistical_results.loc[
        statistical_results["significant_unadjusted"], "population"
    ].tolist()
    fdr_significant = statistical_results.loc[
        statistical_results["significant_after_fdr"], "population"
    ].tolist()
    st.write(
        "Significant difference in relative frequencies at unadjusted α = 0.05:",
        ", ".join(raw_significant) if raw_significant else "None",
    )
    st.write(
        "Significant difference in relative frequencies after Benjamini–Hochberg correction:",
        ", ".join(fdr_significant) if fdr_significant else "None",
    )
else:
    st.warning("Run `make pipeline` to generate the response analysis.")

st.divider()
st.header("Baseline subset analysis")
st.caption(
    "This required subset contains melanoma PBMC samples collected at "
    "baseline (day 0) from subjects treated with miraclib."
)

subset_samples = sample_metadata[
    (sample_metadata["condition"] == "melanoma")
    & (sample_metadata["treatment"] == "miraclib")
    & (sample_metadata["sample_type"] == "PBMC")
    & (sample_metadata["time_from_treatment_start"] == 0)
]

subset_sample_count, subset_subject_count = st.columns(2)
subset_sample_count.metric("Matching samples", subset_samples["sample"].nunique())
subset_subject_count.metric(
    "Matching subjects", subset_samples["subject"].nunique()
)

project_counts = (
    subset_samples.groupby("project", as_index=False)
    .agg(sample_count=("sample", "nunique"))
)
response_counts = (
    subset_samples.groupby("response", as_index=False)
    .agg(subject_count=("subject", "nunique"))
)
sex_counts = (
    subset_samples.groupby("sex", as_index=False)
    .agg(subject_count=("subject", "nunique"))
)
project_column, response_column, sex_column = st.columns(3)
with project_column:
    st.subheader("Samples by project")
    st.dataframe(project_counts, hide_index=True, width="stretch")
with response_column:
    st.subheader("Subjects by response")
    st.dataframe(response_counts, hide_index=True, width="stretch")
with sex_column:
    st.subheader("Subjects by sex")
    st.dataframe(sex_counts, hide_index=True, width="stretch")

with st.expander("View matching samples"):
    st.dataframe(
        subset_samples,
        hide_index=True,
        width="stretch",
    )
