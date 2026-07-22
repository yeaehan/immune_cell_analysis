from pathlib import Path
import sqlite3

import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns
from scipy.stats import false_discovery_control, mannwhitneyu


ROOT = Path(__file__).resolve().parent
DB_PATH = ROOT / "clinical_trial_patients.db"
RESULTS_DIR = ROOT / "results"
FREQUENCY_OUTPUT_PATH = RESULTS_DIR / "frequency_summary.csv"
STATISTICS_OUTPUT_PATH = RESULTS_DIR / "statistical_results.csv"
BOXPLOT_OUTPUT_PATH = RESULTS_DIR / "response_boxplot.png"
BASELINE_SAMPLES_OUTPUT_PATH = RESULTS_DIR / "baseline_samples.csv"
PROJECT_COUNTS_OUTPUT_PATH = RESULTS_DIR / "baseline_project_counts.csv"
RESPONSE_COUNTS_OUTPUT_PATH = RESULTS_DIR / "baseline_response_counts.csv"
SEX_COUNTS_OUTPUT_PATH = RESULTS_DIR / "baseline_sex_counts.csv"
AVERAGE_B_CELL_OUTPUT_PATH = RESULTS_DIR / "average_b_cell_count.csv"


def calculate_frequency(connection: sqlite3.Connection) -> pd.DataFrame:
    """Calculate each cell population's relative frequency per sample."""
    query = """
    WITH sample_counts AS (
        SELECT
            sample_id,
            SUM(count) AS total_count
        FROM cell_counts
        GROUP BY sample_id
    )
    SELECT
        cc.sample_id AS sample,
        totals.total_count,
        cc.population,
        cc.count,
        100.0 * cc.count / NULLIF(totals.total_count, 0) AS percentage
    FROM cell_counts AS cc
    JOIN sample_counts AS totals
        ON totals.sample_id = cc.sample_id
    ORDER BY cc.sample_id, cc.population;
    """
    return pd.read_sql_query(query, connection)


def analyze_treatment_response(
    connection: sqlite3.Connection,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Return response frequencies and population-level statistical tests."""
    query = """
    WITH sample_counts AS (
        SELECT
            sample_id,
            SUM(count) AS total_count
        FROM cell_counts
        GROUP BY sample_id
    )
    SELECT
        s.sample_id AS sample,
        sub.subject_id AS subject,
        sub.response,
        cc.population,
        cc.count,
        totals.total_count,
        100.0 * cc.count / NULLIF(totals.total_count, 0) AS percentage
    FROM samples AS s
    JOIN subjects AS sub
        ON sub.subject_id = s.subject_id
    JOIN cell_counts AS cc
        ON cc.sample_id = s.sample_id
    JOIN sample_counts AS totals
        ON totals.sample_id = s.sample_id
    WHERE LOWER(sub.condition) = 'melanoma'
      AND LOWER(sub.treatment) = 'miraclib'
      AND UPPER(s.sample_type) = 'PBMC'
      AND sub.response IN ('yes', 'no')
    ORDER BY cc.population, sub.response, s.sample_id;
    """
    response_frequencies = pd.read_sql_query(query, connection)

    # Give each subject equal weight by averaging their eligible time points.
    subject_means = (
        response_frequencies.groupby(
            ["subject", "response", "population"], as_index=False
        )["percentage"]
        .mean()
    )

    results = []
    for population, group in subject_means.groupby("population"):
        responders = group.loc[
            group["response"] == "yes", "percentage"
        ]
        non_responders = group.loc[
            group["response"] == "no", "percentage"
        ]

        if responders.empty or non_responders.empty:
            statistic = float("nan")
            p_value = float("nan")
        else:
            test = mannwhitneyu(
                responders,
                non_responders,
                alternative="two-sided",
            )
            statistic = float(test.statistic)
            p_value = float(test.pvalue)

        results.append(
            {
                "population": population,
                "responder_subjects": len(responders),
                "non_responder_subjects": len(non_responders),
                "responder_mean_percentage": responders.mean(),
                "non_responder_mean_percentage": non_responders.mean(),
                "mann_whitney_u": statistic,
                "p_value": p_value,
            }
        )

    statistics = pd.DataFrame(results).sort_values("p_value")
    statistics["adjusted_p_value"] = false_discovery_control(
        statistics["p_value"], method="bh"
    )
    statistics["significant_unadjusted"] = statistics["p_value"] < 0.05
    statistics["significant_after_fdr"] = (
        statistics["adjusted_p_value"] < 0.05
    )
    return response_frequencies, statistics.reset_index(drop=True)


def create_response_boxplot(
    response_frequencies: pd.DataFrame,
    statistical_results: pd.DataFrame,
    output_path: Path = BOXPLOT_OUTPUT_PATH,
) -> None:
    """Save responder versus non-responder boxplots for all populations."""
    plot_data = response_frequencies.copy()
    plot_data["Response"] = plot_data["response"].map(
        {"yes": "Responder", "no": "Non-responder"}
    )

    population_order = statistical_results["population"].tolist()
    p_values = statistical_results.set_index("population")[
        "adjusted_p_value"
    ]

    figure, axis = plt.subplots(figsize=(13, 7))
    sns.boxplot(
        data=plot_data,
        x="population",
        y="percentage",
        hue="Response",
        order=population_order,
        hue_order=["Non-responder", "Responder"],
        palette={
            "Non-responder": "#9E9E9E",
            "Responder": "#D98282",
        },
        showfliers=False,
        ax=axis,
    )
    axis.set(
        title="Immune Cell Frequencies: Miraclib Responders vs Non-responders",
        xlabel="Cell population",
        ylabel="Relative frequency (%)",
    )
    axis.set_xticks(axis.get_xticks(), population_order)
    annotation_heights = []
    for position, population in enumerate(population_order):
        population_data = plot_data.loc[
            plot_data["population"] == population
        ]
        upper_whiskers = []
        for response in ["Non-responder", "Responder"]:
            values = population_data.loc[
                population_data["Response"] == response, "percentage"
            ]
            first_quartile = values.quantile(0.25)
            third_quartile = values.quantile(0.75)
            upper_limit = third_quartile + 1.5 * (
                third_quartile - first_quartile
            )
            upper_whiskers.append(values[values <= upper_limit].max())

        annotation_height = max(upper_whiskers) + 0.8
        annotation_heights.append(annotation_height)
        axis.text(
            position,
            annotation_height,
            f"adj. p = {p_values[population]:.3f}",
            ha="center",
            va="bottom",
            fontsize=10,
            fontweight="bold",
        )
    axis.set_ylim(
        top=max(axis.get_ylim()[1], max(annotation_heights) + 1.0)
    )
    figure.tight_layout()
    figure.savefig(output_path, dpi=300, bbox_inches="tight")
    plt.close(figure)


def analyze_baseline_subset(
    connection: sqlite3.Connection,
) -> tuple[
    pd.DataFrame,
    pd.DataFrame,
    pd.DataFrame,
    pd.DataFrame,
    float,
]:
    """Return the requested miraclib/melanoma/PBMC baseline summaries."""
    query = """
    SELECT
        s.sample_id AS sample,
        sub.subject_id AS subject,
        sub.project_id AS project,
        sub.response,
        sub.sex
    FROM samples AS s
    JOIN subjects AS sub
        ON sub.subject_id = s.subject_id
    WHERE LOWER(sub.condition) = 'melanoma'
      AND LOWER(sub.treatment) = 'miraclib'
      AND UPPER(s.sample_type) = 'PBMC'
      AND s.time_from_treatment_start = 0
    ORDER BY sub.project_id, s.sample_id;
    """
    baseline_samples = pd.read_sql_query(query, connection)

    project_counts = (
        baseline_samples.groupby("project", as_index=False)
        .agg(sample_count=("sample", "nunique"))
    )
    response_counts = (
        baseline_samples.groupby("response", as_index=False)
        .agg(subject_count=("subject", "nunique"))
    )
    sex_counts = (
        baseline_samples.groupby("sex", as_index=False)
        .agg(subject_count=("subject", "nunique"))
    )
    average_query = """
    SELECT AVG(cc.count)
    FROM cell_counts AS cc
    JOIN samples AS s
        ON s.sample_id = cc.sample_id
    JOIN subjects AS sub
        ON sub.subject_id = s.subject_id
    WHERE cc.population = 'b_cell'
      AND LOWER(sub.condition) = 'melanoma'
      AND sub.sex = 'M'
      AND sub.response = 'yes'
      AND s.time_from_treatment_start = 0;
    """
    result = connection.execute(average_query).fetchone()[0]
    if result is None:
        raise ValueError("No samples matched the average B-cell query.")
    average_b_cells = round(float(result), 2)

    return (
        baseline_samples,
        project_counts,
        response_counts,
        sex_counts,
        average_b_cells,
    )


def main() -> None:
    """Generate all required analysis tables and plots."""
    if not DB_PATH.exists():
        raise FileNotFoundError(
            f"Database not found at {DB_PATH}. Run `python load_data.py` first."
        )

    with sqlite3.connect(DB_PATH) as connection:
        frequency_summary = calculate_frequency(connection)
        response_frequencies, statistical_results = (
            analyze_treatment_response(connection)
        )
        (
            baseline_samples,
            project_counts,
            response_counts,
            sex_counts,
            average_b_cells,
        ) = analyze_baseline_subset(connection)

    RESULTS_DIR.mkdir(parents=True, exist_ok=True)
    frequency_summary.to_csv(FREQUENCY_OUTPUT_PATH, index=False)
    statistical_results.to_csv(STATISTICS_OUTPUT_PATH, index=False)
    baseline_samples.to_csv(BASELINE_SAMPLES_OUTPUT_PATH, index=False)
    project_counts.to_csv(PROJECT_COUNTS_OUTPUT_PATH, index=False)
    response_counts.to_csv(RESPONSE_COUNTS_OUTPUT_PATH, index=False)
    sex_counts.to_csv(SEX_COUNTS_OUTPUT_PATH, index=False)
    pd.DataFrame(
        {"average_b_cell_count": [average_b_cells]}
    ).to_csv(AVERAGE_B_CELL_OUTPUT_PATH, index=False)
    create_response_boxplot(response_frequencies, statistical_results)

    print(
        f"Saved analysis outputs to {RESULTS_DIR}.\n"
        "Average B-cell count for male melanoma responders at baseline "
        "across all treatments and sample types: "
        f"{average_b_cells:.2f}"
    )


if __name__ == "__main__":
    main()
