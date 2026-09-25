from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.optimize import curve_fit
from scipy.stats import chi2


PROJECT_ROOT = Path(__file__).resolve().parents[1]

FIGURE_DIR = PROJECT_ROOT / "figures"
FIGURE_DIR.mkdir(parents=True, exist_ok=True)

PANTHEON_URL = (
    "https://raw.githubusercontent.com/dscolnic/Pantheon/"
    "master/lcparam_full_long.txt"
)

PANTHEON_COLUMNS = [
    "name", "zcmb", "zhel", "dz", "mb", "dmb", "x1", "dx1",
    "color", "dcolor", "3rdvar", "d3rdvar", "cov_m_s",
    "cov_m_c", "cov_s_c", "set", "ra", "dec",
]

Z_MIN = 0.01
Z_MAX = 0.30



# Models


def hubble_law_magnitude(redshift, intercept):
    return intercept + 5 * np.log10(redshift)


def cosmographic_magnitude(redshift, intercept, q0):
    distance_shape = redshift * (
        1 + 0.5 * (1 - q0) * redshift
    )

    return intercept + 5 * np.log10(distance_shape)



# Diagnostics


def normalized_residuals(data, model_function, *parameters):

    prediction = model_function(
        data["zcmb"],
        *parameters
    )

    return (
        data["mb"] - prediction
    ) / data["dmb"]


def goodness_of_fit(data, model_function, *parameters):

    residuals = normalized_residuals(
        data,
        model_function,
        *parameters
    )

    chi2_value = np.sum(residuals**2)

    ndf = len(data) - len(parameters)

    reduced_chi2 = chi2_value / ndf

    p_value = chi2.sf(
        chi2_value,
        ndf
    )

    return {
        "chi2": chi2_value,
        "ndf": ndf,
        "reduced_chi2": reduced_chi2,
        "p_value": p_value,
    }



# Load data


def load_data():

    raw_supernovae = pd.read_csv(
        PANTHEON_URL,
        sep=r"\s+",
        skiprows=1,
        names=PANTHEON_COLUMNS,
    )

    mask = (
        raw_supernovae["zcmb"].between(Z_MIN, Z_MAX)
        & np.isfinite(raw_supernovae["mb"])
        & np.isfinite(raw_supernovae["dmb"])
        & (raw_supernovae["dmb"] > 0)
    )

    analysis_data = raw_supernovae.loc[mask].copy()

    analysis_data = (
        analysis_data
        .sort_values("zcmb")
        .reset_index(drop=True)
    )

    return analysis_data



# Fit models


def fit_models(analysis_data):

    x = analysis_data["zcmb"].to_numpy()
    y = analysis_data["mb"].to_numpy()
    sigma_y = analysis_data["dmb"].to_numpy()


    hubble_popt, hubble_pcov = curve_fit(
        hubble_law_magnitude,
        x,
        y,
        p0=[25.0],
        sigma=sigma_y,
        absolute_sigma=True,
    )


    cosmo_popt, cosmo_pcov = curve_fit(
        cosmographic_magnitude,
        x,
        y,
        p0=[25.0, -0.5],
        sigma=sigma_y,
        absolute_sigma=True,
        bounds=(
            [0.0, -10.0],
            [50.0, 3.0],
        ),
    )

    return (
        hubble_popt,
        hubble_pcov,
        cosmo_popt,
        cosmo_pcov,
    )



# Hubble diagram


def plot_hubble_diagram(
    analysis_data,
    hubble_popt,
    cosmo_popt,
):

    x = analysis_data["zcmb"].to_numpy()
    y = analysis_data["mb"].to_numpy()
    sigma_y = analysis_data["dmb"].to_numpy()

    z_grid = np.linspace(
        x.min(),
        x.max(),
        500,
    )

    hubble = hubble_law_magnitude(
        z_grid,
        *hubble_popt
    )

    cosmo = cosmographic_magnitude(
        z_grid,
        *cosmo_popt
    )

    fig, axes = plt.subplots(
        2,
        1,
        figsize=(10, 7),
        sharex=True,
    )

    # Hubble model
    axes[0].errorbar(
        x,
        y,
        yerr=sigma_y,
        fmt=".",
        capsize=2,
        alpha=0.5,
        label="Supernova data",
    )

    axes[0].plot(
        z_grid,
        hubble,
        linewidth=2,
        label="Hubble model",
    )

    axes[0].set_title("Hubble Model")
    axes[0].set_ylabel("Magnitude")
    axes[0].set_xscale("log")
    axes[0].legend()

    # Cosmographic model
    axes[1].errorbar(
        x,
        y,
        yerr=sigma_y,
        fmt=".",
        capsize=2,
        alpha=0.5,
        label="Supernova data",
    )

    axes[1].plot(
        z_grid,
        cosmo,
        linewidth=2,
        label="Cosmographic model",
    )

    axes[1].set_title("Cosmographic Model")
    axes[1].set_xlabel("Redshift")
    axes[1].set_ylabel("Magnitude")
    axes[1].set_xscale("log")
    axes[1].legend()

    fig.tight_layout()

    output_path = (
        FIGURE_DIR / "supernova_hubble_fit.png"
    )

    fig.savefig(
        output_path,
        dpi=200,
        bbox_inches="tight",
    )

    plt.close(fig)

    print(f"Saved {output_path}")



# Residual plots


def plot_residuals(
    analysis_data,
    hubble_popt,
    cosmo_popt,
):

    data = analysis_data.copy()

    data["hubble_residual"] = normalized_residuals(
        data,
        hubble_law_magnitude,
        *hubble_popt
    )

    data["cosmo_residual"] = normalized_residuals(
        data,
        cosmographic_magnitude,
        *cosmo_popt
    )

    fig, axes = plt.subplots(
        2,
        1,
        figsize=(10, 7),
        sharex=True,
        sharey=True,
    )

    # Hubble residuals
    axes[0].axhline(
        0,
        linestyle="--",
        linewidth=1,
    )

    axes[0].scatter(
        data["zcmb"],
        data["hubble_residual"],
        s=15,
    )

    axes[0].set_title("Hubble Model Residuals")
    axes[0].set_ylabel("Normalized residual")
    axes[0].set_xscale("log")

    # Cosmographic residuals
    axes[1].axhline(
        0,
        linestyle="--",
        linewidth=1,
    )

    axes[1].scatter(
        data["zcmb"],
        data["cosmo_residual"],
        s=15,
    )

    axes[1].set_title(
        "Cosmographic Model Residuals"
    )

    axes[1].set_xlabel("Redshift")
    axes[1].set_ylabel("Normalized residual")
    axes[1].set_xscale("log")

    fig.tight_layout()

    output_path = (
        FIGURE_DIR / "supernova_residuals.png"
    )

    fig.savefig(
        output_path,
        dpi=200,
        bbox_inches="tight",
    )

    plt.close(fig)

    print(f"Saved {output_path}")



# Summary table


def build_summary_table(
    analysis_data,
    hubble_popt,
    hubble_pcov,
    cosmo_popt,
    cosmo_pcov,
):

    hubble_gof = goodness_of_fit(
        analysis_data,
        hubble_law_magnitude,
        *hubble_popt
    )

    cosmo_gof = goodness_of_fit(
        analysis_data,
        cosmographic_magnitude,
        *cosmo_popt
    )


    hubble_uncertainties = np.sqrt(
        np.diag(hubble_pcov)
    )

    cosmo_uncertainties = np.sqrt(
        np.diag(cosmo_pcov)
    )


    delta_chi2 = (
        hubble_gof["chi2"]
        - cosmo_gof["chi2"]
    )

    delta_parameters = (
        len(cosmo_popt)
        - len(hubble_popt)
    )

    delta_chi2_p_value = chi2.sf(
        delta_chi2,
        delta_parameters,
    )

    summary = pd.DataFrame({
        "model": [
            "Hubble law",
            "Cosmographic",
        ],

        "intercept": [
            hubble_popt[0],
            cosmo_popt[0],
        ],

        "intercept_uncertainty": [
            hubble_uncertainties[0],
            cosmo_uncertainties[0],
        ],

        "q0": [
            np.nan,
            cosmo_popt[1],
        ],

        "q0_uncertainty": [
            np.nan,
            cosmo_uncertainties[1],
        ],

        "chi2": [
            hubble_gof["chi2"],
            cosmo_gof["chi2"],
        ],

        "ndf": [
            hubble_gof["ndf"],
            cosmo_gof["ndf"],
        ],

        "reduced_chi2": [
            hubble_gof["reduced_chi2"],
            cosmo_gof["reduced_chi2"],
        ],

        "p_value": [
            hubble_gof["p_value"],
            cosmo_gof["p_value"],
        ],

        "delta_chi2": [
            delta_chi2,
            np.nan,
        ],

        "delta_parameters": [
            delta_parameters,
            np.nan,
        ],

        "delta_chi2_p_value": [
            delta_chi2_p_value,
            np.nan,
        ],
    })

    return summary



# Main analysis


def main():

    analysis_data = load_data()


    (
        hubble_popt,
        hubble_pcov,
        cosmo_popt,
        cosmo_pcov,
    ) = fit_models(analysis_data)

    # ----------------------------------------
    # Display fitted parameters
    # ----------------------------------------

    hubble_uncertainty = np.sqrt(
        hubble_pcov[0, 0]
    )

    cosmo_uncertainties = np.sqrt(
        np.diag(cosmo_pcov)
    )

    print("\nBest-fit parameters")

    print(
        f"Hubble intercept: "
        f"{hubble_popt[0]:.6f} "
        f"+/- {hubble_uncertainty:.6f}"
    )

    print(
        f"Cosmographic intercept: "
        f"{cosmo_popt[0]:.6f} "
        f"+/- {cosmo_uncertainties[0]:.6f}"
    )

    print(
        f"Cosmographic q0: "
        f"{cosmo_popt[1]:.6f} "
        f"+/- {cosmo_uncertainties[1]:.6f}"
    )

    # ----------------------------------------
    # Figures
    # ----------------------------------------

    plot_hubble_diagram(
        analysis_data,
        hubble_popt,
        cosmo_popt,
    )

    plot_residuals(
        analysis_data,
        hubble_popt,
        cosmo_popt,
    )

    # ----------------------------------------
    # Summary
    # ----------------------------------------

    summary = build_summary_table(
        analysis_data,
        hubble_popt,
        hubble_pcov,
        cosmo_popt,
        cosmo_pcov,
    )

    print("\nFit summary:")
    print(summary.to_string(index=False))

    summary_path = (
        FIGURE_DIR / "supernova_fit_summary.csv"
    )

    summary.to_csv(
        summary_path,
        index=False,
    )

    print(f"\nSaved {summary_path}")


# ============================================================
# Command-line entry point
# ============================================================

if __name__ == "__main__":
    main()