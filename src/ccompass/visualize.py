"""Various visualization functions."""

from io import BytesIO

import matplotlib.colors as mcolors
import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns
from scipy.cluster.hierarchy import leaves_list, linkage

__all__ = [
    "fract_heatmap",
    "marker_correlation_heatmap",
    "fig_to_bytes",
    "plot_marker_profiles",
]


def fig_to_bytes(
    fig: plt.Figure = None, format: str = "png", close: bool = True
) -> bytes:
    """Convert a Matplotlib figure to bytes."""
    if fig is None:
        fig = plt.gcf()

    bio = BytesIO()
    fig.savefig(bio, format=format)

    if close:
        plt.close(fig)

    return bio.getvalue()


def fract_heatmap(df: pd.DataFrame, title: str = None):
    """Hierarchical clustering of fractionation data for a single condition."""
    # Perform hierarchical clustering on the rows
    linkage_matrix = linkage(df.values, method="ward")
    clustered_rows = leaves_list(
        linkage_matrix
    )  # Order of rows after clustering

    # Reorder the DataFrame rows based on hierarchical clustering
    df_clustered = df.iloc[clustered_rows, :]

    # Custom colormap: from #f2f2f2 (for value 0) to #6d6e71 (for value 1)
    cmap = mcolors.LinearSegmentedColormap.from_list(
        "custom_gray", ["#f2f2f2", "#6d6e71"], N=256
    )

    # Plot the heatmap using seaborn with the custom color gradient
    plt.figure(figsize=(8, 6))
    sns.heatmap(
        df_clustered,
        cmap=cmap,
        cbar=True,
        xticklabels=False,
        yticklabels=False,
        vmin=0,
        vmax=1,
    )

    if title:
        plt.title(title, fontsize=16)

    plt.tight_layout()


def marker_correlation_heatmap(
    dataframe: pd.DataFrame, title: str = None
) -> plt.Figure:
    """Create marker correlation heatmap."""
    fig, ax = plt.subplots(figsize=(8, 8))
    cax = ax.matshow(dataframe, cmap="coolwarm", vmin=-1, vmax=1)
    fig.colorbar(cax)
    ax.set_xticks(range(len(dataframe.columns)))
    ax.set_xticklabels(dataframe.columns, rotation=90, fontsize=8)
    # Rotate x-axis labels 90 degrees
    ax.set_yticks(range(len(dataframe.index)))
    ax.set_yticklabels(dataframe.index, fontsize=8)
    # Adjust the top and bottom margins
    plt.subplots_adjust(top=0.8, bottom=0.1, left=0.2)
    if title:
        plt.title(title)
    return fig


def plot_marker_profiles(df: pd.DataFrame, title=None) -> plt.Figure:
    """Plot marker profiles.

    Line plot with one line per compartment, fractions on the x-axis, and
    normalized intensity on the y-axis.
    """
    fig, ax = plt.subplots(figsize=(8, 6))

    for column in df.columns:
        ax.plot(df.index, df[column], label=column)

    ax.legend(loc="upper center", bbox_to_anchor=(0.5, -0.15), ncol=3)
    ax.set_xlabel("fractions")
    ax.set_ylabel("normalized intensity")
    ax.set_xticks([])  # Remove x-tick labels
    ax.set_yticks([0, 1])
    ax.set_yticklabels(["0", "1"])
    ax.set_xlim(0, len(df.index) - 1)  # Set x-axis limits
    ax.set_ylim(0, 1)

    if title:
        ax.set_title(title)

    # Adjust layout to make room for the legend
    fig.tight_layout(rect=(0, 0, 1, 0.95))
    return fig
