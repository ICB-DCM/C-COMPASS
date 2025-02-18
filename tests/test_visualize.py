"""Visualization tests."""

import matplotlib.pyplot as plt
import pandas as pd

from ccompass.visualize import (
    fract_heatmap,
    marker_correlation_heatmap,
    plot_marker_profiles,
)


def test_fract_heatmap():
    """Test fract_heatmap."""
    df = pd.DataFrame(
        {
            "A": [0.1, 0.2, 0.3, 0.4],
            "B": [0.2, 0.3, 0.4, 0.5],
            "C": [0.3, 0.4, 0.5, 0.6],
            "D": [0.4, 0.5, 0.6, 0.7],
        }
    )
    fract_heatmap(df)
    # plt.show()
    plt.close()


def test_marker_correlation_heatmap():
    """Test marker_correlation_heatmap."""
    df = pd.DataFrame(
        {
            "A": [0.1, 0.2, 0.3, 0.4],
            "B": [0.2, 0.3, 0.4, 0.5],
            "C": [0.3, 0.4, 0.5, 0.6],
            "D": [0.4, 0.5, 0.6, 0.7],
        }
    )

    marker_correlation_heatmap(df)
    # plt.show()
    plt.close()


def test_plot_marker_profiles():
    """Test plot_marker_profiles."""
    df = pd.DataFrame(
        {
            "A": [0.1, 0.2, 0.3, 0.4, 0.5],
            "B": [0.2, 0.3, 0.4, 0.5, 0.6],
            "C": [0.3, 0.4, 0.5, 0.6, 0.7],
            "D": [0.4, 0.5, 0.6, 0.7, 0.8],
        }
    )
    fig = plot_marker_profiles(df)
    # check there is one line for each column
    assert len(fig.get_axes()[0].lines) == len(df.columns)
    # plt.show()
    plt.close()
