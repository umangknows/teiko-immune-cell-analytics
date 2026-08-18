from __future__ import annotations

import os
from pathlib import Path

from .config import OUTPUT_DIR

OUTPUT_DIR.mkdir(exist_ok=True)
os.environ.setdefault("MPLCONFIGDIR", str(OUTPUT_DIR / ".matplotlib"))

import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns


def responder_boxplot(freq: pd.DataFrame, output_path: Path, title: str) -> None:
    sns.set_theme(style="whitegrid")
    fig, ax = plt.subplots(figsize=(11, 6))
    sns.boxplot(
        data=freq,
        x="population",
        y="percentage",
        hue="response",
        order=sorted(freq["population"].unique()),
        ax=ax,
        palette={"yes": "#2b8cbe", "no": "#f03b20"},
    )
    ax.set_title(title)
    ax.set_xlabel("Cell population")
    ax.set_ylabel("Relative frequency (%)")
    ax.legend(title="Response")
    fig.tight_layout()
    fig.savefig(output_path, dpi=180)
    plt.close(fig)
