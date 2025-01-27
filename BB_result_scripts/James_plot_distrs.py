#!/usr/bin/env python

""" Plot log data.

Example:
    ./plot.py -o output.pdf path/to/log_dir1 path/to/log_dir2 ...

The log directory should be of the form "{model_name}_sf/{part}_logs".
"""

from   pathlib import Path
import altair as alt
import argparse
import pandas as pd
import re
import numpy as np

face_pattern = re.compile(r"f(\d+)")
ipc_pattern = re.compile(r"IPC_prob\s?:\s+(-?\d+(\.\d+)?(e-\d+)?)")
bullet_pattern = re.compile(r"bullet_prob\s?:\s+(-?\d+(\.\d+)?(e-\d+)?)")
ours_pattern = re.compile(r"my_prob\s?:\s+(-?\d+(\.\d+)?(e-\d+)?)")


def parse_args():
    parser = argparse.ArgumentParser(description="Plot log data")
    parser.add_argument("log_dir", type=str, nargs="*", help="Log directory")
    parser.add_argument("-o", "--output", type=str, help="Output file")
    return parser.parse_args()


def parse_ipc_log(ipc_log):
    assert ipc_log.exists(), f"{ipc_log} does not exist"
    faces = []
    ipc_prob = []
    ours_prob = []
    with open(ipc_log, "r") as f:
        for line in f:
            r = re.search(face_pattern, line)
            if r is not None:
                face = int(r.group(1))
                faces.append(face)
            r = re.search(ipc_pattern, line)
            if r is not None:
                ipc = float(r.group(1))
                ipc_prob.append(ipc)
            r = re.search(ours_pattern, line)
            if r is not None:
                ours = float(r.group(1))
                ours_prob.append(ours)

    return faces, ipc_prob, ours_prob


def parse_bullet_log(bullet_log):
    assert bullet_log.exists(), f"{bullet_log} does not exist"
    faces = []
    bullet_prob = []
    with open(bullet_log, "r") as f:
        for line in f:
            r = re.search(face_pattern, line)
            if r is not None:
                face = int(r.group(1))
                faces.append(face)
            r = re.search(bullet_pattern, line)
            if r is not None:
                bullet = float(r.group(1))
                bullet_prob.append(bullet)

    return faces, bullet_prob


def extract_logs(log_dir):
    basename = log_dir.stem.replace("_logs", "")
    model_name = log_dir.parent.stem.replace("_sf", "")
    model_name = int(model_name)
    ipc_log = log_dir / f"{basename}_IPC.txt"
    bullet_log = log_dir / f"{basename}_bullet.txt"

    if not ipc_log.exists() or not bullet_log.exists():
        return None

    faces, ipc_prob, ours_prob = parse_ipc_log(ipc_log)
    faces2, bullet_prob = parse_bullet_log(bullet_log)
    assert faces == faces2, "IPC and Bullet have different faces"
    assert len(ipc_prob) == len(ours_prob), "IPC and Ours have different number of data"
    assert len(ipc_prob) == len(
        bullet_prob
    ), "IPC and Bullet have different number of data"

    # Reorder data based on ours prob.
    # order = np.argsort(ours_prob)[::-1]
    # faces = [faces[i] for i in order]
    # bullet_prob = [bullet_prob[i] for i in order]
    # ipc_prob = [ipc_prob[i] for i in order]
    # ours_prob = [ours_prob[i] for i in order]

    data = pd.DataFrame(
        {
            "model": [model_name] * len(faces),
            "part": [basename] * len(faces),
            "face": faces,
            "Bullet": bullet_prob,
            "IPC": ipc_prob,
            "Ours": ours_prob,
        }
    )
    if data.Ours.min() < -0.001:
        print(f"{model_name} {basename} has negative probability")
        return None
    data = data.sort_values("Ours", ascending=False).reset_index(drop=True)
    data["index"] = np.arange(len(data), dtype=int)

    return data


def main():
    args = parse_args()
    log_data = []
    for log_dir in args.log_dir:
        log = extract_logs(Path(log_dir))
        if log is not None:
            log_data.append(log)
    data = pd.concat(log_data)

    chart = (
        alt.Chart(data)
        .mark_bar()
        .encode(
            x=alt.X("value:Q").title(None).scale(domain=[0, 1]),
            y=alt.Y("variable:N")
            .title(None)
            .axis(labelColor="gray", labelPadding=2, labelLimit=30, ticks=False, domain=False),
            color=alt.Color("index:N").legend(None),
        )
        .transform_fold(
            ["Bullet", "IPC", "Ours"],
            as_=["variable", "value"],
        )
        .properties(
            width=150,
            height=30,  # title=f"Model {model_name} - part {basename}"
        )
        .facet(
            row=alt.Row("model:N")
            .header(
                # labelOrient="top",
                labelAngle=0,
                labelAlign="right",
                #labelAnchor="end",
                labelColor="gray",
                labelFontSize=16,
                labelPadding=0,
                titlePadding=0,
                #labelLimit=50,
                #labelFontWeight="bold",
            )
            .title(None),
            column=alt.Column("part:N")
            .header(
                # labelOrient="top",
                labelAnchor="end",
                labelColor="gray",
                labelFontSize=12,
                labelPadding=0,
            )
            .title(None),
            spacing=5,
        )
        .configure_title(fontSize=12, anchor="end", color="gray")
    )
    chart.save(args.output)


if __name__ == "__main__":
    main()