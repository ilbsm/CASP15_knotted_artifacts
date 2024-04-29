#! /usr/bin/env python
# -*- coding: utf-8 -*-
import csv
import glob
import os
import re
from collections import defaultdict

INTERSECTION_POINT = re.compile(
    r'"IP\d> 0: [ACGU](\d+?), (B|C4\'|P) - 0: [ACGU](\d+?), (B|C4\'|P)"'
)
SINGLE_STRAND = re.compile(r'"Single strand \d \(2D: \d+?-\d+?; 3D: 0(\d+?)-0(\d+?)\)"')
PATH = re.compile(r"(R\d{4}v?2?)/(R\d{4}v?2?TS\d{3}_\dx?)_report_normalized.csv")

ends = []
mids = []
by_depth = []

for path in glob.iglob("R*/*_report_normalized.csv"):
    with open(path) as f:
        reader = csv.reader(f)
        rows = [row for row in reader][1:]

    match = PATH.match(path)
    if match is None:
        raise RuntimeError(f"Cannot parse {path}")
    target = match.group(1)
    model = match.group(2)
    residues = set()

    model_paths = glob.glob(f"../pdb_models/{model}*.pdb")
    if len(model_paths) != 1:
        raise RuntimeError(
            f"Ambiguous model selection for model {model}: {model_paths}"
        )

    with open(model_paths[0]) as f:
        for line in f:
            if line.startswith("ATOM"):
                residues.add(int(line[22:26].strip()))

    i = 0
    while i < len(rows):
        entanglement = rows[i : i + 6]

        if entanglement[0][1] in ("L(S)"):
            match = INTERSECTION_POINT.match(entanglement[0][4])
            if match is None:
                raise RuntimeError(f"Cannot parse {entanglement[0][4]}")
            nts = sorted([int(match.group(1)), int(match.group(3))])

            match = SINGLE_STRAND.match(entanglement[3][3])
            if match is None:
                raise RuntimeError(f"Cannot parse {entanglement[3][3]}")
            single_strand = (int(match.group(1)), int(match.group(2)))
            end_5p = nts[0] - min(residues)
            end_3p = max(residues) - nts[1]

            if end_5p <= end_3p:
                ends.append(
                    [
                        model,
                        entanglement[0][1],
                        "5'",
                        end_5p,
                        nts[0],
                        min(residues),
                        max(residues),
                    ]
                )
                depth = end_5p
            else:
                ends.append(
                    [
                        model,
                        entanglement[0][1],
                        "3'",
                        end_3p,
                        nts[1],
                        min(residues),
                        max(residues),
                    ]
                )
                depth = end_3p

            by_depth.append(
                [
                    model,
                    entanglement[0][1],
                    depth,
                    "shallow" if depth <= 5 else "deep",
                ]
            )
        elif entanglement[0][1] in ("L(S.)", "L(L)"):
            match = INTERSECTION_POINT.match(entanglement[0][4])
            if match is None:
                raise RuntimeError(f"Cannot parse {entanglement[0][4]}")
            nts1 = sorted([int(match.group(1)), int(match.group(3))])

            match = INTERSECTION_POINT.match(entanglement[1][4])
            if match is None:
                raise RuntimeError(f"Cannot parse {entanglement[0][4]}")
            nts2 = sorted([int(match.group(1)), int(match.group(3))])

            mids.append(
                [model, entanglement[0][1], nts2[1] - nts1[0], nts1[0], nts2[1]]
            )

            depth = nts2[1] - nts1[0]
            by_depth.append(
                [
                    model,
                    entanglement[0][1],
                    depth,
                    "shallow" if depth <= 5 else "deep",
                ]
            )
        else:
            print(f"Warning: not calculating depth for {entanglement[0][1]} type")
            by_depth.append(
                [
                    model,
                    entanglement[0][1],
                    "n/a",
                    "n/a",
                ]
            )

        i += 7

with open("ends.csv", "w") as f:
    writer = csv.writer(f)
    writer.writerow(
        ["Model", "Type", "Which end", "Length", "IP", "First index", "Last index"]
    )
    writer.writerows(ends)

with open("mids.csv", "w") as f:
    writer = csv.writer(f)
    writer.writerow(["Model", "Type", "Length", "First IP", "Last IP"])
    writer.writerows(mids)

with open("by_depth.csv", "w") as f:
    writer = csv.writer(f)
    writer.writerow(["Model", "Type", "Depth", "Classification"])
    writer.writerows(by_depth)

interlaces = defaultdict(set)
shallows = defaultdict(set)
deeps = defaultdict(set)

is_ml = {
    "TS232": False,
    "TS287": False,
    "TS081": False,
    "TS128": False,
    "TS119": True,
    "TS416": True,
    "TS054": True,
    "TS125": True,
    "TS229": True,
    "TS110": True,
    "TS076": True,
    "TS035": False,
    "TS439": True,
    "TS147": True,
    "TS239": True,
    "TS248": False,
    "TS470": True,
    "TS227": False,
    "TS489": True,
    "TS392": False,
    "TS434": False,
    "TS235": False,
    "TS444": False,
    "TS325": False,
    "TS347": False,
    "TS456": False,
    "TS131": False,
    "TS185": True,
    "TS494": False,
    "TS285": False,
    "TS238": True,
    "TS029": True,
    "TS257": False,
    "TS163": False,
    "TS046": False,
    "TS097": True,
    "TS490": False,
    "TS245": False,
    "TS177": True,
    "TS385": False,
    "TS091": False,
}
is_ml = {k: "ML" if v else "nonML" for k, v in is_ml.items()}

for model, type_, depth, classification in by_depth:
    for subtype in type_.split("+"):
        if "(" in subtype:
            match = re.match(r"([DL]\([DLS])\.*(\))", subtype)
            assert match is not None
            subtype = match.group(1) + match.group(2)
        else:
            match = re.match(r"([DL])\.*(&[DL])\.*", subtype)
            assert match is not None
            subtype = match.group(1) + match.group(2)

        if classification == "shallow":
            shallows[model].add(subtype)
        elif classification == "deep":
            deeps[model].add(subtype)
        else:
            if "(" in subtype:
                deeps[model].add(subtype)
            else:
                interlaces[model].add(subtype)

with open("by_depth_alt.csv", "w") as f:
    reorganized = []

    for model in sorted([x[0] for x in by_depth]):
        reorganized.append(
            (
                model,
                "",
                is_ml[re.match(r".+(TS\d{3}).+", model).group(1)],
                ",".join(sorted(interlaces[model])),
                ",".join(sorted(shallows[model])),
                ",".join(sorted(deeps[model])),
                "",
            )
        )

    writer = csv.writer(f)
    writer.writerow(
        [
            "Model",
            "Human/Webserver",
            "ML/nonML",
            "Interlaces",
            "Shallow",
            "Deep",
            "Knots",
        ]
    )
    writer.writerows(reorganized)
