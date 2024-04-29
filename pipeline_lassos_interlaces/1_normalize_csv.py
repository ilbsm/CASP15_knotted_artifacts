#! /usr/bin/env python
import glob

import pandas as pd

for path in glob.iglob("R*/*report.csv"):
    df = pd.read_csv(path, sep=";")
    # write csv to file
    path = path.replace(".csv", "_normalized.csv")
    df.to_csv(path, index=False)
