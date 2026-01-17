#!/usr/bin/python

import os
import sys
import json
import random
import numpy as np
import pandas as pd

## Set parameters --------------------------------------------------------------

NUM_SAMPLES = 10000
SAMPLE_SIZE = 45
RANDOM_SEED = random.randint(0, 10000)

## Check arguments -------------------------------------------------------------

if len(sys.argv) < 2:
    print("Usage: python Make_sampled_data.py <subject_id>")
    sys.exit(1)
else:
    sid = sys.argv[1]

## Set paths -------------------------------------------------------------------

author = ["Chang_et_al", "Tse_et_al"][0]
task = ["Naming", "LD"][0]
data_folder = os.path.join("Data_Linguistic", author, task)
out_folder = os.path.join(data_folder, "derivatives")
if not os.path.exists(out_folder):
    os.mkdir(out_folder)

## Load data -------------------------------------------------------------------

zscored = ["", "zscored_"][0]
data_path = os.path.join(data_folder, f"{zscored}sub_{sid}.xlsx")
print(f"Loading data from: {data_path}")

DF = pd.read_excel(data_path)
all_indices = list(DF.index)

## Start sampling --------------------------------------------------------------

print(f"Randomly sampling {NUM_SAMPLES} combinations (n={SAMPLE_SIZE}) ...")
random.seed(RANDOM_SEED)
sampled_indices = set()
save_indices = {}
count = 0

while count < NUM_SAMPLES:
    candidate = tuple(sorted(random.sample(all_indices, SAMPLE_SIZE)))

    if candidate not in sampled_indices:
        count += 1
        sampled_indices.add(candidate)
        sub_DF = DF.loc[candidate, :]
        sub_DF.to_excel(
            os.path.join(out_folder, f"{zscored}sub-{sid}_seed={RANDOM_SEED}_{count:05}.xlsx"), 
            index=False
        )
        save_indices[count] = list(sub_DF.index)

fn2 = f"sub-{sid}_seed={RANDOM_SEED}_indices.json"
with open(os.path.join(out_folder, fn2), "w") as f:
    json.dump(save_indices, f)

print("Done!\n")