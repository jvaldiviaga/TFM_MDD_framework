import pandas as pd
import os
import subprocess

deg_file = snakemake.input.deg
script_path = snakemake.input.script
output_file = snakemake.output.go_csv
log_file = snakemake.log[0]
dataset = snakemake.params.dataset

df = pd.read_csv(deg_file)
if df.empty or df["adj.P.Val"].dropna().lt(0.05).sum() == 0:
    print("❌ No DEGs significativos. Se omite enriquecimiento.")
    os.makedirs(os.path.dirname(output_file), exist_ok=True)
    with open(output_file, "w") as f:
        f.write("")  # archivo vacío
else:
    with open(log_file, "w") as logfile:
        subprocess.run(
            ["Rscript", script_path, dataset],
            stdout=logfile, stderr=logfile, check=True
        )
