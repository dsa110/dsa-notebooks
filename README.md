# dsa-notebooks
Jupyter notebooks for sharing DSA array status and data quality

# Contents
- calibration_template.ipynb: Example notebook for inspecting calibration solutions for a single visibility data set.
- T2_template.ipynb: Example notebook for displaying heimdall output and buffer triggers.

# Observer on Duty Instructions
1) log in to dsa-storage using standard ssh config.
2) run `jupyter notebook list` to find server for OoD.
3) Paste notebook server URL into browser.
4) Run template notebook for new day of data
5) Save with name with date and push to github.

# Notebook developer Instructions
1) Create new notebook named "*_template.ipynb".
2) Run analysis on a day of data
3) Save, commit, and push template to github.