# 2024-m1-nemo
This project aims to determine the influence of two squirmers on each other.

# Objectives
- Understand the "Squirmers" model
- Study the dynamics between two squirmers
- Develop a code to calculate the movement of two squirmers
- Verify if our results align with previous studies

# Compilation
To compile the project you need to follow these steps:
```
python3 -m pip install --upgrade pip
pip install -r requirements.txt
sudo apt-get update
sudo apt-get install ffmpeg
```

# Reproducibility
To run the eight simulations done in the report you need to run:
```
python3 Code/run.py sim_choice filename
```
with **sim_choice** a parameter used to choose which simulation is done.

**sim_choice** can be:
- **vid_interact_sq**, a simulation between two interacting is done and a video is stored in the **videos** directory
- **interact_sq**, eight simulations of two interacting squirmers are done and the resulting graphs are stored in the **graphs** directory
- **sq_border**, two simulations of two squirmers interacting with a border are done, the resulting graph is stored in the **graphs** directory

**filename** a name file used with **vid_interact_sq** and **sq_border**.

You can modify the simulation parameters by editing the file **Code/run.py**.

# Contributors
- Justine Ehret
- Robin Roth
