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
python3 Code/run.py sim_choice N filename
```
with **N** the number of squirmers you want to simulate, 
**sim_choice** a parameter used to choose which simulation is done.

**sim_choice** can be:
- **video**, a simulation is done with **N** squirmers and the resulting video is stored in the **videos** directory
- **plot**, a simulation is done with **N** and the resulting graph is stored in the **plot** directory
- **Eo_sim**, simulates the behavior of two interacting squirmers close to each other for 9 different Eo parameter and for each Eo, 8 simulations are done with 8 different initial orientation
- **border**, 4 simulations are done with one squirmer and one border with 4 different initial orientation, the resulting graph is stored in the **graphs/border** directory
- **vicsek**, a simulation with **N** particles is done with 5 time-step the 6 resulting graphs are stored in the **graphs/vicsek** directory

**filename** a name file used with **video**, **plot** and **border**.

You can modify the simulation parameters by editing the file **Code/run.py** and **Code/simulation**.

# Contributors
- Justine Ehret
- Robin Roth
