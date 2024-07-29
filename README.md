# 2024-m1-nemo
This project aims to determine the influence of two squirmers on each other.

# Objectives
- Understand the "Squirmers" model
- Study the dynamics between two squirmers
- Develop a code to calculate the movement of two squirmers
- Verify if our results align with previous studies
- Study the impact of the **D** parameter on the groups behaviors

# Compilation
To compile the project you need to follow these steps:
```
python3 -m pip install --upgrade pip
pip install -r requirements.txt
sudo apt-get update
sudo apt-get install ffmpeg
```

# Reproducibility
```
python3 Code/run.py sim_choice N filename
```
with **N** the number of squirmers you want to simulate, 
**sim_choice** a parameter used to choose which simulation is done.

**sim_choice** can be:
- **video**, a simulation is done with **N** squirmers and the resulting video is stored in the **videos** directory
- **plot**, a simulation is done with **N** and the resulting graph is stored in the **graphs** directory
- **Eo_sim**, simulates the behavior of two interacting squirmers close to each other for 9 different Eo parameter and for each Eo, 8 simulations are done with 8 different initial orientation, not used anymore
- **border**, 20 simulations are done with one squirmer and one border with 4 different initial orientation and 5 beta values, the resulting graphs are stored in the **graphs/simulations/border** directory
- **sim_2_sq**, 40 simulations are done between two squirmers near each other with 8 different initial orientations and 5 different beta values, the resulting graphs or videos are stored in the **graphs/simulations/sim_sq_sq** or **videos/simulations/sim_sq_sq**
- **sim_D**, 3 simulations with **N** squirmers are made, used to study the impact of the **D** parameter on the behaviors of the squirmers
- **vicsek**, a simulation with **N** particles is done with 5 time-step the 6 resulting graphs are stored in the **graphs/vicsek** directory

**filename** a name file used with **video** and **plot**.

You can modify the simulation parameters by editing the file **Code/run.py** and the function **run** in **Code/interactingsquirmers**.

# Contributors
- Justine Ehret
- Robin Roth
