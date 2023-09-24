# Thermodynamic Snookered

Thermodynamic Snookered is a project that gives a visual simulation of the behaviour of ideal gases in a circular container, as well as verifying the ideal gas equation by varying the temperature, pressure and velocity distribution.

The project consists of three Python scripts:
1. `Ball.py` - This script defines the Ball class where each ball represents a gas molecule with mass *m*, radius *R*, position *r* and velocity *v*.
2. `Simulation.py` - This script defines the Simulation class which initialises the initial state of the ideal gas, as well as calculates the post-collision state of the gas using the law of energy and momentum conservation.
3. `Task.py` - This script runs the simulation (without the visualisation) with different parameters (temperature, pressure and velocity distribution) to analyse the property of the ideal gas and to verify the ideal gas equation.
4. `run.py` - This just runs the visual simulation with the default parameters.

## Requirements
- Python 3.7 or later
- Libraries: numpy, scipy, pylab, matplotlib

## Usage

Before running the scripts, install the required Python libraries:

```
pip install -r requirements.txt
```

Then, you can run the scripts:

```
python Task.py
```
When a simulation is run, the frame number is printed to keep track of the progress of the simulation. 


# Example of the simulation


https://github.com/ansonpoon166/Thermodynamics-Snookered/assets/79363169/08e32835-da16-4cd9-aea2-a0f54cbfebd2
