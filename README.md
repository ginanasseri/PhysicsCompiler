# PhysicsCompiler

## A Physics Algorithm Generator

**Author:** Gina Nasseri

---
Compiles user input into physics/stats solutions. User input is translated into Python scripts which are then executed to produce the generated plots and critical values of solution. The values are printed into a table, while the Python script and generated plots will be saved to the current directory. Example problems are provided and detailed in the User Manual. The resulting Python scripts are provided in the `py_scripts` directory.

Input format: create a file containing the following four lines. 

```
subject:
graph:
equation:
data: {}
```

### User Manual (user\_manual.pdf)
- Provides instructions on how to fill in each line to choose equations and specify parameter values/data.
- Includes sample problems. Corresponding input files are available in the `demo/` directory, including cases with grammar/syntax errors and insufficient data to demonstrate error handling.

### System Requirements
- Must have `flex` installed (Linux recommended) 

## Usage 
1. Run `make clean` then `make` to generate the `physwiz` executable. 
2. See User Manual for input format instructions. 
3. Execute `./physwiz your_problem.phys`, replacing `your_problem.phys` with your input file.

A Python script, named after your input file, is generated and executed. The numerical solution is printed, and any requested graphical solutions are saved as images.

### Requirements 
Python libraries:`scipy`, `numpy`, `pandas`, `matplotlib`

**Note:** Developed in Linux environment. 
If encountereing flex-related errors [`ld: library 'fl' not found`], a Linux environment is recommended (e.g., Vagrant Ubuntu; multipass if on MacOS, then set up a Python virtual environment to avoid Python library installation).

---

## Directory Structure
- demo/
  - `bad_equation_num.phys`
  - `bad_grammar.phys`
  - `boeing.phys` *free-fall example*
  - `carronade_cannon.phys` *projectile motion example*
  - `stats_example.phys` *descriptive stats example*
  - `sample_data.csv` *data used in stats example - make sure in directory*
- include/ – *Header files and Python scripts*
- obj/ – *Object files*
- physwiz – *The compiler*
- src/  
  - `lexer.cpp`
  - `lexer.l`
  - `main.cpp`
  - `parser.cpp` *main control flow* 
  - `writer.cpp` *script generator*
- `py\_scripts/`
  - `boeing.py` *script generated from running boeing.phys*
  - `stats_example.py` *script generated from running stats_example.phys*
- `sample_plots/`
  - `free_fall.png` *plots generated from running boeing.phys*
  - `scatter.png` *generated from running stats_example.phys, so are the following two:* 
  - `x_hist.py`
  - `y_hist.py`
