# PhysicsWizard

Author: Gina Nasseri

---
## A Physics Algorithm Generator

Compiles user input into physics/stats solutions. Input format:

```
subject:
graph:
equation:
data: {}
```

### User Manual
- Input values for specifying a physics/stats problem are detailed here. 
- Sample problems are provided with input files provided in the demo directory along with input files with grammar, syntax errors, and insufficient data to demonstrate error output.  

## Usage 
1. Activate Virtual environment (must be using Linux/MacOS): `source venv/bin/activate`
2. Use `make` to generate the `physwiz` executable. 
3. See Physics Solution Generator User Manual for instructions on input format  
4. Use `./physwiz your_problem.phys` to generate solution, replacing `your_problem.phys` with your input file. 

A Python script with the same name as your input file will be generated and executed, and the solution is printed to output, any requested graphical solutions are saved as images.

**Note:** If using MacOS and Flex error encountered, it is recommended to use a Linux environment. A Vagrantfile is included in the repo which can be used (requires Vagrant and virtual machine software such as VirtualBox).

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
- venv/ - *Python environment*
