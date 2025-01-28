# PhysicsCompiler (The Physics Wizard)

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

### User Manual (user\_manual.pdf)
- Provides instructions on how to write input file and fill in each line to define your problem, choose equations, and specify parameter values/data.
- Includes sample problems. Corresponding input files are available in the `demo/` directory, including cases with grammar/syntax errors and insufficient data to demonstrate error handling.


## Usage 
1. Run `make` to generate the `physwiz` executable. 
2. See User Manual for input format instructions. 
3. Execute `./physwiz your_problem.phys`, replacing `your_problem.phys` with your input file.

A Python script, named after your input file, is generated and executed. The numerical solution is printed, and any requested graphical solutions are saved as images.

### Requirements
Python libraries:`scipy`, `numpy`, `pandas`, `matplotlib`

**Note:** If using macOS and encountering Flex-related errors, it is recommended to use a Linux environment. A `Vagrantfile` is included in the repo for this purpose (requires Vagrant and a virtual machine host such as VirtualBox).

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
