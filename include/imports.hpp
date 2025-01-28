#ifndef IMPORTS_HPP
#define IMPORTS_HPP


const char* mechanics_imports = R"(import numpy as np
from scipy.integrate import solve_ivp
import scipy.constants as const
import sys
)";

const char* stats_imports = R"(import os
import numpy as np
import pandas as pd 
from scipy import stats
)";

const char* graph_imports = R"(import matplotlib.pyplot as plt)";

#endif // IMPORTS_HPP
