import sys
import numpy as np
sys.path.append("3rdParty/statelatticeplanner/")
import state_lattice_planner as sp
print("Only Python3 are supported, you are using: ")
print(sys.version)

def callLatticePlanner(LatticeParameter):
    lccs = sp.biased_terminal_state_sampling_test2(LatticeParameter)
    return lccs

def call_numberofline():
    return sp.number_of_curves()

if __name__ == '__main__':
    call_numberofline()
