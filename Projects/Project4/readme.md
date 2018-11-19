## How to run codebase

This program allows for the simulation of the Ising Model using the Monte Carlo method coupled with the Metropolis algorithm.

### User Input for T, L, # of Monte Carlo cycles, and spin orientation to solve their choice of system pertaining to this situation.
```python
python main.py
```

### To run selected test cases (as requested by prompt)
All of the actual computation results are included in the Github Repo and by default the script reads previously computed simulations from /results .txt files. It is highly reccomended that you do not rerun the unit cases due to computational time, but if you would like to actually run the cases and produce new data / .txt files, read below. 

```python
python testing.py
```

All plots will be available in /results

### Running testing.py to generate new data

By default at the top of testing.py, the following states are declared:

```python
run_txt_production_b   = False
run_txt_production_c_1 = False
run_txt_production_c_2 = False
run_txt_production_d   = False
run_txt_production_e   = False
```
To run one of the sections listed in the prompt, change the boolean state to:
```python
True
```
This will allow the algorithm to execute that particular code block and generate new results for analysis.
In solver.py on line 172 you will find the default values outputted to the .txt file. If you are creating new datasets for Part B or Part E, you need to comment out the other two blocks of .txt file writing and uncomment the correct block. ***This means that you cannot run new cases of c_1, c_2, and d at the same time as b or e.***