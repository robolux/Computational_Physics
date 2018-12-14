## How to run codebase

This program allows for the simulation of financial transactions among financial agents using Monte Carlo methods.

### User Input for m0, N (number of agents), transaction count, simulation count, lambda, alpha, and gamma  to solve their choice of system pertaining to this situation.
```python
python main.py
```

### To run selected test cases (as requested by prompt)
All of the actual computation results are included in the Github Repo and by default the script reads previously computed simulations from /results .txt files. It is highly reccomended that you do not rerun the unit cases due to computational time.

```python
python testing.py
```
All .txt files will be in /results

### To generate plots from testing.py results
```
python plotting.py
```

All plots will be in /results