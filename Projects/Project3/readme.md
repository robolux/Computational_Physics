## How to run codebase

This program allows for the simulation of the solar system.


### To run selected test cases (as requested by prompt)
```
python testing.py
```

All data and plots will be available in /results

### Perihelion of Mercury Case

All test cases requested by the prompt are run with the exception of the perihelion of mercury case. This is due to the extremely long runtime of over 2.6 hours to complete. To get past the user having to run this part of the program each instance, there is a bool in testing.py called:

```
run_perihelion = False
```
 
By default this is False and instead of running the code section, it takes the results from a previous run found in a text file and reads them into the program and produces the results of the simulation as expected. It is highly advised to not change this bool but if you need to verify a change, just switching the state will allow for the computations to be completed.

```
run_perihelion = True
```
