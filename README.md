# PAT: finite element solver for the heat transfer analysis of infrastructures subjected to environmental actions

Noemi Schclar Leitão 
Laboratório Nacional de Engenharia Civil (LNEC), Av. do Brasil 101, 1700-066 Lisbon, Portugal
nschclar@lnec.pt

Eloísa Castilho 
Instituto Superior Técnico, Universidade de Lisboa, Av. Rovisco Pais, 1049-001 Lisbon, Portugal
eloisa.castilho@tecnico.ulisboa.pt

## How to use
This program is written in Fortran 95. It is composed of a main program, PAT.f95, and two routine modules, general_routines.f95 and environmental_routines.f95.
To compile them we have used gfortran version 8.1.0 on Windows using the command:
```
gfortran -o PAT.exe -Wextra -Wall PAT.f95 general_routines.f95 environmental_routines.f95
```

See Input_and_Examples.pdf for instructions on how to use the program and how to test the example data provided.
