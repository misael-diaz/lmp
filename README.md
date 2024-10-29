# lmp
repository to learn about lammps

# compile

```sh
make clean && make
```

# run

LAMMPS reads the input data from standard input and so you will need to use redirection
if you wish to supply a file instead (recommended):

```sh
./src/lmp < in.data
```
