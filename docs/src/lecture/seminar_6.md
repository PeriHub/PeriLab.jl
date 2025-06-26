## Seminar 6: Parallelization and contact

- Concept of parallelization

- GPU based
    - avoided, because the whole paradigm is different


## CPU based

$$\begin{bmatrix}1 \\ 2 \\ 3 \\ 4 \\ 5 \\ 6 \\ 7 \end{bmatrix}_{all}\rightarrow\begin{matrix} \begin{bmatrix} 4 \\ 6  \end{bmatrix}_{C1-global}\\
\begin{bmatrix}1 \\ 2  \\ 7 \end{bmatrix}_{C2-global}\\
\begin{bmatrix}3\\5\end{bmatrix}_{C3-global}\end{matrix}\leftarrow \rightarrow\begin{matrix} \begin{bmatrix} 1 \\ 2  \end{bmatrix}_{C1-local}\\
\begin{bmatrix}1 \\ 2  \\ 3 \end{bmatrix}_{C2-local}\\
\begin{bmatrix}1\\2\end{bmatrix}_{C3-local}\end{matrix}$$

- restructuring from fields

!!! info "Data access"
    Not all data is avaiable at each core.


```sh
$ julia
julia> using MPI
julia> MPI.install_mpiexecjl()
```

Run PeriLab with two processors:
```sh
$ mpiexecjl -n 2 julia -e 'using PeriLab; PeriLab.main("examples/DCB/DCBmodel.yaml")'
```
Run PeriLab with one processors:
```sh
julia -e 'using PeriLab; PeriLab.main("examples/DCB/DCBmodel.yaml")'
```

## Multi-core search strategy
In Seminar
