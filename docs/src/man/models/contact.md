## Contact
## Contact Search

| Parameter | Unit | Description |
|---|---|---|
|Contact Radius | $[m]$| Radius to get a list of potential contact pairs.|
|Master| $[-]$| Block ID of the master nodes. |
|Slave | $[-]$| Block ID of the slave nodes.|
**Initialization**

- Step 1 -
Identification of all surface nodes of contact blocks.
- Step 2 -
All surface nodes are stored in a list. This list exists at all cores with the current position of the surface nodes.
- Step 3 -
Connect each surface node to a geometrical surface.

The standard way for points in the parallelization is given here.

$$\begin{bmatrix}1 \\ 2 \\ 3 \\ 4 \\ 5 \\ 6 \\ 7 \end{bmatrix}_{all}\rightarrow\begin{matrix} \begin{bmatrix} 4 \\ 6  \end{bmatrix}_{C1-global}\\
\begin{bmatrix}1 \\ 2  \\ 7 \end{bmatrix}_{C2-global}\\
\begin{bmatrix}3\\5\end{bmatrix}_{C3-global}\end{matrix}\leftarrow \rightarrow\begin{matrix} \begin{bmatrix} 1 \\ 2  \end{bmatrix}_{C1-local}\\
\begin{bmatrix}1 \\ 2  \\ 3 \end{bmatrix}_{C2-local}\\
\begin{bmatrix}1\\2\end{bmatrix}_{C3-local}\end{matrix}$$

A new mapping is needed to the contact surface exchange list. The local numbering (right above) has to be mapped to the sublist and vice, versa.

$$\begin{bmatrix}1 \\ 2   \\ 6 \\ 7 \end{bmatrix}_{surface}\leftarrow\rightarrow \begin{bmatrix}1 \\ 2  \\ 3 \\ 4  \end{bmatrix}_{surface}\rightarrow \begin{matrix} \begin{bmatrix}  6  \end{bmatrix}_{C1}\\
\begin{bmatrix}1 \\ 2  \\ 7 \end{bmatrix}_{C2}\\
\begin{bmatrix}\,\,\,\,\end{bmatrix}_{C3}\end{matrix}\leftarrow \rightarrow\begin{matrix} \begin{bmatrix} 2   \end{bmatrix}_{C1-local}\\
\begin{bmatrix}1 \\ 2  \\ 3 \end{bmatrix}_{C2-local}\\
\begin{bmatrix}\,\,\,\,\end{bmatrix}_{C3-local}\end{matrix}$$


**List mapping**
To fill the postion vector with the needed values mappings are needed.

local point ids -> global point id -> reduced contact block point id
The mapping has to be done in both directions, because the contact forces have to be applied to the local core point id.



---

**Computation**

- Step 1 -
Perform a nearest neighbor search with a user defined Contact Radius. The result is a list of potential contact pairs.
- Step 2 -
A sub list nearest neighbor search is performed.
- Step 3 -
All points found are than in contact.

The quasi code is given here

```julia
# synchronize all contact nodes
synchronize_contact_points("Deformed Coordinates", "NP1",
                                               all_positions)

## loop over all contact pairs
for (cm, block_contact_params) in pairs(contact_params)
    contact_search(...)
end

function contact_search(...)
    if n < search_frequency
        potential_contact_list = global_search(...) # finds all points within the search radius
        local_search(...)
        return
    else
        synchronize_contact_points("Deformed Coordinates", "NP1",
                                               all_positions[potential_contact_list])
        local_search(...)   # finds all points within the contact radius
        return
    end
end

```
The structure tries to minimize the search and the communication between cores.


**Synchronization**
All cores (except core 1) send it's displacement values to core 1. Core 1 sends them back.
!!! note "Efficiency"
    Smaller contact areas are more efficient in numerical analysis, because less synchronisation have to occur.

----

## Contact Models

| Contact Model           | Penalty Contact |
|------------------------|:---------------:|
| Contact Stiffness         | ✔️|

## Penalty Contact

A bond-based contact formulation is computed. The contact stiffness $ c_{contact}$ is defined in YAML input file. Based on that utilzing the horizon $\delta$ the penalty stiffness shown in the table is computed. These parameter are comparable to the [bond-based formulation](@ref "Bond-based Peridynamics"). The algorithm was taken from [Peridigm](https://github.com/peridigm/peridigm/blob/master/src/contact/Peridigm_ShortRangeForceContactModel.cpp)

| Dimension | Penalty Stiffness |
|---|---|
|plane strain:| $c_{Penalty} = \frac{48 c_{contact}}{\pi 5 \delta_{slave}^3}$ |
|plane stress:| $c_{Penalty} = \frac{9 c_{contact}}{\pi \delta_{slave}^3}$ |
|3D:| $c_{Penalty} = \frac{12 c_{contact}}{\pi  \delta_{slave}^4}$ |


### Normal Contact
The distance is currently computed as

$$d = |\underline{\mathbf{Y}}_{master}-\underline{\mathbf{Y}}_{slave}|$$

The normal is

$$\mathbf{n} = \frac{\underline{\mathbf{Y}}_{master}-\underline{\mathbf{Y}}_{slave}}{d}$$

Both can be adopted if needed. Utilizing the contact radius $r_{contact}$ the contact force densities are

$$
\mathbf{f}_{master} = c_{Penalty}\frac{(r_{contact}-d)}{\delta_{slave}}\mathbf{n}V_{slave}$$

$$\mathbf{f}_{slave} = -c_{Penalty}\frac{(r_{contact}-d)}{\delta_{slave}}\mathbf{n}V_{master}$$


Within PeriLab the functions compute_master_force_density and compute_slave_force_density are used to apply the force to the correct postions. This is needed, becaus if a parallel process is run, master and slave nodes not necessarily are on the same core.

### Friction Contact
| Contact Model           | Penalty Contact |
|------------------------|:---------------:|
| Friction Coeffient         | ✔️|

The model is computed if the ""Friction Coeffient"" is greater than zero.
For friction the normal forces
$$\mathbf{f}_n=c_{Penalty}\frac{(r_{contact}-d)}{\delta_{slave}}\mathbf{n}$$
are used. The tangential direction must be computed. Following the [Peridigm](https://github.com/peridigm/peridigm/blob/master/src/contact/Peridigm_ShortRangeForceContactModel.cpp) algorithm, the normal velocity of the master node $m$ and the slave node $s$ are

$$\begin{matrix}
v_{m} = \mathbf{v}_m\mathbf{n}\\
v_{s} = \mathbf{v}_s\mathbf{n}
\end{matrix}$$

To obtain the tangential part $tan$ this value hmust be substrated from the velocity.
$$\begin{matrix}
\mathbf{v}_{tan-m} = \mathbf{v}_m - v_{m}\mathbf{n}\\
\mathbf{v}_{tan-s} = \mathbf{v}_s - v_{s}\mathbf{n}
\end{matrix}$$

Tanking the average
$$\mathbf{v}_{avg} = 0.5(\mathbf{v}_{tan-m}+\mathbf{v}_{tan-s})$$
the relative velocity of each point can be determined.
$$\begin{matrix}
\mathbf{v}_{tan-m} = \mathbf{v}_{tan-m} - \mathbf{v}_{avg}\\
\mathbf{v}_{tan-s} = \mathbf{v}_{tan-s} - \mathbf{v}_{avg}
\end{matrix}$$
If the length of $\mathbf{v}_{tan-m}$ and $\mathbf{v}_{tan-s}$ is greater than zero, respectively, friction occurs. The friction force is computed as using the friction coefficient $\mu$.

$$\begin{matrix}
\mathbf{f}_{fric-m} = -\mu |\mathbf{f}_n|\frac{\mathbf{v}_{tan-m}}{|\mathbf{v}_{tan-m}|}\\
\mathbf{f}_{fric-s} = -\mu |\mathbf{f}_n|\frac{\mathbf{v}_{tan-s}}{|\mathbf{v}_{tan-s}|}
\end{matrix}$$
