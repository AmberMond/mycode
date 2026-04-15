# Unstructured 2D Euler Bump Solver (Implicit FVM) — Theory \& Implementation Guide



## 1. Governing equations (2D Euler)

### 1.1 Conservative variables

$$
U = \begin{bmatrix}
\rho \\
\rho u \\
\rho v \\
\rho E
\end{bmatrix}
$$

### 1.2 Euler equations (conservative form)

$$
\frac{\partial U}{\partial t} + \nabla \cdot \mathbf{F}(U) = 0
$$

where $\mathbf{F}=(F, G)$ and

$$
F(U)=\begin{bmatrix}
\rho u\\
\rho u^2 + p\\
\rho u v\\
u(\rho E + p)
\end{bmatrix},\quad
G(U)=\begin{bmatrix}
\rho v\\
\rho u v\\
\rho v^2 + p\\
v(\rho E + p)
\end{bmatrix}
$$

### 1.3 Ideal gas closure

$$
p = (\gamma - 1)\left(\rho E - \frac{1}{2}\rho(u^2+v^2)\right)
$$

Sound speed:

$$
c = \sqrt{\gamma p/\rho}
$$

---

## 2. Unstructured FVM discretization (cell-centered)

### 2.1 Mesh entities

We store:

- **Nodes**: coordinates $(x,y)$
- **Cells** (triangles/quads): each cell has an area $V_i$ and centroid
- **Faces/edges**: each face has
    - `owner` cell id
    - `neighbor` cell id (or `-1` if boundary)
    - unit normal $\mathbf{n}_f$ pointing outward from owner
    - length $A_f$
    - boundary tag (inlet/outlet/top/bottom)


### 2.2 Semi-discrete finite volume form

Integrate Euler over a control volume (cell) $i$:

$$
\frac{d}{dt}(U_i V_i) + \sum_{f \in \partial i} \hat{F}_f A_f = 0
$$

Define residual:

$$
R_i(U) = \sum_{f \in \partial i} \hat{F}(U_L, U_R, \mathbf{n}_f)\,A_f
$$

Steady solution satisfies:

$$
R(U)=0
$$

---

## 3. Numerical flux

### 3.1 Normal flux for one state

For a face with unit normal $\mathbf{n}=(n_x,n_y)$, define

$$
u_n = u n_x + v n_y
$$

Normal physical flux:

$$
F_n(U)=\begin{bmatrix}
\rho u_n\\
\rho u u_n + p n_x\\
\rho v u_n + p n_y\\
(\rho E + p)u_n
\end{bmatrix}
$$

### 3.2 LLF

$$
\hat{F}=\frac{1}{2}\left(F_n(U_L)+F_n(U_R)\right)
-\frac{1}{2}\,a_{\max}(U_R-U_L)
$$

with

$$
a_{\max}=\max(|u_{n,L}|+c_L,\ |u_{n,R}|+c_R)
$$

---

## 4. Boundary conditions using ghost states (practical FV approach)

For a boundary face, set:

- $U_L$ = interior (owner cell) state
- Construct **ghost** $U_R$ based on boundary type
- Compute flux $\hat{F}(U_L,U_R)$


### 4.1 Slip wall (top and bottom)

Slip wall condition is **no penetration**:

$$
u_n = 0
$$

Ghost construction by reflecting the normal velocity:

1. From interior state, compute $(\rho,u,v,p)$
2. Compute $u_n = u n_x + v n_y$
3. Reflect velocity:

$$
\mathbf{u}^{ghost} = \mathbf{u} - 2 u_n \mathbf{n}
$$
4. Set $\rho^{ghost}=\rho$, $p^{ghost}=p$
5. Convert ghost primitives back to conservative to get $U_R$

This yields near-zero mass flux through the wall (up to numerical dissipation).

---

### 4.2 Subsonic outlet with specified static pressure `p_out`

A standard stable option:

- extrapolate $\rho,u,v$ from interior
- set $p = p_{out}$
- recompute energy:

$$
\rho E = \frac{p}{\gamma-1} + \frac{1}{2}\rho(u^2+v^2)
$$

Ghost state uses the updated $\rho E$.

---

### 4.3 Subsonic inlet with specified total `p0_in`, `T0_in`, given `M=0.5`, `AoA=0`

A first working inlet model: enforce a fixed inflow state derived from total conditions.

Isentropic relations:

$$
T = \frac{T_0}{1+\frac{\gamma-1}{2}M^2}
$$

$$
p = \frac{p_0}{\left(1+\frac{\gamma-1}{2}M^2\right)^{\gamma/(\gamma-1)}}
$$

$$
\rho = \frac{p}{R T}
$$

$$
a=\sqrt{\gamma R T},\quad V=M a
$$

AoA $=0\Rightarrow u=V,\ v=0$.

Then convert $(\rho,u,v,p)$ to conservative ghost state $U_R$.

> Note: a fully characteristic subsonic inlet can be added later; fixed inflow is fine for a first verification solver.

---

## 5. Implicit steady solver (pseudo-time + Newton)

### 5.1 Pseudo-time backward Euler to reach steady solution

Introduce pseudo-time $\tau$:

$$
\frac{V_i}{\Delta\tau}\left(U_i^{k+1}-U_i^{k}\right) + R_i(U^{k+1}) = 0
$$

Define:

$$
G(U) = \frac{V}{\Delta\tau}(U-U^{k}) + R(U)
$$

Solve:

$$
G(U)=0
$$

### 5.2 Newton linearization

At Newton iteration $m$:

$$
J(U^{(m)})\delta U = -G(U^{(m)})
$$

Update:

$$
U^{(m+1)} = U^{(m)} + \alpha \delta U
$$

Use damping $0<\alpha\le 1$ for robustness (line search).

---

## 6. Jacobian-Free Newton–Krylov (JFNK)

Analytic Jacobian assembly on unstructured meshes is complex. Instead:

- Use GMRES to solve the Newton system
- Provide only matrix-vector product $Jv$ using finite differences:

$$
Jv \approx \frac{G(U+\epsilon v)-G(U)}{\epsilon}
$$

Practical epsilon:

$$
\epsilon = 10^{-6}\frac{1+\|U\|}{\|v\|}
$$

This gives a true implicit solver without forming Jacobian matrices.

---

## 7. Linear solver: GMRES

GMRES solves:

$$
J\delta U = -G
$$

Matrix-free GMRES calls `matvec(v)` which returns the finite-difference approximation above.

Preconditioning is optional initially (small meshes) but strongly recommended later.

---

## 8. Local pseudo-time step (CFL-like) for faster convergence

Even for implicit, a local pseudo-time helps:

$$
\Delta\tau_i = \text{CFL}\,\frac{V_i}{\sum_{f\in i} (|u_n|+c)_f A_f}
$$

Typical starting CFL for implicit pseudo-time:

- 5–50 (increase gradually)
- use damping if Newton becomes unstable

---

## 9. Face-based residual assembly (conservative)

Algorithm to compute residual $R(U)$:

For each face `f`:

1. owner cell `i` → `UL = U[i]`
2. if interior face: neighbor `j` → `UR = U[j]`
3. if boundary face: construct ghost `UR` from BC type
4. compute flux `F = flux(UL, UR, n_f)`
5. accumulate:
    - `R[i] += F * A_f`
    - if neighbor exists: `R[j] -= F * A_f`

This ensures conservation and correct coupling.

---

## 10. Mesh processing 

Given cell connectivities (tri/quad):

1. Loop edges of each cell (node pairs)
2. Use a hash map with key `(min(n0,n1), max(n0,n1))`
3. If edge first seen: create face with `owner = cell`
4. If edge seen again: set `neighbor = cell`
5. Boundary faces have `neighbor = -1`
6. Compute geometry:
    - face endpoints $P_0,P_1$
    - length $A = \|P_1-P_0\|$
    - candidate normal from edge vector (rotate by 90°)
    - orient normal outward: check dot(normal, face_center - cell_center) > 0

Boundary tags come from mesh physical groups or entity tags.

---

## 11. Convergence monitoring and stopping criteria

Residual norms:

$$
\|R\|_2 = \sqrt{\sum_i \|R_i\|^2}
$$

Stop when:

- residual reduced by $10^{-6}$ to $10^{-8}$ relative
- and solution changes are small

Also check physical sanity:

- wall-normal velocity near zero on slip walls
- outlet pressure near `p_out`
- symmetry for AoA=0
