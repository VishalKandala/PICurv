# Grid, Cell, and Variable Architecture Guide

## 1. Overview

This document provides a detailed explanation of the solver's grid and variable architecture. A clear understanding of this architecture is critical for correctly implementing new physics, boundary conditions, or post-processing routines.

The solver employs a structured curvilinear grid system managed by PETSc's `DMDA`. The core variables are defined on a staggered grid, with face-centered fluxes (`ucont`) and cell-centered quantities (`ucat`, `P`).

The most important architectural feature to understand is the **Shifted Index Architecture** for cell-centered data. In this design, the physical value for a geometric cell at index `i` is stored in the corresponding array at index `i+1`. This seemingly simple offset creates a clean, symmetric, and robust system for handling boundary conditions, resolving what might otherwise appear to be bugs or asymmetries in the code.

---

## 2. The Geometric Foundation: Nodes and Cells

The grid is defined by a set of nodes with physical coordinates. The space between nodes forms the computational cells (or control volumes).

*   **Nodes:** A grid with `IM` nodes in the i-direction is indexed from `0` to `IM-1`. Node coordinates are stored in the `coor` array, where `coor[i]` holds the `(x,y,z)` position of the i-th node.
*   **Cells:** `IM` nodes form `IM-1` physical cells, indexed `0` to `IM-2`. `Cell i` is the physical volume spanning between `Node i` and `Node i+1`.

**Example in 1D (i-direction):**

```
Nodes:     |-----------|-----------|-----------| ..... |-------------|
Index:     0           1           2           3       IM-2         IM-1

Cells:       [ Cell 0 ]  [ Cell 1 ]  [ Cell 2 ]          [ Cell IM-2 ]
Index:          0           1           2                   IM-2
```
*   **Geometric Reality:** The simulation space contains `(IM-1) x (JM-1) x (KM-1)` physical cells.

---

## 3. The Primary Variable: Face-Centered Flux (`ucont`)

The solver's fundamental fluid dynamics variable is `ucont`, the contravariant flux, which is physically located on the faces of the cells. Its mapping is direct and intuitive.

*   **Location:** `ucont.x` is on i-faces, `ucont.y` is on j-faces, etc.
*   **Indexing:** The flux through the face located at `Node i` is stored at index `i`.
    *   `ucont.x[0]`: Flux through the `-Xi` boundary face.
    *   `ucont.x[IM-1]`: Flux through the `+Xi` boundary face.
    *   The solver computes interior fluxes, while `FormBCS` sets boundary fluxes. This is a standard, symmetric approach.

**Example in 1D (i-direction):**
```
Faces:     |           |           |           | ..... |             |
ucont.x: ucont[0]   ucont[1]    ucont[2]    ucont[3]      ucont[IM-1]

Cells:     [ Cell 0 ]  [ Cell 1 ]  [ Cell 2 ]          [ Cell IM-2 ]
```

---

## 4. The Shifted Index Architecture for Cell-Centered Variables (`ucat`, `P`)

The Cartesian velocity (`ucat`) and pressure (`P`) are cell-centered variables. The key to this architecture is understanding their mapping:

> **The physical value for `Cell i` is stored in the array at index `i+1`.**

This means there is a one-to-one correspondence between the `IM-1` physical cells and `IM-1` physical `ucat` values, but with a simple offset. The first and last elements of the array (`ucat[0]` and `ucat[IM]`) are reserved for ghost values.

### 4.1. How it Works at the Boundaries

This shifted index system creates a clean and symmetric ghost-cell framework.

**At the "-Xi" Boundary (Min Boundary):**
*   The first physical cell is **`Cell 0`**. Its velocity is stored at **`ucat[1]`**.
*   The array element **`ucat[0]`** is now cleanly a **ghost value** for `Cell 0`.
*   The `FormBCS` function correctly sets this ghost value (e.g., `ucat[0] = -ucat[1]` for a wall) based on the first physical value.

**Diagram for the `-Xi` Boundary:**
```
           Face 0 (Wall)
              |
  [  Cell 0  ]-----[  Cell 1  ]-----[  Cell 2  ]
       ^                 ^                 ^
       |                 |                 |
   ucat[1]           ucat[2]           ucat[3]   <-- Physical Values
ucat[0]                                            <-- Ghost Value for Cell 0
```

**At the "+Xi" Boundary (Max Boundary):**
*   The last physical cell is **`Cell IM-2`**. Its velocity is stored at **`ucat[IM-1]`**.
*   The `Contra2Cart` kernel correctly computes up to `ucat[IM-1]`, which is the value for the last physical cell. There is no "leaked physics."
*   The array element **`ucat[IM]`** is now cleanly a **ghost value** for `Cell IM-2`.
*   `FormBCS` correctly sets this ghost value (e.g., `ucat[IM] = -ucat[IM-1]` for a wall) based on the last physical value.

**Diagram for the `+Xi` Boundary:**
```
...[ Cell IM-3 ]-----[ Cell IM-2 ]-----| Face IM-1 (Wall)
         ^                 ^           |
         |                 |           |
     ucat[IM-2]        ucat[IM-1]      | <-- Physical Values
                                   ucat[IM]  <-- Ghost Value for Cell IM-2
```

---

## 5. Effective Computational Domain and Resolution

This architecture means there is **no loss of resolution**. The effective computational domain for cell-centered variables matches the geometric domain.

*   **Geometric Cell Count (i-dir):** `IM - 1` (indices `0` to `IM-2`)
*   **Physical `ucat` Value Count (i-dir):** `IM - 1` (stored at indices `1` to `IM-1`)

For a grid defined with `IM x JM x KM` nodes, the resolution is:
*   **Geometric and Effective `ucat`/`P` Resolution:** `(IM-1) x (JM-1) x (KM-1)` cells.

---

## 6. Implications for Post-Processing (`ComputeNodalAverage`)

This understanding is critical for validating post-processing functions. The `ComputeNodalAverage` function interpolates cell values to the nodes.

*   **Stencil:** To compute the value at `Node i`, the function averages the 8 cells that meet at that node.
*   **At the `+Xi` Boundary:** To compute the value at `Node IM-1` (the physical boundary), the stencil correctly averages the surrounding cells. The key inputs in the i-direction are `Cell IM-2` and its ghost representation.
    *   The value for `Cell IM-2` is stored at `ucat[IM-1]`.
    *   The ghost value for `Cell IM-2` is stored at `ucat[IM]`.
*   **Validation:** The `ComputeNodalAverage` function, by sampling `ucat[IM-1]` and `ucat[IM]`, is correctly averaging the last physical cell with its proper ghost value. For a wall where `ucat[IM] = -ucat[IM-1]`, the average will be near zero, correctly representing the no-slip condition at the node. **The post-processing code is therefore correctly implemented for this architecture.**

---

## 7. Summary Table of `ucat` Anatomy (i-direction)

This table provides a quick reference for the mapping.

| Array Index `i` | Geometric Association | Role in Solver |
| :--- | :--- | :--- |
| `0` | Ghost value for `Cell 0` | **Boundary Condition Tool** |
| `1` | **Center of `Cell 0`** | **First True Physical Value** |
| `...` | ... | ... |
| `i+1` | **Center of `Cell i`** | Interior Physical Value |
| `...` | ... | ... |
| `IM-1` | **Center of `Cell IM-2`** | **Last True Physical Value** |
| `IM` | Ghost value for `Cell IM-2` | **Boundary Condition Tool** |