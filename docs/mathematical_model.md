# Mathematical Model of 2D Two-Body Orbital Dynamics

This document details the mathematical formulation of the two-dimensional two-body gravitational model implemented in OrbitSim (Milestone M0).

---

## 1. Two-Body Problem

The classical two-body problem describes the gravitational interaction between two point masses. Under the standard astrodynamics assumptions for a lightweight satellite orbiting Earth:

1. The Earth is treated as an ideal central body located at the coordinate origin $(0, 0)$ and assumed stationary.
2. The satellite is modeled as a point mass without spatial extent or internal degrees of freedom.
3. Only the central Newtonian gravitational attraction of Earth is considered.
4. No atmospheric drag, solar radiation pressure, third-body gravity (Sun, Moon), Earth oblateness ($J_2$), or thrust forces are included at this stage.
5. The motion is constrained to a 2D orbital plane.

---

## 2. Coordinate System

We define a 2D Geocentric Cartesian coordinate system:

* **Origin $(0, 0)$**: Center of mass of the Earth.
* **$X$-axis**: In-plane reference axis (horizontal).
* **$Y$-axis**: Orthogonal in-plane axis (vertical).

---

## 3. State Vector

The instantaneous kinematic state of the satellite is represented by the 4-dimensional state vector $\mathbf{s}(t)$:

$$\mathbf{s}(t) = \begin{bmatrix} x(t) \\ y(t) \\ v_x(t) \\ v_y(t) \end{bmatrix}$$

where:
* $x, y$ are the position coordinates.
* $v_x = \dot{x} = \frac{dx}{dt}$ is the velocity component along the $X$-axis.
* $v_y = \dot{y} = \frac{dy}{dt}$ is the velocity component along the $Y$-axis.

---

## 4. Gravitational Parameter

The standard gravitational parameter of Earth is denoted by $\mu$ ($GM$):

$$\mu = 398600.435 \text{ km}^3/\text{s}^2$$

The equatorial radius of the Earth is:

$$R_E = 6378.137 \text{ km}$$

---

## 5. Distance Calculation

The scalar distance $r$ from the coordinate origin (Earth's center of mass) to the satellite is the Euclidean norm of the position vector $\mathbf{r} = [x, y]^T$:

$$r = \|\mathbf{r}\| = \sqrt{x^2 + y^2}$$

For physical validity, $r > 0$. At $r = 0$, a physical and mathematical singularity occurs.

---

## 6. Gravitational Acceleration

According to Newton's law of universal gravitation, the gravitational force acts towards the center of mass. The vector acceleration $\mathbf{a} = [a_x, a_y]^T$ is:

$$\mathbf{a} = -\frac{\mu}{r^3} \mathbf{r}$$

In Cartesian components:

$$a_x = -\mu \frac{x}{r^3}$$

$$a_y = -\mu \frac{y}{r^3}$$

Notice that the acceleration is strictly antiparallel to the position vector:

$$\mathbf{r} \times \mathbf{a} = x a_y - y a_x = x \left(-\mu \frac{y}{r^3}\right) - y \left(-\mu \frac{x}{r^3}\right) = 0$$

---

## 7. State Derivative

The system of first-order ordinary differential equations (ODEs) governing orbital motion is:

$$\frac{d\mathbf{s}}{dt} = \mathbf{f}(t, \mathbf{s}) = \begin{bmatrix} \dot{x} \\ \dot{y} \\ \dot{v}_x \\ \dot{v}_y \end{bmatrix} = \begin{bmatrix} v_x \\ v_y \\ a_x \\ a_y \end{bmatrix} = \begin{bmatrix} v_x \\ v_y \\ -\mu \frac{x}{r^3} \\ -\mu \frac{y}{r^3} \end{bmatrix}$$

---

## 8. Initial Conditions

For our initial baseline orbit, we simulate a low Earth circular orbit at an altitude of $400\text{ km}$:

$$r_0 = R_E + 400.0\text{ km} = 6378.137 + 400.0 = 6778.137\text{ km}$$

The initial state is positioned along the positive $X$-axis with velocity purely in the positive $Y$-direction:

$$\mathbf{s}_0 = \begin{bmatrix} x_0 \\ y_0 \\ v_{x0} \\ v_{y0} \end{bmatrix} = \begin{bmatrix} r_0 \\ 0 \\ 0 \\ v_{\text{circ}} \end{bmatrix} = \begin{bmatrix} 6778.137\text{ km} \\ 0\text{ km} \\ 0\text{ km/s} \\ \sqrt{\frac{\mu}{r_0}} \end{bmatrix}$$

---

## 9. Circular Orbit Velocity

A circular orbit requires the gravitational attraction to equal the required centripetal acceleration:

$$\frac{v_{\text{circ}}^2}{r} = \frac{\mu}{r^2}$$

Solving for speed yields:

$$v_{\text{circ}} = \sqrt{\frac{\mu}{r}}$$

For $r_0 = 6778.137\text{ km}$:

$$v_{\text{circ}} = \sqrt{\frac{398600.435}{6778.137}} \approx 7.67264\text{ km/s}$$

---

## 10. Units of Measurement

The simulation adheres to a consistent set of astrodynamic units:

| Quantity | Unit | Symbol |
| :--- | :--- | :--- |
| Length / Distance | Kilometer | $\text{km}$ |
| Time | Second | $\text{s}$ |
| Velocity | Kilometer per second | $\text{km/s}$ |
| Acceleration | Kilometer per second squared | $\text{km/s}^2$ |
| Gravitational Parameter $\mu$ | Cubic kilometers per second squared | $\text{km}^3/\text{s}^2$ |

