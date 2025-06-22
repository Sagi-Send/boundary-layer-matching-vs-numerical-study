# Boundary-Layer BVP — Uniform Approximation vs. Numerical Solution

This mini-repo demonstrates how a **matched-asymptotic (composite) uniform
approximation** to a boundary layer compares with a fully numerical
solution of the underlying singularly-perturbed ODE.

---

## Problem statement

We consider  

\[
\varepsilon\,y''(x)+\frac{1}{x}\,y'(x)+y(x)=0,
\qquad
x\in[\varepsilon,1],\qquad
y(\varepsilon)=0,\;
y(1)=e^{-1/2},
\]

where \(0<\varepsilon\ll1\) creates a thin boundary layer near
\(x=\varepsilon\).

* **Uniform approximation**  
  Built with matched-asymptotic expansions:  
  outer solution + inner layer → matched and combined.

* **Numerical reference**  
  Reformulated as a first-order system and solved with MATLAB `bvp4c`.

We sweep several values of ε and track the **L₁-error** between the two
solutions.

---

## Repository layout

| Path               | Purpose                                    |
|--------------------|--------------------------------------------|
| `compare_uniform_vs_numerical.m` | Main script—runs everything and saves the figures |
| `fig/`             | Output images are written here (auto-created) |
| `README.md`        | What you are reading |

---

## Quick start

```bash
git clone https://github.com/<your-user>/boundary-layer-bvp.git
cd boundary-layer-bvp
matlab -batch "compare_uniform_vs_numerical"