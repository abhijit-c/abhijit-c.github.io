---
layout: post
title: "Sketching Regularized Least Squares: Primal ≠ Dual"
date: 2026-04-15 12:00:00
description: "Sketching the forward problem and converting via Woodbury is not the same as sketching the dual—here's a proof and the intuition behind it."
tags: numerical linear algebra, randomized algorithms, least squares, sketching
categories: math
published: false
giscus_comments: true
toc:
  sidebar: left
---

Here's something surprising I encountered recently.
Consider the standard Tikhonov-regularized least squares problem

$$
\min_{x \in \mathbb{R}^n} 
\frac{1}{2}\|Ax - b\|_2^2 + \frac{\lambda}{2} \|x\|_2^2, 
\qquad A \in \mathbb{R}^{m \times n},\ b \in \mathbb{R}^m,\ \lambda > 0. \tag{P}
$$

When $m$ (or $n$) is relatively large, your matrix naturally contains many redundancies.
Hence, a natural idea is to *sketch* the problem. 
Apply a matrix $S \in \mathbb{R}^{k \times m}$ with $k \ll m$ to compress the
$m$ rows down to $k$ and solve the cheaper surrogate.

Moreover, from classic optimization, the original system can be viewed in both
its primal and dual formulation.
We're considering the over-determined case, but in the under-determined setting
that can computationally be the right perspective.
Because our system is convex, the duality gap is zero; that is, they should
give me the same solution.
A natural question, therefore, is if that's preserved under sketching. 
That is, do the following agree:

1. **Sketch the primal** directly, then rewrite via the Woodbury identity.
2. **Sketch the dual** system directly.

They don't. That was super weird to me!

---

## Primal and dual forms

The primal is what we saw

$$\min_{x,\, r} \frac{1}{2}\|r\|_2^2 + \frac{\lambda}{2}\|x\|_2^2 \quad \text{subject to} \quad Ax - r = b, \tag{P}$$

where we have introduced the residual $r = Ax - b$ as a free variable to make
the constraint explicit. 
The normal equations give the familiar solution

$$x^* = (A^T A + \lambda I_n)^{-1} A^T b.$$

We can obtain the dual via the Lagrangian with multiplier $\alpha \in
\mathbb{R}^m$:

$$\mathcal{L}(x, r, \alpha) = \frac{1}{2}\|r\|_2^2 + \frac{\lambda}{2}\|x\|_2^2 + \alpha^T(b - Ax + r).$$

Differentiating for $(x, r)$:

$$\nabla_x \mathcal{L} = \lambda x - A^T\alpha = 0 \implies x = \tfrac{1}{\lambda} A^T\alpha,$$

$$\nabla_r \mathcal{L} = r + \alpha = 0 \implies r = -\alpha.$$

Substituting back, the dual function is

$$
\begin{aligned}
g(\alpha) 
&= \frac{1}{2}\|\alpha\|_2^2 + \frac{1}{2\lambda}\|A^T\alpha\|_2^2 + \alpha^T b - \frac{1}{\lambda}\alpha^T A A^T \alpha - \|\alpha\|_2^2  \\
&= \alpha^T b - \frac{1}{2\lambda}\alpha^T(A A^T + \lambda I_m)\alpha.
\end{aligned}
$$

Since $(\mathrm{P})$ is convex and the equality constraint is affine, strong duality holds. The dual problem itself is

$$\max_{\alpha \in \mathbb{R}^m}\; \alpha^T b - \frac{1}{2\lambda}\,\alpha^T(A A^T + \lambda I_m)\,\alpha, \tag{D}$$

and the primal solution is recovered via $x^* = \frac{1}{\lambda}A^T\alpha^*$. Setting $\nabla g = 0$ gives

$$(A A^T + \lambda I_m)\,\alpha^* = \lambda b, \qquad x^* = A^T(A A^T + \lambda I_m)^{-1} b.$$

---

## Sketching the primal

Apply $S \in \mathbb{R}^{k\times m}$ to the data in $(\mathrm{P})$:

$$\min_x \frac{1}{2}\|SAx - Sb\|_2^2 + \frac{\lambda}{2}\|x\|_2^2.$$

The normal equations for this sketched problem are $(A^T S^T S A + \lambda I_n)\tilde{x} = A^T S^T S b$, giving

$$\tilde{x}_P = (A^T S^T S A + \lambda I_n)^{-1} A^T S^T S b.$$

Applying the Woodbury identity $(A^T S^T S A + \lambda I_n)^{-1} A^T S^T = A^T S^T (S A A^T S^T + \lambda I_k)^{-1}$ to this solution,

$$\boxed{\tilde{x}_P = A^T S^T \bigl(S A A^T S^T + \lambda I_k\bigr)^{-1} S b.}$$

The regularizer $\frac{\lambda}{2}\|x\|_2^2$ was never touched by $S$, so the identity in the $k\times k$ inverse is exactly $\lambda I_k$.

---

## Sketching the dual

Restrict the dual variable to the row-space of $S$ by setting $\alpha = S^T\beta$ for $\beta\in\mathbb{R}^k$ and substituting into $(\mathrm{D})$:

$$\max_{\beta \in \mathbb{R}^k}\; \beta^T S b - \frac{1}{2\lambda}\,\beta^T S(A A^T + \lambda I_m) S^T \beta.$$

Expanding, this is a quadratic in $\beta$ with coefficient matrix $S(AA^T + \lambda I_m)S^T = SAA^TS^T + \lambda SS^T$. Setting the gradient to zero,

$$Sb = \frac{1}{\lambda}(S A A^T S^T + \lambda S S^T)\,\beta^*,$$

so $\beta^* = \lambda(S A A^T S^T + \lambda S S^T)^{-1} S b$ and $\tilde{\alpha} = S^T \beta^*$. Recovering the primal variable via $\tilde{x} = \frac{1}{\lambda}A^T\tilde{\alpha}$,

$$\boxed{\tilde{x}_D = A^T S^T \bigl(S A A^T S^T + \lambda S S^T\bigr)^{-1} S b.}$$

Now the $k\times k$ inverse contains $\lambda SS^T$, not $\lambda I_k$.

---

## Why do they differ?

These are only equal if $S^T S = I_k$.
That might be true in expectation (for certain clases of sketches), but that's not true in general!

As far as I understand, the intuition here is that the regularizer
$\lambda\|x\|_2^2$ lives in primal (parameter) space, not data space, and the
two sketching strategies treat it differently.

Sketching the primal applies $S$ only to the *data*,the rows of $A$ and the entries of $b$. 
It never touches $\frac{\lambda}{2}\|x\|_2^2$. 
The regularizer remains perfectly isotropic after the sketch. 
When we then use Woodbury to rewrite the solution in terms of a $k\times k$
system, the regularization shows up as $\lambda I_k$: a scalar multiple of the
identity in the compressed space, penalizing all directions equally.

Sketching the dual applies $S$ to the entire $m\times m$ system $(\mathrm{D})$, including the $\lambda I_m$ term that encodes the primal prior. 
Under the substitution $\alpha = S^T\beta$, that term becomes $\lambda SS^T$ in the $k\times k$ system. 
Unless $SS^T = I_k$, this compressed regularizer is *anisotropic*: some directions in $\mathbb{R}^k$ are penalized more than others.

Put differently: the Woodbury identity is an exact algebraic fact about the *full* problem. 
Sketching is not an algebraic operation, it changes the optimization problem itself. 
The sketched primal and the sketched dual are two genuinely different approximate problems, and their solutions need not agree.

---

## Which is better?

This is context-dependent, and the answer turns on what sketch you're using and which regime you're in.

### The primal sketch is the cleaner approximation

From a numerical standpoint, $\tilde{x}_P$ has a decisive structural advantage: it solves a genuine regularized least squares problem in its own right,

$$\tilde{x}_P = \arg\min_x \tfrac{1}{2}\|SAx - Sb\|_2^2 + \tfrac{\lambda}{2}\|x\|_2^2.$$

The regularizer is never touched by $S$, so the $\lambda I_k$ in the Woodbury-rewritten inverse is genuinely isotropic — every direction in the compressed space is penalized equally, exactly as in the original problem.
Error analysis is correspondingly clean: approximation quality depends on how faithfully $S^TS$ preserves the spectrum of $A$ restricted to its row space, and for standard sketching families (Gaussian, subsampled randomized Hadamard, sparse sign) this is a textbook guarantee that holds with high probability once $k$ is large enough relative to $\text{rank}(A)$ and $1/\lambda$.

### The dual sketch is the Nyström approximation

The dual sketch is not just an alternative formula — it corresponds to a well-known and heavily studied approximation in machine learning. 
If $S$ is a *column-selection* sketch (i.e., $S$ picks $k$ rows out of the $m$-row identity), then $SAA^TS^T$ is a $k\times k$ principal submatrix of the Gram matrix $K = AA^T$, and $\tilde{x}_D$ is exactly the **Nyström approximation to kernel ridge regression** (KRR).
In that literature the $\lambda SS^T$ term is a feature: it reflects how the regularization interacts with the geometry of the selected subspace.
When $S$ is chosen to capture most of the spectral mass of $K$ (e.g., by selecting the columns corresponding to the $k$ largest diagonal entries or via leverage-score sampling), the Nyström approximation enjoys sharp bounds that decay with the tail of the spectrum of $K$. That's a stronger, problem-adapted guarantee than what generic sketching gives you.

### Summary

| Regime | Prefer |
|---|---|
| $m \gg n$ (tall, overdetermined) | **Primal** — cheap, clean, textbook guarantees |
| $m \approx n$ or $m \ll n$ (wide/kernel regime) | **Dual** — Nyström-KRR literature applies |
| Generic unstructured sketch ($S$ Gaussian, SRHT, …) | **Primal** — isotropic $\lambda I_k$ is well-behaved |
| Column-selection or leverage-score sketch | **Dual** — designed for spectral approximation of $K$ |
| $SS^T = I_k$ (orthogonal rows) | **Either** — the two solutions coincide |

The practical upshot: if you want a general-purpose, principled approximation and do not have strong prior knowledge about the spectrum of $AA^T$, use the primal sketch.
Its approximation error is easy to control and the regularizer remains well-behaved in the compressed space.

If you are explicitly working in the kernel or wide-matrix regime and can afford to design $S$ to capture the leading spectral structure of $K = AA^T$, use the dual sketch.
The anisotropy of $\lambda SS^T$ is then not a defect but an adaptation to the geometry of the problem, and the Nyström KRR literature gives you sharp, spectral-tail-dependent error bounds rather than generic oblivious-sketch guarantees.