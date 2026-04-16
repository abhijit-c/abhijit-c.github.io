---
layout: post
title: "Sketching Regularized Least Squares: Primal ≠ Dual"
date: 2026-04-15 12:00:00
description: "Sketching the forward problem and converting via Woodbury is not the same as sketching the dual—here's a proof and the intuition behind it."
tags: numerical linear algebra, randomized algorithms, least squares, sketching
categories: math
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