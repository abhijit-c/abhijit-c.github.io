---
layout: post
title: "Sketching Tikhonov in the Underdetermined Case"
date: 2026-04-15 12:00:00
description: "The usual left-sketching story is built for overdetermined problems. For underdetermined Tikhonov problems, dualizing makes the matrix tall again, but it is not obvious what the right sketch should be."
tags: numerical linear algebra, randomized algorithms, least squares, sketching
categories: math
published: False
giscus_comments: true
toc:
  sidebar: left
---

Consider the Tikhonov-regularized least squares problem

$$
\min_{x \in \mathbb{R}^n}
\frac{1}{2}\|Ax - b\|_2^2 + \frac{\lambda}{2}\|x\|_2^2,
\qquad A \in \mathbb{R}^{m \times n},\ b \in \mathbb{R}^m,\ \lambda > 0,
\tag{P}
$$

and suppose that the system is **underdetermined**, so $m < n$.

The regularization makes the optimization problem perfectly well posed. 
But from the point of view of sketching, this regime is kind of awkward.
Most sketching theory is formulated in an overdetermined setting, i.e. $m \gg
n$ and $A$ is tall.
Then one introduces a sketch

$$
S \in \mathbb{R}^{k \times m},
\qquad k \ll m,
$$

and replaces the original problem by

$$
\min_x \frac{1}{2}\|S(Ax-b)\|_2^2 + \frac{\lambda}{2}\|x\|_2^2.
\tag{S}
$$

This is natural as the sketch compresses the large dimension: the number of
rows.
But, in the underdetermined case $m < n$, the rows are not the problematic
dimension.
So the standard recipe does not really line up with the geometry of the regime
I am interested in.

If the primal matrix $A$ is wide, then its transpose $A^T$ is tall.
So a natural idea is to simply convert to the dual problem, which is formulated
in terms of $A^T$, and simply sketch there.
At least at the level of dimensions, this is appealing. 

---

## Deriving the dual via the Lagrangian

Introduce the residual explicitly:

$$
\min_{x,r} \frac{1}{2}\|r\|_2^2 + \frac{\lambda}{2}\|x\|_2^2
\quad \text{subject to} \quad Ax - r = b.
$$

Let $\alpha \in \mathbb{R}^m$ be the Lagrange multiplier for the constraint.
The Lagrangian is

$$
\mathcal{L}(x,r,\alpha)
= \frac{1}{2}\|r\|_2^2 + \frac{\lambda}{2}\|x\|_2^2
+ \alpha^T(b - Ax + r).
$$

The stationarity conditions are

$$
\nabla_x \mathcal{L} = \lambda x - A^T\alpha = 0
\qquad \Longrightarrow \qquad
x = \frac{1}{\lambda}A^T\alpha,
$$

and

$$
\nabla_r \mathcal{L} = r + \alpha = 0
\qquad \Longrightarrow \qquad
r = -\alpha.
$$

Substituting these back into the Lagrangian gives the dual objective

$$
g(\alpha)
= \alpha^T b - \frac{1}{2\lambda}\|A^T\alpha\|_2^2 - \frac{1}{2}\|\alpha\|_2^2.
$$

So the dual problem is

$$
\max_{\alpha \in \mathbb{R}^m}
\alpha^T b - \frac{1}{2\lambda}\|A^T\alpha\|_2^2 - \frac{1}{2}\|\alpha\|_2^2.
\tag{D}
$$

Equivalently,

$$
\max_{\alpha \in \mathbb{R}^m}
\alpha^T b - \frac{1}{2\lambda}\alpha^T(AA^T + \lambda I_m)\alpha.
$$

The first-order condition is

$$
(AA^T + \lambda I_m)\alpha^* = \lambda b,
$$

and the primal solution is recovered from

$$
x^* = \frac{1}{\lambda}A^T\alpha^*.
$$

This is the point I wanted to make explicit. In the dual, the quadratic term is
written through $A^T\alpha$, and $A^T \in \mathbb{R}^{n \times m}$ is now a
tall matrix when $m < n$.

So at least formally, dualizing converts the underdetermined primal problem
into an overdetermined-looking dual one.

---

## But what exactly should be sketched?

This is where the issue becomes less clear.

Once the dual is written as

$$
\max_\alpha \;
\alpha^T b - \frac{1}{2\lambda}\|A^T\alpha\|_2^2 - \frac{1}{2}\|\alpha\|_2^2,
$$

it is tempting to say: now sketch the tall matrix $A^T$ from the left.
For example, with a sketch

$$
S \in \mathbb{R}^{k \times n},
\qquad k \ll n,
$$

one might try to replace $A^T$ by $SA^T$.

But that is exactly the point where I am no longer sure what the most natural
formulation is.

What is the right sketched dual model?

- Should one simply replace $\|A^T\alpha\|_2^2$ by $\|SA^T\alpha\|_2^2$ and
  leave the $\|\alpha\|_2^2$ term untouched?
- Should one instead derive a sketched dual from some sketched constrained
  formulation?
- If one sketches in the dual and then maps back to the primal, what primal
  optimization problem has actually been approximated?

For overdetermined least squares, sketching the primal from the left is the
canonical move. For underdetermined least squares, dualizing first is an
obvious idea because it turns the wide matrix $A$ into the tall matrix $A^T$.
But after doing that, the notion of the "right" sketch is not automatic in the
same way.

That is the discussion I want to pin down: in the underdetermined case, what is
the most natural way to formulate sketching at all?
