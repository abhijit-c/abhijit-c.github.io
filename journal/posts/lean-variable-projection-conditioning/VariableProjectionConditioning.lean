import Mathlib.LinearAlgebra.Matrix.Hermitian
import Mathlib.Tactic

/-!
# Variable projection does not worsen quadratic conditioning

For a positive-definite block Hessian

    H = [A  B]
        [Bᵀ D],

eliminating the first variable gives the Schur complement

    S = D - Bᵀ A⁻¹ B.

This file separates the proof into two layers.

1. A norm-agnostic variational theorem: if `xStar y` minimizes the full
   quadratic at fixed `y`, every lower and upper curvature bound for the full
   problem is inherited by the reduced problem.
2. The matrix identity in mathlib that identifies the reduced quadratic with
   the Schur-complement quadratic and proves the minimizing property.
-/

open Matrix

namespace VariableProjectionConditioning

/-! ## The variational conditioning argument -/

section Variational

variable {X Y : Type*}

/-- The reduced quadratic obtained by eliminating `x`. -/
def reducedQuadratic (fullQ : X → Y → ℝ) (xStar : Y → X) (y : Y) : ℝ :=
  fullQ (xStar y) y

/-- A lower curvature bound for a quadratic on the product space. -/
def JointLowerBound
    (fullQ : X → Y → ℝ) (sqX : X → ℝ) (sqY : Y → ℝ) (μ : ℝ) : Prop :=
  ∀ x y, μ * (sqX x + sqY y) ≤ fullQ x y

/-- An upper curvature bound for a quadratic on the product space. -/
def JointUpperBound
    (fullQ : X → Y → ℝ) (sqX : X → ℝ) (sqY : Y → ℝ) (L : ℝ) : Prop :=
  ∀ x y, fullQ x y ≤ L * (sqX x + sqY y)

/-- A lower curvature bound for a quadratic on the retained variable. -/
def LowerBound (q : Y → ℝ) (sqY : Y → ℝ) (μ : ℝ) : Prop :=
  ∀ y, μ * sqY y ≤ q y

/-- An upper curvature bound for a quadratic on the retained variable. -/
def UpperBound (q : Y → ℝ) (sqY : Y → ℝ) (L : ℝ) : Prop :=
  ∀ y, q y ≤ L * sqY y

/-- `μ` is the largest lower curvature bound of `q`. -/
def IsBestLowerBound (q : Y → ℝ) (sqY : Y → ℝ) (μ : ℝ) : Prop :=
  LowerBound q sqY μ ∧ ∀ c, LowerBound q sqY c → c ≤ μ

/-- `L` is the smallest upper curvature bound of `q`. -/
def IsBestUpperBound (q : Y → ℝ) (sqY : Y → ℝ) (L : ℝ) : Prop :=
  UpperBound q sqY L ∧ ∀ c, UpperBound q sqY c → L ≤ c

/-- The spectral condition number represented by sharp curvature bounds. -/
noncomputable def conditionNumber (μ L : ℝ) : ℝ :=
  L / μ

/--
The lower curvature of the full quadratic survives elimination.

The key point is that the norm of `(xStar y, y)` is at least the norm of `y`.
-/
theorem reduced_lower_bound
    (fullQ : X → Y → ℝ) (xStar : Y → X)
    (sqX : X → ℝ) (sqY : Y → ℝ) {μ : ℝ}
    (hμ : 0 ≤ μ) (hsqX : ∀ x, 0 ≤ sqX x)
    (hLower : JointLowerBound fullQ sqX sqY μ) :
    LowerBound (reducedQuadratic fullQ xStar) sqY μ := by
  intro y
  calc
    μ * sqY y ≤ μ * (sqX (xStar y) + sqY y) := by
      nlinarith [hsqX (xStar y)]
    _ ≤ fullQ (xStar y) y := hLower (xStar y) y
    _ = reducedQuadratic fullQ xStar y := rfl

/--
The upper curvature of the full quadratic survives elimination.

Here we use exact minimization: compare `xStar y` with the admissible point
`x = xZero`.
-/
theorem reduced_upper_bound
    (fullQ : X → Y → ℝ) (xStar : Y → X) (xZero : X)
    (sqX : X → ℝ) (sqY : Y → ℝ) {L : ℝ}
    (hsqXZero : sqX xZero = 0)
    (hUpper : JointUpperBound fullQ sqX sqY L)
    (hMin : ∀ y x, fullQ (xStar y) y ≤ fullQ x y) :
    UpperBound (reducedQuadratic fullQ xStar) sqY L := by
  intro y
  calc
    reducedQuadratic fullQ xStar y
        = fullQ (xStar y) y := rfl
    _ ≤ fullQ xZero y := hMin y xZero
    _ ≤ L * (sqX xZero + sqY y) := hUpper xZero y
    _ = L * sqY y := by rw [hsqXZero, zero_add]

/--
If the lower endpoint moves up and the upper endpoint moves down, their ratio
cannot increase.
-/
theorem conditionNumber_mono
    {μFull μReduced LReduced LFull : ℝ}
    (hμFull : 0 < μFull)
    (hμ : μFull ≤ μReduced)
    (hLReduced : 0 ≤ LReduced)
    (hL : LReduced ≤ LFull) :
    conditionNumber μReduced LReduced ≤ conditionNumber μFull LFull := by
  have hμReduced : 0 < μReduced := lt_of_lt_of_le hμFull hμ
  have h₁ : LReduced * μFull ≤ LReduced * μReduced :=
    mul_le_mul_of_nonneg_left hμ hLReduced
  have h₂ : LReduced * μReduced ≤ LFull * μReduced :=
    mul_le_mul_of_nonneg_right hL hμReduced.le
  unfold conditionNumber
  exact (div_le_div_iff₀ hμReduced hμFull).2 (h₁.trans h₂)

/--
The complete abstract theorem.

For Euclidean Hessians, the best lower and upper curvature bounds are the
smallest and largest eigenvalues. Thus the conclusion is exactly

    κ(reduced Hessian) ≤ κ(full Hessian).
-/
theorem variable_projection_conditioning
    (fullQ : X → Y → ℝ) (xStar : Y → X) (xZero : X)
    (sqX : X → ℝ) (sqY : Y → ℝ)
    {μFull LFull μReduced LReduced : ℝ}
    (hμFull : 0 < μFull)
    (hLReduced : 0 ≤ LReduced)
    (hsqX : ∀ x, 0 ≤ sqX x)
    (hsqXZero : sqX xZero = 0)
    (hLower : JointLowerBound fullQ sqX sqY μFull)
    (hUpper : JointUpperBound fullQ sqX sqY LFull)
    (hMin : ∀ y x, fullQ (xStar y) y ≤ fullQ x y)
    (hBestLower :
      IsBestLowerBound (reducedQuadratic fullQ xStar) sqY μReduced)
    (hBestUpper :
      IsBestUpperBound (reducedQuadratic fullQ xStar) sqY LReduced) :
    conditionNumber μReduced LReduced ≤ conditionNumber μFull LFull := by
  have hReducedLower :
      LowerBound (reducedQuadratic fullQ xStar) sqY μFull :=
    reduced_lower_bound fullQ xStar sqX sqY hμFull.le hsqX hLower
  have hReducedUpper :
      UpperBound (reducedQuadratic fullQ xStar) sqY LFull :=
    reduced_upper_bound fullQ xStar xZero sqX sqY hsqXZero hUpper hMin
  have hμ : μFull ≤ μReduced := hBestLower.2 μFull hReducedLower
  have hL : LReduced ≤ LFull := hBestUpper.2 LFull hReducedUpper
  exact conditionNumber_mono hμFull hμ hLReduced hL

end Variational

/-! ## The Schur-complement bridge -/

section SchurComplement

variable {m n : Type*}
variable [Fintype m] [DecidableEq m] [Fintype n]

/-- The full block quadratic `[x; y]ᵀ [A B; Bᵀ D] [x; y]`. -/
noncomputable def fullBlockQuadratic
    (A : Matrix m m ℝ) (B : Matrix m n ℝ) (D : Matrix n n ℝ)
    (x : m → ℝ) (y : n → ℝ) : ℝ :=
  (star (Sum.elim x y) ᵥ*
      Matrix.fromBlocks A B B.conjTranspose D) ⬝ᵥ
    Sum.elim x y

/-- The quadratic form of the Schur complement `D - Bᵀ A⁻¹ B`. -/
noncomputable def schurQuadratic
    (A : Matrix m m ℝ) (B : Matrix m n ℝ) (D : Matrix n n ℝ)
    [Invertible A] (y : n → ℝ) : ℝ :=
  (star y ᵥ* (D - B.conjTranspose * A⁻¹ * B)) ⬝ᵥ y

/-- The eliminated variable `x*(y) = -A⁻¹ B y`. -/
noncomputable def xStar
    (A : Matrix m m ℝ) (B : Matrix m n ℝ) [Invertible A]
    (y : n → ℝ) : m → ℝ :=
  -((A⁻¹ * B).mulVec y)

/--
Mathlib's completing-the-square identity for a Hermitian block matrix.
-/
theorem schur_complement_identity
    (A : Matrix m m ℝ) (B : Matrix m n ℝ) (D : Matrix n n ℝ)
    [Invertible A] (hA : A.IsHermitian)
    (x : m → ℝ) (y : n → ℝ) :
    fullBlockQuadratic A B D x y =
      (star (x + (A⁻¹ * B).mulVec y) ᵥ* A) ⬝ᵥ
          (x + (A⁻¹ * B).mulVec y) +
        schurQuadratic A B D y := by
  exact Matrix.schur_complement_eq₁₁ B D x y hA

/-- At `x*(y)`, the square term vanishes and only the Schur complement remains. -/
theorem fullBlockQuadratic_at_xStar
    (A : Matrix m m ℝ) (B : Matrix m n ℝ) (D : Matrix n n ℝ)
    [Invertible A] (hA : A.IsHermitian) (y : n → ℝ) :
    fullBlockQuadratic A B D (xStar A B y) y =
      schurQuadratic A B D y := by
  rw [schur_complement_identity A B D hA]
  simp [xStar]

/--
If the leading-block quadratic is nonnegative, `x*(y)` minimizes the full
quadratic over `x`.
-/
theorem xStar_minimizes
    (A : Matrix m m ℝ) (B : Matrix m n ℝ) (D : Matrix n n ℝ)
    [Invertible A] (hA : A.IsHermitian)
    (hAPos : ∀ u : m → ℝ, 0 ≤ (star u ᵥ* A) ⬝ᵥ u) :
    ∀ y x,
      fullBlockQuadratic A B D (xStar A B y) y ≤
        fullBlockQuadratic A B D x y := by
  intro y x
  rw [fullBlockQuadratic_at_xStar A B D hA y]
  rw [schur_complement_identity A B D hA x y]
  exact le_add_of_nonneg_left (hAPos _)

/--
The matrix-specialized conclusion.

When `sqX` and `sqY` are squared Euclidean norms, the sharp constants in
`hBestLower` and `hBestUpper` are the extreme eigenvalues of the Schur
complement. If `μFull` and `LFull` are the extreme eigenvalues of the full
block matrix, the conclusion is the usual spectral condition-number bound.
-/
theorem schur_complement_conditioning
    (A : Matrix m m ℝ) (B : Matrix m n ℝ) (D : Matrix n n ℝ)
    [Invertible A] (hA : A.IsHermitian)
    (hAPos : ∀ u : m → ℝ, 0 ≤ (star u ᵥ* A) ⬝ᵥ u)
    (sqX : (m → ℝ) → ℝ) (sqY : (n → ℝ) → ℝ)
    {μFull LFull μReduced LReduced : ℝ}
    (hμFull : 0 < μFull)
    (hLReduced : 0 ≤ LReduced)
    (hsqX : ∀ x, 0 ≤ sqX x)
    (hsqXZero : sqX 0 = 0)
    (hLower :
      JointLowerBound (fullBlockQuadratic A B D) sqX sqY μFull)
    (hUpper :
      JointUpperBound (fullBlockQuadratic A B D) sqX sqY LFull)
    (hBestLower :
      IsBestLowerBound (schurQuadratic A B D) sqY μReduced)
    (hBestUpper :
      IsBestUpperBound (schurQuadratic A B D) sqY LReduced) :
    conditionNumber μReduced LReduced ≤ conditionNumber μFull LFull := by
  have hReducedEq :
      reducedQuadratic (fullBlockQuadratic A B D) (xStar A B) =
        schurQuadratic A B D := by
    funext y
    exact fullBlockQuadratic_at_xStar A B D hA y
  rw [← hReducedEq] at hBestLower hBestUpper
  exact variable_projection_conditioning
    (fullBlockQuadratic A B D) (xStar A B) (0 : m → ℝ)
    sqX sqY hμFull hLReduced hsqX hsqXZero hLower hUpper
    (xStar_minimizes A B D hA hAPos) hBestLower hBestUpper

end SchurComplement

end VariableProjectionConditioning
