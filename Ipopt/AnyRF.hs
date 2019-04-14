{-# LANGUAGE ConstraintKinds, RankNTypes, TypeFamilies, FlexibleInstances,
      UndecidableInstances #-}
{- | Description: representing functions that can be differentiated

The 'AnyRF' wrapper holds functions that can be used
for the objective (`f`) or for constraints (`g`). Many functions
in the instances provided are partial: this seems to be unavoidable
because the input variables haven't been decided yet, so you should
not be allowed to use 'compare' on these. But for now just use the
standard Prelude classes, and unimplementable functions (which
would not produce an 'AnyRF') are calls to 'error'.

Values of type @AnyRF Identity@ can be generated using functions 
defined in "Ipopt.NLP" (also exported by "Ipopt"). Directly using the
constructor is another option: @AnyRF $ Identity . V.sum@, calculates
the sum of all variables in the problem.

convergence can sometimes be improved by exposing additional
variables/derivatives to the solver.  IE. instead of maximizing f(u) where
we internally calculate x(u), maximize f(u,y) with another constraint that
y = x(u). Albersmeyer 2010 SIAM. Also generally known in process optimization
(ie. "equation-oriented" which has many variables works, while "sequential-modular"
modes that have very few variables left does not work)

AnyRF could be used to generate some C and feed it to something like casadi,
though this has to happen a runtime unless there's a clean way to handle constant
parameters such as (fromIntegral 3)

-}
module Ipopt.AnyRF where

import Data.Sequence (Seq)
import Data.Vector (Vector)
import Data.Monoid
import Control.Monad.Identity
import qualified Data.VectorSpace as VectorSpace
import Data.VectorSpace (VectorSpace, Scalar)

import qualified Numeric.AD as AD
import qualified Numeric.AD.Mode as AD
import qualified Numeric.AD.Internal.Forward as AD
import qualified Numeric.AD.Internal.Identity as AD
import qualified Numeric.AD.Internal.Kahn as AD
import qualified Numeric.AD.Internal.On as AD
import qualified Numeric.AD.Internal.Reverse as AD
import qualified Numeric.AD.Internal.Sparse as AD


-- | @AnyRF cb@ is a function that uses variables from the nonlinear
-- program in a way supported by 'AnyRFCxt'. The @cb@ is
-- usually 'Identity'
data AnyRF cb = AnyRF (forall a. AnyRFCxt a => Vector a -> cb a)

-- | RealFloat gives most numerical operations,
-- 'VectorSpace' is involved to allow using definitions from the
-- <http://hackage.haskell.org/package/splines splines> package
type AnyRFCxt a = (VectorSpace a, RealFloat a, VectorSpace.Scalar a ~ a)

-- *** helpers for defining instances
liftOp0 :: (forall a. AnyRFCxt a => a) -> AnyRF Identity
liftOp0 op = AnyRF $ \x -> Identity op

liftOp1 :: (forall a. AnyRFCxt a => a -> a) -> AnyRF Identity -> AnyRF Identity
liftOp1 op (AnyRF a) = AnyRF $ \x -> Identity (op (runIdentity (a x)))

liftOp2 :: (forall a. AnyRFCxt a => a -> a -> a) -> AnyRF Identity -> AnyRF Identity -> AnyRF Identity
liftOp2 op (AnyRF a) (AnyRF b) = AnyRF $ \x -> Identity (runIdentity (a x) `op` runIdentity (b x))

instance Num (AnyRF Identity) where
    (+) = liftOp2 (+)
    (*) = liftOp2 (*)
    (-) = liftOp2 (-)
    abs = liftOp1 abs
    signum = liftOp1 signum
    fromInteger n = liftOp0 (fromInteger n)

instance Fractional (AnyRF Identity) where
    (/) = liftOp2 (/)
    recip = liftOp1 recip
    fromRational n = liftOp0 (fromRational n)

instance Floating (AnyRF Identity) where
    pi = liftOp0 pi
    exp = liftOp1 exp
    sqrt = liftOp1 sqrt
    log = liftOp1 log
    sin = liftOp1 sin
    tan = liftOp1 tan
    cos = liftOp1 cos
    asin = liftOp1 asin
    atan = liftOp1 atan
    acos = liftOp1 acos
    sinh = liftOp1 sinh
    tanh = liftOp1 tanh
    cosh = liftOp1 cosh
    asinh = liftOp1 asinh
    atanh = liftOp1 atanh
    acosh = liftOp1 acosh
    (**) = liftOp2 (**)
    logBase = liftOp2 logBase

instance Real (AnyRF Identity) where
    toRational _ = error "Real AnyRF Identity"

instance Ord (AnyRF Identity) where
    compare _ = error "anyRF compare"
    max = liftOp2 max
    min = liftOp2 min
instance Eq (AnyRF Identity) where
        (==) = error "anyRF =="
instance RealFrac (AnyRF Identity) where
    properFraction = error "properFraction AnyRF"
instance RealFloat (AnyRF Identity) where
    isInfinite = error "isInfinite AnyRF"
    isNaN = error "isNaN AnyRF"
    decodeFloat = error "decodeFloat AnyRF"
    floatRange = error "floatRange AnyRF"
    isNegativeZero = error "isNegativeZero AnyRF"
    isIEEE = error "isIEEE AnyRF"
    isDenormalized = error "isDenormalized AnyRF"
    floatDigits _ = floatDigits (error "RealFrac AnyRF Identity floatDigits" :: Double)
    floatRadix _ = floatRadix   (error "RealFrac AnyRF Identity floatRadix" :: Double)
    atan2 = liftOp2 atan2
    significand = liftOp1 significand
    scaleFloat n = liftOp1 (scaleFloat n)
    encodeFloat a b = liftOp0 (encodeFloat a b)

instance Semigroup (AnyRF Seq) where
    AnyRF f <> AnyRF g = AnyRF (f `mappend` g)

instance Monoid (AnyRF Seq) where
    AnyRF f `mappend` AnyRF g = AnyRF (f `mappend` g)
    mempty = AnyRF mempty

instance VectorSpace.VectorSpace (AnyRF Identity) where
        type Scalar (AnyRF Identity) = Double
        x *^ v = realToFrac x*v

instance VectorSpace.AdditiveGroup (AnyRF Identity) where
        zeroV = liftOp0 0
        (^+^) = (+)
        negateV = negate

-- * orphan instances 

-- $orphans
-- these belong somewhere between the @ad@ package and @vector-space@


instance (AD.Mode a) => VectorSpace.AdditiveGroup (AD.AD s a) where
  zeroV = AD.zero
  (^+^) = (+)
  negateV = negate

instance (AD.Mode a) => VectorSpace.VectorSpace (AD.AD s a) where
  type Scalar (AD.AD s a) = AD.AD s a
  (*^) = (*)

instance (AD.Mode (AD.On a)) => VectorSpace.AdditiveGroup (AD.On a) where
  zeroV = AD.zero
  (^+^) = (+)
  negateV = negate

instance (AD.Mode (AD.On a)) => VectorSpace.VectorSpace (AD.On a) where
  type Scalar (AD.On a) = AD.On a
  (*^) = (*)



instance (Num a, AD.Mode (AD.Reverse s a)) => VectorSpace.AdditiveGroup (AD.Reverse s a) where
  zeroV = AD.zero
  (^+^) = (+)
  negateV = negate

instance (Num a, AD.Mode (AD.Reverse s a)) => VectorSpace.VectorSpace (AD.Reverse s a) where
  type Scalar (AD.Reverse s a) = AD.Reverse s a
  (*^) = (*)



instance Num a => VectorSpace.VectorSpace (AD.Sparse a) where
  type Scalar (AD.Sparse a) = AD.Sparse a
  (*^) = (*)

instance (Num a, AD.Mode (AD.Sparse a)) => VectorSpace.AdditiveGroup (AD.Sparse a) where
  zeroV = AD.zero
  (^+^) = (+)
  negateV = negate
