{-# LANGUAGE TypeFamilies,DeriveDataTypeable, RankNTypes, ConstraintKinds, FlexibleContexts #-}
module Nlopt.Raw where

import Ipopt.AnyRF

import Foreign.C.Types
import Foreign.Ptr
import Foreign.Marshal
import Foreign.Storable
import Data.Int
import C2HS

import qualified Data.Vector.Storable.Mutable as VM
import qualified Data.Vector as V
import qualified Data.Vector.Storable as VS
import qualified Data.Vector.Generic as VG

import Control.Exception
import Data.Typeable
import Numeric.AD (grad, jacobian, hessian)
import Control.Monad

#include "nlopt.h"

-- * converting haskell functions

{- $conventions

NLOpt has three different types for functions 'FunPtrFunc' 'FunPtrMFunc' and 'FunPtrPrecond'.

[@AD@] means derivatives are calculated, and the haskell function does no IO.

[@G@] mean you are providing the derivatives

[@M@] means m functions are calculated at a time

-}

-- | an exact hessian calculated with AD. See 'toPrecondG'
-- XXX BFGS approx could also be done...
toPrecondAD :: (forall a. AnyRFCxt a => V.Vector a -> a) -> Precond
toPrecondAD f n xs vs r _usrData =
    toPrecondG (\x v -> return (hessian f x `mXv` v)) n xs vs r _usrData

-- | see <http://ab-initio.mit.edu/wiki/index.php/NLopt_Reference#Preconditioning_with_approximate_Hessians>
-- only applies to 'NLOPT_LD_CCSAQ'
toPrecondG :: (VG.Vector vec Double, v ~ vec Double)
    => (v -> v -> IO v) -- ^ given @x,v@ calculate  @H(x)v@ where H the hessian
    -> Precond
toPrecondG f n xs vs r _usrData = do
    xs' <- ptrToV n xs
    vs' <- ptrToV n vs
    copyInto n r =<< f xs' vs'

toFuncM :: VG.Vector v Double => (v Double -> IO (v Double) ) -> MFunc
toFuncM f m r n xs _g _usrData = do
    copyInto m r =<< f =<< ptrToV n xs

-- | @n@ and @m@ type variables indicate the vector size as
-- number of inputs and number of outputs respectively
toFuncMG :: (VG.Vector n Double, VG.Vector n (m Double), n ~ m)
    => (n Double -> IO (m Double) ) -- @f@
    -> (n Double -> IO (n (m Double))) -- @grad f@
    -> MFunc
toFuncMG f g m r n xs g' _usrData = do
    xs' <- ptrToV n xs
    copyInto m r =<< f xs'
    -- probably should do this better...
    copyInto (n*m) g' . VG.concat . VG.toList =<< g xs'

toFuncMAD :: (forall a. AnyRFCxt a => V.Vector a -> V.Vector a) -> MFunc
toFuncMAD f m r n xs g _usrData = do
    xs' <- ptrToV n xs
    copyInto m r (f xs')
    -- grad[i*n + j] should be satisfied...
    copyInto (n*m) g (V.concat (V.toList (jacobian f xs')))

toFunc :: (VG.Vector v Double) => (v Double -> IO Double) -> Func
toFunc f n xs _g _usrData = fmap CDouble $ f =<< ptrToV n xs

-- | where the gradient happens via AD
toFuncAD :: (forall a. AnyRFCxt a => V.Vector a -> a) -> Func
toFuncAD f n xs g _usrData = do
    xs' <- ptrToV n xs
    copyInto n g (grad f xs')
    return (CDouble (f xs'))

toFuncG :: (VG.Vector v Double) => (v Double -> IO Double) -- ^ @f@
    -> (v Double -> IO (v Double)) -- ^ @grad(f)@
    -> Func
toFuncG f g n xs g' _usrData = do
    xs' <- ptrToV n xs
    copyInto n g' =<< g xs'
    CDouble `fmap` f xs'

type family UnFunPtr a
type instance UnFunPtr (FunPtr a) = a

type Func = UnFunPtr FunPtrFunc
type MFunc = UnFunPtr FunPtrMFunc
type Precond = UnFunPtr FunPtrPrecond

type FunPtrFunc    = {# type nlopt_func #}
type FunPtrMFunc   = {# type nlopt_mfunc #}
type FunPtrPrecond = {# type nlopt_precond #}

{#pointer *nlopt_opt as NLOpt foreign newtype #}
{#enum nlopt_algorithm as ^ {} deriving (Show,Bounded,Ord,Eq) #}

-- | negative (above NLOPT_SUCCESS) values of these are thrown as exceptions. The positive ones are
-- return values.
{#enum nlopt_result as ^ {} deriving (Show,Typeable) #}

instance Exception NloptResult

checkEC :: CInt -> IO NloptResult
checkEC n | n `elem` [-1,-2,-3] = throwIO e
          | otherwise = return e
    where e = fromCInt n :: NloptResult

{#fun nlopt_srand as ^ { `Int' } -> `()' #}
{#fun nlopt_srand_time as ^ { } -> `()' #}
{#fun nlopt_version as ^ {
        alloca- `Int' peekInt*,
        alloca- `Int' peekInt*,
        alloca- `Int' peekInt*} -> `()' #}
{#fun nlopt_create as ^ { toCInt `NloptAlgorithm', `Int' } -> `NLOpt' ptrToNLOpt * #}
{#fun nlopt_copy as ^ { withNLOpt_* `NLOpt' } -> `NLOpt' ptrToNLOpt * #}

-- | should not need to be called manually
{#fun nlopt_destroy as ^ { id `Ptr ()' } -> `()' #}

{#fun nlopt_optimize as ^ {
    withNLOpt_* `NLOpt', vmUnsafeWith* `Vec', alloca- `Double' peekDouble* }
        -> `NloptResult' checkEC* #}

{#fun nlopt_set_min_objective as ^ {
    withNLOpt_* `NLOpt', withFunc* `Func', withNull- `()' }
        -> `NloptResult' checkEC* #}
{#fun nlopt_set_max_objective as ^ {
    withNLOpt_* `NLOpt', withFunc* `Func', withNull- `()' }
        -> `NloptResult' checkEC* #}

{#fun nlopt_set_precond_min_objective as ^ {
    withNLOpt_* `NLOpt',
    withFunc* `Func',
    withPrecond* `Precond',
    withNull- `()' }
    -> `NloptResult' checkEC*  #}

{#fun nlopt_set_precond_max_objective as ^ {
    withNLOpt_* `NLOpt',
    withFunc* `Func',
    withPrecond* `Precond',
    withNull- `()' } -> `NloptResult' checkEC* #}

{#fun nlopt_get_algorithm as ^ { withNLOpt_* `NLOpt' } -> `NloptAlgorithm' fromCInt #}
{#fun nlopt_get_dimension as ^ { withNLOpt_* `NLOpt' } -> `Int' #}

-- * constraints
{#fun nlopt_set_lower_bounds as ^ { withNLOpt_* `NLOpt', vmUnsafeWith* `Vec' }
    -> `NloptResult' checkEC* #}
{#fun nlopt_set_upper_bounds as ^ { withNLOpt_* `NLOpt', vmUnsafeWith* `Vec' }
    -> `NloptResult' checkEC* #}

{#fun nlopt_get_lower_bounds as ^ { withNLOpt_* `NLOpt', vmUnsafeWith* `Vec' }
    -> `NloptResult' checkEC* #}
{#fun nlopt_get_upper_bounds as ^ { withNLOpt_* `NLOpt', vmUnsafeWith* `Vec' }
    -> `NloptResult' checkEC* #}

{#fun nlopt_set_lower_bounds1 as ^ { withNLOpt_* `NLOpt', `Double' }
    -> `NloptResult' checkEC* #}
{#fun nlopt_set_upper_bounds1 as ^ { withNLOpt_* `NLOpt', `Double' }
    -> `NloptResult' checkEC* #}

{#fun nlopt_remove_inequality_constraints as ^ { withNLOpt_* `NLOpt' }
    -> `NloptResult' checkEC* #}

{#fun nlopt_add_inequality_constraint as ^
    { withNLOpt_* `NLOpt',
      withFunc* `Func',
      withNull- `()',
      `Double' } -> `NloptResult' checkEC* #}

{#fun nlopt_add_precond_inequality_constraint as ^
    { withNLOpt_* `NLOpt',
      withFunc* `Func',
      withPrecond* `Precond',
      withNull- `()',
      `Double' } -> `NloptResult' checkEC* #}

{#fun nlopt_add_inequality_mconstraint as ^
    { withNLOpt_* `NLOpt',
      `Int',
      withMFunc* `MFunc',
      withNull- `()',
      vmUnsafeWith* `Vec' } -> `NloptResult' checkEC* #}

{#fun nlopt_remove_equality_constraints as ^ { withNLOpt_* `NLOpt' }
    -> `NloptResult' checkEC* #}

{#fun nlopt_add_equality_constraint as ^
    { withNLOpt_* `NLOpt',
      withFunc* `Func',
      withNull- `()',
      `Double' } -> `NloptResult' checkEC* #}

{#fun nlopt_add_precond_equality_constraint as ^
    { withNLOpt_* `NLOpt',
      withFunc* `Func',
      withPrecond* `Precond',
      withNull- `()',
      `Double' } -> `NloptResult' checkEC* #}

{#fun nlopt_add_equality_mconstraint as ^
    { withNLOpt_* `NLOpt',
      `Int',
      withMFunc* `MFunc',
      withNull- `()',
      vmUnsafeWith* `Vec' } -> `NloptResult' checkEC* #}

-- * stopping criteria
{#fun nlopt_set_stopval as ^ { withNLOpt_* `NLOpt', `Double' }
    -> `NloptResult' checkEC* #}

{#fun nlopt_get_stopval as ^ { withNLOpt_* `NLOpt' } -> `Double' #}

{#fun nlopt_set_ftol_rel as ^ { withNLOpt_* `NLOpt', `Double' }
    -> `NloptResult' checkEC* #}

{#fun nlopt_get_ftol_rel as ^ { withNLOpt_* `NLOpt' } -> `Double' #}

{#fun nlopt_set_ftol_abs as ^ { withNLOpt_* `NLOpt', `Double' }
    -> `NloptResult' checkEC* #}

{#fun nlopt_get_ftol_abs as ^ { withNLOpt_* `NLOpt' } -> `Double' #}

{#fun nlopt_set_xtol_rel as ^ { withNLOpt_* `NLOpt', `Double' }
    -> `NloptResult' checkEC* #}

{#fun nlopt_get_xtol_rel as ^ { withNLOpt_* `NLOpt' } -> `Double' #}

{#fun nlopt_set_xtol_abs1 as ^ { withNLOpt_* `NLOpt', `Double' }
    -> `NloptResult' checkEC* #}

{#fun nlopt_set_xtol_abs as ^ { withNLOpt_* `NLOpt', vmUnsafeWith* `Vec' }
    -> `NloptResult' checkEC* #}

{#fun nlopt_get_xtol_abs as ^ { withNLOpt_* `NLOpt', vmUnsafeWith* `Vec' }
    -> `NloptResult' checkEC* #}

{#fun nlopt_set_maxeval as ^ { withNLOpt_* `NLOpt', `Int' }
    -> `NloptResult' checkEC* #}

{#fun nlopt_get_maxeval as ^ { withNLOpt_* `NLOpt' } -> `Int'  #}

{#fun nlopt_set_maxtime as ^ { withNLOpt_* `NLOpt', `Double' }
    -> `NloptResult' checkEC* #}

{#fun nlopt_get_maxtime as ^ { withNLOpt_* `NLOpt' } -> `Double' #}

{#fun nlopt_force_stop as ^ { withNLOpt_* `NLOpt' } -> `NloptResult' checkEC* #}
{#fun nlopt_set_force_stop as ^ { withNLOpt_* `NLOpt', `Int' } -> `NloptResult' checkEC* #}
{#fun nlopt_get_force_stop as ^ { withNLOpt_* `NLOpt' } -> `Int' #}


-- * more algorithm-specific parameters
{#fun nlopt_set_local_optimizer as ^
    { withNLOpt_* `NLOpt', withNLOpt_* `NLOpt' } -> `NloptResult' checkEC* #}

{#fun nlopt_set_population as ^
    { withNLOpt_* `NLOpt', `Int' } -> `NloptResult' checkEC* #}

{#fun nlopt_get_population as ^ { withNLOpt_* `NLOpt' } -> `Int' #}

{#fun nlopt_set_vector_storage as ^
    { withNLOpt_* `NLOpt', `Int' } -> `NloptResult' checkEC* #}

{#fun nlopt_get_vector_storage as ^ { withNLOpt_* `NLOpt' } -> `Int' #}

{#fun nlopt_set_initial_step as ^ { withNLOpt_* `NLOpt',
    vmUnsafeWith* `Vec' } -> `NloptResult' checkEC* #}

{#fun nlopt_get_initial_step as ^ { withNLOpt_* `NLOpt',
    vmUnsafeWith* `Vec',
    vmUnsafeWith* `Vec' } -> `NloptResult' checkEC* #}

{#fun nlopt_set_initial_step1 as ^ { withNLOpt_* `NLOpt',
    `Double' } -> `NloptResult' checkEC* #}

-- * utils for ffi

-- $note
-- much is copied from Ipopt.Raw

withNLOpt_ x f = withNLOpt x (f . castPtr)

vmUnsafeWith :: VM.IOVector Double -> (Ptr CDouble -> IO b) -> IO b
vmUnsafeWith v f = VM.unsafeWith v (f . castPtr)

ptrToVS :: Integral n => n -> Ptr CDouble -> IO Vec
ptrToVS n p = do
    fp <- newForeignPtr_ (castPtr p)
    return (VM.unsafeFromForeignPtr0 fp (fromIntegral n))

ptrToV :: (VG.Vector v Double, Integral n) => n -> Ptr CDouble -> IO (v Double)
ptrToV n p = fmap V.convert $ VS.unsafeFreeze =<< ptrToVS n p

copyInto _ p _ | p == nullPtr = return ()
copyInto n dest src = do
    to <- ptrToVS n dest
    VM.copy to =<< VS.unsafeThaw (V.convert src)

type Vec = VM.IOVector Double

toCInt x = fromIntegral (fromEnum x)
fromCInt x = toEnum (fromIntegral x)
peekInt ptr = fmap fromIntegral (peek ptr)


ptrToNLOpt :: Ptr () -> IO NLOpt
ptrToNLOpt p = do
    fp <- newForeignPtr nloptDestroyFP p
    return (NLOpt (castForeignPtr fp))

foreign import ccall "wrapper" mkNloptFinalizer :: (Ptr () -> IO ()) -> IO (FunPtr (Ptr () -> IO ()))

foreign import ccall unsafe "& nlopt_destroy" nloptDestroyFP :: FunPtr (Ptr () -> IO ())

foreign import ccall "wrapper" mkFunc :: Func -> IO (FunPtr Func)
foreign import ccall "wrapper" mkPrecond :: Precond -> IO (FunPtr Precond)
foreign import ccall "wrapper" mkMFunc :: MFunc -> IO (FunPtr MFunc)

withFunc :: Func -> (FunPtr Func -> IO b) -> IO b
withFunc f g = do
    f' <- mkFunc f
    r <- g f'
    -- freeHaskellFunPtr f'
    return r

withPrecond :: Precond -> (FunPtr Precond -> IO b) -> IO b
withPrecond f g = do
    f' <- mkPrecond f
    r <- g f'
    -- freeHaskellFunPtr f'
    return r

withMFunc :: MFunc -> (FunPtr MFunc -> IO b) -> IO b
withMFunc f g = do
    f' <- mkMFunc f
    r <- g f'
    -- freeHaskellFunPtr f'
    return r

withNull f = f nullPtr

-- | c2hs generates CDouble peek a Double instead
peekDouble :: Ptr CDouble -> IO Double
peekDouble p = peek (castPtr p)


-- | naive matrix × vector
mXv :: Num a => V.Vector (V.Vector a) -> V.Vector a -> V.Vector a
mXv m v = V.map (V.sum . V.zipWith (*) v) m

-- * c2hs-generated
