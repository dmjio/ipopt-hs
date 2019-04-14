{-#  LANGUAGE ForeignFunctionInterface #-}
{- |

Copyright: (C) 2013 Adam Vogt
Maintainer: Adam Vogt <vogt.adam@gmail.com>
Stability: unstable
Description: lowest-level parts of the binding


-}
module Ipopt.Raw (
    -- * specifying problem
    createIpoptProblemAD,
    createIpoptProblemADSparse,

    -- ** solve
    ipoptSolve,
    IpOptSolved(..),

    -- ** solver options
    addIpoptNumOption,
    addIpoptStrOption,
    addIpoptIntOption,
    openIpoptOutputFile,

    -- * types
    Vec,

    IpNumber(..),
    IpIndex(..),
    IpInt(..),
    IpBool(..),

    IpF(..),
    IpGradF(..),
    IpG(..),
    IpJacG(..),
    IpH(..),

    IpProblem(..),

    IntermediateCB,

    ApplicationReturnStatus(..),

    -- * lower-level parts of the binding
    createIpoptProblem,

    freeIpoptProblem,
    setIpoptProblemScaling,
    setIntermediateCallback,


    -- ** marshalling functions
    wrapIpF,
    wrapIpGradF,
    wrapIpG,
    wrapIpJacG,
    wrapIpH,

    wrapIpF1,
    wrapIpGradF1,
    wrapIpG1,
    wrapIpJacG1,
    wrapIpH1,

    wrapIpF2,
    wrapIpGradF2,
    wrapIpG2,
    wrapIpJacG2,
    wrapIpH2,

    ) where

import C2HS
import Control.Exception
import Control.Monad
import Data.IORef
import Foreign.C
import Foreign.ForeignPtr
import Foreign.Ptr
import Foreign.Storable
import Numeric.AD
import qualified Data.Vector as V
import qualified Data.Vector.Generic as VG
import qualified Data.Vector.Storable as VS
import qualified Data.Vector.Storable.Mutable as VM

import qualified Numeric.AD.Internal.Identity as AD

import qualified Numeric.AD.Rank1.Sparse as Sparse
import qualified Numeric.AD.Internal.Sparse as Sparse
import Numeric.AD.Rank1.Sparse (Sparse)

import Ipopt.AnyRF
import Data.VectorSpace (VectorSpace,Scalar)

#include "IpStdCInterface.h"

type IpNumber = {# type Number #}
type IpIndex = {# type Index #}
type IpInt = {# type Int #}
type IpBool = {# type Bool #}

type IpF = {# type Eval_F_CB #}
type IpGradF = {# type Eval_Grad_F_CB #}
type IpG = {# type Eval_G_CB #}
type IpJacG = {# type Eval_Jac_G_CB #}
type IpH = {# type Eval_H_CB #}

type IpIntermediateCB = {# type Intermediate_CB #}

type IntermediateCB = CInt -- ^ alg_mod (0 regular, 1 is resto)
  -> CInt -- ^ iter count
  -> Double -- ^ obj value
  -> Double -- ^ inf_pr
  -> Double -- ^ inf_du
  -> Double -- ^ mu
  -> Double -- ^ d_norm
  -> Double -- ^ regularization_size
  -> Double -- ^ alpha_du
  -> Double -- ^ alpha_pr
  -> CInt -- ^ ls_trials
  -> Ptr () -- ^ user_data (usually null)
  -> IO IpBool

{#enum ApplicationReturnStatus as ^ {underscoreToCase} deriving (Show) #}
{#enum AlgorithmMode as ^ {underscoreToCase} deriving (Show) #}


{#pointer IpoptProblem as IpProblem foreign newtype #}
ipp x = withIpProblem x

type family UnFunPtr a
type instance UnFunPtr (FunPtr a) = a

{#enum define IpBoolT { TRUE as IpTrue, FALSE as IpFalse } deriving (Eq,Ord,Show) #}

ipTrue = toEnum (fromEnum IpTrue) :: IpBool
ipFalse = toEnum (fromEnum IpFalse) :: IpBool


-- | likely an unsafe method for getting a "Data.Vector.Storable.Mutable" out of a 'Ptr'
ptrToVS n p = do
    fp <- newForeignPtr_ p
    return (VM.unsafeFromForeignPtr0 fp (fromIntegral n))

foreign import ccall "wrapper" wrapIpF1 :: UnFunPtr IpF -> IO IpF
foreign import ccall "wrapper" wrapIpG1 :: UnFunPtr IpG -> IO IpG
foreign import ccall "wrapper" wrapIpGradF1 :: UnFunPtr IpGradF -> IO IpGradF
foreign import ccall "wrapper" wrapIpJacG1 :: UnFunPtr IpJacG -> IO IpJacG
foreign import ccall "wrapper" wrapIpH1 :: UnFunPtr IpH -> IO IpH

foreign import ccall "wrapper" wrapIpIntermediateCB :: IntermediateCB -> IO IpIntermediateCB


toB x =  either (\ e@SomeException {} -> print e >> return ipFalse)
                (\ _ -> return ipTrue ) =<< try x

wrapIpF2' fun n xin new_x obj_val _userData = do
    toB $ poke obj_val =<< fun =<< ptrToVS n xin

wrapIpF2 fun n xin new_x obj_val _userData = do
    toB $ poke obj_val =<< fun =<< ptrToVS n xin

wrapIpG2 fun n xin new_x m gout _userData = do
    toB $ join $ liftM2 VM.copy (ptrToVS m gout) (fun =<< ptrToVS n xin)

wrapIpGradF2 fun n x new_x grad_f _userData = do
    toB $ join $ liftM2 VM.copy (ptrToVS n grad_f) (fun =<< ptrToVS n x)

wrapIpJacG2 fun1 fun2 n x new_x m nj iRow jCol jacs _userData
    | jacs == nullPtr = do
            toB $ join $ liftM2 fun1 (ptrToVS nj iRow) (ptrToVS nj jCol)
    | otherwise = do
            toB $ join $ liftM2 fun2 (ptrToVS n x) (ptrToVS nj jacs)

wrapIpH2 funSparsity funEval n x new_x obj_factor m lambda new_lambda nHess iRow jCol values _userData
    | iRow == nullPtr = do
            toB $ join $ liftM3 (funEval obj_factor)
                        (ptrToVS m lambda)
                        (ptrToVS n x)
                        (ptrToVS nHess values)
    | otherwise = do
            toB $ join $ liftM2 funSparsity
                        (ptrToVS nHess iRow)
                        (ptrToVS nHess jCol)

wrapIpF f = wrapIpF1 (wrapIpF2 f)
wrapIpG f = wrapIpG1 (wrapIpG2 f)
wrapIpGradF f = wrapIpGradF1 (wrapIpGradF2 f)
wrapIpJacG f1 f2 = wrapIpJacG1 (wrapIpJacG2 f1 f2)
wrapIpH fSparsity fEval = wrapIpH1 (wrapIpH2 fSparsity fEval)


vmUnsafeWith :: Vec -> (Ptr CDouble -> IO r) -> IO r
vmUnsafeWith v f = VM.unsafeWith v (f . castPtr)

-- | Vector of numbers
type Vec = VM.IOVector Double

-- depend on CDouble being just a newtype on Double which
-- is the same as the IpNumber defined in the ipopt header, so this
-- should cause a compile failure if that's not the case...
_ = CDouble (5 :: Double) :: IpNumber

fromVec :: VG.Vector v Double => Vec -> IO (v Double)
fromVec mv = do
    v <- VS.freeze mv 
    return (VG.convert v)


createIpoptProblem :: Vec -> Vec -> Vec -> Vec
    -> Int -> Int -> IpF -> IpG -> IpGradF -> IpJacG -> IpH -> IO IpProblem
createIpoptProblem xL xU gL gU nJac nHess f g gradF jacG hess
    | lx <- VM.length xL,
      lx == VM.length xU,
      lg <- VM.length gL,
      lg == VM.length gU = do
        p <- createIpoptProblem3 lx xL xU lg gL gU nJac nHess 0 f g gradF jacG hess
        p' <- newForeignPtr freeIpoptProblem (castPtr p)
        return (IpProblem (castForeignPtr p'))
   | otherwise = error "dimensions wrong!"


{#fun CreateIpoptProblem as createIpoptProblem3
        { `Int', vmUnsafeWith* `Vec', vmUnsafeWith* `Vec',
          `Int', vmUnsafeWith* `Vec', vmUnsafeWith* `Vec',
          `Int', `Int', `Int', id `IpF', id `IpG', id `IpGradF',
          id `IpJacG', id `IpH' } -> `Ptr IpProblem' id #}

{#fun AddIpoptNumOption as ^ { ipp* `IpProblem', `String', `Double' } -> `Bool' #}

{#fun AddIpoptStrOption as ^ { ipp* `IpProblem', `String', `String' } -> `Bool' #}

{#fun AddIpoptIntOption as ^ { ipp* `IpProblem', `String', `Int' } -> `Bool' #}

foreign import ccall unsafe "&FreeIpoptProblem"
    freeIpoptProblem :: FunPtr (Ptr () -> IO ())

{#fun OpenIpoptOutputFile as ^ { ipp* `IpProblem', `String', `Int' } -> `Bool' #}

{#fun SetIpoptProblemScaling as ^
    { ipp* `IpProblem',
    `Double',
     vmUnsafeWith* `Vec',
     vmUnsafeWith* `Vec'
     } -> `Bool' #}

{#fun SetIntermediateCallback as setIntermediateCallback1
  { ipp* `IpProblem',
    id `IpIntermediateCB' } -> `Bool' #}

setIntermediateCallback pp cb = do
  cb' <- wrapIpIntermediateCB cb
  setIntermediateCallback1 pp cb'

-- | lenses are in "Ipopt.PP"
data IpOptSolved v = IpOptSolved
  { _ipOptSolved_status :: ApplicationReturnStatus,
    _ipOptSolved_objective :: Double,
    _ipOptSolved_x,
    _ipOptSolved_g,
    _ipOptSolved_mult_g,
    _ipOptSolved_mult_x_L,
    _ipOptSolved_mult_x_U :: v Double }

ipoptSolve :: VG.Vector v Double => IpProblem
    -> Vec -- ^ starting point @x@. Note that the value is overwritten with the final @x@.
    -> IO (IpOptSolved v)
ipoptSolve problem x = do
    g <- VM.new (VM.length x)
    mult_g <- VM.new (VM.length x)
    mult_x_L <- VM.new (VM.length x)
    mult_x_U <- VM.new (VM.length x)

    out <- ipoptSolve2
     problem
     x
     g
     mult_g
     mult_x_L
     mult_x_U
     nullPtr
    
    x'        <- fromVec x
    g'        <- fromVec g
    mult_g'   <- fromVec mult_g 
    mult_x_L' <- fromVec mult_x_L
    mult_x_U' <- fromVec mult_x_U

    return $ IpOptSolved
      (fst out)
      (snd out)
      x'
      g'
      mult_g'
      mult_x_L'
      mult_x_U'

{#fun IpoptSolve as ipoptSolve2
        { ipp* `IpProblem',
        vmUnsafeWith* `Vec',
        vmUnsafeWith* `Vec',
        alloca- `Double' peekFloatConv*,
        vmUnsafeWith* `Vec',
        vmUnsafeWith* `Vec',
        vmUnsafeWith* `Vec',
        id `Ptr ()' } -> `ApplicationReturnStatus' cToEnum #}


{- | Set-up an 'IpProblem' to be solved later. Only objective function (@f@)
and constraint functions (@g@) need to be specified. Derivatives needed by ipopt
are computed by "Numeric.AD".

To solve the optimization problem:

>              min f(x)
>   such that
>           xL <=  x     <= xU
>           gL <=  g(x)  <= gU

First create an opaque 'IpProblem' object (nlp):

> nlp <- createIpOptProblemAD xL xU gL gU f g

Then pass it off to 'ipoptSolve'.

> ipoptSolve nlp x0

Refer to @examples/HS71ad.hs@ for details of setting up the vectors supplied.
-}
createIpoptProblemAD
    :: Vec -- ^ @xL@ 'VM.Vector' of lower bounds for decision variables with length @n@
    -> Vec -- ^ @xU@ 'VM.Vector' of upper bounds for decision variables
    -> Vec -- ^ @gL@ 'VM.Vector' of lower bounds for constraint functions @g(x)@ with length @m@
    -> Vec -- ^ @gU@ 'VM.Vector' of upper bounds for constraint functions @g(x)@
    -> (forall a. AnyRFCxt a => V.Vector a -> a) -- ^ objective function @f : R^n -> R@
    -> (forall a. AnyRFCxt a => V.Vector a -> V.Vector a) -- ^ constraint functions @g : R^n -> R^m@
    -> IO IpProblem
createIpoptProblemAD xL xU gL gU f g
    | n <- VM.length xL,
      n == VM.length xU,
      m <- VM.length gL,
      m == VM.length gU = do
    (eval_f, eval_grad_f, eval_g, eval_jac_g, eval_h) <- mkFs n m f g
    createIpoptProblem xL xU gL gU (n*m) (((n+1)*n) `div` 2)
            eval_f eval_g eval_grad_f eval_jac_g eval_h


{- | this is 50% slower than 'createIpoptProblemAD' in one instance
<http://code.haskell.org/~aavogt/ipopt-hs/examples/bench.html#williams-otto-process>).
But the benefit is that no RankN types are used (so it is possible to implement
more functions without having to modify 'AnyRFCxt')
-}
createIpoptProblemADSparse

    :: Vec -- ^ @xL@ 'VM.Vector' of lower bounds for decision variables with length @n@
    -> Vec -- ^ @xU@ 'VM.Vector' of upper bounds for decision variables
    -> Vec -- ^ @gL@ 'VM.Vector' of lower bounds for constraint functions @g(x)@ with length @m@
    -> Vec -- ^ @gU@ 'VM.Vector' of upper bounds for constraint functions @g(x)@
    -> (V.Vector (Sparse CDouble) -> Sparse CDouble) -- ^ objective function @f : R^n -> R@
    -> (V.Vector (Sparse CDouble) -> V.Vector (Sparse CDouble)) -- ^ constraint functions @g : R^n -> R^m@
    -> IO IpProblem
createIpoptProblemADSparse xL xU gL gU f g
    | n <- VM.length xL,
      n == VM.length xU,
      m <- VM.length gL,
      m == VM.length gU = do
    (eval_f, eval_grad_f, eval_g, eval_jac_g, eval_h) <- mkFsSparse n m f g
    createIpoptProblem xL xU gL gU (n*m) (((n+1)*n) `div` 2)
            eval_f eval_g eval_grad_f eval_jac_g eval_h


mkFs ::
       Int -- ^ @n@ number of variables
    -> Int -- ^ @m@ number of constraints
    -> (forall a. AnyRFCxt a => V.Vector a -> a) -- ^ objective function @R^n -> R@
    -> (forall a. AnyRFCxt a => V.Vector a -> V.Vector a) -- ^ constraint functions @R^n -> R^m@
    -> IO (IpF, IpGradF, IpG, IpJacG, IpH)
mkFs n m f g = do
  ipF <- wrapIpF $ \x -> do
    x <- VG.convert `fmap` VS.unsafeFreeze x
    return $ f x

  ipGradF <- wrapIpGradF $ \x -> do
    x <- VG.convert `fmap` VS.unsafeFreeze x
    VS.unsafeThaw $ VG.convert (grad f x)

  ipG <- wrapIpG $ \x -> do
    x <- VG.convert `fmap` VS.unsafeFreeze x
    VS.unsafeThaw $ VG.convert (g x)

  ipJacG <- wrapIpJacG (denseIJ n m) $ \x y -> do
    x <- VG.convert `fmap` VS.unsafeFreeze x
    jac <- VS.unsafeThaw $ VG.convert $ VG.concat $ VG.toList $ jacobian g x
    VM.copy y jac

  ipH <- wrapIpH (denseIJh n m)
     ( \ obj_factor lambda x values -> do

        x <- VG.convert `fmap` VS.unsafeFreeze x
        lambda <- VG.convert `fmap` VS.unsafeFreeze lambda

        let tri = VG.concat . VG.toList . V.imap (\n -> V.take (n+1))
            obj = V.map (*obj_factor) $ tri $ hessian f x
            gj = V.zipWith (\l v -> V.map (l*) v) lambda (V.map tri (hessianF g x))
            lagrangian = V.foldl (V.zipWith (+)) obj gj

        VM.copy values =<< VS.unsafeThaw (VG.convert lagrangian)
        )

  return (ipF, ipGradF, ipG, ipJacG, ipH)

sparsePrimal :: Num a => Sparse a -> a
sparsePrimal (Sparse.Sparse a _) = a
sparsePrimal Sparse.Zero = 0

mkFsSparse ::
       Int -- ^ @n@ number of variables
    -> Int -- ^ @m@ number of constraints
    -> (V.Vector (Sparse CDouble) -> Sparse CDouble) -- ^ objective function @R^n -> R@
    -> (V.Vector (Sparse CDouble) -> V.Vector (Sparse CDouble)) -- ^ constraint functions @R^n -> R^m@
    -> IO (IpF, IpGradF, IpG, IpJacG, IpH)
mkFsSparse n m f g = do
  ipF <- wrapIpF $ \x -> do
    x <- VG.convert `fmap` VS.unsafeFreeze x
    return $ sparsePrimal $ f (Sparse.vars x)

  ipGradF <- wrapIpGradF $ \x -> do
    x <- VG.convert `fmap` VS.unsafeFreeze x
    VS.unsafeThaw $ VG.convert (Sparse.grad f x)

  ipG <- wrapIpG $ \x -> do
    x <- VG.convert `fmap` VS.unsafeFreeze x
    VS.unsafeThaw $ VG.convert $ V.map sparsePrimal $ g $ Sparse.vars x

  ipJacG <- wrapIpJacG (denseIJ n m) $ \x y -> do
    x <- VG.convert `fmap` VS.unsafeFreeze x
    jac <- VS.unsafeThaw $ VG.convert $ VG.concat $ VG.toList $ Sparse.jacobian g x
    VM.copy y jac

  ipH <- wrapIpH (denseIJh n m)
     ( \ obj_factor lambda x values -> do

        x <- VG.convert `fmap` VS.unsafeFreeze x
        lambda <- VG.convert `fmap` VS.unsafeFreeze lambda

        let tri = VG.concat . VG.toList . V.imap (\n -> V.take (n+1))
            obj = V.map (*obj_factor) $ tri $ Sparse.hessian f x
            gj = V.zipWith (\l v -> V.map (l*) v) lambda (V.map tri (Sparse.hessianF g x))
            lagrangian = V.foldl (V.zipWith (+)) obj gj

        VM.copy values =<< VS.unsafeThaw (VG.convert lagrangian)
        )

  return (ipF, ipGradF, ipG, ipJacG, ipH)

-- | indexes the same as http://www.coin-or.org/Ipopt/documentation/node40.html
denseIJ n m iRow jCol = do
    VM.copy iRow =<< VS.unsafeThaw (VS.generate (n*m) (\x -> fromIntegral $ x `div` n))
    VM.copy jCol =<< VS.unsafeThaw (VS.generate (n*m) (\x -> fromIntegral $ x `mod` n))

-- | indexes the same as http://www.coin-or.org/Ipopt/documentation/node41.html
denseIJh n m iRow jCol = do
    i <- newIORef 0
    forM_ [0 .. fromIntegral n-1] $ \ row ->
        forM_ [ 0 .. row ] $ \col -> do
            ii <- readIORef i
            VM.write iRow ii row
            VM.write jCol ii col
            writeIORef i (ii+1)

