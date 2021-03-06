{-# LANGUAGE TemplateHaskell #-}
-- | Description : a EDSL for describing nonlinear programs
--
-- see usage in @examples/Test3.hs@ (and other examples)
--
-- IPOPT does support naming variables if you use c++
-- (by overriding a @virtual void finalize_metadata@), but
-- it's not clear that we can set that from c/haskell
module Ipopt.NLP where

import Ipopt.AnyRF
import Control.Monad.Fail
import Control.Applicative
import Control.Lens
import Control.Monad
import Control.Monad.Identity
import Control.Monad.State
import Data.Foldable (toList)
import Data.IntMap (IntMap)
import Data.List
import Data.Monoid
import Data.Monoid (First)
import Data.Sequence (Seq)
import Data.Vector (Vector)
import Foreign.C.Types (CDouble(..))
import qualified Data.Foldable as F
import qualified Data.IntMap as IM
import qualified Data.Map as M
import qualified Data.Sequence as Seq
import qualified Data.Set as S
import qualified Data.Vector as V
import qualified Data.Vector.Generic as VG
import Data.Vector ((!))
import qualified Data.Vector.Storable as VS
import Data.Maybe

import qualified Data.VectorSpace as VectorSpace

import Ipopt.Raw

-- * state
data NLPFun = NLPFun
    { _funF, _funG :: AnyRF Seq,
      _boundX, _boundG :: Seq (Double,Double) }
instance Show NLPFun where
    show (NLPFun f g a b) = "NLPFun <f> <g> {" ++ show a ++ "}{" ++ show b ++ "}"

data NLPState = NLPState
    { -- | current maximum index
      _nMax :: Ix,
      -- | what namespace are we currently in (see 'inEnv')
      _currentEnv :: [String],
      -- | fully qualified (see 'inEnv') name
      _variables :: M.Map String Ix,
      -- | invert _variables
      _variablesInv :: IxMap String,
      -- | human-readable descriptions for the constraint, objective and
      -- variables
      _constraintLabels, _objLabels :: IntMap String,
      _varLabels :: IxMap String,
      -- | in what environments is a given var used?
      _varEnv :: IxMap (S.Set [String]),
      _constraintEnv, _objEnv :: IntMap [String],
      _nlpfun :: NLPFun,
      -- | the default @(xL,xU)@ for @xL < x < xU@
      _defaultBounds :: (Double,Double),
      -- | for nlopt (lower/upper)
      _defaultConstraintTol :: (Double,Double),
      _constraintTol :: Seq (Double,Double),
      -- | inital state variable for the solver
      _initX :: Vector Double }
    deriving (Show)

newtype IxMap a = IxMap (IntMap a)
    deriving (Show, Monoid, Functor, Semigroup)

-- | the solver deals with arrays. This type is for indexes into the array
-- for the current variables that the solver is trying to find.
newtype Ix = Ix { _varIx :: Int } deriving Show


-- | the initial state to use when you actually have to get to IO
-- with the solution
nlpstate0 = NLPState (Ix (-1)) mempty mempty
    mempty mempty mempty -- labels
    mempty
    mempty mempty mempty -- env
    (mempty :: NLPFun)
    (-1/0, 1/0) -- bounds at infinity
    (1e-6,1e-6) -- arbitrary constraint tol
    mempty
    V.empty

type NLPT = StateT NLPState
type NLP = NLPT IO

instance Semigroup NLPFun where
    (<>) = mappend

instance Monoid NLPFun where
    NLPFun f g bx bg `mappend` NLPFun f' g' bx' bg' = NLPFun (f <> f') (g <> g') (bx <> bx') (bg <> bg')
    mempty = NLPFun mempty mempty mempty mempty

-- ** low level lenses to NLPState
makeLenses ''NLPFun
makeLenses ''Ix
makeLenses ''NLPState
ixMap x = iso IxMap (\(IxMap a) -> a) x -- could makeLenses ''IxMap

-- | should be a way to write an instance of At that'll make the normal
-- at/ix work with the IxMap / Ix (as opposed to IntMap/Int)
ix_ (Ix i) = from ixMap . ix i
at_ (Ix i) = from ixMap . at i

-- | @to@ is one of 'varEnv', 'constraintEnv', 'objEnv'
copyEnv to n = do
    ev <- use currentEnv
    cloneLens to %= IM.insert n ev

-- | @to@ should be 'constraintLabels', 'objLabels', 'varLabels'
addDesc to Nothing n = return ()
addDesc to (Just x) n = do
    to %= IM.insert n x

-- * high-level functions

-- | calls 'createIpoptProblemAD' and 'ipoptSolve'. To be used at the
-- end of a do-block.
solveNLP' :: (VG.Vector v Double, MonadIO m) =>
    (IpProblem -> IO ()) -- ^ set ipopt options (using functions from "Ipopt.Raw") or the 'ipopts' quasiquoter
    -> NLPT m (IpOptSolved v)
solveNLP' setOpts = do
    (xl,xu) <- join $ uses (nlpfun . boundX) seqToVecs
    (gl,gu) <- join $ uses (nlpfun . boundG) seqToVecs
    AnyRF fs <- use (nlpfun . funF)
    AnyRF gs <- use (nlpfun . funG)

    p <- liftIO (createIpoptProblemAD xl xu gl gu (F.sum . fs) (V.fromList . toList . gs))
    liftIO (setOpts p)
    x0 <- uses initX V.convert
    r <- liftIO (ipoptSolve p =<< VS.thaw x0)
    return r

-- | a slower version of 'solveNLP'' that uses 'createIpoptProblemADSparse'
solveNLP'sparse :: (VG.Vector v Double, MonadIO m) =>
    (IpProblem -> IO ()) -- ^ set ipopt options (using functions from "Ipopt.Raw") or the 'ipopts' quasiquoter
    -> NLPT m (IpOptSolved v)
solveNLP'sparse setOpts = do
    (xl,xu) <- join $ uses (nlpfun . boundX) seqToVecs
    (gl,gu) <- join $ uses (nlpfun . boundG) seqToVecs
    AnyRF fs <- use (nlpfun . funF)
    AnyRF gs <- use (nlpfun . funG)

    p <- liftIO (createIpoptProblemADSparse xl xu gl gu (F.sum . fs) (V.fromList . toList . gs))
    liftIO (setOpts p)
    x0 <- uses initX V.convert
    r <- liftIO (ipoptSolve p =<< VS.thaw x0)
    return r


-- | add a constraint
addG :: Monad m
    => Maybe String -- ^ optional description
    -> (Double,Double) -- ^ bounds @(gl,gu)@ for the single inequality @gl_i <= g_i(x) <= gu_i@
    -> AnyRF Identity -- ^ @g_i(x)@
    -> NLPT m ()
addG d b (AnyRF f) = do
    nlpfun . boundG %= (Seq.|> b)
    nlpfun . funG %= \(AnyRF fs) -> AnyRF $ \x -> fs x Seq.|> runIdentity (f x)
    n <- use (nlpfun . boundG . to Seq.length)
    copyEnv constraintEnv n
    join $ uses defaultConstraintTol $ \t -> constraintTol %= (<> Seq.singleton t)
    addDesc constraintLabels d n

{- | add a piece of the objective function, which is added in the form
`f_1 + f_2 + ...`, to make it easier to understand (at some point)
which components are responsible for the majority of the cost, and
which are irrelevant.
-}
addF :: Monad m
    => Maybe String -- ^ description
    -> AnyRF Identity -- ^ `f_i(x)`
    -> NLPT m ()
addF d (AnyRF f) = do
    nlpfun . funF %= \(AnyRF fs) -> AnyRF $ \x -> fs x Seq.|> runIdentity (f x)
    n <- use (objEnv . to ((+1) . IM.size))
    copyEnv objEnv n
    addDesc objLabels d n

-- | add a variable, or get a reference to the the same variable if it has
-- already been used
var' :: (Monad m, Functor m)
    => Maybe (Double,Double) -- ^ bounds @(xl,xu)@ to request that @xl <= x <= xu@.
                             -- if Nothing, you get whatever is in 'defaultBounds'
    -> Maybe String -- ^ optional longer description
    -> String -- ^ variable name (namespace from the 'pushEnv' / 'popEnv' can
              -- make an @"x"@ you request here different from one you
              -- previously requested
    -> NLPT m Ix -- ^ the index (into the rawvector of variables that the solver sees)
var' bs longDescription s = do
    ev <- use currentEnv
    m <- use variables
    let s' = intercalate "." (reverse (s:ev))
    n <- case M.lookup s' m of
        Nothing -> do
            nMax %= over varIx (+1)
            db <- use defaultBounds
            nlpfun . boundX %= (Seq.|> db)
            n' <- use nMax
            variables %= M.insert s' n'
            variablesInv . at_ n' .= Just s'
            F.for_ longDescription $ \d -> varLabels . at_ n' %= (<> Just d)

            -- try to find a sane initial value
            initX %= let x0 = fromMaybe 0 $
                            bs >>= \(a,b) -> find valid [(a+b)/2, a, b]
                         valid x = not $ isInfinite x || isNaN x
                    in (`V.snoc` x0)
            return n'
        Just n -> return n
    varEnv . at_ n %= (<> Just (S.singleton ev))
    F.for_ longDescription $ \str -> varLabels . at_ n %= (<> Just str)

    F.traverse_ (narrowBounds n) bs
    return n

-- | a combination of 'var'' and 'ixToVar'
var bs s = ixToVar <$> var' bs Nothing s

{- | 'var', except this causes the solver to get a new variable,
so that you can use:

> [a,b,c,d,e] <- replicateM 5 (varFresh (Just (0, 10)) "x")

and the different letters can take different values (between 0 and 10)
in the optimal solution (depending on what you do with @a@ and similar
in the objective function and other constraints).
-}
varFresh' :: (MonadFail m, Monad m, Functor m) => Maybe (Double,Double) -> String -> NLPT m Ix
varFresh' bs s = do
    existing <- gets (^? variables . ix s)
    case existing of
        Just _ -> do
            m <- use variables
            let n = M.size m + 1
            -- get the first of "x", "x1", "x1_", "x1__", "x1___"
            Just sUniq <- return $ find (`M.notMember` m) $ s : iterate (++"_") (s ++ show n)
            var' bs Nothing sUniq
        Nothing -> var' bs Nothing s

-- | see 'varFresh''
varFresh bs s = fmap ixToVar $ varFresh' bs s

-- *** namespace
{- $namespace

When you build up an optimization problem, it may be composed of pieces.
Functions in this section help to ease the pain of making unique variables.
To illustrate:

> m <- inEnv "A" (var b "x")
> n <- var b "A.x" 

@m@ and @n@ above should refer to the same variable. In some sense this
is \"better\" that using 'varFresh' all the time, since perhaps you would
like to add dependencies between components (say the size of a header pipe,
refridgeration unit, foundation etc. has to satisfy sizes of individual
components).

-}

-- | combination of 'pushEnv' and 'popEnv'
inEnv :: (MonadFail m, Monad m) => String -> NLPT m a -> NLPT m a
inEnv s action = do
    pushEnv s
    r <- action
    popEnv
    return r

pushEnv :: Monad m => String -> NLPT m ()
pushEnv s = currentEnv %= (s:)

popEnv :: (Monad m, MonadFail m) => NLPT m String
popEnv = do
    n : ns <- use currentEnv
    currentEnv .= ns
    return n

-- *** piecewise linear

-- $piecewise
-- see for example chapter 20 of <http://www.ampl.com/BOOK/download.html>
-- and use of the splines package in @examples\/Test4@ and @examples\/Test5@

{- | splits a variable @x@ into two positive variables such that
@x = x^+ - x^-@ the new variables represent the positive and negative
parts of @x - b@

> (xMinus, xPlus) <- splitVar b x

Using @max (x-b) 0@ instead of xPlus (ie. not telling the solver that @b@ is
a special point) seems to work just as well: additional special treatment
is needed. For example see chapter 11 of 

> Nonlinear Programming: Concepts, Algorithms, and Applications to Chemical Processes
> Lorenz T. Biegler 
> SIAM 2010

which discusses several ways to reformulate the problem so that
an ordinary NLP solver will not have trouble with the fact that one of
the pair of constraints (@x+ = 0  | x- = 0@) is tight at an optimum.

-}
splitVar :: (Monad m, Functor m, MonadFail m)
    => Double -- ^ @b@
    -> Ix -- ^ index for @x@
    -> NLPT m (AnyRF Identity, AnyRF Identity) -- ^ @(b-x)_+, (x-b)_+@
splitVar b i = do
    -- need to have a variable name... 
    Just s <- gets (^? variablesInv . ix_ i)
    let x = ixToVar i
    xPlus  <- varFresh' (Just (0, 1/0)) (s ++ "+")
    xMinus <- varFresh' (Just (0, 1/0)) (s ++ "-")
    x0 <- uses initX (! view varIx i)
    initX %= (V.// [(view varIx xPlus, max 0 (x0 - b)), (view varIx xMinus, max 0 ( b - x0) )])
    addG (Just ("splitVar: " ++ s)) (b, b) (x + ixToVar xMinus - ixToVar xPlus)
    return (ixToVar xMinus, ixToVar xPlus)

ixToVar :: Ix -> AnyRF Identity
ixToVar (Ix i) = AnyRF (\v -> Identity (v V.! i))

-- *** bounds

-- | override bounds. Should be unnecessary given 'var' takes bounds.
setBounds :: Monad m => Ix -> (Double,Double) -> NLPT m ()
setBounds (Ix i) bs = nlpfun . boundX %= Seq.update i bs

-- | shrink the interval in which that variable is allowed.
narrowBounds :: Monad m => Ix -> (Double,Double) -> NLPT m ()
narrowBounds (Ix i) (a,b) = nlpfun . boundX . ix i %= \(a',b') -> (max a a', min b b')

-- * internal
seqToVecs :: MonadIO m => Seq (Double,Double) -> m (Vec,Vec)
seqToVecs x = let (a,b) = unzip (toList x) in liftIO $ do
    a' <- VS.thaw (VS.fromList a)
    b' <- VS.thaw (VS.fromList b)
    return (a',b')
