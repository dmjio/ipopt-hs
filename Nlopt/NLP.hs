module Nlopt.NLP where

import Nlopt.Raw
import Ipopt.NLP
import Ipopt.AnyRF
import qualified Data.Sequence as Seq
import qualified Data.Vector.Generic as VG
import qualified Data.Vector.Storable as VS
import qualified Data.Vector.Storable.Mutable as VM
import qualified Data.Vector as V
import Control.Monad
import Control.Monad.State
import Control.Lens
import Foreign.C
import Control.Monad.Trans
import qualified Data.Foldable as F
import Control.Lens

solveNlopt :: (VG.Vector v Double, MonadIO m)
    => NloptAlgorithm
    -> (NLOpt -> IO t) -- ^ set additional options
    -> NLPT m (v Double, Double, NloptResult) -- ^ @x,objective,exitCode@
solveNlopt alg moreOpts = do
    (xl,xu) <- join $ uses (nlpfun . boundX) seqToVecs
    gSeq <- use (nlpfun . boundG)
    let ng = Seq.length gSeq

    AnyRF fs <- use (nlpfun . funF)
    AnyRF gs <- use (nlpfun . funG)

    m <- liftIO . nloptCreate alg . (+1) =<< use (nMax . varIx)

    x0 <- join $ uses initX (liftIO . VS.thaw . VG.convert)

    gTol <- liftIO . VS.thaw . VS.fromList
         =<< uses constraintTol (F.concatMap (\(a,b) -> [a,b]))
    
    (status,obj) <- liftIO $ do
        nloptSetMinObjective m (toFuncAD (F.sum . fs))
        nloptSetLowerBounds m xl
        nloptSetUpperBounds m xu

        let toLE :: (AnyRFCxt a) => Seq.Seq a -> V.Vector a
            toLE xs = V.fromList $ concatMap (\(x,(l,u)) -> [realToFrac l-x,x- realToFrac u])
                        (F.toList xs `zip` F.toList gSeq)
        -- XXX drop constraints with bounds more than 1e20 or so? as ipopt does
        nloptAddInequalityMconstraint m (2*ng) (toFuncMAD (toLE . gs)) gTol
        moreOpts m
        nloptOptimize m x0

    x <- liftIO (VS.freeze x0)
    return (VG.convert x, obj, status)
