{-# LANGUAGE TemplateHaskell, NoMonomorphismRestriction, OverloadedStrings #-}
module Ipopt.PP where

import Data.List
import Data.Ord
import Ipopt.Raw
import Ipopt.NLP
import Ipopt.AnyRF
import Text.PrettyPrint.ANSI.Leijen hiding ((<>), (<$>))
import qualified Text.PrettyPrint.ANSI.Leijen as PP

import Text.Printf
import qualified Data.IntMap as IM
import qualified Data.Vector.Storable as VS
import qualified Data.Vector.Generic as VG
import Data.Vector (Vector)
import qualified Data.Vector as V
import Control.Monad.State
import Control.Monad.Writer
import Foreign.C.Types(CDouble(..))
import Control.Lens
import Data.List
import Data.Foldable (toList,for_, foldMap)

-- * lenses for IpOptSolved
status f s    = (\b -> s { _ipOptSolved_status    = b }) `fmap` f (_ipOptSolved_status s)
objective f s = (\b -> s { _ipOptSolved_objective = b }) `fmap` f (_ipOptSolved_objective s)
x f s         = (\b -> s { _ipOptSolved_x         = b }) `fmap` f (_ipOptSolved_x s)
g f s         = (\b -> s { _ipOptSolved_g         = b }) `fmap` f (_ipOptSolved_g s)
mult_g f s    = (\b -> s { _ipOptSolved_mult_g    = b }) `fmap` f (_ipOptSolved_mult_g s)
mult_x_L f s  = (\b -> s { _ipOptSolved_mult_x_L  = b }) `fmap` f (_ipOptSolved_mult_x_L s)
mult_x_U f s  = (\b -> s { _ipOptSolved_mult_x_U  = b }) `fmap` f (_ipOptSolved_mult_x_U s)

-- * pretty printing
ppSoln state0 problem = flip evalStateT state0 $ runWriterT $ do
    s <- lift problem
    st <- get

    tell (dullgreen "status: " <> (s^.status.to colorStatus))
    br
    tell $ "obj_tot" <> double (s^.objective)
    join $ uses (nlpfun . funF) $ \(AnyRF f) -> case sortBy (comparing fst) $
                                toList (f (s^.x&VG.convert)) `zip` [1 .. ] of
        [_] -> return ()
        [] -> return ()
        xs -> for_ xs $ \(x,i) -> tell $ "obj" <> int i <> colon PP.<$> double x

    for_ (st ^. variablesInv . from ixMap . to IM.toList) $ \(k,desc) -> do
        br
        tell $ string desc <> "(" <> int k <> ")" <> "=" <>
            string (printf "%.3g" (s ^?! x . ix k))
    br
    tell $ "g" PP.<$> foldMap (\e -> mempty <$$> double e) (s ^. g & VG.convert :: V.Vector Double)

    return s

-- * internal

statusOk :: ApplicationReturnStatus -> Bool
statusOk x = case x of
  SolveSucceeded -> True
  SolvedToAcceptableLevel -> True
  UserRequestedStop -> True
  FeasiblePointFound -> True
  _ -> False

colorStatus x = (if statusOk x then id else black . onred) (string (show x))

br = tell (mempty <$$> mempty)


