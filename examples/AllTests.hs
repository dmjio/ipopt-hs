import Ipopt.PP
import HS71manual
import HS71ad
import HS71nlpMonad
import HS71adN
import HS71nlpMonadN
-- import HS71c
import Spline1
import System.Environment

import Text.Printf
import qualified Data.Vector.Storable as V
import Text.PrettyPrint.ANSI.Leijen (putDoc)
import Text.Read
import Control.Lens
import Nlopt.Raw

ppError v x = do
    let x_official = V.fromList [1, 4.74299964, 3.82114998, 1.37940829]

    let sse :: Double
        sse = V.sum $ V.zipWith (\a b -> (a-b)^2) (x^.v) x_official
    printf "\n||x - x_official|| = %f\n" sse

main = do
    as <- getArgs
    case as of
        [] -> error "run like: ipopt-hs_Tests 1\nipopt-hs_Tests all"
        [a] | Just n <- readMaybe a,
              n < length examples -> examples !! n
        ["all"] -> sequence_ examples

-- probably should get autogenerated if lots of examples get generated...
examples = 
 [ppError x =<< hs71manual
 ,ppError x =<< hs71ad
 ,ppError (_1.x) =<< hs71nlpMonad
 ,ppError _1 =<< hs71adN NLOPT_LD_SLSQP
 ,ppError (_1.to V.convert) =<< hs71nlpMonadN NLOPT_LD_SLSQP
 ,spline1
 ]
