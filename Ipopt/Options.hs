{-# LANGUAGE TemplateHaskell #-}
-- | Description: a quasiquote for the solver's options
module Ipopt.Options where

import Text.ParserCombinators.UU
import Text.ParserCombinators.UU.Utils
import Text.ParserCombinators.UU.BasicInstances hiding (Parser)

import Ipopt.Raw
import Language.Haskell.TH.Quote
import Language.Haskell.TH
import Language.Haskell.TH.Syntax
import Control.Monad
import Data.Char
import qualified Data.Map as M

{- | an expression-only quasiquote intended to be the second argument for
'Ipopt.NLP.solveNLP'' (so @solveNLP [ipopts| tol = 1e-3 |]@). This is a shortcut
for calling 'addIpoptNumOption' 'addIpoptStrOption' or 'addIpoptIntOption'. 

Refer to <http://www.coin-or.org/Ipopt/documentation/node39.html ipopt's options reference>,
the syntax is like @option_name = value;option2 = value2@. The semicolons are optional and
whitespace alone can separate the option name from the value. A few examples of the parser:

>>> :set -XQuasiQuotes
>>> let f = [ipopts| tol = 1e-3; print_level = 0 |]
>>> :t f
f :: IpProblem -> IO ()

>>> let p x = uuParseTest parseOpts x
>>> p "tol = 3"
([("tol",ANum 3.0)],[])

>>> p "tol = 3 tol = 4" -- the last one wins. No warnings done (yet).
([("tol",ANum 3.0),("tol",ANum 4.0)],[])

>>> p "tol = 3\n\ntol = 4"
([("tol",ANum 3.0),("tol",ANum 4.0)],[])

>>> p "acceptable_iter = 25; output_file = foobar" -- quotation marks optional
([("acceptable_iter",AInt 25),("output_file",AStr "foobar")],[])

>>> p "output_file = \"foo bar\"" -- but needed here
([("output_file",AStr "foo bar")],[])

>>> putStrLn $ (++"...") $ take 100 $ show $ p "atol = 1.8" -- typo gets corrected
([("tol",ANum 1.8)],[--    Deleted   'a' at position LineColPos 0 0 0 expecting one of [Whitespace, ...

>>> p "tol = $xY" -- interpolating haskell variables
([("tol",AVar OptionNum "xY")],[])

-}
ipopts :: QuasiQuoter
ipopts = QuasiQuoter { quoteExp = \s -> do
    Loc _filename _pkg _mod (lStart,cStart) (lEnd,cEnd) <- location
    let (settings, errs) = parse ((,) <$> parseOpts <*> pEnd)
                        (createStr (LineColPos 0 lStart cStart) s)
    unless (null errs) $ reportWarning (unlines (map show errs))
    [| \c -> $( foldr
                (\(n,v) next -> [| do $(addIpoptOpt 'c n v); $next |])
                [| return () |]
                settings)
        |],
 quotePat = error "ipopts",
 quoteDec = error "ipopts",
 quoteType = error "ipopts"
 }

uuParseTest p x = parse ((,) <$> p <*> pEnd)
                        (createStr (LineColPos 0 0 0) x)

addIpoptOpt :: Name -> String -> OptionVal -> ExpQ
addIpoptOpt c n (AStr s) = [| addIpoptStrOption $(varE c) n s |]
addIpoptOpt c n (AInt i) = [| addIpoptIntOption $(varE c) n (i :: Int) |]
addIpoptOpt c n (ANum i) = [| addIpoptNumOption $(varE c) n $(liftDouble i) |]
addIpoptOpt c n (AVar OptionNum  x) = [| addIpoptNumOption $(varE c) n ($(dyn x)::Double) |]
addIpoptOpt c n (AVar OptionInt  x) = [| addIpoptIntOption $(varE c) n ($(dyn x)::Int) |]
addIpoptOpt c n (AVar OptionBool x) = [| addIpoptStrOption $(varE c) n
                                            (if $(dyn x) then "yes" else "no") |]
addIpoptOpt c n (AVar OptionStr x) = [| addIpoptStrOption $(varE c) n $(dyn x) |]

-- | what 'lift' should do for Double
liftDouble :: Double -> ExpQ
liftDouble n | isNaN n = [| 0/0 :: Double |]
    | isInfinite n, n > 0 = [| 1/0 :: Double |]
    | isInfinite n, n < 0 = [| - 1 / 0 :: Double |]
    | otherwise = [| $(litE (rationalL (toRational n))) :: Double |]

parseOpts = pSpaces *> pListSep (optional (pSymbol ";") *> pSpaces) (parseOpt <* pSpaces)

parseOpt =
    pAny (\(a,b) -> (,) <$>
                pSymbol a <*
                pSymbol "=" <*>
                (parseVar b <<|> parseLit b ))
        (M.toList ipoptOptions)
 where
  parseVar b = AVar b <$> (pSym '$' *> hsIdent)

  -- [_a-z][A-Za-z_']*
  hsIdent = (:) <$> (pRange ('a','z') <|> pSym '_')
                <*> pMany (pRange ('A','Z') <|>
                            pRange ('a','z') <|>
                            pSym '_' <|>
                            pSym '\'')
  parseLit b = case b of
    OptionNum  -> ANum <$> pDouble
    OptionBool -> AStr <$> pAny pSymbol ["yes","no"]
    OptionStr  -> AStr <$> (strLit <<|> pMany notSpaceOrSemicolonAscii)
    OptionInt  -> AInt <$> pInteger

  notSpaceOrSemicolonAscii = pRange ('!',':') <|> pRange ('<','~')
  
  strLit = pSym '"' *>
      pMany
          (pSatisfy (/= '"')
                      (Insertion "\\\"" '"' 0))
      <* pSym '"'

data OptionVal = ANum Double | AStr String | AInt Int
               | AVar OptionType String -- ^ @$x@
    deriving (Show)
data OptionType = OptionNum
                | OptionStr
                | OptionInt
                | OptionBool -- ^ actually string yes or string no
    deriving (Show,Eq)

-- | a list of all the options in
-- <http://www.coin-or.org/Ipopt/documentation/node39.html>
ipoptOptions :: M.Map String OptionType
ipoptOptions = M.fromList [
  ("print_level",                            OptionInt),
  ("print_user_options",                     OptionBool),
  ("print_options_documentation",            OptionBool),
  ("print_frequency_iter",                   OptionInt),
  ("print_frequency_time",                   OptionNum),
  ("output_file",                            OptionStr),
  ("file_print_level",                       OptionInt),
  ("option_file_name",                       OptionStr),
  ("print_info_string",                      OptionBool),
  ("inf_pr_output",                          OptionStr),
  ("print_timing_statistics" ,               OptionBool),
  -- tolerances
  ("tol",                                    OptionNum),
  ("max_iter",                               OptionInt),
  ("max_cpu_time",                           OptionNum),
  ("dual_inf_tol",                           OptionNum),
  ("constr_viol_tol",                        OptionNum),
  ("compl_inf_tol",                          OptionNum),
  ("acceptable_tol",                         OptionNum),
  ("acceptable_iter",                        OptionInt),
  ("acceptable_constr_viol_tol",             OptionNum),
  ("acceptable_dual_inf_tol",                OptionNum),
  ("acceptable_compl_inf_tol",               OptionNum),
  ("acceptable_obj_change_tol",              OptionNum),
  ("diverging_iterates_tol",                 OptionNum),

  -- scaling
  ("obj_scaling_factor",                     OptionNum),
  ("nlp_scaling_method",                     OptionStr),
  ("nlp_scaling_max_gradient",               OptionNum),
  ("nlp_scaling_min_value",                  OptionNum),

  -- NLP
  ("bound_relax_factor",                     OptionNum),
  ("honor_original_bounds",                  OptionBool),
  ("check_derivatives_for_naninf",           OptionBool),
  ("nlp_lower_bound_inf",                    OptionNum),
  ("nlp_upper_bound_inf",                    OptionNum),
  ("fixed_variable_treatment",               OptionStr),
  ("jac_c_constant",                         OptionBool),
  ("jac_d_constant",                         OptionBool),
  ("hessian_constant",                       OptionBool),


  -- initialization
  ("bound_frac",                             OptionNum),
  ("bound_push",                             OptionNum),
  ("slack_bound_frac",                       OptionNum),
  ("slack_bound_push",                       OptionNum),
  ("bound_mult_init_val",                    OptionNum),
  ("constr_mult_init_max",                   OptionNum),
  ("bound_mult_init_method",                 OptionStr),

  -- barrier parameter
  ("mehrotra_algorithm",                     OptionBool),
  ("mu_strategy",                            OptionStr),
  ("mu_oracle",                              OptionStr),
  ("quality_function_max_section_steps",     OptionInt),
  ("fixed_mu_oracle",                        OptionStr),
  ("adaptive_mu_globalization",              OptionStr),
  ("mu_init",                                OptionNum),
  ("mu_max_fact",                            OptionNum),
  ("mu_max",                                 OptionNum),
  ("mu_min",                                 OptionNum),
  ("mu_target",                              OptionNum),
  ("barrier_tol_factor",                     OptionNum),
  ("mu_linear_decrease_factor",              OptionNum),
  ("mu_superlinear_decrease_power",          OptionNum),


  -- multiplier updates
  ("alpha_for_y",                            OptionStr),
  ("alpha_for_y_tol",                        OptionNum),
  ("recalc_y",                               OptionBool),
  ("recalc_y_feas_tol",                      OptionNum),

  -- line search
  ("max_soc",                                OptionInt),
  ("watchdog_shortened_iter_trigger",        OptionInt),
  ("watchdog_trial_iter_max",                OptionInt),
  ("accept_every_trial_step",                OptionBool),
  ("corrector_type",                         OptionStr),

  -- warm start
  ("warm_start_init_point",                  OptionBool),
  ("warm_start_bound_push",                  OptionNum),
  ("warm_start_bound_frac",                  OptionNum),
  ("warm_start_slack_bound_frac",            OptionNum),
  ("warm_start_slack_bound_push",            OptionNum),
  ("warm_start_mult_bound_push",             OptionNum),
  ("warm_start_mult_init_max",               OptionNum),

  -- restoration phase
  ("expect_infeasible_problem",              OptionBool),
  ("expect_infeasible_problem_ctol",         OptionNum),
  ("expect_infeasible_problem_ytol",         OptionNum),
  ("start_with_resto",                       OptionBool),
  ("soft_resto_pderror_reduction_factor",    OptionNum),
  ("required_infeasibility_reduction",       OptionNum),
  ("bound_mult_reset_threshold",             OptionNum),
  ("constr_mult_reset_threshold",            OptionNum),
  ("evaluate_orig_obj_at_resto_trial",       OptionBool),

  -- linear solver
  ("linear_solver",                          OptionStr),
  ("linear_system_scaling",                  OptionStr),
  ("linear_scaling_on_demand",               OptionBool),
  ("max_refinement_steps",                   OptionInt),
  ("min_refinement_steps",                   OptionInt),

  -- hessian perturbation
  ("max_hessian_perturbation",               OptionNum),
  ("min_hessian_perturbation",               OptionNum),
  ("first_hessian_perturbation",             OptionNum),
  ("perturb_inc_fact_first",                 OptionNum),
  ("perturb_inc_fact",                       OptionNum),
  ("perturb_dec_fact",                       OptionNum),
  ("jacobian_regularization_value",          OptionNum),

  -- quasi newton
  ("hessian_approximation",                  OptionStr),
  ("limited_memory_update_type",             OptionStr),
  ("limited_memory_max_history",             OptionInt),
  ("limited_memory_max_skipping",            OptionInt),
  ("limited_memory_initialization",          OptionStr),
  ("limited_memory_init_val",                OptionNum),
  ("limited_memory_init_val_max",            OptionNum),
  ("limited_memory_init_val_min",            OptionNum),
  ("limited_memory_special_for_resto",       OptionBool),

  -- derivative test
  ("derivative_test",                        OptionStr),
  ("derivative_test_perturbation",           OptionNum),
  ("derivative_test_tol",                    OptionNum),
  ("derivative_test_print_all",              OptionBool),
  ("derivative_test_first_index",            OptionInt),
  ("point_perturbation_radius",              OptionNum),

  -- ma27
  ("ma27_pivtol",                            OptionNum),
  ("ma27_pivtolmax",                         OptionNum),
  ("ma27_liw_init_factor",                   OptionNum),
  ("ma27_la_init_factor",                    OptionNum),
  ("ma27_meminc_factor",                     OptionNum),

  -- ma57
  ("ma57_pivtol",                            OptionNum),
  ("ma57_pivtolmax",                         OptionNum),
  ("ma57_pre_alloc",                         OptionNum),
  ("ma57_pivot_order",                       OptionInt),
  ("ma57_automatic_scaling",                 OptionBool),
  ("ma57_block_size",                        OptionInt),
  ("ma57_node_amalgamation",                 OptionInt),
  ("ma57_small_pivot_flag",                  OptionInt),

  -- ma77
  ("ma77_print_level",                       OptionInt),
  ("ma77_buffer_lpage",                      OptionInt),
  ("ma77_buffer_npage",                      OptionInt),
  ("ma77_file_size",                         OptionInt),
  ("ma77_maxstore",                          OptionInt),
  ("ma77_nemin",                             OptionInt),
  ("ma77_order",                             OptionStr),
  ("ma77_small",                             OptionNum),
  ("ma77_static",                            OptionNum),
  ("ma77_u",                                 OptionNum),
  ("ma77_umax",                              OptionNum),

  -- ma86
  ("ma86_print_level",                       OptionInt),
  ("ma86_nemin",                             OptionInt),
  ("ma86_order",                             OptionStr),
  ("ma86_scaling",                           OptionStr),
  ("ma86_small",                             OptionNum),
  ("ma86_static",                            OptionNum),
  ("ma86_u",                                 OptionNum),
  ("ma86_umax",                              OptionNum),

  -- ma97
  ("ma97_print_level",                       OptionInt),
  ("ma97_nemin",                             OptionInt),
  ("ma97_order",                             OptionStr),
  ("ma97_scaling",                           OptionStr),
  ("ma97_scaling1",                          OptionStr),
  ("ma97_scaling2",                          OptionStr),
  ("ma97_scaling3",                          OptionStr),
  ("ma97_small",                             OptionNum),
  ("ma97_solve_blas3",                       OptionBool),
  ("ma97_switch1",                           OptionStr),
  ("ma97_switch2",                           OptionStr),
  ("ma97_switch3",                           OptionStr),
  ("ma97_u",                                 OptionNum),
  ("ma97_umax",                              OptionNum),

  -- mumps
  ("mumps_pivtol",                           OptionNum),
  ("mumps_pivtolmax",                        OptionNum),
  ("mumps_mem_percent",                      OptionInt),
  ("mumps_permuting_scaling",                OptionInt),
  ("mumps_pivot_order",                      OptionInt),
  ("mumps_scaling",                          OptionInt),

  -- paradiso
  ("pardiso_matching_strategy",              OptionStr),
  ("pardiso_max_iterative_refinement_steps", OptionInt),
  ("pardiso_msglvl",                         OptionInt),
  ("pardiso_order",                          OptionStr),

  -- wsmp
  ("wsmp_num_threads",                       OptionInt),
  ("wsmp_ordering_option",                   OptionInt),
  ("wsmp_pivtol",                            OptionNum),
  ("wsmp_pivtolmax",                         OptionNum),
  ("wsmp_scaling",                           OptionInt),
  ("wsmp_singularity_threshold:",            OptionNum)]
