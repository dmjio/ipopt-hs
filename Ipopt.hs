-- | Description: exports most things you should need
--
-- This module exports most things you should need.
-- Also take a look at "Ipopt.NLP" and "Ipopt.Raw" and @examples/@
module Ipopt (
  -- * high-level
  -- ** variables
  var', var, varFresh', varFresh,
  AnyRF(..), Identity(..),
  -- ** functions
  addG, addF,
  -- ** running the solver
  ppSoln,
  NLPT, nlpstate0,
  module Control.Monad.State,
  solveNLP',
  -- *** solver options
  ipopts,
  setIpoptProblemScaling,
  openIpoptOutputFile,
  setIntermediateCallback,
  IntermediateCB,

  -- * low-level bits still needed
  IpOptSolved(..),
  ApplicationReturnStatus(..),
  -- ** types
  Vec, IpNumber,
) where

import Ipopt.PP
import Ipopt.Options
import Ipopt.NLP
import Ipopt.Raw
import Ipopt.AnyRF
import Control.Monad.State
import Control.Monad.Identity (Identity(..))

