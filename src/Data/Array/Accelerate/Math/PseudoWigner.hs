{-# LANGUAGE FlexibleContexts #-}
{-# LANGUAGE TypeFamilies #-}

-- |
-- Module      : Data.Array.Accelerate.Math.PseudoWigner
-- Copyright   : [2017] Rinat Stryungis
-- License     : BSD3
--
-- Maintainer  : Rinat Stryungis <lazybonesxp@gmail.com>
-- Stability   : experimental
-- Portability : non-portable (GHC extensions)
--
-- Computation of a PsudoWigner transform using the accelerate-fft library.
--
-- This module uses the accelerate-fft library. And the base implementation of fft 
-- uses a naive divide-and-conquer fft implementation
-- whose absolute performance is appalling. It also requires that you know on
-- the Haskell side the size of the data being transformed, and that this is
-- a power-of-two in each dimension.
--
-- For performance, compile accelerate-fft against the foreign library bindings (using any
-- number of '-fllvm-ptx', and '-fllvm-cpu' for the accelerate-llvm-ptx, and
-- accelerate-llvm-native backends, respectively), which have none of the above
-- restrictions.
-- Both this flags are enabled by default. 

module Data.Array.Accelerate.Math.PseudoWigner(pWignerVille) where

import Data.Array.Accelerate.Math.Hilbert
import Data.Array.Accelerate.Math.WindowFunc
import qualified Data.Array.Accelerate as A
import Data.Array.Accelerate.Array.Sugar as S 
import qualified Data.Array.Accelerate.Math.FFT as AMF
import qualified Data.Array.Accelerate.Data.Complex as ADC



pWignerVille :: (A.RealFloat e, A.IsFloating e, A.FromIntegral Int e, Elt e, sh ~ DIM1) => 
  sh -> A.Acc (A.Array A.DIM1 e) -> A.Acc (A.Array A.DIM1 (ADC.Complex e)) -> A.Acc (A.Array A.DIM2 e)
pWignerVille sh window arr = 
  let times = A.enumFromN (A.index1 leng) 0 :: A.Acc (Array DIM1 Int)
      leng = A.length arr
      taumx = taumaxs times window
      lims = limits taumx
  in A.map ADC.real $ A.transpose $ AMF.fft1D_2r' AMF.Forward sh $ createMatrix arr window taumx lims

taumax :: A.Exp Int -> A.Exp Int -> A.Exp Int -> A.Exp Int
taumax leng lh t = min (min (min t (leng - t - 1) ) (A.round (((A.fromIntegral leng :: A.Exp Double)/2.0) - 1))) lh

taumaxs :: (A.RealFloat e, Elt e) => 
  A.Acc (A.Array A.DIM1 Int) -> A.Acc (A.Array A.DIM1 e) -> A.Acc (A.Array A.DIM1 Int)
taumaxs times window = 
  let leng = A.length times
      lh = (A.length window - 1) `div` 2
  in A.map (taumax leng lh) times                  

times :: Elt a => A.Acc (A.Array A.DIM1 a) -> A.Acc (A.Array A.DIM1 Int)
times arr = 
  let leng = A.length arr 
  in A.enumFromN (A.index1 leng) 0 :: A.Acc (Array DIM1 Int)

limits :: A.Acc (A.Array A.DIM1 Int) -> A.Acc (A.Array A.DIM1 Int)
limits taumaxs = 
  let funk = (\x -> 2*x + 1)
  in A.map funk taumaxs

moveUp ::  A.Acc (A.Array A.DIM1 Int) -> A.Exp Int -> A.Exp DIM2 -> A.Exp DIM2
moveUp taumaxs leng sh = 
  let taum t = taumaxs A.!! t 
  in (\(x,t) -> A.index2 ((x+(taum t)) `A.mod` leng) t) $ A.unlift $ A.unindex2 sh

generateValue :: (A.RealFloat e, Elt e) => 
  A.Acc (A.Array A.DIM1 (ADC.Complex e)) -> A.Exp Int -> A.Exp Int -> A.Exp e -> A.Exp (ADC.Complex e)
generateValue arr time tau h = (makeComplex h) * (arr A.!! (time + tau)) * (ADC.conjugate $ arr A.!! (time - tau))


createMatrix :: (A.RealFloat e, Elt e) => 
  A.Acc (A.Array A.DIM1 (ADC.Complex e)) -> A.Acc (A.Array A.DIM1 e) -> A.Acc (A.Array A.DIM1 Int) -> A.Acc (A.Array A.DIM1 Int) -> A.Acc (A.Array A.DIM2 (ADC.Complex e)) 
createMatrix arr window taumaxs lims = A.backpermute (A.index2 leng leng) (moveUp taumaxs leng) raw 
  where
    raw = A.generate (A.index2 leng leng) (\sh -> let (A.Z A.:.x A.:. t) = A.unlift sh
                                                      lim = lims A.!! t
                                                      taum = taumaxs A.!! t
                                                      h = window A.!! (lh + (x - taum))
                                                  in gen x t lim taum h)
    leng = A.length arr
    lh = (A.length window - 1) `div` 2 
    gen x t lim taum h = A.cond (x A.< lim) (generateValue arr t (x - taum) h) 0

sinc :: (Floating (A.Exp e), Elt e, A.Ord e) => A.Exp e -> A.Exp e
sinc x = 
  A.cond (ax A.< eps_0) 1 (A.cond (ax A.< eps_2) (1 - x2/6) (A.cond (ax A.< eps_4) (1 - x2/6 + x2*x2/120) ((A.sin x)/x)))
  where 
    ax = A.abs x
    x2 = x*x
    eps_0 = 1.8250120749944284e-8 -- sqrt (6ε/4)
    eps_2 = 1.4284346431400855e-4 --   (30ε)**(1/4) / 2
    eps_4 = 4.043633626430947e-3  -- (1206ε)**(1/6) / 2