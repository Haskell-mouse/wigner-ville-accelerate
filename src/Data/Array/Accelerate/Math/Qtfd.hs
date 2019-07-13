{-# LANGUAGE FlexibleContexts #-}
-- |
-- Module      : Data.Array.Accelerate.Math.Qtfd
-- Copyright   : [2017] Rinat Stryungis
-- License     : BSD3
--
-- Maintainer  : Rinat Stryungis <lazybonesxp@gmail.com>
-- Stability   : experimental
-- Portability : non-portable (GHC extensions)
--
-- Computation of quadratic time-frequency distributions using the accelerate-fft library.
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
-- Both of this flags are enabled by default.

module Data.Array.Accelerate.Math.Qtfd where

import Data.Array.Accelerate as A
import qualified Data.Array.Accelerate.Data.Complex as ADC
import qualified Data.Array.Accelerate.Math.BornJordan as BJ
import qualified Data.Array.Accelerate.Math.ChoiWilliams as CW
import qualified Data.Array.Accelerate.Math.FFT as AMF
import qualified Data.Array.Accelerate.Math.ModifiedB as EMB
import qualified Data.Array.Accelerate.Math.PseudoWigner as P
import Data.Array.Accelerate.Math.Wigner'

import Data.Complex
import Data.List
import Math.Gamma
-- | Wigner-ville distribution. It takes 1D array of complex floating numbers and returns 2D array of real numbers.
--  Columns of result array represents time and rows - frequency. Frequency range is from 0 to n/4, where n is a sampling frequency.

wignerVille :: (A.RealFloat e, A.IsFloating e, A.FromIntegral Int e, Elt e, Elt (ADC.Complex e), AMF.Numeric e)
  => Acc (Array DIM1 (ADC.Complex e))  -- ^ Data array
  -> Acc (Array DIM2 e)
wignerVille arr =
  let times = A.enumFromN (A.index1 leng) 0 :: Acc (Array DIM1 Int)
      leng = A.length arr
      taumx = taumaxs times
      lims = limits taumx
  in A.map ADC.real . A.transpose . AMF.fft AMF.Forward $ createMatrix arr taumx lims

wignerVilleNew :: (A.RealFloat e, A.IsFloating e, A.FromIntegral Int e, Elt e, Elt (ADC.Complex e), AMF.Numeric e)
  => Acc (Array DIM1 (ADC.Complex e))  -- ^ Data array
  -> Acc (Array DIM2 e)
wignerVilleNew arr =
  let times = A.enumFromN (A.index1 leng) 0 :: Acc (Array DIM1 Int)
      leng = A.length arr
      taumx = taumaxs times
      lims = limits taumx
  in A.map ADC.real . A.transpose . AMF.fft AMF.Forward $ createMatrix arr taumx lims

-- | Pseudo Wigner-ville distribution.
-- It takes 1D array of complex floating numbers, window and returns 2D array of real numbers.
-- Columns of result array represents time and rows - frequency. Frequency range is from 0 to n/4, where n is a sampling frequency.

pWignerVille :: (A.RealFloat e, A.IsFloating e, A.FromIntegral Int e, Elt e, Elt (ADC.Complex e), AMF.Numeric e)
  => Acc (Array DIM1 e)                -- ^ Smoothing window. Length of it must be odd.
  -> Acc (Array DIM1 (ADC.Complex e))  -- ^ Data array
  -> Acc (Array DIM2 e)
pWignerVille window arr =
  let times = A.enumFromN (A.index1 leng) 0 :: Acc (Array DIM1 Int)
      leng = A.length arr
      taumx = P.taumaxs times window
      lims = P.limits taumx
  in A.map ADC.real $ A.transpose $ AMF.fft AMF.Forward $ P.createMatrix arr window taumx lims

-- | Choi-Williams distribution. It takes 1D array of complex floating numbers,
-- and returns 2D array of real numbers.
-- Columns of result array represents time and rows - frequency. Frequency range is from 0 to n/4, where n is a sampling frequency.

choiWilliams :: (A.RealFloat e, A.IsFloating e, A.FromIntegral Int e, Elt e, Elt (ADC.Complex e), AMF.Numeric e)
  => A.Exp e                           -- ^ sigma
  -> Maybe (Acc (Array DIM1 e), Acc (Array DIM1 e))
  -> A.Exp Int                         -- ^ Smoothing window over mu, must be odd and symmetrical
  -> A.Exp Int                         -- ^ Smoothing window over tau, must be odd and symmetrical
  -> Bool                              -- ^ If normalise
  -> Acc (Array DIM1 (ADC.Complex e))  -- ^ Data array
  -> Acc (Array DIM2 e)
choiWilliams sigma mWindowArrays uWindow nWindow normalise arr =
  let times = A.enumFromN (A.index1 leng) 0 :: Acc (Array DIM1 Int)
      leng = A.length arr
      taumx = taumaxs times
      lims = limits taumx
  in A.transpose $ A.map ((*4) . ADC.real) $
    AMF.fft AMF.Forward $
      CW.summedOverMu arr taumx lims sigma mWindowArrays uWindow nWindow normalise

bornJordan :: (A.RealFloat e, A.IsFloating e, A.FromIntegral Int e, Elt e, Elt (ADC.Complex e), AMF.Numeric e)
  => A.Exp e                           -- ^ sigma
  -> Maybe (Acc (Array DIM1 e), Acc (Array DIM1 e))
  -> A.Exp Int                         -- ^ Smoothing window over mu, must be odd and symmetrical
  -> A.Exp Int                         -- ^ Smoothing window over tau, must be odd and symmetrical
  -> Acc (Array DIM1 (ADC.Complex e))  -- ^ Data array
  -> Acc (Array DIM2 e)
bornJordan sigma mWindowArrays uWindow nWindow arr =
  let times = A.enumFromN (A.index1 leng) 0 :: Acc (Array DIM1 Int)
      leng = A.length arr
      taumx = taumaxs times
      lims = limits taumx
  in A.transpose $ A.map ((*2) . ADC.real) $
    AMF.fft AMF.Forward $
      BJ.summedOverMu arr taumx lims sigma mWindowArrays uWindow nWindow


eModifiedB :: (A.RealFloat e, A.IsFloating e, A.FromIntegral Int e, Elt e, Elt (ADC.Complex e), AMF.Numeric e)
  => A.Exp e                           -- ^ alpha
  -> (Acc (Array DIM1 e))              -- ^ gammas
  -> Maybe (Acc (Array DIM1 e), Acc (Array DIM1 e))
  -> A.Exp Int
  -> A.Exp Int
  -> Bool 
  -> Acc (Array DIM1 (ADC.Complex e))  -- ^ Data array
  -> Acc (Array DIM2 e)
eModifiedB alpha gammas mWindowArrays uWindow nWindow normalise arr =
  let times = A.enumFromN (A.index1 leng) 0 :: Acc (Array DIM1 Int)
      leng = A.length arr
      taumx = taumaxs times
      lims = limits taumx
  in A.transpose $ A.map (*2) $ A.map ADC.real $ AMF.fft AMF.Forward $ A.transpose $ EMB.summedOverMu arr taumx lims alpha gammas mWindowArrays uWindow nWindow normalise

makeGammaArray2 :: Double -> Int -> [Double]
makeGammaArray2 beta wLength =
  let dwLength = Prelude.fromIntegral wLength :: Double
      v = [(-0.5), ((-0.5) + 1.0/dwLength)..0.5]
      s = Prelude.map (\x -> beta :+ pi*x) v
      f = Prelude.map (\x -> abs $ ((magnitude $ gamma x) Prelude.^ 2)/((gamma beta) Prelude.^ 2)) s
      n = Prelude.floor $ (Prelude.fromIntegral $ Prelude.length f)/2.0
      frst = Data.List.take (n-1) f
      secnd  = Data.List.drop n f
  in secnd Prelude.++ frst
