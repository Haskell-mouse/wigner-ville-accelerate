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

import Data.Array.Accelerate.Math.Wigner'
import qualified Data.Array.Accelerate.Math.PseudoWigner as P
import qualified Data.Array.Accelerate.Math.ChoiWilliams as CW
import qualified Data.Array.Accelerate.Math.BornJordan as BJ
import Data.Array.Accelerate as A
import qualified Data.Array.Accelerate.Math.FFT as AMF
import qualified Data.Array.Accelerate.Data.Complex as ADC

-- | Wigner-ville distribution. It takes 1D array of complex floating numbers and returns 2D array of real numbers. 
--  Columns of result array represents time and rows - frequency. Frequency range is from 0 to n/4, where n is a sampling frequency.

wignerVille :: (A.RealFloat e, A.IsFloating e, A.FromIntegral Int e, Elt e)
  => Acc (Array DIM1 (ADC.Complex e))  -- ^ Data array
  -> Acc (Array DIM2 e) 
wignerVille arr = 
  let times = A.enumFromN (A.index1 leng) 0 :: Acc (Array DIM1 Int)
      leng = A.length arr 
      taumx = taumaxs times
      lims = limits taumx
  in A.map ADC.real $ A.transpose $ AMF.fft AMF.Forward $ createMatrix arr taumx lims 

-- | Pseudo Wigner-ville distribution. 
-- It takes 1D array of complex floating numbers, window and returns 2D array of real numbers. 
-- Columns of result array represents time and rows - frequency. Frequency range is from 0 to n/4, where n is a sampling frequency.

pWignerVille :: (A.RealFloat e, A.IsFloating e, A.FromIntegral Int e, Elt e)                                   -- ^ shape of the data array. It is ignored, when compiled with Native or PTX backend. 
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

choiWilliams :: (A.RealFloat e, A.IsFloating e, A.FromIntegral Int e, Elt e)
  => Acc (Array DIM1 (ADC.Complex e))  -- ^ Data array
  -> A.Exp e                           -- ^ sigma
  -> Acc (Array DIM2 e) 
choiWilliams arr sigma = 
  let times = A.enumFromN (A.index1 leng) 0 :: Acc (Array DIM1 Int)
      leng = A.length arr 
      taumx = taumaxs times
      lims = limits taumx
  in A.transpose $ A.map (*2) $ A.map ADC.real $ AMF.fft AMF.Forward $ A.transpose $ CW.sFunc (CW.coreFunction leng sigma) (CW.amatrix arr taumx lims)

-- | Choi-Willams with smoothing window in frequency domain

choiWilliams_w :: (A.RealFloat e, A.IsFloating e, A.FromIntegral Int e, Elt e)
  => Acc (Array DIM1 e)                -- ^ Smoothing window. Length of it must be odd.
  -> A.Exp e                           -- ^ sigma
  -> Acc (Array DIM1 (ADC.Complex e))  -- ^ Data array
  -> Acc (Array DIM2 e) 
choiWilliams_w window sigma arr = 
  let times = A.enumFromN (A.index1 leng) 0 :: Acc (Array DIM1 Int)
      leng = A.length arr 
      taumx = P.taumaxs times window
      lims = P.limits taumx
  in A.transpose $ A.map (*2) $ A.map ADC.real $ AMF.fft AMF.Forward $ A.transpose $ CW.sFunc (CW.coreFunction leng sigma) (CW.amatrix_w arr taumx lims window)

bornJordan :: (A.RealFloat e, A.IsFloating e, A.FromIntegral Int e, Elt e)
  => Acc (Array DIM1 (ADC.Complex e))  -- ^ Data array
  -> A.Exp e                           -- ^ alpha
  -> Acc (Array DIM2 e) 
bornJordan arr alpha = 
  let times = A.enumFromN (A.index1 leng) 0 :: Acc (Array DIM1 Int)
      leng = A.length arr 
      taumx = taumaxs times
      lims = limits taumx
  in A.transpose $ A.map (*2) $ A.map ADC.real $ AMF.fft AMF.Forward $ A.transpose $ BJ.sFunc (BJ.coreFunction leng alpha) (CW.amatrix arr taumx lims)

bornJordan_test :: (A.RealFloat e, A.IsFloating e, A.FromIntegral Int e, Elt e)
  => Acc (Array DIM1 (ADC.Complex e))  -- ^ Data array
  -> A.Exp e                           -- ^ alpha
  -> Acc (Array DIM2 (ADC.Complex e)) 
bornJordan_test arr alpha = 
  let times = A.enumFromN (A.index1 leng) 0 :: Acc (Array DIM1 Int)
      leng = A.length arr 
      taumx = taumaxs times
      lims = limits taumx
  in BJ.sFunc (BJ.coreFunction leng alpha) (BJ.amatrix arr taumx lims)  

bornJordan_matrix :: (A.RealFloat e, A.IsFloating e, A.FromIntegral Int e, Elt e)
  => Acc (Array DIM1 (ADC.Complex e))  -- ^ Data array
  -> A.Exp e                           -- ^ alpha
  -> Acc (Array DIM3 (ADC.Complex e)) 
bornJordan_matrix arr alpha = 
  let times = A.enumFromN (A.index1 leng) 0 :: Acc (Array DIM1 Int)
      leng = A.length arr 
      taumx = taumaxs times
      lims = limits taumx
  in  BJ.amatrix arr taumx lims

bornJordan_matrix2 :: (A.RealFloat e, A.IsFloating e, A.FromIntegral Int e, Elt e)
  => Acc (Array DIM1 (ADC.Complex e))  -- ^ Data array
  -> A.Exp e                           -- ^ alpha
  -> Acc (Array DIM3 (ADC.Complex e)) 
bornJordan_matrix2 arr alpha = 
  let times = A.enumFromN (A.index1 leng) 0 :: Acc (Array DIM1 Int)
      leng = A.length arr 
      taumx = taumaxs times
      lims = limits taumx
  in  BJ.amatrix arr taumx lims