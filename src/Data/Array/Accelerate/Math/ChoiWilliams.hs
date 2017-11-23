{-# LANGUAGE FlexibleContexts #-}
-- |
-- Module      : Data.Array.Accelerate.Math.Choi-Williams
-- Copyright   : [2017] Rinat Stryungis
-- License     : BSD3
--
-- Maintainer  : Rinat Stryungis <lazybonesxp@gmail.com>
-- Stability   : experimental
-- Portability : non-portable (GHC extensions)
--
-- Computation of Choi-Williams transform using the accelerate-fft library.
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

module Data.Array.Accelerate.Math.ChoiWilliams where

import Data.Array.Accelerate.Math.Wigner'
import Data.Array.Accelerate.Math.Hilbert
import qualified Data.Array.Accelerate.Math.PseudoWigner as P
import Data.Array.Accelerate as A
import qualified Data.Array.Accelerate.Math.FFT as AMF
import qualified Data.Array.Accelerate.Data.Complex as ADC

choiWilliams :: (A.RealFloat e, A.IsFloating e, A.FromIntegral Int e, Elt e)
  => Acc (Array DIM1 (ADC.Complex e))  -- ^ Data array
  -> A.Exp e
  -> Acc (Array DIM2 e) 
choiWilliams arr sigma = 
  let times = A.enumFromN (A.index1 leng) 0 :: Acc (Array DIM1 Int)
      leng = A.length arr 
      taumx = taumaxs times
      lims = limits taumx
  in A.transpose $ A.map (*2) $ A.map ADC.real $ AMF.fft AMF.Forward $ A.transpose $ sFunc (coreFunction leng sigma) (amatrix arr taumx lims)

amatrix :: (A.RealFloat e, A.IsFloating e, A.FromIntegral Int e, Elt e)
  => Acc (Array DIM1 (ADC.Complex e))
  -> Acc (Array DIM1 Int)
  -> Acc (Array DIM1 Int)
  -> Acc (Array DIM3 (ADC.Complex e))
amatrix arr taumx lims = 
  let a = A.transpose $ createMatrix arr taumx lims
      leng = A.length arr  
  in A.replicate (A.lift $ A.Z A.:. All A.:. leng A.:. All) a

coreFunction :: (A.RealFloat e, A.IsFloating e, A.FromIntegral Int e, Elt e)
  => Exp Int
  -> Exp e 
  -> A.Acc (Array DIM3 e)
coreFunction leng sigma =  
  A.generate (A.index3 leng leng leng) (\sh -> 
                let (A.Z A.:.u A.:. l A.:. n) = A.unlift sh 
                in genCore u l n leng sigma)

genCore :: (A.RealFloat e, A.IsFloating e, A.FromIntegral Int e, Elt e) 
  => Exp Int 
  -> Exp Int 
  -> Exp Int
  -> Exp Int 
  -> Exp e
  -> Exp e
genCore n l u leng sigma = 
  let u1 = ((A.fromIntegral u) - h )
      l1 = (A.fromIntegral l) - h
      n1 = cond ((A.fromIntegral n) A.< h) (A.fromIntegral n) (A.fromIntegral (n - leng))
      h = A.fromIntegral (leng `div` 2)
  in cond (n1 A.== 0) (cond (u1 A.== l1) 1.0 0.0) ((1.0 / sqrt ((4.0*pi*n1*n1/sigma))) * exp ((-1.0)*(((u1 - l1)*(u1 - l1))/(4.0*n1*n1/sigma)) ) )

sFunc :: (A.RealFloat e, A.IsFloating e, A.FromIntegral Int e, Elt e)
  => A.Acc (Array DIM3 e) 
  -> A.Acc (Array DIM3 (ADC.Complex e))
  -> A.Acc (Array DIM2 (ADC.Complex e))
sFunc core a =  A.fold (+) 0 $ A.zipWith (*) (A.map makeComplex core) a

divp ::  
  A.Exp Int 
  -> A.Exp Int 
  -> A.Exp Int 
divp a b = cond (a `A.mod` b A.== 0) (a `A.div` b) ((a `A.div` b) + 1)