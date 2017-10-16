{-# LANGUAGE FlexibleContexts#-}
{-# LANGUAGE TypeFamilies #-}

-- |
-- Module      : Data.Array.Accelerate.Math.Hylbert
-- Copyright   : [2017] Rinat Stryungis
-- License     : BSD3
--
-- Maintainer  : Rinat Stryungis <lazybonesxp@gmail.com>
-- Stability   : experimental
-- Portability : non-portable (GHC extensions)
--
-- Computation of a Hilbert Transform using the accelerate-fft library.
-- It just makes fft transform, remove signal with negative frequencies and makes inverse fft.  
-- The time complexity is O(n log n) in the size of the input.
--
-- The base implementation of fft uses a naïve divide-and-conquer fft implementation
-- whose absolute performance is appalling. It also requires that you know on
-- the Haskell side the size of the data being transformed, and that this is
-- a power-of-two in each dimension.
--
-- For performance, compile accelerate-fft against the foreign library bindings (using any
-- number of '-fllvm-ptx', and '-fllvm-cpu' for the accelerate-llvm-ptx, and
-- accelerate-llvm-native backends, respectively), which have none of the above
-- restrictions.
--

module Data.Array.Accelerate.Math.Hilbert(hilbert, makeComplex) where

import qualified Data.Array.Accelerate as A
import Data.Array.Accelerate.Array.Sugar as S 
import qualified Data.Array.Accelerate.Math.FFT as AMF
import qualified Data.Array.Accelerate.Data.Complex as ADC

-- | hylbert transform. It removes a negative frequencies from the signal. 
-- The default implementation requires the array dimension to be a power of two
-- (else error).
-- The FFI-backed implementations ignore the Haskell-side size parameter (first
-- argument).

hilbert :: (A.RealFloat e, Fractional (A.Exp e), Floating (A.Exp e), A.IsFloating e, A.FromIntegral Int e, Elt e, sh ~ DIM1) => 
  sh -> A.Acc (A.Array A.DIM1 e) -> A.Acc (A.Array A.DIM1 (ADC.Complex e))
hilbert sh arr = 
  let leng = A.length arr 
      hVect = h (A.unit leng)
  in  inverseFFT sh (applyFFt sh arr) hVect

-- | Scalar myltiplies our vector with h vector and make inverse FFT 

inverseFFT :: (A.RealFloat e, Fractional (A.Exp e), Floating (A.Exp e), A.IsFloating e, A.FromIntegral Int e, Elt e, sh ~ DIM1) => 
  sh -> A.Acc (A.Array A.DIM1 (ADC.Complex e)) -> A.Acc (A.Array A.DIM1 e) -> A.Acc (A.Array A.DIM1 (ADC.Complex e))
inverseFFT sh arr h = AMF.fft1D' AMF.Inverse sh $ A.zipWith (A.*) arr (A.map makeComplex h)

-- | Load vector to GPU, make it complex and apply FFT

applyFFt :: (A.RealFloat e, Fractional (A.Exp e), Floating (A.Exp e), A.IsFloating e, A.FromIntegral Int e, Elt e, sh ~ DIM1) => 
  sh -> A.Acc (A.Array A.DIM1 e) -> A.Acc (A.Array A.DIM1 (ADC.Complex e))
applyFFt sh arr = AMF.fft1D' AMF.Forward sh $ A.map makeComplex $ arr

-- | Form Vector that will be scalar multiplied with our spectre vector.

h :: (A.RealFloat e, Fractional (A.Exp e), Floating (A.Exp e), A.IsFloating e, A.FromIntegral Int e, Elt e) =>  
  A.Acc (A.Scalar Int) -> A.Acc (A.Array A.DIM1 e)
h s = A.generate (A.index1 size) (\ix -> let A.Z A.:. x = A.unlift ix in def (A.fromIntegral x :: A.Exp Float))
  where 
    size = A.the s 
    dsize = A.fromIntegral size
    def x = A.ifThenElse (A.even size) (defEven x) (defOdd x) 
    defEven x = 
      A.caseof x [
      (\y ->(y A.== 0) A.|| (y A.== (dsize/2.0)), 1),
  	  (\y -> (y A.<= (dsize/2.0)), 2)] 0
    defOdd x = A.caseof x [((\y -> (y A.== 0)), 1),
      (\y -> (y A.< ((dsize A.+ 1.0)/2.0)), 2)] 0

makeComplex :: (Floating (A.Exp e), Elt e) => A.Exp e -> A.Exp (ADC.Complex e)
makeComplex = (flip ADC.mkPolar) 0.0