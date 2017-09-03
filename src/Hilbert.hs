{-# LANGUAGE BangPatterns #-}
{-# LANGUAGE ScopedTypeVariables #-}
{-# LANGUAGE TypeFamilies #-}
{-# LANGUAGE TypeOperators #-}
{-# LANGUAGE FlexibleContexts#-}
{-# LANGUAGE ConstraintKinds #-}
{-# LANGUAGE EmptyDataDecls #-}

module Hilbert where

import qualified Data.Array.Accelerate as A
import Data.Array.Accelerate.Array.Sugar as S 
import qualified Data.Array.Accelerate.Math.FFT as AMF
import qualified Data.Array.Accelerate.Data.Complex as ADC


hilbert :: A.Acc (A.Array A.DIM1 Float) -> A.Acc (A.Array A.DIM1 (ADC.Complex Float))
hilbert arr = 
  let leng = A.length arr 
      hVect = h (A.unit leng)
  in  inverseFFT (applyFFt arr) hVect

-- | Scalar myltiplies our vector with h vector and make inverse FFT 

inverseFFT :: A.Acc (A.Array A.DIM1 (ADC.Complex Float)) -> A.Acc (A.Array A.DIM1 Float) -> A.Acc (A.Array A.DIM1 (ADC.Complex Float))
inverseFFT arr h = AMF.fft1D' AMF.Inverse (A.Z A.:. 16) $ A.zipWith (A.*) arr (A.map makeComplex h)

-- | Load vector to GPU, make it complex and apply FFT

applyFFt :: A.Acc (A.Array A.DIM1 Float) -> A.Acc (A.Array A.DIM1 (ADC.Complex Float))
applyFFt arr = AMF.fft1D' AMF.Forward (A.Z A.:. 16) $ A.map makeComplex $ arr

-- | Form Vector that will be scalar multiplied with our spectre vector.

h :: A.Acc (A.Scalar Int) -> A.Acc (A.Array A.DIM1 Float)
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
