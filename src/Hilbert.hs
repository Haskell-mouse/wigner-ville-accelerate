{-# LANGUAGE FlexibleContexts #-}

module Hilbert where

import qualified Data.Array.Accelerate as A
import Data.Array.Accelerate.Array.Sugar as S 
import qualified Data.Array.Accelerate.Math.FFT as AMF
import qualified Data.Array.Accelerate.Data.Complex as ADC

{- hylbert :: A.Acc (A.Array A.DIM1 Double) -> A.Acc (A.Array A.DIM1 (ADC.Complex Double))
hylbert arr = 
  let hVect = h $ A.unit $ A.length arr
  in  inverseFFT (loadArrayFFt arr) hVect

-- | Scalar myltiplies our vector with h vector and make inverse FFT 

inverseFFT :: A.Acc (A.Array A.DIM1 (ADC.Complex Double)) -> A.Acc (A.Array A.DIM1 Double) -> A.Acc (A.Array A.DIM1 (ADC.Complex Double))
inverseFFT arr h = AMF.fft1D' AMF.Inverse (A.Z A.:. 10) $ A.zipWith (A.*) arr (makeComplex h)

-- | Load vector to GPU, make it complex and apply FFT

loadArrayFFt :: A.Acc (A.Array A.DIM1 Double) -> A.Acc (A.Array A.DIM1 (ADC.Complex Double))
loadArrayFFt arr = AMF.fft1D' AMF.Forward (A.Z A.:. 10) $ makeComplex arr

-- | Form Vector that will be scalar multiplied with our spectre vector.

h :: A.Acc (Scalar Int) -> A.Acc (A.Array A.DIM1 Double)
h size | (even size) && (size > 0) = A.generate (A.index1 size) (\ix -> let Z :. x = A.unlift ix in defEven x)
  where 
  	defEven x = A.caseof x [(\y ->(y A.== (A.constant 0)) A.|| (y A.== (A.constant (size `A.div` 2))), A.constant 1),
  	  (\y -> (x A.<= (A.constant (size `div` 2))), A.constant 2)] (A.constant 0)
h size | otherwise                 = A.generate (A.index1 (A.constant size)) (\ix -> let Z :. x = A.unlift ix in defOdd x)
  where
    defOdd x = A.caseof x [((\y -> (y A.== (A.constant 0))), A.constant 1),
      (\y -> (y A.>= (A.constant ((size+1) `div` 2))), A.constant 0)] (A.constant 2) -}

-- | Make our array complex with zero imaginary part 

makeComplex :: (A.Shape sh, Floating e, Floating (A.Exp e), Elt e) => A.Acc (A.Array sh e) -> A.Acc (A.Array sh (ADC.Complex e))
makeComplex = A.map ((flip ADC.mkPolar) 0.0)
