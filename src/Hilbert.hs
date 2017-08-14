module Hilbert where

import qualified Data.Array.Accelerate as A
import Data.Array.Accelerate.Array.Sugar as S 
import qualified Data.Array.Accelerate.Math.FFT as AMF
import qualified Data.Array.Accelerate.Data.Complex as ADC

ezero :: A.Exp Double
ezero = A.constant 0

eone :: A.Exp Double
eone = A.constant 1 

etwo :: A.Exp Double 
etwo = A.constant 2

hilbert :: A.Acc (A.Array A.DIM1 Double) -> A.Acc (A.Array A.DIM1 (ADC.Complex Double))
hilbert arr = 
  let leng = A.length arr 
      hVect = h (A.unit leng)
  in  inverseFFT (loadArrayFFt arr) hVect

-- | Scalar myltiplies our vector with h vector and make inverse FFT 

inverseFFT :: A.Acc (A.Array A.DIM1 (ADC.Complex Double)) -> A.Acc (A.Array A.DIM1 Double) -> A.Acc (A.Array A.DIM1 (ADC.Complex Double))
inverseFFT arr h = AMF.fft1D' AMF.Inverse (A.Z A.:. 0) $ A.zipWith (A.*) arr (makeComplex h)

-- | Load vector to GPU, make it complex and apply FFT

applyFFt :: A.Acc (A.Array A.DIM1 Double) -> A.Acc (A.Array A.DIM1 (ADC.Complex Double))
applyFFt arr = AMF.fft1D' AMF.Forward (A.Z A.:. 0) $ makeComplex arr

-- | Form Vector that will be scalar multiplied with our spectre vector.

h :: A.Acc (A.Scalar Int) -> A.Acc (A.Array A.DIM1 Double)
h s = A.generate (A.index1 size) (\ix -> let Z :. x = A.unlift ix in def (A.fromIntegral x))
  where 
    size = A.the s 
    dsize = A.fromIntegral size
    def x = A.ifThenElse (A.even size) (defEven x) (defOdd x)
    defEven x = 
      A.caseof x [
      (\y ->(y A.== ezero) A.|| (y A.== (dsize/etwo)), eone),
  	  (\y -> (x A.<= (dsize/etwo)), etwo)] ezero
    defOdd x = A.caseof x [((\y -> (y A.== ezero)), eone),
      (\y -> (y A.>= ((dsize A.+ eone)/etwo)), eone)] etwo

-- | Make our array complex with zero imaginary part 

makeComplex :: A.Acc (A.Array A.DIM1 Double) -> A.Acc (A.Array A.DIM1 (ADC.Complex Double))
makeComplex = A.map ((flip ADC.mkPolar) 0.0)