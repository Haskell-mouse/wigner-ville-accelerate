{-# LANGUAGE FlexibleContexts#-}
{-# LANGUAGE TypeFamilies #-}

module Data.Array.Accelerate.Math.Wigner'(wignerVille) where

import Data.Array.Accelerate.Math.Hilbert
import qualified Data.Array.Accelerate as A
import Data.Array.Accelerate.Array.Sugar as S 
import qualified Data.Array.Accelerate.Math.FFT as AMF
import qualified Data.Array.Accelerate.Data.Complex as ADC

-- | Wigner-ville distribution. It takes 1D array of complex floating numbers and returns 2D array of real numbers. 
-- | Columns represents time and rows - frequency. Frequency range is from 0 to n/4, where n is a sampling frequency frequancy 

wignerVille :: (A.RealFloat e, Fractional (A.Exp e), Floating (A.Exp e), A.IsFloating e, A.FromIntegral Int e, Elt e, sh ~ DIM2) => 
  sh -> A.Acc (A.Array A.DIM1 (ADC.Complex e)) -> A.Acc (A.Array A.DIM2 e)
wignerVille sh arr = 
  let times = A.enumFromN (A.index1 leng) 0 :: A.Acc (Array DIM1 Int)
      leng = A.length arr 
      taumx = taumaxs times
      lims = limits taumx
  in A.map ADC.real $ A.transpose $ AMF.fft1D_2r' AMF.Forward sh $ createMatrix arr taumx lims 

taumax :: A.Exp Int -> A.Exp Int -> A.Exp Int
taumax leng t = min (min t (leng - t - 1) ) (A.round (((A.fromIntegral leng)/2.0) - 1.0 :: A.Exp Double))

taumaxs :: A.Acc (A.Array A.DIM1 Int) -> A.Acc (A.Array A.DIM1 Int)
taumaxs times = 
  let leng = A.length times
  in A.map (taumax leng) times                  

times :: Elt a => A.Acc (A.Array A.DIM1 a) -> A.Acc (A.Array A.DIM1 Int)
times arr = 
  let leng = A.length arr 
  in A.enumFromN (A.index1 leng) 0 :: A.Acc (Array DIM1 Int)

limits :: A.Acc (A.Array A.DIM1 Int) -> A.Acc (A.Array A.DIM1 Int)
limits taumaxs = 
  let funk = (\x -> 2*x + 1)
  in A.map funk taumaxs

moveUp :: A.Acc (A.Array A.DIM1 Int) -> A.Exp Int -> A.Exp DIM2 -> A.Exp DIM2
moveUp taumaxs leng sh = 
  let taum t = taumaxs A.!! t 
  in (\(x,t) -> A.index2 ((x+(taum t)) `A.mod` leng) t) $ A.unlift $ A.unindex2 sh

generateValue :: (A.RealFloat e, Fractional (A.Exp e), Floating (A.Exp e), A.IsFloating e, A.FromIntegral Int e, Elt e) => A.Acc (A.Array A.DIM1 (ADC.Complex e)) -> A.Exp Int -> A.Exp Int -> A.Exp (ADC.Complex e)
generateValue arr time tau = (arr A.!! (time + tau)) * (ADC.conjugate $ arr A.!! (time - tau))


createMatrix :: (A.RealFloat e, Fractional (A.Exp e), Floating (A.Exp e), A.IsFloating e, A.FromIntegral Int e, Elt e) => A.Acc (A.Array A.DIM1 (ADC.Complex e)) -> A.Acc (A.Array A.DIM1 Int) -> A.Acc (A.Array A.DIM1 Int) -> A.Acc (A.Array A.DIM2 (ADC.Complex e)) 
createMatrix arr taumaxs lims = A.transpose $ A.backpermute (A.index2 leng leng) (moveUp taumaxs leng) raw 
  where
    raw = A.generate (A.index2 leng leng) (\sh -> let (A.Z A.:.x A.:. t) = A.unlift sh
                                                      lim = lims A.!! t
                                                      taum = taumaxs A.!! t
                                                  in gen x t lim taum)
    leng = A.length arr
    gen x t lim taum = A.cond (x A.< lim) (generateValue arr t (x - taum)) 0