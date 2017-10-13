{-# LANGUAGE FlexibleContexts#-}
{-# LANGUAGE DeriveDataTypeable #-}

module PseudoWigner where

import Hilbert
import qualified Data.Array.Accelerate as A
import Data.Array.Accelerate.Array.Sugar as S 
--import qualified Data.Array.Accelerate.Data.Monoid as AM
import qualified Data.Array.Accelerate.Math.FFT as AMF
import qualified Data.Array.Accelerate.Data.Complex as ADC
import Data.Data
import Data.Typeable


data WindowFunc = Rect | Sin | Lanczos | Hanning | Hamming | Bartlett 
  deriving (Read, Show, Data, Typeable)
 
makeWindow :: (A.RealFloat e, Fractional (A.Exp e), Floating (A.Exp e), A.IsFloating e, A.FromIntegral Int e, Ord e) => 
  WindowFunc -> A.Acc (Scalar Int) -> A.Acc (A.Array A.DIM1 e)
makeWindow func leng = 
  let gen = A.generate (A.index1 $ A.the leng)
  in case func of 
       Rect -> A.fill (A.index1 $ A.the leng) 1.0
       Sin  -> gen (\sh -> let (A.Z A.:.x) = A.unlift sh in sin (pi*(A.fromIntegral x)/(A.fromIntegral $ A.the leng - 1)))
       Lanczos -> gen (\sh -> let (A.Z A.:.x) = A.unlift sh in sinc ((2*(A.fromIntegral x)/(A.fromIntegral $ A.the leng - 1)) - 1.0)) 
       Hanning -> gen (\sh -> let (A.Z A.:.x) = A.unlift sh in 0.5 - (0.5 * (cos (2*pi*(A.fromIntegral (x + 1))/(A.fromIntegral $ A.the leng + 1)))))
       Hamming -> gen (\sh -> let (A.Z A.:.x) = A.unlift sh in 0.54 - (0.46 * (cos (2*pi*(A.fromIntegral (x + 1))/(A.fromIntegral $ A.the leng + 1)))))
       Bartlett -> gen (\sh -> let (A.Z A.:.x) = A.unlift sh in 1.0 - A.abs (((A.fromIntegral x)/(A.fromIntegral (A.the leng - 1)/2.0)) - 1.0))
 
pWignerVille :: (A.RealFloat e, A.IsFloating e, A.FromIntegral Int e, Elt e) => 
  (A.Acc (A.Array A.DIM1 e), A.Acc (A.Array A.DIM1 (ADC.Complex e))) -> A.Acc (A.Array A.DIM2 e)
pWignerVille (window, arr) = 
  let times = A.enumFromN (A.index1 leng) 0 :: A.Acc (Array DIM1 Int)
      leng = A.length arr
      taumx = taumaxs times window
      lims = limits taumx
  in A.map ADC.real $ A.transpose $ AMF.fft1D_m' AMF.Forward (A.Z A.:. 16 A.:. 16) $ createMatrix arr window taumx lims

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
{-
-- | Compute sinc function @sin(x)\/x@
sinc :: (Ord e, Floating e) => e -> e
sinc x
  | ax < eps_0 = 1
  | ax < eps_2 = 1 - x2/6
  | ax < eps_4 = 1 - x2/6 + x2*x2/120
  | otherwise  = sin x / x
  where
    ax    = abs x
    x2    = x*x
    -- For explanation of choice see `doc/sinc.hs'
    eps_0 = 1.8250120749944284e-8 -- sqrt (6ε/4)
    eps_2 = 1.4284346431400855e-4 --   (30ε)**(1/4) / 2
    eps_4 = 4.043633626430947e-3  -- (1206ε)**(1/6) / 2 -}