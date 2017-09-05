{-# LANGUAGE BangPatterns #-}
{-# LANGUAGE ScopedTypeVariables #-}
{-# LANGUAGE TypeFamilies #-}
{-# LANGUAGE TypeOperators #-}
{-# LANGUAGE FlexibleContexts#-}
{-# LANGUAGE ConstraintKinds #-}
{-# LANGUAGE EmptyDataDecls #-}

module Smoothed(WindowFunc(..), spWignerVille) where

import Hilbert
import PseudoWigner
import qualified Data.Array.Accelerate as A
import Data.Array.Accelerate.Array.Sugar as S 
import qualified Data.Array.Accelerate.Data.Fold as AF
import qualified Data.Array.Accelerate.Data.Monoid as AM
import qualified Data.Array.Accelerate.Interpreter as ALI
import qualified Data.Array.Accelerate.LLVM.Native as ALN
import qualified Data.Array.Accelerate.LLVM.PTX as ALP
import qualified Data.Array.Accelerate.Math.FFT as AMF
import qualified Data.Array.Accelerate.Data.Complex as ADC
import qualified Data.Vector.Storable as VS
import qualified Data.Array.Accelerate.IO as AI
import qualified Data.Vector as V 
import Debug.Trace

spWignerVille :: ((A.Acc (A.Array A.DIM1 Float), A.Acc (A.Array A.DIM1 Float)), A.Acc (A.Array A.DIM1 (ADC.Complex Float))) -> A.Acc (A.Array A.DIM2 Float)
spWignerVille ((hWindow, fWindow), arr) = 
  let times = A.enumFromN (A.index1 leng) 0 :: A.Acc (Array DIM1 Int)
      leng = A.length arr
      lh = (A.length hWindow - 1) `div` 2
      lg = (A.length fWindow - 1) `div` 2
      taumxs = A.map (taumax leng lh lg)  times
      lims = limits taumxs
  in A.map ADC.real $ A.transpose $ AMF.fft1D_m' AMF.Forward (A.Z A.:. 16 A.:. 16) $ createMatrix arr hWindow fWindow taumxs lims

taumax :: A.Exp Int -> A.Exp Int -> A.Exp Int -> A.Exp Int -> A.Exp Int
taumax leng lh lg t = min (min (min (t + lg) (leng - t - 1 + lg) ) (A.round (((A.fromIntegral leng :: A.Exp Double)/2.0) - 1))) lh               

times :: Elt a => A.Acc (A.Array A.DIM1 a) -> A.Acc (A.Array A.DIM1 Int)
times arr = 
  let leng = A.length arr 
  in A.enumFromN (A.index1 leng) 0 :: A.Acc (Array DIM1 Int)

limits :: A.Acc (A.Array A.DIM1 Int) -> A.Acc (A.Array A.DIM1 Int)
limits taumaxs = 
  let funk = (\taumax -> 2*taumax + 1)
  in A.map funk taumaxs

moveUp ::  A.Acc (A.Array A.DIM1 Int) -> A.Exp Int -> A.Exp DIM2 -> A.Exp DIM2
moveUp taumaxs leng sh = 
  let taum t = taumaxs A.!! t 
  in (\(x,t) -> A.index2 ((x+(taum t)) `A.mod` leng) t) $ A.unlift $ A.unindex2 sh

generateValue :: A.Acc (A.Array A.DIM1 (ADC.Complex Float)) -> A.Exp Int -> A.Exp Int -> A.Exp Float ->  A.Acc (A.Array DIM1 Int) -> A.Acc (A.Array DIM1 Float) -> A.Exp (ADC.Complex Float)
generateValue arr time tau h points g2 = (makeComplex h) * ( A.the $ (A.sum (A.zipWith (*) (A.map makeComplex g2) (A.map (\p -> (arr A.!! (time + tau -p)) * (ADC.conjugate $ arr A.!! (time - tau - p))) points))))


createMatrix :: A.Acc (A.Array A.DIM1 (ADC.Complex Float)) -> A.Acc (A.Array A.DIM1 Float) -> A.Acc (A.Array A.DIM1 Float)  -> A.Acc (A.Array A.DIM1 Int) -> A.Acc (A.Array A.DIM1 Int) -> A.Acc (A.Array A.DIM2 (ADC.Complex Float)) 
createMatrix arr hWindow gWindow taumaxs lims = A.transpose $ A.backpermute (A.index2 leng leng) (moveUp taumaxs leng) raw 
  where
    raw = A.generate (A.index2 leng leng) (\sh -> let (A.Z A.:.x A.:. t) = A.unlift sh
                                                      lim = lims A.!! t
                                                      taum = taumaxs A.!! t
                                                      h = hWindow A.!! (lh + (x - taum))
                                                      g = A.map (\p -> gWindow A.!! (lg + p)) (points (x - taum) t)
                                                      g2 = A.map (\g1 -> g1/(A.the $ A.sum g)) g
                                                  in gen x t lim taum h g2)
    leng = A.length arr
    lh = (A.length hWindow - 1) `div` 2 
    lg = (A.length gWindow - 1) `div` 2
    points x t = point lg leng t x
    gen x t lim taum h g2 = A.cond (x A.< lim) (generateValue arr t (x - taum) h (points (x - taum) t) g2) 0

point :: A.Exp Int -> A.Exp Int -> A.Exp Int -> A.Exp Int -> A.Acc (A.Array DIM1 Int)
point lg leng t tau = 
  let startVal = min lg (leng - t - tau)
      pLength = startVal + (min lg (t - tau -1))
  in A.enumFromN (A.index1 pLength) startVal

epsilon :: RealFloat a => a
epsilon = encodeFloat 1 (fromIntegral $ 1-floatDigits epsilon)

sinc :: (A.RealFloat a) => A.Exp a -> A.Exp  a
sinc x =
   if abs x >= taylor_n_bound
     then sin x / x
     else 1 - x^2/6 + x^4/120
 where
  taylor_n_bound = sqrt $ sqrt epsilon